using LinearAlgebra
using SparseArrays
using ProgressMeter
using Parameters

# SuiteSparse's sparse `lu` (UMFPACK) is not safe to call concurrently from
# multiple Julia threads (mirrors the GMSH_LOCK precedent in io/gmsh_utils.jl).
# Serialize factorization + solve so parallel task workers don't segfault.
const SPARSE_LU_LOCK = ReentrantLock()

"""
    assemble_system_A(mdl::Stokes, cache::BasisFunctionCache)

Assembles the finite element system matrix for the velocity field using pre-computed basis functions.

# Arguments:
- `mdl::Stokes` : Material model
- `cache::BasisFunctionCache` : Pre-computed basis functions cache

# Returns:
- `K::SparseMatrixCSC{Float64,Int64}` : Sparse stiffness matrix
"""
function assemble_system_A(mdl::Stokes, cache::BasisFunctionCache)::SparseMatrixCSC{Float64,Int64}
    # unpack the model parameters to local variables
    # this is done to avoid the need to pass the model object around
    @unpack ne, ndim, nDof_u = mdl
    @unpack NodeList, IEN, ID, volume_element_shape, basis_order = mdl.mesh_u

    IEN_u_cached::Matrix{Int} = IEN
    ID_cached::Matrix{Int} = ID
    ne_cached::Int = ne
    ndim_cached::Int = ndim
    nDof_u_cached::Int = nDof_u
    element_shape_u_cached, basis_order_u_cached = volume_element_shape, basis_order

    @unpack NodeList, IEN, volume_element_shape, basis_order = mdl.mesh_x

    NodeList_x_cached::Matrix{Float64} = NodeList
    IEN_x_cached::Matrix{Int} = IEN
    element_shape_x_cached, basis_order_x_cached = volume_element_shape, basis_order

    C::Matrix{Float64} = get_c_mat(1.0,0.0,type="standard")
    IEN_u_rows::Int = size(IEN_u_cached,1)
  
    # (I,J,V) vectors for COO sparse matrix
    if nDof_u_cached == 1
        E = zeros(  Int64, ne_cached*IEN_u_rows^2)
        J = zeros(  Int64, ne_cached*IEN_u_rows^2)
        V = zeros(Float64, ne_cached*IEN_u_rows^2)
    else
        ID_rows::Int = size(ID_cached,1)
        E = zeros(  Int64, ne_cached*((ID_rows*IEN_u_rows)^2))
        J = zeros(  Int64, ne_cached*((ID_rows*IEN_u_rows)^2))
        V = zeros(Float64, ne_cached*((ID_rows*IEN_u_rows)^2))  
    end

    # Unpack pre-computed basis functions from cache
    N_u_gp, ΔN_u_gp, N_p_gp, ΔN_p_gp, N_x_gp, ΔN_x_gp, wpoints = cache.bf_volume
    
    # initialize the Jacobian matrix
    Jac = zeros(Float64, ndim_cached, ndim_cached)  # Jacobian matrix
    vol::Float64 = 0.0
    e_iter = 1:ne_cached # iterator for elements loop 
    gpiter = 1:length(wpoints) # iterator for integration loop
    
    # integration loop - use pre-computed basis functions
    for gp::Int in gpiter
        N_u = N_u_gp[gp]
        ΔN_u = ΔN_u_gp[gp]
        N_x, ΔN_x = N_x_gp[gp], ΔN_x_gp[gp]

        dNdX_u = zeros(Float64, size(N_u, 1), ndim_cached) # Gradient of basis functions
        if nDof_u_cached == 2
            B = zeros(Float64, ndim_cached*length(N_u), 3)
        elseif nDof_u_cached == 3
            B = zeros(Float64, ndim_cached*length(N_u), 6)
        end

        szB::Int = size(B, 1)  # Number of basis functions
        Ke = zeros(Float64, szB, szB)  # Element stiffness matrix
        Ke_row::Int, Ke_col::Int = szB, szB
        Ke_len::Int = length(Ke) 

        # element loop
        for e::Int in e_iter
            coords::Matrix{Float64} = NodeList_x_cached[:, IEN_x_cached[:, e]]

            mul!(Jac, coords, ΔN_x)  # Jacobian matrix [dx/dxi dx/deta; dy/dxi dy/deta]
            w::Float64 = wpoints[gp] * abs(det(Jac))
            vol = vol + w
            dNdX_u .= ΔN_u / Jac 

            if nDof_u_cached == 1
                szN::Int = size(N_u, 1)  # Number of basis functions
                # Loop between basis functions of the element
                @inbounds @fastmath for i::Int in 1:szN
                    for j::Int in 1:szN
                        inz::Int = (szN)^2 * (e - 1) + szN * (i - 1) + j  # Index for the COO sparse matrix
                        E[inz] = IEN_u_cached[i, e]  # Row index
                        J[inz] = IEN_u_cached[j, e]  # Column index
                        V[inz] += w * dot(dNdX_u[i, :], dNdX_u[j, :])  # Inner product of the gradient of the basis functions
                    end
                end
            else
                if nDof_u_cached == 2
                    B .= 0.0
                    B[1:nDof_u_cached:end, 1] = dNdX_u[:, 1]
                    B[2:nDof_u_cached:end, 2] = dNdX_u[:, 2]
                    B[1:nDof_u_cached:end, 3] = dNdX_u[:, 2]
                    B[2:nDof_u_cached:end, 3] = dNdX_u[:, 1]
                elseif nDof_u_cached == 3
                    B .= 0.0
                    B[1:nDof_u_cached:end, 1] = dNdX_u[:, 1]
                    B[2:nDof_u_cached:end, 2] = dNdX_u[:, 2]
                    B[3:nDof_u_cached:end, 3] = dNdX_u[:, 3]
                    B[2:nDof_u_cached:end, 4] = dNdX_u[:, 3]
                    B[3:nDof_u_cached:end, 4] = dNdX_u[:, 2]
                    B[1:nDof_u_cached:end, 5] = dNdX_u[:, 3]
                    B[3:nDof_u_cached:end, 5] = dNdX_u[:, 1]
                    B[1:nDof_u_cached:end, 6] = dNdX_u[:, 2]
                    B[2:nDof_u_cached:end, 6] = dNdX_u[:, 1]
                end
                mul!(Ke, 2*w*B, C*B')  # Element stiffness matrix

                # Loop between basis functions of the element
                iNodes = 1:div(Ke_row, nDof_u_cached)
                jNodes = 1:div(Ke_col, nDof_u_cached)
                iDofs = 1:ID_rows
                jDofs = 1:ID_rows
                
                @inbounds @fastmath for iNode::Int in iNodes
                    for jNode::Int in jNodes
                        for iDof::Int in iDofs
                            for jDof::Int in jDofs
                                i::Int = (iNode - 1) * nDof_u_cached + iDof
                                j::Int = (jNode - 1) * nDof_u_cached + jDof
                                inz::Int = Ke_len * (e - 1) + (iNode - 1) * nDof_u_cached * Ke_col + (jNode - 1) * nDof_u_cached^2 + (iDof - 1) * nDof_u_cached + jDof  # Index for the COO sparse matrix
                                E[inz] = ID_cached[iDof, IEN_u_cached[iNode, e]]  # Row index
                                J[inz] = ID_cached[jDof, IEN_u_cached[jNode, e]]  # Column index
                                V[inz] += Ke[i, j]
                            end
                        end
                    end
                end
            end
        end
    end
    K = sparse(E,J,V)
    return K
end

"""
    assemble_system_B(mdl::Stokes, cache::BasisFunctionCache)

Assembles the pressure-velocity coupling matrix (B matrix) for the Stokes system using pre-computed basis functions.

# Arguments:
- `mdl::Stokes` : Material model
- `cache::BasisFunctionCache` : Pre-computed basis functions cache

# Returns:
- `K::SparseMatrixCSC{Float64,Int64}` : Sparse coupling matrix relating velocity and pressure degrees of freedom
"""
function assemble_system_B(mdl::Stokes, cache::BasisFunctionCache)::SparseMatrixCSC{Float64,Int64}
    # unpack the model parameters to local variables
    # this is done to avoid the need to pass the model object around
    @unpack ne, ndim, nDof_u = mdl
    @unpack IEN, volume_element_shape, basis_order = mdl.mesh_p

    IEN_p_cached::Matrix{Int} = IEN
    element_shape_p_cached, basis_order_p_cached = volume_element_shape, basis_order

    @unpack NodeList, IEN, ID, volume_element_shape, basis_order = mdl.mesh_u

    IEN_u_cached::Matrix{Int} = IEN
    nodeList_cached::Matrix{Float64} = NodeList
    ID_cached::Matrix{Int} = ID
    ne_cached::Int = ne
    ndim_cached::Int = ndim
    nDof_u_cached::Int = nDof_u
    element_shape_u_cached, basis_order_u_cached = volume_element_shape, basis_order

    @unpack NodeList, IEN, volume_element_shape, basis_order = mdl.mesh_x

    NodeList_x_cached::Matrix{Float64} = NodeList
    IEN_x_cached::Matrix{Int} = IEN
    element_shape_x_cached, basis_order_x_cached = volume_element_shape, basis_order

    # (I,J,V) vectors for COO sparse matrix
    IEN_u_rows::Int = size(IEN_u_cached,1)
    IEN_p_rows::Int = size(IEN_p_cached,1)
    ID_u_rows::Int = size(ID_cached,1)

    if nDof_u_cached == 1
        E = zeros(  Int64, ne_cached*IEN_u_rows*IEN_p_rows)
        J = zeros(  Int64, ne_cached*IEN_u_rows*IEN_p_rows)
        V = zeros(Float64, ne_cached*IEN_u_rows*IEN_p_rows)
    else
        E = zeros(  Int64, ne_cached*(ID_u_rows*IEN_u_rows)*IEN_p_rows)
        J = zeros(  Int64, ne_cached*(ID_u_rows*IEN_u_rows)*IEN_p_rows)
        V = zeros(Float64, ne_cached*(ID_u_rows*IEN_u_rows)*IEN_p_rows)  
    end

    # Unpack pre-computed basis functions
    N_u_gp, ΔN_u_gp, N_p_gp, ΔN_p_gp, N_x_gp, ΔN_x_gp, wpoints = cache.bf_volume
    
    # initialize the Jacobian matrix
    Jac = zeros(Float64, ndim_cached, ndim_cached)  # Jacobian matrix

    e_iter = 1:ne_cached    # iterator for elements loop
    gpiter = 1:length(wpoints)  # iterator for integration loop
    
    vol = 0.0
    # integration loop - use pre-computed basis functions
    for gp::Int in gpiter
        N_u = N_u_gp[gp]
        ΔN_u = ΔN_u_gp[gp]
        N_p = N_p_gp[gp]
        ΔN_p = ΔN_p_gp[gp]
        N_x, ΔN_x = N_x_gp[gp], ΔN_x_gp[gp]
        
        dNdX_u = zeros(Float64, size(ΔN_u, 1), ndim_cached) # Gradient of basis functions
        Be = zeros(Float64, 1, 3*size(dNdX_u, 1)) # Gradient of basis functions
        colB::Int = size(Be, 2)  # Number of basis functions
        rowNp::Int = size(N_p, 1)  # Number of basis functions
        Ke = zeros(Float64, rowNp, colB)  # Element stiffness matrix
        Ke_row::Int, Ke_col::Int = rowNp, colB

        # element loop
        for e::Int in e_iter

            coords::Matrix{Float64} = NodeList_x_cached[:, IEN_x_cached[:, e]]

            mul!(Jac, coords, ΔN_x)  # Jacobian matrix [dx/dxi dx/deta; dy/dxi dy/deta]
            w::Float64 = wpoints[gp] * abs(det(Jac))
            vol = vol + w
            dNdX_u .= ΔN_u / Jac  # Compute ΔN_u * inv(Jac) via right division

            if nDof_u_cached == 1
                szN = size(N_u,1) # number of basis functions
                # loop between basis functions of the element
                @inbounds @fastmath for i::Int in 1:szN
                    for j::Int in 1:szN
                        inz::Int = (szN)^2*(e-1) + szN*(i-1) + j # index for the COO sparse matrix
                        E[inz] = IEN_u_cached[i,e] # row index 
                        J[inz] = IEN_u_cached[j,e] # column index
                        V[inz] += w*dot(dNdX_u[i,:],dNdX_u[j,:])# inner product of the gradient of the basis functions
                    end
                end
            else  
                Be .= 0.0
                Be[1:nDof_u_cached:end] = dNdX_u[:,1]
                Be[2:nDof_u_cached:end] = dNdX_u[:,2]
                Be[3:nDof_u_cached:end] = dNdX_u[:,3]
                # Ke = -w*N_p*Be # element stiffness matrix
                mul!(Ke, -w*N_p, Be) # element stiffness matrix

                # loop between basis functions of the element
                iNodes = 1:div(Ke_col,nDof_u_cached)  # column node index
                iDofs = 1:ID_u_rows # column dof index
                jNodes = 1:Ke_row # row index

                @inbounds @fastmath for iNode::Int in iNodes 
                    for jNode::Int in jNodes 
                        for iDof::Int in iDofs
                            i::Int = (iNode-1)*nDof_u_cached + iDof
                            j::Int = jNode
                            inz::Int = length(Ke)*(e-1) + (jNode-1)*Ke_col + (iNode-1)*nDof_u_cached + iDof # index for the COO sparse matrix
                            E[inz] = ID_cached[iDof,IEN_u_cached[iNode,e]] # row index 
                            J[inz] = IEN_p_cached[jNode,e] # column index
                            V[inz] += Ke[j,i] 
                        end
                    end
                end
            end
        end
    end
    K = sparse(E,J,V)
    return K
end

"""
    apply_boundary_conditions(mdl::Stokes, cache::BasisFunctionCache)

Apply the Neumann slip boundary conditions to the global stiffness matrix.

# Arguments:
- `mdl::Stokes` : Material model
- `cache::BasisFunctionCache` : Pre-computed basis functions cache

# Returns:
- `K::SparseMatrixCSC{Float64,Int64}` : Sparse stiffness matrix with boundary conditions applied
"""
function apply_boundary_conditions(mdl::Stokes, cache::BasisFunctionCache)::SparseMatrixCSC{Float64,Int64}
    @unpack ne, ndim, nDof_u = mdl
    @unpack NodeList, IEN_top, IEN_bottom, ID, volume_element_shape, basis_order = mdl.mesh_u

    NodeList_u_cached::Matrix{Float64} = NodeList
    IEN_u_top_cached::Matrix{Int} = IEN_top
    IEN_u_btm_cached::Matrix{Int} = IEN_bottom
    ID_cached::Matrix{Int} = ID
    ne_cached::Int = ne
    ndim_cached::Int = ndim
    nDof_u_cached::Int = nDof_u
    element_shape_u_cached, basis_order_u_cached = volume_element_shape, basis_order

    @unpack NodeList, IEN_top, IEN_bottom, ID, volume_element_shape, basis_order = mdl.mesh_x

    NodeList_x_cached::Matrix{Float64} = NodeList
    IEN_x_top_cached::Matrix{Int} = IEN_top
    IEN_x_btm_cached::Matrix{Int} = IEN_bottom
    element_shape_x_cached, basis_order_x_cached = volume_element_shape, basis_order

    ID_u_rows::Int = size(ID_cached,1)
    IEN_btm_rows::Int = size(IEN_u_btm_cached,1)

    ne_surface = size(IEN_u_top_cached, 2) # number of elements on the surface
    E = zeros(  Int64, ne_surface*(ID_u_rows*IEN_btm_rows)^2*2) # *2 because we have two surfaces
    J = zeros(  Int64, ne_surface*(ID_u_rows*IEN_btm_rows)^2*2) # *2 because we have two surfaces
    V = zeros(Float64, ne_surface*(ID_u_rows*IEN_btm_rows)^2*2) # *2 because we have two surfaces
    sz = size(NodeList_x_cached,2)*ID_u_rows # size of the global stiffness matrix = number of nodes x number of dofs

    # Unpack pre-computed surface basis functions
    N_u_surf_gp, ΔN_u_surf_gp, N_p_surf_gp, ΔN_p_surf_gp, N_x_surf_gp, ΔN_x_surf_gp, wpoints = cache.bf_surface
    
    e_iter = 1:size(IEN_u_top_cached, 2)   # iterator for elements loop
    gpiter = 1:length(wpoints) # iterator for integration loop

    A_btm = 0
    A_top = 0
    # integration loop
    for gp::Int in gpiter
        # Use pre-computed basis functions
        N_u_top = N_u_surf_gp[gp]
        ΔN_u_top = ΔN_u_surf_gp[gp]
        N_u_btm = N_u_top
        ΔN_u_btm = ΔN_u_top
        N_x_top = N_u_top
        ΔN_x_top = ΔN_u_top
        N_x_btm = N_u_btm
        ΔN_x_btm = ΔN_u_btm
            
        M_top = zeros(Float64, 3, ndim_cached*length(N_u_top))
        M_btm = zeros(Float64, 3, ndim_cached*length(N_u_top))
        rowM::Int = size(M_top,2)
        be_row::Int, be_col::Int = rowM, rowM
        be_top = zeros(Float64, be_col, be_row)
        be_btm = zeros(Float64, be_col, be_row)
        len_be::Int = length(be_top)

        for e::Int in e_iter
            coords_top::Matrix{Float64} = NodeList_x_cached[:,IEN_x_top_cached[:,e]]
            coords_btm::Matrix{Float64} = NodeList_x_cached[:,IEN_x_btm_cached[:,e]]

            dxdξ_top = coords_top*ΔN_x_top         # Jacobian matrix [dx/dxi dx/deta; dy/dxi dy/deta; dz/dxi dz/deta]
            dxdξ_btm = coords_btm*ΔN_x_btm         # Jacobian matrix [dx/dxi dx/deta; dy/dxi dy/deta; dz/dxi dz/deta]

            w_top::Float64 = wpoints[gp]*norm(cross(dxdξ_top[:,1],dxdξ_top[:,2]))     # weight of the quadrature point top surface
            w_btm::Float64 = wpoints[gp]*norm(cross(dxdξ_btm[:,1],dxdξ_btm[:,2]))     # weight of the quadrature point bottom surface
            
            A_top += w_top
            A_btm += w_btm

            M_top .= 0.0
            M_top[1,1:nDof_u_cached:end] = N_u_top 
            M_top[2,2:nDof_u_cached:end] = N_u_top
            M_top[3,3:nDof_u_cached:end] = zeros(Float64, size(N_u_top)) # slip boundary condition only affects the tangential components of the velocity   

            M_btm .= 0.0
            M_btm[1,1:nDof_u_cached:end] = N_u_btm
            M_btm[2,2:nDof_u_cached:end] = N_u_btm
            M_btm[3,3:nDof_u_cached:end] = zeros(Float64, size(N_u_btm)) # slip boundary condition only affects the tangential components of the velocity 

            mul!(be_btm,M_btm',M_btm)
            mul!(be_top,M_top',M_top)

            # loop between basis functions of the element
            iNodes = 1:be_row÷nDof_u_cached
            jNodes = 1:be_col÷nDof_u_cached
            iDofs = 1:ID_u_rows
            jDofs = 1:ID_u_rows

            @inbounds @fastmath for iNode::Int in iNodes
                for jNode::Int in jNodes
                    for iDof::Int in iDofs
                        for jDof::Int in jDofs
                            i::Int = (iNode-1)*nDof_u_cached + iDof
                            j::Int = (jNode-1)*nDof_u_cached + jDof
                            
                            inz_btm::Int = len_be*(e-1) + (iNode-1)*nDof_u_cached*be_col + (jNode-1)*nDof_u_cached^2 + (iDof-1)*nDof_u_cached + jDof # index for the COO sparse matrix
                            inz_top::Int = len_be*ne_surface + len_be*((e-1)) + (iNode-1)*nDof_u_cached*be_col + (jNode-1)*nDof_u_cached^2 + (iDof-1)*nDof_u_cached + jDof # index for the COO sparse matrix
                            
                            E[inz_top] = ID_cached[iDof,IEN_u_top_cached[iNode,e]]    # row index 
                            J[inz_top] = ID_cached[jDof,IEN_u_top_cached[jNode,e]]   # column index
                            V[inz_top] += w_top*be_top[i,j] 

                            E[inz_btm] = ID_cached[iDof,IEN_u_btm_cached[iNode,e]]    # row index 
                            J[inz_btm] = ID_cached[jDof,IEN_u_btm_cached[jNode,e]]   # column index
                            V[inz_btm] += w_btm*be_btm[i,j] 

                        end
                    end
                end
            end
        end
    end
    K = sparse(E,J,V,sz,sz)
    return  K
end

"""
    apply_boundary_conditions_dense(mdl::Stokes, cache::BasisFunctionCache)

Applies the Neumann slip boundary conditions to the global stiffness matrix (dense format).

# Arguments:
- `mdl::Stokes` : Material model
- `cache::BasisFunctionCache` : Pre-computed basis functions cache

# Returns:
- `K::Matrix{Float64}` : Dense stiffness matrix with boundary conditions applied
"""
function apply_boundary_conditions_dense(mdl::Stokes, cache::BasisFunctionCache)::Matrix{Float64}
    @unpack ne, ndim, nDof_u = mdl
    @unpack NodeList, IEN_top, IEN_bottom, ID, volume_element_shape, basis_order = mdl.mesh_u

    NodeList_u_cached::Matrix{Float64} = NodeList
    IEN_u_top_cached::Matrix{Int} = IEN_top
    IEN_u_btm_cached::Matrix{Int} = IEN_bottom
    ID_cached::Matrix{Int} = ID
    ne_cached::Int = ne
    ndim_cached::Int = ndim
    nDof_u_cached::Int = nDof_u
    element_shape_u_cached, basis_order_u_cached = volume_element_shape, basis_order

    @unpack NodeList, IEN_top, IEN_bottom, ID, volume_element_shape, basis_order = mdl.mesh_x

    NodeList_x_cached::Matrix{Float64} = NodeList
    IEN_x_top_cached::Matrix{Int} = IEN_top
    IEN_x_btm_cached::Matrix{Int} = IEN_bottom
    element_shape_x_cached, basis_order_x_cached = volume_element_shape, basis_order

    ID_u_rows::Int = size(ID_cached,1)
    IEN_btm_rows::Int = size(IEN_u_btm_cached,1)
    ne_surface = size(IEN_u_top_cached, 2) # number of elements on the surface
    nNodes_cached = size(NodeList_x_cached, 2)

    sz = nDof_u_cached*nNodes_cached
    K = zeros(Float64, sz,sz)

    # Unpack pre-computed surface basis functions
    N_u_surf_gp, ΔN_u_surf_gp, N_p_surf_gp, ΔN_p_surf_gp, N_x_surf_gp, ΔN_x_surf_gp, wpoints = cache.bf_surface
    
    e_iter = 1:size(IEN_u_top_cached, 2)   # iterator for elements loop
    gpiter = 1:length(wpoints) # iterator for integration loop
    A_btm = 0
    A_top = 0
    # integration loop
    for gp::Int in gpiter
        # Use pre-computed basis functions
        N_u_top = N_u_surf_gp[gp]
        ΔN_u_top = ΔN_u_surf_gp[gp]
        N_u_btm = N_u_top
        ΔN_u_btm = ΔN_u_top
        N_x_top = N_u_top
        ΔN_x_top = ΔN_u_top
        N_x_btm = N_u_btm
        ΔN_x_btm = ΔN_u_btm
        # element loop
        for e::Int in e_iter
            coords_top::Matrix{Float64} = NodeList_x_cached[:,IEN_x_top_cached[:,e]]
            coords_btm::Matrix{Float64} = NodeList_x_cached[:,IEN_x_btm_cached[:,e]]

            M_top = zeros(Float64, 3, ndim_cached*length(N_u_top))
            M_btm = zeros(Float64, 3, ndim_cached*length(N_u_top))
            rowM::Int = size(M_top,2)
            be_row::Int, be_col::Int = rowM, rowM
            be_top = zeros(Float64, be_col, be_row)
            be_btm = zeros(Float64, be_col, be_row)
            len_be::Int = length(be_top)

            dxdξ_top = coords_top*ΔN_x_top         # Jacobian matrix [dx/dxi dx/deta; dy/dxi dy/deta; dz/dxi dz/deta]
            dxdξ_btm = coords_btm*ΔN_x_btm         # Jacobian matrix [dx/dxi dx/deta; dy/dxi dy/deta; dz/dxi dz/deta]

            w_top::Float64 = wpoints[gp]*norm(cross(dxdξ_top[:,1],dxdξ_top[:,2]))     # weight of the quadrature point top surface
            w_btm::Float64 = wpoints[gp]*norm(cross(dxdξ_btm[:,1],dxdξ_btm[:,2]))     # weight of the quadrature point bottom surface
            
            A_top += w_top
            A_btm += w_btm

            M_top .= 0.0
            M_top[1,1:nDof_u_cached:end] = N_u_top
            M_top[2,2:nDof_u_cached:end] = N_u_top
            M_top[3,3:nDof_u_cached:end] = N_u_top

            M_btm .= 0.0
            M_btm[1,1:nDof_u_cached:end] = N_u_btm
            M_btm[2,2:nDof_u_cached:end] = N_u_btm
            M_btm[3,3:nDof_u_cached:end] = N_u_btm

            # be = M'*M
            mul!(be_btm,M_btm',M_btm)
            mul!(be_top,M_top',M_top)

            # loop between basis functions of the element
            iNodes = 1:be_row÷nDof_u_cached
            jNodes = 1:be_col÷nDof_u_cached
            iDofs = 1:ID_u_rows
            jDofs = 1:ID_u_rows
            for iNode in iNodes
                for jNode in jNodes
                    for iDof in iDofs
                        for jDof in jDofs
                            i = (iNode-1)*nDof_u_cached + iDof
                            j = (jNode-1)*nDof_u_cached + jDof
                            
                            K[ID_cached[iDof,IEN_u_top_cached[iNode,e]] ,ID_cached[jDof,IEN_u_top_cached[jNode,e]]] += w_top*be_top[i,j] 
                            K[ID_cached[iDof,IEN_u_btm_cached[iNode,e]] ,ID_cached[jDof,IEN_u_btm_cached[jNode,e]]] += w_btm*be_btm[i,j] 
                        end
                    end
                end
            end
        end
    end
    return  K
end

"""
    set_boundary_cond(mdl::Stokes; DENSE::Bool=false)

Set the Dirichlet boundary conditions for the problem.

# Arguments:
- `mdl::Stokes` : Material model
- `DENSE::Bool` : Flag to use dense matrices or sparse matrices

# Returns:
- `q_upper::Vector{Float64}` : Dirichlet boundary conditions upper surface
- `q_side::Vector{Float64}` : Dirichlet boundary conditions side surface
- `q_lower::Vector{Float64}` : Dirichlet boundary conditions lower surface
- `C_uc::SparseMatrixCSC{Float64,Int64}` : Constraint matrix
"""
function set_boundary_cond(mdl::Stokes; DENSE::Bool=false)

    @unpack ndim, nDof_u = mdl
    @unpack NodeList, nNodes, ID = mdl.mesh_u

    nodeList_cached::Matrix{Float64} = NodeList
    ID_cached::Matrix{Int} = ID
    ndim_cached::Int = ndim
    nDof_u_cached::Int = nDof_u
    nNodes_cached::Int = nNodes    
    rCol = Array{Int}(undef,0)

    if DENSE                # initialize the vector of the Neumann boundary conditions (for ndof = 1) / Dirichlet boundary conditions lower surface (for ndof > 1)
        C = Matrix{Float64}(I,nDof_u_cached*nNodes_cached,nDof_u_cached*nNodes_cached)
        
        q_upper = zeros(Float64, nDof_u_cached*nNodes_cached,1)                  # initialize the vector of the Dirichlet boundary conditions (for ndof = 1) / Dirichlet boundary conditions upper surface (for ndof > 1)
        q_lower = zeros(Float64, nDof_u_cached*nNodes_cached,1)                  # initialize the vector of the Neumann boundary conditions (for ndof = 1) / Dirichlet boundary conditions lower surface (for ndof > 1)
        q_side = zeros(Float64, nDof_u_cached*nNodes_cached,1) 
    else                 # initialize the vector of the Neumann boundary conditions (for ndof = 1) / Dirichlet boundary conditions lower surface (for ndof > 1)
        C = sparse(I,nDof_u_cached*nNodes_cached,nDof_u_cached*nNodes_cached)
        
        q_upper = spzeros(nDof_u_cached*nNodes_cached,1)                  # initialize the vector of the Dirichlet boundary conditions (for ndof = 1) / Dirichlet boundary conditions upper surface (for ndof > 1)
        q_lower = spzeros(nDof_u_cached*nNodes_cached,1)                  # initialize the vector of the Neumann boundary conditions (for ndof = 1) / Dirichlet boundary conditions lower surface (for ndof > 1)
        q_side = spzeros(nDof_u_cached*nNodes_cached,1) 
    end

    if nDof_u_cached == 1
        if ndim_cached == 3
            Dbound1 = 0
            Dbound2 = 1
            iter = 1:size(nodeList_cached,2)
            for n::Int in iter
                coord = nodeList_cached[:,n] # get the coordinates of the node
                if coord[3] == Dbound1 # bottom boundary
                    q_upper[n] = 0
                    push!(rCol, n)
                elseif coord[3] == Dbound2 # top boundary
                    q_upper[n] = -1
                    push!(rCol, n)
                end
            end
        elseif ndim_cached == 2
            Dbound1 = 0
            Dbound2 = 1
            iter = 1:size(nodeList_cached,2)
            for n::Int in iter
                coord = nodeList_cached[:,n]
                if coord[2] == Dbound1
                    q_upper[n] = 0
                    push!(rCol, n)
                elseif coord[2] == Dbound2
                    q_upper[n] = -1
                    push!(rCol, n)
                end
            end
        elseif ndim_cached == 1
            Dbound1 = 0
            Dbound2 = 1
            iter = 1:size(nodeList_cached,2)
            for n::Int in iter
                coord = nodeList_cached[:,n]
                if coord[1] == Dbound1
                    q_upper[n] = 0
                    push!(rCol, n)
                elseif coord[1] == Dbound2
                    q_upper[n] = -1
                    push!(rCol, n)
                end
            end
        end

        C_uc = C[:,setdiff(1:size(C,2),rCol)]

    else
        z0Bound = 0
        z1Bound = maximum(mdl.mesh_u.NodeList[ndim_cached, :])

        iter = 1:size(nodeList_cached,2)
        for nNode::Int in iter
            coord = nodeList_cached[:,nNode]
            if coord[ndim_cached] == z1Bound
                q_upper[ID_cached[ndim_cached,nNode]] = 1
                push!(rCol,ID_cached[ndim_cached,nNode])
            elseif coord[ndim_cached] == z0Bound
                q_lower[ID_cached[ndim_cached,nNode]] = 1
                push!(rCol,ID_cached[ndim_cached,nNode])
            end
        end

        C_uc = C[:,setdiff(1:size(C,2),rCol)]

    end
    return q_upper, q_side, q_lower, C_uc
end

"""
    simulate(mdl, scene, conditions; sensitivities=true)

Time-step a squeeze-flow simulation, returning the deformed state and the analytic parameter
sensitivities at every step. The sensitivities `∂q/∂η` and `∂q/∂β` are carried through the solve
rather than finite-differenced afterwards, which is what lets `fit_model` run Gauss-Newton
directly on the forward solve.

Two driving modes are selected by `scene.control`. Under `"force"` the prescribed force history
`scene.cParam` is imposed and the top-plate velocity is an unknown. Under `"velocity"`,
`scene.cParam` is the plate velocity itself. Any other value throws an `ArgumentError`.

Each step assembles the saddle-point system into a fixed sparsity pattern and solves it by
iterative refinement against a factorization reused across steps (refactorized once it has gone
stale), starting from an extrapolation of the previous steps. Results agree with a direct solve
to a relative residual of 1e-9 (state) and 1e-8 (sensitivities).

Output written to `conditions.filepath` is controlled by the flags on `conditions`: `WRITEVTK`
for ParaView scenes, `WRITECONTOUR` for projected 2D contours, `ANIMATE` for GIFs.

# Arguments
- `mdl::Stokes`: Model holding the velocity, pressure and geometry meshes.
- `scene::SqueezeFlow`: Scenario supplying `β`, the control history, the viscosity law and the
  time grid.
- `conditions::Conditions`: Camera model and output flags.
- `μu_btm::Real`, `μu_side::Real`: prescribed velocity multipliers for the bottom plate and the
  lateral surface. Only the default `0.0` — a stationary bottom plate and a free side, the
  configuration every result to date was produced with — is supported; anything else throws
  an `ArgumentError`.
- `sensitivities::Bool`: Compute the parameter sensitivities. `false` skips their solve and
  projection — for prediction only — and leaves `gradList` and `gradList_3d` holding only
  their initial entry.

# Returns
- `output::Vector{Float64}`: Top-plate displacement per time step.
- `gradList`: Contour sensitivities in 2D per time step, one entry per view.
- `borderPts2DList`: Projected border points per time step, one entry per view.
- `displacement`: Nodal displacements per time step.
- `surface_pts_3D`: Side-surface node positions per time step.
- `pos2D`: Projected surface points per time step.
- `pos3D`: Nodal positions per time step.
- `velocity`: Nodal velocities per time step.
- `pressure`: Nodal pressures per time step.
- `gradList_3d`: Nodal sensitivities in 3D per time step.
"""
function simulate(mdl::Stokes, scene::SqueezeFlow, conditions::Conditions;
                  μu_btm::Real=0.0, μu_side::Real=0.0, sensitivities::Bool=true)

    reset_model!(mdl)
    
    h_cached::Float64 = hasproperty(mdl.mesh_x, :h) ? mdl.mesh_x.h : mdl.mesh_x.lz

    IEN_u_cached::Matrix{Int} = mdl.mesh_u.IEN
    side_node_list_cached::Vector{Int} = mdl.mesh_u.side_nodes

    NodeList_p_cached::Matrix{Float64} = mdl.mesh_p.NodeList

    @unpack η = mdl
    η_cached::Any = η

    @unpack β, viscosity_type, t_steps, sim_time, control = scene
    β_cached::Any = β
    viscosity_type_cached::String = viscosity_type
    t_steps_cached::Float64 = t_steps
    sim_time_cached::Float64 = sim_time
    control_cached::String = control
    cParam_cached::Vector{Float64} = scene.cParam

    @unpack camera_matrix, obj_pose, viewing_angles = conditions    
    rot_angle_cached::Vector{Float64} = viewing_angles
    camera_matrix_cached::Matrix{Float64} = camera_matrix
    obj_pose_cached::AbstractArray{Float64} = obj_pose
    
    nodeList_cached::Matrix{Float64} = mdl.mesh_x.NodeList
    time = collect(Float64, range(start=t_steps_cached, stop=sim_time_cached, step=t_steps_cached))
    len_t = length(time)

    if conditions.filepath != ""
        isnothing(conditions.filepath) && throw(AssertionError("Please provide a filepath to write the data"))
        set_file(conditions.filepath)
    end

    
    BorderPts2D, surface_pts_2d, obs_border_pts = _get_2D_data(nodeList_cached, camera_matrix_cached, obj_pose_cached, h_cached, BorderNodesList=side_node_list_cached, angles=rot_angle_cached, FIT=conditions.WRITECONTOUR)

    velocity = AbstractArray[zeros(Float64,size(nodeList_cached,1),size(nodeList_cached,2))]
    pressure = AbstractArray[zeros(Float64,size(nodeList_cached,1),1)]
    displacement = AbstractArray[zeros(Float64,size(nodeList_cached,1),size(nodeList_cached,2))]
    surface_pts_3D = AbstractArray[nodeList_cached[:,side_node_list_cached]']
    gradList = AbstractArray[[zeros(Float64, 3, size(b,2), 2) for b in BorderPts2D]] # store the solution fields of the border nodes in 2D, one entry per view; the initial contour is θ-independent
    gradList_3d = AbstractArray[zeros(Float64, size(nodeList_cached,1),size(nodeList_cached,2),2)] # store the solution fields of the border nodes in 3D
    pos3D = AbstractArray[nodeList_cached]
    pos2D = AbstractArray[surface_pts_2d]
    borderPts2DList = AbstractArray[BorderPts2D] # store the border points in 2D for each angle at each time step
    output = Float64[]
    writeborderList = [obs_border_pts]

    dNodeList_dη = zeros(Float64, size(nodeList_cached))
    dNodeList_dβ = zeros(Float64, size(nodeList_cached))

    iter::Int = 1
    pr = progress_guard(len_t; desc= "Simulating with prescribed $(control_cached) ...", showspeed=true)
    
    # Pre-compute basis functions before time loop and create cache
    cache = BasisFunctionCache(mdl)

    control_cached in ("force", "velocity") || throw(ArgumentError("Control type is unknown"))
    (μu_btm == 0 && μu_side == 0) ||
        throw(ArgumentError("the solver assumes an undriven bottom plate and sides (μu_btm = μu_side = 0)"))
    velocity_control = control_cached == "velocity"

    sys = _cached_squeeze_system(mdl, scene)
    nf, np = sys.nf, sys.np
    n = size(sys.K, 1)
    # Solve for the response to a unit plate force (or unit plate velocity) and scale by the
    # control value: the unit response depends on the geometry alone, so it extrapolates well
    # even when a measured control signal is noisy.
    rhs_z = zeros(Float64, n)
    rhs_z[end] = 1.0
    z = zeros(Float64, n)
    w = zeros(Float64, n)                       # ∂z/∂β
    # M ∂z/∂β = -(∂M/∂β) z. Under velocity control the plate column of ∂M/∂β times z_μ = 1 is
    # the right-hand side's own derivative, and the plate row stays fixed.
    rhs_w!(r, z) = (mul!(r, sys.dKdβ, z, -1.0, 0.0); velocity_control && (r[end] = 0.0); r)
    z_hist = [zeros(Float64, n) for _ in 1:4]
    w_hist = [zeros(Float64, n) for _ in 1:4]
    sol = zeros(Float64, n)
    dsoldη = zeros(Float64, n)
    dsoldβ = zeros(Float64, n)
    dvdη = zeros(Float64, size(sys.node_dofs))
    dvdβ = zeros(Float64, size(sys.node_dofs))
    # z's velocity part scales as η^(a-1) and its pressure as η^a, w = ∂z/∂β one power of η
    # lower. The history is kept at unit viscosity, where it varies only with the geometry and
    # β/η — bulk viscosity changes η every step.
    a = velocity_control ? 1 : 0

    for t in time
        η_t = viscosity_type_cached == "bulk_viscosity" ? η_cached[iter] : η_cached[1]
        β_t = viscosity_type_cached == "bulk_viscosity" && length(β_cached) > 1 ? β_cached[iter] : β_cached[1]
        _assemble_squeeze_system!(sys, mdl, cache, η_t, β_t)
        velocity_control && _prescribe_plate_velocity!(sys, rhs_z)
        # a factorization left over from an earlier call would make results depend on it
        iter == 1 && _refactor!(sys)

        n_hist = min(iter - 1, 4)
        _scale_unknowns!(_extrapolate!(z, z_hist, n_hist), sys, η_t^(a - 1), η_t^a)
        if sensitivities
            _scale_unknowns!(_extrapolate!(w, w_hist, n_hist), sys, η_t^(a - 2), η_t^(a - 1))
            _lagged_solve!(sys, z, rhs_z, w, rhs_w!)
            _scale_unknowns!(_push_history!(w_hist, w)[1], sys, η_t^(2 - a), η_t^(1 - a))
        else
            _lagged_solve!(sys, z, rhs_z)
        end
        _scale_unknowns!(_push_history!(z_hist, z)[1], sys, η_t^(1 - a), η_t^(-a))

        c_t = cParam_cached[iter]               # plate force, or plate velocity
        @. sol = c_t * z
        μ_tp = sol[end]
        p = sol[(nf + 1):(nf + np)]'   # a copy: sol is reused next step
        velocity_field = _nodal_velocity(sol, sys)
        motion = velocity_field*t_steps_cached
        nodeList_cached = nodeList_cached + motion

        if sensitivities
            @. dsoldβ = c_t * w
            # M = ηĀ + βb + (θ-free divergence blocks), so scaling (η, β) by s scales the
            # velocity by 1/s under force control and the pressure by s under velocity control:
            # η ∂sol/∂η + β ∂sol/∂β = -[u; 0] or [0; p]. That gives ∂/∂η without a second solve.
            if velocity_control
                @. dsoldη = -β_t * dsoldβ / η_t
                @views @. dsoldη[(nf + 1):(nf + np)] += sol[(nf + 1):(nf + np)] / η_t
            else
                @. dsoldη = -(sol + β_t * dsoldβ) / η_t
                @views @. dsoldη[(nf + 1):(nf + np)] = -β_t * dsoldβ[(nf + 1):(nf + np)] / η_t
            end

            _nodal_velocity!(dvdη, dsoldη, sys)
            _nodal_velocity!(dvdβ, dsoldβ, sys)
            mat_nan_inf_check(dvdη)
            mat_nan_inf_check(dvdβ)
            @. dNodeList_dη += dvdη * t_steps_cached
            @. dNodeList_dβ += dvdβ * t_steps_cached
            dmdθ_out = @views cat(dNodeList_dη,dNodeList_dβ,dims=3)

            BorderPts2D, dudθ, surface_pts_2d, _, obs_border_pts = _get_2D_data(nodeList_cached, camera_matrix_cached, obj_pose_cached, h_cached, BorderNodesList=side_node_list_cached, GRAD=true, dqdθ=dmdθ_out, angles=rot_angle_cached, FIT=conditions.WRITECONTOUR)
            push!(gradList, dudθ)
            push!(gradList_3d, dmdθ_out)
        else
            BorderPts2D, surface_pts_2d, obs_border_pts = _get_2D_data(nodeList_cached, camera_matrix_cached, obj_pose_cached, h_cached, BorderNodesList=side_node_list_cached, angles=rot_angle_cached, FIT=conditions.WRITECONTOUR)
        end
        mdl.mesh_x.NodeList = nodeList_cached

        push!(output, μ_tp*t_steps_cached)
        push!(velocity, velocity_field)
        push!(pressure, p)
        push!(displacement, motion)
        push!(surface_pts_3D, nodeList_cached[:,side_node_list_cached]')
        push!(pos2D, surface_pts_2d)
        push!(pos3D, nodeList_cached)
        push!(borderPts2DList, BorderPts2D)
        push!(writeborderList, obs_border_pts)

        iter += 1
        next!(pr, showvalues = [(:iterations,iter),(:time,string(t," seconds"))])
    end

    if conditions.WRITEVTK
        # write_scene(string(conditions.filepath,"/data"), NodeList_p_cached, mdl.mesh_p.IEN, mdl.ne, mdl.ndim, pressure, ID=ID_cached, FunctionClass=mdl.mesh_p.FunctionClass)
        write_stokes_scene(string(conditions.filepath,"/data"), mdl.mesh_u.NodeList, mdl.mesh_u.IEN, NodeList_p_cached, mdl.mesh_p.IEN, mdl.ne, mdl.ndim, velocity, pressure, pos3D=pos3D,
            element_shape_u=mdl.mesh_u.volume_element_shape, basis_order_u=mdl.mesh_u.basis_order,
            element_shape_p=mdl.mesh_p.volume_element_shape, basis_order_p=mdl.mesh_p.basis_order)
    end

    if conditions.WRITECONTOUR
        write_2d_data(string(conditions.filepath,"/data/sim_data/contour_data"), writeborderList)
    end
    
    # write the data to a file
    if conditions.ANIMATE
        animate_fields(filepath = string(conditions.filepath,"/Results/images/"), Nodes=pos3D , IEN=IEN_u_cached, sim_border_nodes_2d=borderPts2DList, sim_pts_2d=pos2D)
        animate_fields(filepath = string(conditions.filepath,"/Results/images/surface"), Nodes=surface_pts_3D)
    end
    
    return output, gradList, borderPts2DList, displacement, surface_pts_3D, pos2D, pos3D, velocity, pressure, gradList_3d
end
