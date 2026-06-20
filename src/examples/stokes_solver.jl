using LinearAlgebra
using SparseArrays
using ProgressMeter
using Parameters

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

    C::Matrix{Float64} = get_cMat(1.0,0.0,type="standard")
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
            coords::Matrix{Float64} = NodeList_x_cached[:, IEN_x_cached[:, e]]  # Get the coordinates of the nodes of the element

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
    # println("volume A:", vol)
    return K
end

"""
    assemble_system_A_dense(mdl::Stokes)

Assembles the finite element system matrix (dense format) for the velocity field.

# Arguments:
- `mdl::Stokes` : Material model

# Returns:
- `K::Matrix{Float64}` : Dense stiffness matrix of size [ndof, ndof]
"""
function assemble_system_A_dense(mdl::Stokes)::Matrix{Float64}
    @unpack ne, ndim, nDof_u = mdl
    @unpack NodeList, IEN, ID, volume_element_shape, basis_order = mdl.mesh_u

    C::Matrix{Float64} = get_cMat(1.0,0.0,type="standard")
    IEN_u_rows::Int = size(IEN,1)
    ID_u_rows::Int = size(ID,1)

    if nDof_u == 1
        sz = ne^ndim*IEN_u_rows # no elements x no nodes
        K = zeros(Float64, sz,sz)
    else
        sz = nDof_u*(ne)^ndim # no elements x no nodes x no dofs
        K = zeros(Float64, sz,sz)  
    end

    if ndim == 1
        # gaussian quadrature points for the element [-1,1] 
        ξ, w_ξ = gaussian_quadrature(-1,1)
        
        wpoints =  [w_ξ[1], w_ξ[2]]
        
        x = [ξ[1], ξ[2]]
    elseif ndim == 2
        # gaussian quadrature points for the element [-1,1]x[-1,1] 
        ξ, w_ξ = gaussian_quadrature(-1,1,3)
        η, w_η = gaussian_quadrature(-1,1,3)

        # Pre-allocate arrays for efficiency (faster than push!)
        npts = size(ξ,1) * size(η,1)
        x = Vector{Float64}(undef, npts)
        y = Vector{Float64}(undef, npts)
        wpoints = Vector{Float64}(undef, npts)
        
        idx = 1
        for j::Int in 1:size(η,1)
            for i::Int in 1:size(ξ,1)
                x[idx] = ξ[i]
                y[idx] = η[j]
                wpoints[idx] = w_ξ[i]*w_η[j]
                idx += 1
            end
        end

    elseif ndim == 3
        # gaussian quadrature points for the element [-1,1]x[-1,1]x[-1,1] 
        ξ, w_ξ = gaussian_quadrature(-1,1,3)
        η, w_η = gaussian_quadrature(-1,1,3)
        ζ, w_ζ = gaussian_quadrature(-1,1,3)

        # Pre-allocate arrays for efficiency (faster than push!)
        npts = size(ξ,1) * size(η,1) * size(ζ,1)
        x = Vector{Float64}(undef, npts)
        y = Vector{Float64}(undef, npts)
        z = Vector{Float64}(undef, npts)
        wpoints = Vector{Float64}(undef, npts)
        
        idx = 1
        for k::Int in 1:size(ζ,1)
            for j::Int in 1:size(η,1)
                for i::Int in 1:size(ξ,1)
                    x[idx] = ξ[i]
                    y[idx] = η[j]
                    z[idx] = ζ[k]
                    wpoints[idx] = w_ξ[i]*w_η[j]*w_ζ[k]
                    idx += 1
                end
            end
        end
    end
    # initialize the Jacobian matrix
    Jac = zeros(Float64, ndim_cached, ndim_cached)  # Jacobian matrix

    e_iter = 1:ne^ndim # iterator for elements loop
    gpiter = 1:length(wpoints) # iterator for integration loop
    # integration loop
    for gp::Int in gpiter
        if ndim == 1
            N, ΔN = basis_function(x[gp], nothing, nothing, volume_element_shape, basis_order)
        elseif ndim == 2
            N, ΔN = basis_function(x[gp], y[gp], nothing, volume_element_shape, basis_order)
        elseif ndim == 3
            N, ΔN = basis_function(x[gp], y[gp], z[gp], volume_element_shape, basis_order)
        end

        dNdX = zeros(Float64, size(ΔN, 1), ndim_cached) # Gradient of basis functions
        if nDof_u == 2
            B = zeros(Float64, ndim*length(N), 3)
        elseif nDof_u == 3
            B = zeros(Float64, ndim*length(N), 6)
        end
        szB::Int = size(B, 1)  # Number of basis functions
        Ke = zeros(Float64, szB, szB)  # Element stiffness matrix
        Ke_row::Int, Ke_col::Int = szB, szB
        # element loop
        for e::Int in e_iter
            coords::Matrix{Float64} = nodeList_cached[:, IEN_cached[:, e]]  # Get the coordinates of the nodes of the element

            mul!(Jac, coords, ΔN)  # Jacobian matrix [dx/dxi dx/deta; dy/dxi dy/deta]
            w::Float64 = wpoints[gp] * abs(det(Jac))
            dNdX .= ΔN / Jac  # Solve Jac * X = ΔN' (5-10x faster than inv())

            if nDof_u == 1
                szN = size(N,1) # number of basis functions
                # loop between basis functions of the element
                for i::Int in 1:szN
                    for j::Int in 1:szN
                        inz = (szN)^2*(e-1) + szN*(i-1) + j # index for the COO sparse matrix
                        E[inz] = IEN[i,e] # row index 
                        J[inz] = IEN[j,e] # column index
                        V[inz] += w*dot(dNdX[i,:],dNdX[j,:])# inner product of the gradient of the basis functions
                    end
                end
            else   
                if nDof_u == 2
                    B .= 0.0
                    B[1:nDof_u:end,1] = dNdX[:,1]
                    B[2:nDof_u:end,2] = dNdX[:,2]
                    B[1:nDof_u:end,3] = dNdX[:,2]
                    B[2:nDof_u:end,3] = dNdX[:,1]
                elseif nDof_u == 3
                    B .= 0.0
                    B[1:nDof_u:end,1] = dNdX[:,1]
                    B[2:nDof_u:end,2] = dNdX[:,2]
                    B[3:nDof_u:end,3] = dNdX[:,3]
                    B[2:nDof_u:end,4] = dNdX[:,3]
                    B[3:nDof_u:end,4] = dNdX[:,2]
                    B[1:nDof_u:end,5] = dNdX[:,3]
                    B[3:nDof_u:end,5] = dNdX[:,1]
                    B[1:nDof_u:end,6] = dNdX[:,2]
                    B[2:nDof_u:end,6] = dNdX[:,1]
                end
                # Ke = 2*w*B*C*B' # element stiffness matrix
                mul!(Ke, 2*w*B, C*B') # element stiffness matrix

                # loop between basis functions of the element
                iNodes = 1:div(Ke_row,nDof_u)  # column node index
                jNodes = 1:div(Ke_col,nDof_u)  # row index
                iDofs = 1:ID_u_rows
                jDofs = 1:ID_u_rows
                for iNode::Int in iNodes
                    for jNode::Int in jNodes
                        for iDof::Int in iDofs
                            for jDof in jDofs
                                i::Int = (iNode-1)*nDof_u + iDof
                                j::Int = (jNode-1)*nDof_u + jDof
                                K[ID[iDof,IEN[iNode,e]],ID[jDof,IEN[jNode,e]]] += Ke[i,j] 
                            end
                        end
                    end
                end
            end
        end
    end
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

            coords::Matrix{Float64} = NodeList_x_cached[:, IEN_x_cached[:, e]]  # Get the coordinates of the nodes of the element

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
    # println("volume B:", vol)
    return K
end

"""
    assemble_system_B_dense(mdl::Stokes)

Assembles the pressure-velocity coupling matrix (B matrix) in dense format for the Stokes system.

# Arguments:
- `mdl::Stokes` : Material model

# Returns:
- `K::Matrix{Float64}` : Dense coupling matrix of size [ndof_u, ndof_p]
"""
function assemble_system_B_dense(mdl::Stokes)::Matrix{Float64}

    @unpack ne, ndim, nDof_u = mdl
    @unpack IEN, volume_element_shape, basis_order = mdl.mesh_p

    IEN_p_cached::Matrix{Int} = IEN
    element_shape_p_cached, basis_order_p_cached = volume_element_shape, basis_order

    @unpack NodeList, IEN, ID, volume_element_shape, basis_order = mdl.mesh_p

    IEN_u_cached::Matrix{Int} = IEN
    nodeList_cached::Matrix{Float64} = NodeList
    ID_cached::Matrix{Int} = ID
    ne_cached::Int = ne
    ndim_cached::Int = ndim
    nDof_u_cached::Int = nDof_u
    element_shape_u_cached, basis_order_u_cached = volume_element_shape, basis_order

    IEN_u_rows = size(IEN_u_cached,1)
    ID_u_rows = size(ID_cached,1)
    if nDof_u_cached == 1
        sz = ne_cached*IEN_u_rows # no elements x no nodes
        K = zeros(Float64, sz,sz)
    else
        sz1 = nDof_u_cached*nNodes_cached # no elements x no nodes x no dofs
        sz2 = nNodes_cached
        K = zeros(Float64, sz1, sz2)  
    end

    # element loop
    if ndim_cached == 1
        # gaussian quadrature points for the element [-1,1] 
        ξ, w_ξ = gaussian_quadrature(-1,1)
    
        wpoints = [w_ξ[1], w_ξ[2]]
        
        x = [ξ[1], ξ[2]]
    elseif ndim_cached == 2
        # gaussian quadrature points for the element [-1,1]x[-1,1] 
        ξ, w_ξ = gaussian_quadrature(-1,1,3)
        η, w_η = gaussian_quadrature(-1,1,3)

        # Pre-allocate arrays for efficiency (faster than push!)
        npts = size(ξ,1) * size(η,1)
        x = Vector{Float64}(undef, npts)
        y = Vector{Float64}(undef, npts)
        wpoints = Vector{Float64}(undef, npts)
        
        idx = 1
        for j::Int in 1:size(η,1)
            for i::Int in 1:size(ξ,1)
                x[idx] = ξ[i]
                y[idx] = η[j]
                wpoints[idx] = w_ξ[i]*w_η[j]
                idx += 1
            end
        end

    elseif ndim_cached == 3
        # gaussian quadrature points for the element [-1,1]x[-1,1]x[-1,1] 
        ξ, w_ξ = gaussian_quadrature(-1,1,3)
        η, w_η = gaussian_quadrature(-1,1,3)
        ζ, w_ζ = gaussian_quadrature(-1,1,3)

        # Pre-allocate arrays for efficiency (faster than push!)
        npts = size(ξ,1) * size(η,1) * size(ζ,1)
        x = Vector{Float64}(undef, npts)
        y = Vector{Float64}(undef, npts)
        z = Vector{Float64}(undef, npts)
        wpoints = Vector{Float64}(undef, npts)
        
        idx = 1
        for k::Int in 1:size(ζ,1)
            for j::Int in 1:size(η,1)
                for i::Int in 1:size(ξ,1)
                    x[idx] = ξ[i]
                    y[idx] = η[j]
                    z[idx] = ζ[k]
                    wpoints[idx] = w_ξ[i]*w_η[j]*w_ζ[k]
                    idx += 1
                end
            end
        end
    end
    # initialize the Jacobian matrix
    Jac = zeros(Float64, ndim_cached, ndim_cached)  # Jacobian matrix

    e_iter = 1:ne_cached    # iterator for elements loop
    gpiter = 1:length(wpoints) # iterator for integration loop
    # integration loop
    for gp::Int in gpiter
        if ndim_cached == 1
            N_u, ΔN_u = basis_function(x[gp], nothing, nothing, element_shape_u_cached, basis_order_u_cached)
            N_p, ΔN_p = basis_function(x[gp], nothing, nothing, element_shape_p_cached, basis_order_p_cached)
        elseif ndim_cached == 2
            N_u, ΔN_u = basis_function(x[gp], y[gp], nothing, element_shape_u_cached, basis_order_u_cached)
            N_p, ΔN_p = basis_function(x[gp], y[gp], nothing, element_shape_p_cached, basis_order_p_cached)
        elseif ndim_cached == 3
            N_u, ΔN_u = basis_function(x[gp], y[gp], z[gp], element_shape_u_cached, basis_order_u_cached)
            N_p, ΔN_p = basis_function(x[gp], y[gp], z[gp], element_shape_p_cached, basis_order_p_cached)
        end

        dNdX = zeros(Float64, size(ΔN_u, 1), ndim_cached) # Gradient of basis functions
        Be = zeros(Float64, 1, 3*size(dNdX, 1)) # Gradient of basis functions
        colB::Int = size(Be, 2)  # Number of basis functions
        rowNp::Int = size(N_p, 1)  # Number of basis functions
        Ke = zeros(Float64, rowNp, colB)  # Element stiffness matrix
        Ke_row::Int, Ke_col::Int = rowNp, colB

        # element loop
        for e::Int in e_iter
            coords_u::Matrix{Float64} = nodeList_cached[:, IEN_u_cached[:, e]]  # Get the coordinates of the nodes of the element

            mul!(Jac, coords_u, ΔN_u)  # Jacobian matrix [dx/dxi dx/deta; dy/dxi dy/deta]
            w::Float64 = wpoints[gp] * abs(det(Jac))
            dNdX .= ΔN_u / Jac # Solve Jac * X = ΔN_u' (5-10x faster than inv())

            if nDof_u_cached == 1
                szN = size(N_u,1) # number of basis functions
                # loop between basis functions of the element
                for i in 1:szN
                    for j in 1:szN
                        inz = (szN)^2*(e-1) + szN*(i-1) + j # index for the COO sparse matrix
                        E[inz] = mdl.IEN[i,e] # row index 
                        J[inz] = mdl.IEN[j,e] # column index
                        V[inz] += w*dot(dNdX[i,:],dNdX[j,:])# inner product of the gradient of the basis functions
                    end
                end
            else  
                Be .= 0.0
                Be[1:nDof_u_cached:end] = dNdX[:,1]
                Be[2:nDof_u_cached:end] = dNdX[:,2]
                Be[3:nDof_u_cached:end] = dNdX[:,3]
                # Ke = -w*N_p*Be # element stiffness matrix
                mul!(Ke, -w*N_p, Be) # element stiffness matrix

                # loop between basis functions of the element
                iNodes = 1:Ke_col÷nDof_u_cached  # column node index
                iDofs = 1:ID_u_rows # column dof index
                jNodes = 1:Ke_row # row index

                for iNode in iNodes 
                    for jNode in jNodes 
                        for iDof in iDofs
                            i::Int = (iNode-1)*nDof_u_cached + iDof
                            j::Int = jNode
                            K[ID_cached[iDof,IEN_u_cached[iNode,e]],IEN_p_cached[jNode,e]] += Ke[i,j] 
                        end
                    end
                end
            end
        end
    end
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
            coords_top::Matrix{Float64} = NodeList_x_cached[:,IEN_x_top_cached[:,e]] # get the coordinates of the nodes of the element
            coords_btm::Matrix{Float64} = NodeList_x_cached[:,IEN_x_btm_cached[:,e]] # get the coordinates of the nodes of the element

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

            mul!(be_btm,M_btm',M_btm) # multiply the matrix by itself to get the stiffness matrix
            mul!(be_top,M_top',M_top) # multiply the matrix by itself to get the stiffness matrix

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
    # println("Maximum row :",maximum(E))
    # println("Maximum column :",maximum(J))
    # println("A_top = ", A_top)
    # println("A_btm = ", A_btm)
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
            coords_top::Matrix{Float64} = NodeList_x_cached[:,IEN_x_top_cached[:,e]] # get the coordinates of the nodes of the element
            coords_btm::Matrix{Float64} = NodeList_x_cached[:,IEN_x_btm_cached[:,e]] # get the coordinates of the nodes of the element

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
            mul!(be_btm,M_btm',M_btm) # multiply the matrix by itself to get the stiffness matrix
            mul!(be_top,M_top',M_top) # multiply the matrix by itself to get the stiffness matrix

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
        C = Matrix{Float64}(I,nDof_u_cached*nNodes_cached,nDof_u_cached*nNodes_cached)      # definition of the constraint matrix
        
        q_upper = zeros(Float64, nDof_u_cached*nNodes_cached,1)                  # initialize the vector of the Dirichlet boundary conditions (for ndof = 1) / Dirichlet boundary conditions upper surface (for ndof > 1)
        q_lower = zeros(Float64, nDof_u_cached*nNodes_cached,1)                  # initialize the vector of the Neumann boundary conditions (for ndof = 1) / Dirichlet boundary conditions lower surface (for ndof > 1)
        q_side = zeros(Float64, nDof_u_cached*nNodes_cached,1) 
    else                 # initialize the vector of the Neumann boundary conditions (for ndof = 1) / Dirichlet boundary conditions lower surface (for ndof > 1)
        C = sparse(I,nDof_u_cached*nNodes_cached,nDof_u_cached*nNodes_cached)      # definition of the constraint matrix
        
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

function simulate(mdl::Stokes, scene::SqueezeFlow, conditions::Conditions)

    reset_model!(mdl)
    
    @unpack volume_element_shape, basis_order, IEN, ID, NodeList = mdl.mesh_x
    h_cached::Float64 = hasproperty(mdl.mesh_x, :h) ? mdl.mesh_x.h : mdl.mesh_x.lz
    element_shape_x_cached, basis_order_x_cached = volume_element_shape, basis_order
    NodeList_x_cached::Matrix{Float64} = NodeList
    IEN_x_cached::Matrix{Int} = IEN
    ID_x_cached::Matrix{Int} = ID

    @unpack IEN, ID, volume_element_shape, basis_order, top_nodes, bottom_nodes, side_nodes, nNodes = mdl.mesh_u
    element_shape_u_cached, basis_order_u_cached = volume_element_shape, basis_order
    IEN_u_cached::Matrix{Int} = IEN
    ID_u_cached::Matrix{Int} = ID
    nNodes_u_cached::Int = nNodes
    NodeList_u_cached::Matrix{Float64} = NodeList
    top_node_list_cached::Vector{Int} = top_nodes # top nodes
    bottom_node_list_cached::Vector{Int} = bottom_nodes # bottom nodes
    side_node_list_cached::Vector{Int} = side_nodes

    @unpack nNodes, NodeList = mdl.mesh_p
    nNodes_p_cached::Int = nNodes
    NodeList_p_cached::Matrix{Float64} = NodeList

    @unpack η, nDof_u, nDof_p, ndim = mdl
    η_cached::Any = η
    nDof_u_cached::Int = nDof_u
    nDof_p_cached::Int = nDof_p
    ndim_cached::Int = ndim

    @unpack β, viscosity_type, q_d, C_uc, t_steps, sim_time, control, cParam = scene
    β_cached::Any = β
    viscosity_type_cached::String = viscosity_type
    q_d_cached::Dict{Symbol, Matrix{Float64}} = q_d
    q_d_cached_top::SparseMatrixCSC{Bool, Int64} = q_d_cached[:top]
    q_d_cached_btm::SparseMatrixCSC{Bool, Int64} = q_d_cached[:bottom]
    q_d_cached_brdr::SparseMatrixCSC{Bool, Int64} = q_d_cached[:border]
    C_uc_cached::SparseMatrixCSC{Bool, Int64} = C_uc
    t_steps_cached::Float64 = t_steps
    sim_time_cached::Float64 = sim_time
    control_cached::String = control
    cParam_cached::Vector{Float64} = scene.cParam

    @unpack camera_matrix, obj_pose, viewing_angles = conditions    
    rot_angle_cached::Vector{Float64} = viewing_angles
    camera_matrix_cached::Matrix{Float64} = camera_matrix
    obj_pose_cached::Vector{Float64} = obj_pose
    
    nodeList_cached::Matrix{Float64} = NodeList_u_cached
    ID_cached::Matrix{Int} = ID_u_cached
    time = collect(Float64, range(start=t_steps_cached, stop=sim_time_cached, step=t_steps_cached))
    len_t = length(time)

    C_Tu = transpose(C_uc_cached) # transpose the constraint matrix

    if conditions.filepath != ""
        isnothing(conditions.filepath) && throw(AssertionError("Please provide a filepath to write the data"))
        set_file(conditions.filepath)
    end

    μu_btm = 0  
    μu_side = 0
    
    BorderPts2D, surface_pts_2d, obs_border_pts = _get_2D_data(nodeList_cached, camera_matrix_cached, obj_pose_cached, h_cached, BorderNodesList=side_node_list_cached, angles=rot_angle_cached)
    dqdη = zeros(Float64, size(q_d_cached_top))
    dqdβ = zeros(Float64, size(q_d_cached_top))

    velocity = AbstractArray[zeros(Float64,size(nodeList_cached,1),size(nodeList_cached,2))] # store the velocity of the mesh in 3D
    pressure = AbstractArray[zeros(Float64,size(nodeList_cached,1),1)] # store the pressure of the mesh in 3D
    displacement = AbstractArray[zeros(Float64,size(nodeList_cached,1),size(nodeList_cached,2))] # store the displacement of the mesh in 3D
    surface_fields = AbstractArray[]
    surface_pts_3D = AbstractArray[vcat(nodeList_cached[:,top_node_list_cached]', 
                                        nodeList_cached[:,bottom_node_list_cached]', 
                                        nodeList_cached[:,side_node_list_cached]')'] # store the solution fields of the mesh in 3D
    gradList = AbstractArray[zeros(Float64, size(BorderPts2D,1),size(BorderPts2D,2),2)] # store the solution fields of the border nodes in 2D
    gradList_3d = AbstractArray[zeros(Float64, size(nodeList_cached,1),size(nodeList_cached,2),2)] # store the solution fields of the border nodes in 3D
    pos3D = AbstractArray[nodeList_cached]       # store the solution fields of the mesh in 3D
    pos3D_cp = AbstractArray[nodeList_cached]
    pos2D = AbstractArray[surface_pts_2d]          # store the solution fields of the mesh in 2D
    borderPts2DList = AbstractArray[BorderPts2D] # store the border points in 2D for each angle at each time step
    splinep = AbstractArray[BorderPts2D[1,:]]  # store the x coordinates samples of the spline parameters of the border nodes
    splineq = AbstractArray[BorderPts2D[2,:]]  # store the y coordinates samples of the spline parameters of the border nodes
    output = Float64[] 
    writeborderList = [obs_border_pts]

    dac_list = Float64[]
    bd_list = Float64[]
    dad_list = Float64[]
    A_list = Float64[]

    _A_bar = SparseMatrixCSC{Float64,Int}(I, nDof_u_cached*nNodes_u_cached, nDof_u_cached*nNodes_u_cached)  # initialize the stiffness matrix
    B = SparseMatrixCSC{Float64,Int}(I, nDof_u_cached*nNodes_u_cached, nDof_p_cached*nNodes_p_cached)      # initialize the stiffness matrix
    b = SparseMatrixCSC{Float64,Int}(I, nDof_u_cached*nNodes_u_cached, nDof_u_cached*nNodes_u_cached)      # initialize the stiffness matrix
    q_d = spzeros(nDof_u_cached*nNodes_u_cached,1)                                                                       # initialize the vector of the Dirichlet boundary conditions (for ndof = 1) / Dirichlet boundary conditions upper surface (for ndof > 1)
    A = similar(_A_bar)
    A_bar = similar(_A_bar)

    A_free = SparseMatrixCSC{Float64, Int64}(I, size(C_Tu,1),size(C_uc_cached,2)) # convert to sparse matrix
    B_free = SparseMatrixCSC{Float64, Int64}(I, size(C_Tu,1),size(B,2)) # convert to sparse matrix

    dA_freedη = similar(A_free)                         
    dA_freedβ = similar(A_free)                         
    dB_free = spzeros(size(B_free))         

    dAdη = similar(_A_bar)
    dAdβ = similar(_A_bar)
    dB = spzeros(size(B))

    q = similar(q_d)

    dNodeList_dη = zeros(Float64, size(nodeList_cached))
    dNodeList_dβ = zeros(Float64, size(nodeList_cached))

    iter::Int = 1
    pr = progress_guard(len_t; desc= "Simulating with prescribed $(control_cached) ...", showspeed=true)
    
    # Pre-compute basis functions before time loop and create cache
    cache = BasisFunctionCache(mdl)

    if control_cached == "force"

        M = spzeros((size(A_free,1)+size(B_free,2)+1),(size(A_free,2)+size(B_free,2)+1))
        dMdη = spzeros(size(M))
        dMdβ = spzeros(size(M))
        
        # Pre-allocate zero matrices for sensitivity RHS (Problem 5 optimization)
        zero_matrix_p_q = zeros(Float64, size(B,2), size(q_d,2))
        zero_matrix_np = spzeros(size(B_free,2), size(B_free,2))  # Pre-allocate zero pressure-pressure block
        
        # Pre-allocate RHS vectors
        r = zeros(Float64, size(A_free,1) + size(B_free,2) + 1)
        drdη = zeros(Float64, size(A_free,1) + size(B_free,2) + 1)
        drdβ = zeros(Float64, size(A_free,1) + size(B_free,2) + 1)

        sol = zeros(Float64, size(r))

        for t in time
            _A_bar .= assemble_system_A(mdl, cache)    # assemble the stiffness matrix
            B .= assemble_system_B(mdl, cache)         # assemble the stiffness matrix
            b .= apply_boundary_conditions(mdl, cache) # apply the neumann boundary conditions
            q_d .= (μu_btm*q_d_cached_btm + μu_side*q_d_cached_brdr) # apply the Dirichlet boundary conditions
   
            if viscosity_type_cached == "bulk_viscosity"
                if length(β_cached) == 1
                    A .= η_cached[iter]*_A_bar + β_cached[1]*b
                else
                    A .= η_cached[iter]*_A_bar + β_cached[iter]*b
                end
            else
                A .= η_cached[1]*_A_bar + β_cached[1]*b
            end

            dAdη .= _A_bar
            dAdβ .= b
 
            A_free .= C_Tu*A*C_uc_cached # extract the free part of the stiffness matrix
            B_free .= C_Tu*B             # extract the free part of the stiffness matrix

            dA_freedη .= C_Tu*dAdη*C_uc_cached # extract the free part of the stiffness matrix
            dA_freedβ .= C_Tu*dAdβ*C_uc_cached # extract the free part of the stiffness matrix

            M = [A_free B_free C_Tu*A*q_d_cached_top; 
                B_free' zero_matrix_np B'*q_d_cached_top;
                q_d_cached_top'*A*C_uc_cached q_d_cached_top'*B (q_d_cached_top'*A*q_d_cached_top)[end]]
            
            dMdη = [dA_freedη dB_free C_Tu*dAdη*q_d_cached_top; 
                dB_free' zero_matrix_np dB'*q_d_cached_top;
                q_d_cached_top'*dAdη*C_uc_cached q_d_cached_top'*dB (q_d_cached_top'*dAdη*q_d_cached_top)[end]]

            dMdβ = [dA_freedβ dB_free C_Tu*dAdβ*q_d_cached_top; 
                dB_free' zero_matrix_np dB'*q_d_cached_top;
                q_d_cached_top'*dAdβ*C_uc_cached q_d_cached_top'*dB (q_d_cached_top'*dAdβ*q_d_cached_top)[end]]
                
            r .= [-C_Tu*A*q_d; -B'*q_d; cParam_cached[iter].-q_d_cached_top'*A*q_d]    # assemble the system of equations
            drdη .= -[C_Tu*dAdη*q_d; zero_matrix_p_q; q_d_cached_top'*dAdη*q_d] # solve the system of equations
            drdβ .= -[C_Tu*dAdβ*q_d; zero_matrix_p_q; q_d_cached_top'*dAdβ*q_d] # solve the system of equations           
            
            @debug begin
                println("Time: ", t)
                println("Iteration: ", iter)
                println("η: ", η_cached[1])
                println("β: ", β_cached[1])
                println("Norm of _A_bar: ", maximum(_A_bar))
                println("Norm of A: ", maximum(A))
                println("Norm of b: ", maximum(b))
                println("Norm of A: ", maximum(A))
                println("Norm of B: ", maximum(B))
                println("Norm of M: ", maximum(M))
                push!(dac_list, norm(q_d_cached_top'*A*C_uc_cached))
                push!(bd_list, norm(q_d_cached_top'*B))
                push!(dad_list, norm((q_d_cached_top'*A*q_d_cached_top)[end]))
                push!(A_list, norm(A))
            end

            lum = lu(M) # LU decomposition of the system of equations

            sol = lum\r            # solve the system of equations
            dsoldη = lum\(drdη - dMdη*sol) # solve the system of equations
            dsoldβ = lum\(drdβ - dMdβ*sol) # solve the system of equations

            q_f = view(sol, 1:size(A_free, 1))
            dqfdη = view(dsoldη, 1:size(A_free, 1))
            dqfdβ = view(dsoldβ, 1:size(A_free, 1))

            p_f = view(sol, (size(A_free, 1) + 1):(size(A_free, 1) + size(B_free, 2)))
            dpfdη = view(dsoldη, (size(A_free, 1) + 1):(size(A_free, 1) + size(B_free, 2)))
            dpfdβ = view(dsoldβ, (size(A_free, 1) + 1):(size(A_free, 1) + size(B_free, 2)))

            μ_tp = sol[end]
            dμdη = dsoldη[end]
            dμdβ = dsoldβ[end]
            
            q .= q_d + C_uc_cached*q_f + μ_tp*q_d_cached_top;       # assemble the solution 
            dqdη = C_uc_cached*dqfdη + dμdη*q_d_cached_top; # assemble the solution
            dqdβ = C_uc_cached*dqfdβ + dμdβ*q_d_cached_top; # assemble the solution

            p = p_f'       # assemble the solution;
            dpdη = dpfdη'; # assemble the solution
            dpdβ = dpfdβ'; # assemble the solution

            velocity_field = @views hcat(q[ID_cached[1,:]], q[ID_cached[2,:]], q[ID_cached[3,:]])' # reshape the solution to get the velocity field
            dvdη = @views hcat(dqdη[ID_cached[1,:]], dqdη[ID_cached[2,:]], dqdη[ID_cached[3,:]])'
            dvdβ = @views hcat(dqdβ[ID_cached[1,:]], dqdβ[ID_cached[2,:]], dqdβ[ID_cached[3,:]])'
            
            motion = velocity_field*t_steps_cached # extract the motion of the mesh grid
            dmotiondη = dvdη*t_steps_cached
            dmotiondβ = dvdβ*t_steps_cached

            nodeList_cached = nodeList_cached + motion # update the mesh grid
            mdl.mesh_x.NodeList = nodeList_cached      # update the mesh grid
            dNodeList_dη += dmotiondη
            dNodeList_dβ += dmotiondβ
            
            mat_nan_inf_check(dvdη)
            mat_nan_inf_check(dvdβ)
            
            dmdθ_out = @views cat(dNodeList_dη,dNodeList_dβ,dims=3) # concatenate the gradients in to a tensor

            BorderPts2D, dudθ, surface_pts_2d, _, obs_border_pts = _get_2D_data(nodeList_cached, camera_matrix_cached, obj_pose_cached, h_cached, BorderNodesList=side_node_list_cached, GRAD=true, dqdθ=dmdθ_out, angles=rot_angle_cached)
            
            push!(output, μ_tp*t_steps_cached) # store displacement at the top surface
            push!(velocity, velocity_field) # store the velocity of the mesh in 3D
            push!(pressure, p) # store the pressure of the mesh in 3D
            push!(displacement, motion)
            push!(surface_fields, motion[:,side_node_list_cached])
            push!(surface_pts_3D, vcat(nodeList_cached[:,top_node_list_cached]', nodeList_cached[:,bottom_node_list_cached]', nodeList_cached[:,side_node_list_cached]')')
            push!(gradList,dudθ)
            push!(gradList_3d, dmdθ_out)
            push!(pos2D, surface_pts_2d)
            push!(pos3D, nodeList_cached)
            push!(pos3D_cp, nodeList_cached)
            push!(borderPts2DList, BorderPts2D)
            push!(splinep, BorderPts2D[1,:])
            push!(splineq, BorderPts2D[2,:])
            push!(writeborderList, obs_border_pts)

            iter += 1
            next!(pr, showvalues = [(:iterations,iter),(:time,string(t," seconds"))])
        end
    elseif control_cached == "velocity"
        for t in time
            _A_bar .= assemble_system_A(mdl, cache)    # assemble the stiffness matrix
            B .= assemble_system_B(mdl, cache)         # assemble the stiffness matrix
            b = apply_boundary_conditions(mdl, cache) # apply the neumann boundary conditions
            q_d .= (μu_btm*q_d_cached_btm + cParam_cached[iter]*q_d_cached_top + μu_side*q_d_cached_brdr)      # apply the Dirichlet boundary conditions

            if viscosity_type_cached == "bulk_viscosity"
                if length(β_cached) == 1
                    A .= η_cached[iter]*_A_bar + β_cached[1]*b
                else
                    A .= η_cached[iter]*_A_bar + β_cached[iter]*b
                end
            else
                A .= η_cached[1]*_A_bar + β_cached[1]*b
            end

            dAdη .= _A_bar
            dAdβ .= b
            dB .= zeros(Float64, size(B))

            C_Tu = transpose(C_uc_cached)      # transpose the constraint matrix
        
            A_free .= C_Tu*A*C_uc_cached        # extract the free part of the stiffness matrix
            B_free .= C_Tu*B                    # extract the free part of the stiffness matrix

            dA_freedη .= C_Tu*dAdη*C_uc_cached        # extract the free part of the stiffness matrix
            dA_freedβ .= C_Tu*dAdβ*C_uc_cached        # extract the free part of the stiffness matrix
            dB_free .= zeros(Float64, size(B_free))
        
            K_free = [A_free B_free; B_free' zeros(Float64, size(B_free,2),size(B_free,2))]      # assemble the system of equations
            dKdη = [C_Tu*dAdη*C_uc_cached dB_free; dB_free' zeros(Float64, size(B,2),size(B,2))] # assemble the system of equations
            dKdβ = [C_Tu*dAdβ*C_uc_cached dB_free; dB_free' zeros(Float64, size(B,2),size(B,2))] # assemble the system of equations
            
            luk = lu(K_free) # LU decomposition of the system of equations
        
            r = [C_Tu*A*q_d; B'*q_d]    # assemble the system of equations
            drdη = [C_Tu*dAdη*q_d; zero_matrix_p_q] # solve the system of equations
            drdβ = [C_Tu*dAdβ*q_d; zero_matrix_p_q] # solve the system of equations

            sol = luk\-Matrix(r)                    # solve the system of equations
            dsoldη = luk\-(drdη + dKdη*sol) # solve the system of equations
            dsoldβ = luk\-(drdβ + dKdβ*sol) # solve the system of equations
        
            q_f = sol[1:size(A_free,1)]      # extract the free part of the solution
            dqfdη = dsoldη[1:size(A_free,1)] # extract the free part of the solution
            dqfdβ = dsoldβ[1:size(A_free,1)] # extract the free part of the solution

            p_f = sol[size(A_free,1)+1:end]      # extract the free part of the solution
            dpfdη = dsoldη[size(A_free,1)+1:end] # extract the free part of the solution 
            dpfdβ = dsoldβ[size(A_free,1)+1:end] # extract the free part of the solution
        
            q = q_d + C_uc_cached*q_f;         # assemble the solution 
            dqdη = C_uc_cached*dqfdη;   # assemble the solution
            dqdβ = C_uc_cached*dqfdβ;   # assemble the solution

            p = p_f';       # assemble the solution
            dpdη = dpfdη';  # assemble the solution
            dpdβ = dpfdβ';  # assemble the solution

            velocity_field = @views hcat(q[ID_cached[1,:]], q[ID_cached[2,:]], q[ID_cached[3,:]])' # reshape the solution to get the velocity field
            dmdη = @views hcat(dqdη[ID_cached[1,:]], dqdη[ID_cached[2,:]], dqdη[ID_cached[3,:]])'*t_steps_cached
            dmdβ = @views hcat(dqdβ[ID_cached[1,:]], dqdβ[ID_cached[2,:]], dqdβ[ID_cached[3,:]])'*t_steps_cached
            
            motion_y = velocity_field*t_steps_cached # extract the motion of the mesh grid
            motion =  motion_y # extract the motion of the mesh grid

            nodeList_cached = nodeList_cached + motion # update the mesh grid
            mdl.mesh_x.NodeList = nodeList_cached      # update the mesh grid
            dNodeList_dη += dmdη
            dNodeList_dβ += dmdβ
  
            mat_nan_inf_check(dmdη)
            mat_nan_inf_check(dmdβ)

            dmdθ_out = @views cat(dNodeList_dη,dNodeList_dβ,dims=3) # concatenate the gradients in to a tensor

            BorderPts2D, dudθ, surface_pts_2d, _, obs_border_pts = _get_2D_data(nodeList_cached, camera_matrix_cached, obj_pose_cached, h_cached, BorderNodesList=side_node_list_cached, GRAD=true, dqdθ=dmdθ_out, angles=rot_angle_cached)

            # push!(output, μ_tp*t_steps_cached) # store displacement at the top surface
            push!(velocity, velocity_field) # store the velocity of the mesh in 3D
            push!(pressure, p) # store the pressure of the mesh in 3D
            push!(displacement, motion)
            push!(surface_fields, motion[:,side_node_list_cached])
            push!(surface_pts_3D, vcat(nodeList_cached[:,top_node_list_cached]', nodeList_cached[:,bottom_node_list_cached]', nodeList_cached[:,side_node_list_cached]')')
            push!(gradList,dudθ)
            push!(pos2D, surface_pts_2d)
            push!(pos3D, nodeList_cached)
            push!(pos3D_cp, nodeList_cached)
            push!(borderPts2DList, BorderPts2D)
            push!(splinep, BorderPts2D[1,:])
            push!(splineq, BorderPts2D[2,:])
            push!(writeborderList, obs_border_pts)
        
            iter += 1
            next!(pr, showvalues = [(:iterations,iter),(:time,string(t," seconds"))])
        end
    else
            throw(ArgumentError("Control type is unknown"))
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
        animate_fields(filepath = string(conditions.filepath,"/Results/images/"), Nodes=pos3D , IEN=IEN_u_cached, border_nodes_2d=borderPts2DList, sim_pts_2d=pos2D)
        animate_fields(filepath = string(conditions.filepath,"/Results/images/surface"), Nodes=surface_pts_3D)
    end
    
    return output, gradList, borderPts2DList, displacement, surface_pts_3D, pos2D, pos3D, splinep, splineq, velocity, pressure, gradList_3d
end