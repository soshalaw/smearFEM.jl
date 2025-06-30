using LinearAlgebra
using SparseArrays
using ProgressMeter
using Parameters
""" 
    assemble_system_A(mdl::Stokes)

Assembles the finite element system. # Returns the global stiffness matrix

# Arguments:
- `ne::Interger`: number of elements in each direction
- `NodeList::Matrix{Float64}{ndim,nNodes}` : coordinates of the nodes
- `IEN::Matrix{Int}{nElements,nLocalNodes}` : connectivity matrix
- `ndim::Interger`: number of dimensions
- `nDof::Interger`: number of degree of freedom per node
- `FunctionClass::String`: type of basis functions to be considered (Q1:quadratic or Q2:Lagrange)
- `ID::Matrix{Int}{nDof, nNodes}` : matrix that maps the global degrees of freedom to the local degrees of freedom
- `Young::Float64`: Young's modulus
- `ν::Float64`: Poisson's ratio

# Returns:
- `K::SparseMatrixCSC{Float64,Int64}{ndof,ndof}` : sparse stiffness matrix 
"""
function assemble_system_A(mdl::Stokes)::SparseMatrixCSC{Float64,Int64}
    # unpack the model parameters to local variables
    # this is done to avoid the need to pass the model object around
    @unpack ne, ndim, nDof_u = mdl
    @unpack NodeList, IEN, ID, FunctionClass = mdl.mesh_u

    NodeList_cached::Matrix{Float64} = NodeList
    IEN_cached::Matrix{Int} = IEN
    ID_cached::Matrix{Int} = ID
    ne_cached::Int = ne
    ndim_cached::Int = ndim
    nDof_u_cached::Int = nDof_u
    FunctionClass_cached::String = FunctionClass

    C::Matrix{Float64} = get_cMat(1.0,0.0,type="standard")
    IEN_u_rows::Int = size(IEN_cached,1)
  
    # (I,J,V) vectors for COO sparse matrix
    if nDof_u_cached == 1
        E = zeros(  Int64, ne_cached^ndim*IEN_u_rows^2)
        J = zeros(  Int64, ne_cached^ndim*IEN_u_rows^2)
        V = zeros(Float64, ne_cached^ndim*IEN_u_rows^2)
    else
        ID_u_rows::Int = size(ID_cached,1)
        E = zeros(  Int64, ne_cached^ndim*((ID_u_rows*IEN_u_rows)^2))
        J = zeros(  Int64, ne_cached^ndim*((ID_u_rows*IEN_u_rows)^2))
        V = zeros(Float64, ne_cached^ndim*((ID_u_rows*IEN_u_rows)^2))  
    end

    # element loop
    if ndim_cached == 1
        # gaussian quadrature points for the element [-1,1] 
        ξ, w_ξ = gaussian_quadrature(-1,1)
        
        wpoints =  [w_ξ[1], w_ξ[2]]
        
        x = [ξ[1], ξ[2]]
    elseif ndim_cached == 2
        # gaussian quadrature points for the element [-1,1]x[-1,1] 
        ξ, w_ξ = gaussian_quadrature(-1,1,nGaussPoints=3)
        η, w_η = gaussian_quadrature(-1,1,nGaussPoints=3)

        x = Float64[]
        y = Float64[]
        wpoints =  Float64[]
        
        n = 1:size(ξ,1)
        m = 1:size(η,1)
        for j::Int in m # loop over η
            for i::Int in n # loop over ξ
                push!(x, ξ[i])
                push!(y, η[j])
                push!(wpoints, w_ξ[i]*w_η[j])
            end
        end
    elseif ndim_cached == 3
        # gaussian quadrature points for the element [-1,1]x[-1,1]x[-1,1] 
        ξ, w_ξ = gaussian_quadrature(-1,1,nGaussPoints=3)
        η, w_η = gaussian_quadrature(-1,1,nGaussPoints=3)
        ζ, w_ζ = gaussian_quadrature(-1,1,nGaussPoints=3)

        x = Float64[]
        y = Float64[]
        z = Float64[]
        wpoints = Float64[]
        
        l = 1:size(ζ,1)
        m = 1:size(η,1)
        n = 1:size(ξ,1)
        for k::Int in l # loop over ζ
            for j::Int in m # loop over η
                for i::Int in n # loop over ξ
                    push!(x, ξ[i])
                    push!(y, η[j])
                    push!(z, ζ[k])
                    push!(wpoints, w_ξ[i]*w_η[j]*w_ζ[k])
                end
            end
        end
    end
    # initialize the Jacobian matrix and its inverse
    Jac = zeros(Float64, ndim_cached, ndim_cached)  # Jacobian matrix
    invJ = similar(Jac) # Inverse of the Jacobian matrix

    e_iter = 1:ne_cached^ndim_cached # iterator for elements loop 
    gpiter = 1:length(wpoints) # iterator for integration loop
    # integration loop
    for gp::Int in gpiter
        if ndim_cached == 1
            N, ΔN = basis_function(x[gp], nothing, nothing, FunctionClass_cached) # add function type
        elseif ndim_cached == 2
            N, ΔN = basis_function(x[gp], y[gp], nothing, FunctionClass_cached) 
        elseif ndim_cached == 3
            N, ΔN = basis_function(x[gp], y[gp], z[gp], FunctionClass_cached) 
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
        Ke_len::Int = length(Ke)  

        # element loop
        for e::Int in e_iter
            coords::Matrix{Float64} = NodeList_cached[:, IEN_cached[:, e]]  # Get the coordinates of the nodes of the element

            mul!(Jac, coords, ΔN)  # Jacobian matrix [dx/dxi dx/deta; dy/dxi dy/deta]
            w::Float64 = wpoints[gp] * abs(det(Jac))
            invJ .= inv(Jac)  # Inverse of the Jacobian matrix
            mul!(dNdX, ΔN, invJ)

            if nDof_u_cached == 1
                szN::Int = size(N, 1)  # Number of basis functions
                # Loop between basis functions of the element
                for i::Int in 1:szN
                    for j::Int in 1:szN
                        inz::Int = (szN)^2 * (e - 1) + szN * (i - 1) + j  # Index for the COO sparse matrix
                        E[inz] = IEN_cached[i, e]  # Row index
                        J[inz] = IEN_cached[j, e]  # Column index
                        V[inz] += w * dot(dNdX[i, :], dNdX[j, :])  # Inner product of the gradient of the basis functions
                    end
                end
            else
                if nDof_u_cached == 2
                    B .= 0.0
                    B[1:nDof_u_cached:end, 1] = dNdX[:, 1]
                    B[2:nDof_u_cached:end, 2] = dNdX[:, 2]
                    B[1:nDof_u_cached:end, 3] = dNdX[:, 2]
                    B[2:nDof_u_cached:end, 3] = dNdX[:, 1]
                elseif nDof_u_cached == 3
                    B .= 0.0
                    B[1:nDof_u_cached:end, 1] = dNdX[:, 1]
                    B[2:nDof_u_cached:end, 2] = dNdX[:, 2]
                    B[3:nDof_u_cached:end, 3] = dNdX[:, 3]
                    B[2:nDof_u_cached:end, 4] = dNdX[:, 3]
                    B[3:nDof_u_cached:end, 4] = dNdX[:, 2]
                    B[1:nDof_u_cached:end, 5] = dNdX[:, 3]
                    B[3:nDof_u_cached:end, 5] = dNdX[:, 1]
                    B[1:nDof_u_cached:end, 6] = dNdX[:, 2]
                    B[2:nDof_u_cached:end, 6] = dNdX[:, 1]
                end
                mul!(Ke, 2*w*B, C*B')  # Element stiffness matrix

                # Loop between basis functions of the element
                iNodes = 1:div(Ke_row, nDof_u_cached)
                jNodes = 1:div(Ke_col, nDof_u_cached)
                iDofs = 1:ID_u_rows
                jDofs = 1:ID_u_rows

                for iNode::Int in iNodes
                    for jNode::Int in jNodes
                        for iDof::Int in iDofs
                            for jDof::Int in jDofs
                                i::Int = (iNode - 1) * nDof_u_cached + iDof
                                j::Int = (jNode - 1) * nDof_u_cached + jDof
                                inz::Int = Ke_len * (e - 1) + (iNode - 1) * nDof_u_cached * Ke_col + (jNode - 1) * nDof_u_cached^2 + (iDof - 1) * nDof_u_cached + jDof  # Index for the COO sparse matrix
                                E[inz] = ID_cached[iDof, IEN_cached[iNode, e]]  # Row index
                                J[inz] = ID_cached[jDof, IEN_cached[jNode, e]]  # Column index
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

function assemble_system_A_dense(mdl::Stokes)::Matrix{Float64}
    C::Matrix{Float64} = get_cMat(1.0,0.0,type="standard")
    IEN_u_rows::Int = size(IEN,1)
    ID_u_rows::Int = size(ID,1)

    if nDof_u == 1
        sz = ne^ndim*IEN_u_rows # no elements x no nodes
        K = zeros(Float64, sz,sz)
    else
        sz = nDof_u*(nNodes)^ndim # no elements x no nodes x no dofs
        K = zeros(Float64, sz,sz)  
    end

    if ndim == 1
        # gaussian quadrature points for the element [-1,1] 
        ξ, w_ξ = gaussian_quadrature(-1,1)
        
        wpoints =  [w_ξ[1], w_ξ[2]]
        
        x = [ξ[1], ξ[2]]
    elseif ndim == 2
        # gaussian quadrature points for the element [-1,1]x[-1,1] 
        ξ, w_ξ = gaussian_quadrature(-1,1,nGaussPoints=3)
        η, w_η = gaussian_quadrature(-1,1,nGaussPoints=3)

        x = Float64[]
        y = Float64[]
        wpoints =  Float64[]
        
        n = 1:size(ξ,1)
        m = 1:size(η,1)
        for j::Int in m # loop over η
            for i::Int in n # loop over ξ
                push!(x, ξ[i])
                push!(y, η[j])
                push!(wpoints, w_ξ[i]*w_η[j])
            end
        end

    elseif ndim == 3
        # gaussian quadrature points for the element [-1,1]x[-1,1]x[-1,1] 
        ξ, w_ξ = gaussian_quadrature(-1,1,nGaussPoints=3)
        η, w_η = gaussian_quadrature(-1,1,nGaussPoints=3)
        ζ, w_ζ = gaussian_quadrature(-1,1,nGaussPoints=3)

        x = Float64[]
        y = Float64[]
        z = Float64[]
        wpoints = Float64[]
        
        l = 1:size(ζ,1)
        m = 1:size(η,1)
        n = 1:size(ξ,1)
        for k::Int in l # loop over ζ
            for j::Int in m # loop over η
                for i::Int in n # loop over ξ
                    push!(x, ξ[i])
                    push!(y, η[j])
                    push!(z, ζ[k])
                    push!(wpoints, w_ξ[i]*w_η[j]*w_ζ[k])
                end
            end
        end
    end
    # initialize the Jacobian matrix and its inverse
    Jac = zeros(Float64, ndim_cached, ndim_cached)  # Jacobian matrix
    invJ = similar(Jac) # Inverse of the Jacobian matrix

    e_iter = 1:ne^ndim # iterator for elements loop
    gpiter = 1:length(wpoints) # iterator for integration loop
    # integration loop
    for gp::Int in gpiter
        if ndim == 1
            N, ΔN = basis_function(x[gp], nothing, nothing, FunctionClass)
        elseif ndim == 2
            N, ΔN = basis_function(x[gp], y[gp], nothing, FunctionClass) 
        elseif ndim == 3
            N, ΔN = basis_function(x[gp], y[gp], z[gp], FunctionClass) 
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
            coords::Matrix{Float64} = NodeList_cached[:, IEN_cached[:, e]]  # Get the coordinates of the nodes of the element

            mul!(Jac, coords, ΔN)  # Jacobian matrix [dx/dxi dx/deta; dy/dxi dy/deta]
            w::Float64 = wpoints[gp] * abs(det(Jac))
            invJ .= inv(Jac)  # Inverse of the Jacobian matrix
            mul!(dNdX, ΔN, invJ)

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

function assemble_system_B(mdl::Stokes)::SparseMatrixCSC{Float64,Int64}
    # unpack the model parameters to local variables
    # this is done to avoid the need to pass the model object around
    @unpack ne, ndim, nDof_u = mdl
    @unpack IEN, FunctionClass = mdl.mesh_p

    IEN_p_cached::Matrix{Int} = IEN
    FunctionClass_p_cached::String = FunctionClass

    @unpack NodeList, IEN, ID, FunctionClass = mdl.mesh_u

    IEN_u_cached::Matrix{Int} = IEN
    NodeList_cached::Matrix{Float64} = NodeList
    ID_cached::Matrix{Int} = ID
    ne_cached::Int = ne
    ndim_cached::Int = ndim
    nDof_u_cached::Int = nDof_u
    FunctionClass_u_cached::String = FunctionClass

    # (I,J,V) vectors for COO sparse matrix
    IEN_u_rows::Int = size(IEN_u_cached,1)
    IEN_p_rows::Int = size(IEN_p_cached,1)
    ID_u_rows::Int = size(ID_cached,1)

    if nDof_u_cached == 1
        E = zeros(  Int64, ne_cached^ndim_cached*IEN_u_rows*IEN_p_rows)
        J = zeros(  Int64, ne_cached^ndim_cached*IEN_u_rows*IEN_p_rows)
        V = zeros(Float64, ne_cached^ndim_cached*IEN_u_rows*IEN_p_rows)
    else
        E = zeros(  Int64, ne_cached^ndim_cached*(ID_u_rows*IEN_u_rows)*IEN_p_rows)
        J = zeros(  Int64, ne_cached^ndim_cached*(ID_u_rows*IEN_u_rows)*IEN_p_rows)
        V = zeros(Float64, ne_cached^ndim_cached*(ID_u_rows*IEN_u_rows)*IEN_p_rows)  
    end

    if ndim_cached == 1
        # gaussian quadrature points for the element [-1,1] 
        ξ, w_ξ = gaussian_quadrature(-1,1)
    
        wpoints = [w_ξ[1], w_ξ[2]]
        
        x = [ξ[1], ξ[2]]
    elseif ndim_cached == 2
        # gaussian quadrature points for the element [-1,1]x[-1,1] 
        ξ, w_ξ = gaussian_quadrature(-1,1,nGaussPoints=3)
        η, w_η = gaussian_quadrature(-1,1,nGaussPoints=3)

        x = Float64[]
        y = Float64[]
        wpoints = Float64[]
        
        n = 1:size(ξ,1)
        m = 1:size(η,1)
        for j in m # loop over η
            for i in n # loop over ξ
                push!(x, ξ[i])
                push!(y, η[j])
                push!(wpoints, w_ξ[i]*w_η[j])
            end
        end

    elseif ndim_cached == 3
        # gaussian quadrature points for the element [-1,1]x[-1,1]x[-1,1] 
        ξ, w_ξ = gaussian_quadrature(-1,1,nGaussPoints=3)
        η, w_η = gaussian_quadrature(-1,1,nGaussPoints=3)
        ζ, w_ζ = gaussian_quadrature(-1,1,nGaussPoints=3)

        x = Float64[]
        y = Float64[]
        z = Float64[]
        wpoints = Float64[]
        
        l = 1:size(ζ,1)
        m = 1:size(η,1)
        n = 1:size(ξ,1)
        for k::Int in l # loop over ζ
            for j::Int in m # loop over η
                for i::Int in n # loop over ξ
                    push!(x, ξ[i])
                    push!(y, η[j])
                    push!(z, ζ[k])
                    push!(wpoints, w_ξ[i]*w_η[j]*w_ζ[k])
                end
            end
        end
    end
    # initialize the Jacobian matrix and its inverse
    Jac = zeros(Float64, ndim_cached, ndim_cached)  # Jacobian matrix
    invJ = similar(Jac) # Inverse of the Jacobian matrix

    e_iter = 1:ne_cached^ndim_cached    # iterator for elements loop
    gpiter = 1:length(wpoints)  # iterator for integration loop
    # integration loop
    for gp::Int in gpiter
        if ndim_cached == 1
            N_u, ΔN_u = basis_function(x[gp], nothing, nothing, FunctionClass_u_cached)
            N_p, ΔN_p = basis_function(x[gp], nothing, nothing, FunctionClass_p_cached)
        elseif ndim_cached == 2
            N_u, ΔN_u = basis_function(x[gp], y[gp], nothing, FunctionClass_u_cached)
            N_p, ΔN_p = basis_function(x[gp], y[gp], nothing, FunctionClass_p_cached)
        elseif ndim_cached == 3
            N_u, ΔN_u = basis_function(x[gp], y[gp], z[gp], FunctionClass_u_cached)
            N_p, ΔN_p = basis_function(x[gp], y[gp], z[gp], FunctionClass_p_cached)
        end
        dNdX = zeros(Float64, size(ΔN_u, 1), ndim_cached) # Gradient of basis functions
        Be = zeros(Float64, 1, 3*size(dNdX, 1)) # Gradient of basis functions
        colB::Int = size(Be, 2)  # Number of basis functions
        rowNp::Int = size(N_p, 1)  # Number of basis functions
        Ke = zeros(Float64, rowNp, colB)  # Element stiffness matrix
        Ke_row::Int, Ke_col::Int = rowNp, colB
        # element loop
        for e::Int in e_iter
            coords_u::Matrix{Float64} = NodeList_cached[:, IEN_u_cached[:, e]]  # Get the coordinates of the nodes of the element

            mul!(Jac, coords_u, ΔN_u)  # Jacobian matrix [dx/dxi dx/deta; dy/dxi dy/deta]
            w::Float64 = wpoints[gp] * abs(det(Jac))
            invJ .= inv(Jac)  # Inverse of the Jacobian matrix
            mul!(dNdX, ΔN_u, invJ)

            if nDof_u_cached == 1
                szN = size(N_u,1) # number of basis functions
                # loop between basis functions of the element
                for i::Int in 1:szN
                    for j::Int in 1:szN
                        inz::Int = (szN)^2*(e-1) + szN*(i-1) + j # index for the COO sparse matrix
                        E[inz] = IEN_u_cached[i,e] # row index 
                        J[inz] = IEN_u_cached[j,e] # column index
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
                iNodes = 1:div(Ke_col,nDof_u_cached)  # column node index
                iDofs = 1:ID_u_rows # column dof index
                jNodes = 1:Ke_row # row index

                for iNode::Int in iNodes 
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

function assemble_system_B_dense(mdl::Stokes)::Matrix{Float64}

    @unpack ne, ndim, nDof_u = mdl
    @unpack IEN, FunctionClass = mdl.mesh_p

    IEN_p_cached::Matrix{Int} = IEN
    FunctionClass_p_cached::String = FunctionClass

    @unpack NodeList, IEN, ID, FunctionClass = mdl.mesh_p

    IEN_u_cached::Matrix{Int} = IEN
    NodeList_cached::Matrix{Float64} = NodeList
    ID_cached::Matrix{Int} = ID
    ne_cached::Int = ne
    ndim_cached::Int = ndim
    nDof_u_cached::Int = nDof_u
    FunctionClass_u_cached::String = FunctionClass

    IEN_u_rows = size(IEN_u_cached,1)
    ID_u_rows = size(ID_cached,1)
    if nDof_u_cached == 1
        sz = ne_cached^ndim_cached*IEN_u_rows # no elements x no nodes
        K = zeros(Float64, sz,sz)
    else
        sz1 = nDof_u_cached*(nNodes_cached)^ndim_cached # no elements x no nodes x no dofs
        sz2 = (nNodes_cached)^ndim_cached
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
        ξ, w_ξ = gaussian_quadrature(-1,1,nGaussPoints=3)
        η, w_η = gaussian_quadrature(-1,1,nGaussPoints=3)

        x = Float64[]
        y = Float64[]
        wpoints = Float64[]
        
        n = 1:size(ξ,1)
        m = 1:size(η,1)
        for j::Int in m # loop over η
            for i::Int in n # loop over ξ
                push!(x, ξ[i])
                push!(y, η[j])
                push!(wpoints, w_ξ[i]*w_η[j])
            end
        end

    elseif ndim_cached == 3
        # gaussian quadrature points for the element [-1,1]x[-1,1]x[-1,1] 
        ξ, w_ξ = gaussian_quadrature(-1,1,nGaussPoints=3)
        η, w_η = gaussian_quadrature(-1,1,nGaussPoints=3)
        ζ, w_ζ = gaussian_quadrature(-1,1,nGaussPoints=3)

        x = Float64[]
        y = Float64[]
        z = Float64[]
        wpoints = Float64[]
        
        l = 1:size(ζ,1)
        m = 1:size(η,1)
        n = 1:size(ξ,1)
        for k::Int in l # loop over ζ
            for j::Int in m # loop over η
                for i::Int in n # loop over ξ
                    push!(x, ξ[i])
                    push!(y, η[j])
                    push!(z, ζ[k])
                    push!(wpoints, w_ξ[i]*w_η[j]*w_ζ[k])
                end
            end
        end
    end
    # initialize the Jacobian matrix and its inverse
    Jac = zeros(Float64, ndim_cached, ndim_cached)  # Jacobian matrix
    invJ = similar(Jac) # Inverse of the Jacobian matrix

    e_iter = 1:ne_cached^ndim_cached    # iterator for elements loop
    gpiter = 1:length(wpoints) # iterator for integration loop
    # integration loop
    for gp::Int in gpiter
        if ndim_cached == 1
            N_u, ΔN_u = basis_function(x[gp], nothing, nothing, FunctionClass_u_cached)
            N_p, ΔN_p = basis_function(x[gp], nothing, nothing, FunctionClass_p_cached)
        elseif ndim_cached == 2
            N_u, ΔN_u = basis_function(x[gp], y[gp], nothing, FunctionClass_u_cached)
            N_p, ΔN_p = basis_function(x[gp], y[gp], nothing, FunctionClass_p_cached)
        elseif ndim_cached == 3
            N_u, ΔN_u = basis_function(x[gp], y[gp], z[gp], FunctionClass_u_cached)
            N_p, ΔN_p = basis_function(x[gp], y[gp], z[gp], FunctionClass_p_cached)
        end

        dNdX = zeros(Float64, size(ΔN_u, 1), ndim_cached) # Gradient of basis functions
        Be = zeros(Float64, 1, 3*size(dNdX, 1)) # Gradient of basis functions
        colB::Int = size(Be, 2)  # Number of basis functions
        rowNp::Int = size(N_p, 1)  # Number of basis functions
        Ke = zeros(Float64, rowNp, colB)  # Element stiffness matrix
        Ke_row::Int, Ke_col::Int = rowNp, colB

        # element loop
        for e::Int in e_iter
            coords_u::Matrix{Float64} = NodeList_cached[:, IEN_u_cached[:, e]]  # Get the coordinates of the nodes of the element

            mul!(Jac, coords_u, ΔN_u)  # Jacobian matrix [dx/dxi dx/deta; dy/dxi dy/deta]
            w::Float64 = wpoints[gp] * abs(det(Jac))
            invJ .= inv(Jac)  # Inverse of the Jacobian matrix
            mul!(dNdX, ΔN_u, invJ)

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

""" Apply the Neumann slip boundary conditions to the global stiffness matrix

# Arguments:
K: {[ndof,ndof] SparseMatrixCSC{Float64,Int64}} : sparse stiffness matrix 
ID: {[nDof,nNodes] Matrix{Int}} : matrix that maps the global degrees of freedom to the local degrees of freedom
q_d: {[ndof] Vector{Float64}} : Dirichlet boundary conditions
q_n: {[ndof] Vector{Float64}} : Neumann boundary conditions

# Returns:
K: {[ndof,ndof] SparseMatrixCSC{Float64,Int64}} : sparse stiffness matrix with the boundary conditions applied
F: {[ndof] Vector{Float64}} : force vector
"""
function apply_boundary_conditions(mdl::Stokes)::SparseMatrixCSC{Float64,Int64}

    @unpack ne, ndim, nDof_u = mdl
    @unpack NodeList, IEN_top, IEN_bottom, ID, FunctionClass = mdl.mesh_u

    NodeList_cached::Matrix{Float64} = NodeList
    IEN_top_cached::Matrix{Int} = IEN_top
    IEN_btm_cached::Matrix{Int} = IEN_bottom
    ID_cached::Matrix{Int} = ID
    ne_cached::Int = ne
    ndim_cached::Int = ndim
    nDof_u_cached::Int = nDof_u
    FunctionClass_cached::String = FunctionClass

    ID_u_rows::Int = size(ID_cached,1)
    IEN_btm_rows::Int = size(IEN_btm_cached,1)

    E = zeros(  Int64, ne_cached^(ndim_cached-1)*(ID_u_rows*IEN_btm_rows)^2*2) # *2 because we have two surfaces
    J = zeros(  Int64, ne_cached^(ndim_cached-1)*(ID_u_rows*IEN_btm_rows)^2*2) # *2 because we have two surfaces
    V = zeros(Float64, ne_cached^(ndim_cached-1)*(ID_u_rows*IEN_btm_rows)^2*2) # *2 because we have two surfaces

    if ndim_cached == 2
        # gaussian quadrature points for the element [-1,1] 
        ξ, w_ξ = gaussian_quadrature(-1,1)

        wpoints = [w_ξ[1], w_ξ[2]]
        
        x = [ξ[1], ξ[2]]
    elseif ndim_cached == 3
        # gaussian quadrature points for the element [-1,1]x[-1,1] 
        ξ, w_ξ = gaussian_quadrature(-1,1)
        η, w_η = gaussian_quadrature(-1,1)
        
        wpoints = [w_ξ[1]*w_η[1], w_ξ[2]*w_η[1], w_ξ[2]*w_η[2], w_ξ[1]*w_η[2]]
        
        x = [ξ[1], ξ[2], ξ[2], ξ[1]]
        y = [η[1], η[1], η[2], η[2]]
    end

    e_iter = 1:ne_cached^(ndim_cached-1)   # iterator for elements loop
    gpiter = 1:length(wpoints) # iterator for integration loop
    # integration loop
    for gp::Int in gpiter

        if ndim_cached == 2
            N, ΔN = basis_function(x[gp], nothing, nothing, FunctionClass_cached)
        elseif ndim_cached == 3
            N, ΔN = basis_function(x[gp], y[gp], nothing, FunctionClass_cached) 
        end

        M = zeros(Float64, 3, ndim_cached*length(N))
        rowM::Int = size(M,2)
        be_row::Int, be_col::Int = rowM, rowM
        be = zeros(Float64, be_col, be_row)
        len_be::Int = length(be)

        # element loop
        for e::Int in e_iter
        
            coords_u_top::Matrix{Float64} = NodeList_cached[:,IEN_top_cached[:,e]] # get the coordinates of the nodes of the element
            coords_u_btm::Matrix{Float64} = NodeList_cached[:,IEN_btm_cached[:,e]] # get the coordinates of the nodes of the element

            dxdξ_top = coords_u_top*ΔN         # Jacobian matrix [dx/dxi dx/deta; dy/dxi dy/deta; dz/dxi dz/deta]
            dxdξ_btm = coords_u_btm*ΔN         # Jacobian matrix [dx/dxi dx/deta; dy/dxi dy/deta; dz/dxi dz/deta]

            w_top::Float64 = wpoints[gp]*norm(cross(dxdξ_top[:,1],dxdξ_top[:,2]))     # weight of the quadrature point top surface
            w_btm::Float64 = wpoints[gp]*norm(cross(dxdξ_btm[:,1],dxdξ_btm[:,2]))     # weight of the quadrature point bottom surface
            
            M .= 0.0
            M[1,1:nDof_u_cached:end] = N
            M[2,2:nDof_u_cached:end] = N
            M[3,3:nDof_u_cached:end] = N

            # be = M'*M
            mul!(be,M',M) # multiply the matrix by itself to get the stiffness matrix

            # loop between basis functions of the element
            iNodes = 1:be_row÷nDof_u_cached
            jNodes = 1:be_col÷nDof_u_cached
            iDofs = 1:ID_u_rows
            jDofs = 1:ID_u_rows
            for iNode::Int in iNodes
                for jNode::Int in jNodes
                    for iDof::Int in iDofs
                        for jDof::Int in jDofs
                            i::Int = (iNode-1)*nDof_u_cached + iDof
                            j::Int = (jNode-1)*nDof_u_cached + jDof
                            
                            inz_btm::Int = len_be*(e-1) + (iNode-1)*nDof_u_cached*be_col + (jNode-1)*nDof_u_cached^2 + (iDof-1)*nDof_u_cached + jDof # index for the COO sparse matrix
                            inz_top::Int = len_be*(ne_cached)^(ndim_cached-1)+ len_be*((e-1)) + (iNode-1)*nDof_u_cached*be_col + (jNode-1)*nDof_u_cached^2 + (iDof-1)*nDof_u_cached + jDof # index for the COO sparse matrix
                            
                            E[inz_top] = ID_cached[iDof,IEN_top_cached[iNode,e]]    # row index 
                            J[inz_top] = ID_cached[jDof,IEN_top_cached[jNode,e]]   # column index
                            V[inz_top] += w_top*be[i,j] 

                            E[inz_btm] = ID_cached[iDof,IEN_btm_cached[iNode,e]]    # row index 
                            J[inz_btm] = ID_cached[jDof,IEN_btm_cached[jNode,e]]   # column index
                            V[inz_btm] += w_btm*be[i,j] 
                        end
                    end
                end
            end
            # TODO include in tha assembly function
        end
    end
    K = sparse(E,J,V)
    return  K
end

function apply_boundary_conditions_dense(mdl::Stokes)::Matrix{Float64}
    @unpack ne, ndim, nDof_u = mdl
    @unpack NodeList, IEN_top, IEN_bottom, ID, FunctionClass = mdl.mesh_u

    NodeList_cached::Matrix{Float64} = NodeList
    IEN_top_cached::Matrix{Int} = IEN_top
    IEN_btm_cached::Matrix{Int} = IEN_bottom
    ID_cached::Matrix{Int} = ID
    ne_cached::Int = ne
    ndim_cached::Int = ndim
    nDof_u_cached::Int = nDof_u
    FunctionClass_cached::String = FunctionClass

    IEN_btm_rows = size(IEN_btm_cached,1)
    IEN_rows = size(mdl.IEN,1)
    ID_u_rows = size(ID_cached,1)
    sz =    nDof_u_cached*(nNodes_cached)^ndim_cached
    K = zeros(Float64, sz,sz)

    if ndim_cached == 2
        # gaussian quadrature points for the element [-1,1] 
        ξ, w_ξ = gaussian_quadrature(-1,1)

        wpoints = [w_ξ[1], w_ξ[2]]
        
        x = [ξ[1], ξ[2]]
    elseif ndim_cached == 3
        # gaussian quadrature points for the element [-1,1]x[-1,1] 
        ξ, w_ξ = gaussian_quadrature(-1,1)
        η, w_η = gaussian_quadrature(-1,1)
        
        wpoints = [w_ξ[1]*w_η[1], w_ξ[2]*w_η[1], w_ξ[2]*w_η[2], w_ξ[1]*w_η[2]]
        
        x = [ξ[1], ξ[2], ξ[2], ξ[1]]
        y = [η[1], η[1], η[2], η[2]]
    end 
    e_iter = 1:ne_cached^(ndim_cached-1)  # iterator for elements loop
    gpiter = 1:length(wpoints)  # iterator for integration loop
    # integration loop
    for gp in gpiter

        if ndim_cached == 2
            N, ΔN = basis_function(x[gp], nothing, nothing, FunctionClass_cached)
        elseif ndim_cached == 3
            N, ΔN = basis_function(x[gp], y[gp], nothing, FunctionClass_cached) 
        end

        M = zeros(Float64, 3, ndim_cached*length(N))
        rowM::Int = size(M,1)
        be_row::Int, be_col::Int = rowM, rowM
        be = zeros(Float64, be_col, be_row)
        len_be::Int = length(be)

        # element loop
        for e in e_iter
        
            coords_u_top::Matrix{Float64} = NodeList_cached[:,IEN_top_cached[:,e]] # get the coordinates of the nodes of the element
            coords_u_btm::Matrix{Float64} = NodeList_cached[:,IEN_btm_cached[:,e]] # get the coordinates of the nodes of the element

            dxdξ_top = coords_u_top*ΔN         # Jacobian matrix [dx/dxi dx/deta; dy/dxi dy/deta; dz/dxi dz/deta]
            dxdξ_btm = coords_u_btm*ΔN         # Jacobian matrix [dx/dxi dx/deta; dy/dxi dy/deta; dz/dxi dz/deta]

            w_top::Float64 = wpoints[gp]*norm(cross(dxdξ_top[:,1],dxdξ_top[:,2]))     # weight of the quadrature point top surface
            w_btm::Float64 = wpoints[gp]*norm(cross(dxdξ_btm[:,1],dxdξ_btm[:,2]))     # weight of the quadrature point bottom surface
            
            M .= 0.0
            M[1,1:nDof_u_cached:end] = N
            M[2,2:nDof_u_cached:end] = N
            M[3,3:nDof_u_cached:end] = N

            # be = M'*M
            mul!(be, M', M)

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
                            
                            inz_btm = len_be*(e-1) + (iNode-1)*nDof_u_cached*be_col + (jNode-1)*nDof_u_cached^2 + (iDof-1)*nDof_u_cached + jDof # index for the COO sparse matrix
                            inz_top = len_be*(ne_cached)^(ndim_cached-1)+ len_be*((e-1)) + (iNode-1)*nDof_u_cached*be_col + (jNode-1)*nDof_u_cached^2 + (iDof-1)*nDof_u_cached + jDof # index for the COO sparse matrix
                            
                            E[inz_top] = ID_cached[iDof,IEN_top_cached[iNode,e]]    # row index 
                            J[inz_top] = ID_cached[jDof,IEN_top_cached[jNode,e]]   # column index
                            V[inz_top] += w_top*be[i,j] 

                            E[inz_btm] = ID_cached[iDof,IEN_btm_cached[iNode,e]]    # row index 
                            J[inz_btm] = ID_cached[jDof,IEN_btm_cached[jNode,e]]   # column index
                            V[inz_btm] += w_btm*be[i,j] 
                        end
                    end
                end
            end
            # TODO include in tha assembly function
        end
    end
    K = sparse(E,J,V)
    return  K
end

"""
Set the Dirichlet boundary conditions for the problem

# Arguments:
- `NodeList::Matrix{Float64}{nNodes,ndim}` : array of nodes
- `ne::Int` : number of elements
- `ndim::Int` : number of dimensions
- `FunctionClass::String` : type of basis function

# Returns:
- `q_upper::Vector{Float64}` : vector of the Dirichlet boundary conditions (for ndof = 1) / Dirichlet boundary conditions upper surface (for ndof > 1)
- `q_lower::Vector{Float64}` : vector of the Neumann boundary conditions (for ndof = 1) / Dirichlet boundary conditions lower surface (for ndof > 1)
- `C_uc_cached::SparseMatrixCSC{Float64,Int64}` : onstraint matrix
"""
function set_boundary_cond(mdl::Stokes; DENSE::Bool=false)

    @unpack ndim, nDof_u = mdl
    @unpack NodeList, nNodes = mdl.mesh_u

    NodeList_cached::Matrix{Float64} = NodeList
    ndim_cached::Int = ndim
    nDof_u_cached::Int = nDof_u
    nNodes_cached::Int = nNodes

    q_upper = zeros(Float64, nDof_u_cached*(nNodes_cached)^ndim_cached,1)                  # initialize the vector of the Dirichlet boundary conditions (for ndof = 1) / Dirichlet boundary conditions upper surface (for ndof > 1)
    q_lower = zeros(Float64, nDof_u_cached*(nNodes_cached)^ndim_cached,1)                  # initialize the vector of the Neumann boundary conditions (for ndof = 1) / Dirichlet boundary conditions lower surface (for ndof > 1)
    q_side = zeros(Float64, nDof_u_cached*(nNodes_cached)^ndim_cached,1)   

    if DENSE                # initialize the vector of the Neumann boundary conditions (for ndof = 1) / Dirichlet boundary conditions lower surface (for ndof > 1)
        C = Matrix{Float64}(I,nDof_u_cached*(nNodes_cached)^ndim_cached,nDof_u_cached*(nNodes_cached)^ndim_cached)      # definition of the constraint matrix
    else                 # initialize the vector of the Neumann boundary conditions (for ndof = 1) / Dirichlet boundary conditions lower surface (for ndof > 1)
        C = sparse(I,nDof_u_cached*(nNodes_cached)^ndim_cached,nDof_u_cached*(nNodes_cached)^ndim_cached)      # definition of the constraint matrix
    end

    if nDof_u_cached == 1
        if ndim_cached == 3
            Dbound1 = 0
            Dbound2 = 1
            iter = 1:size(NodeList_cached,2)
            for n::Int in iter
                coord = NodeList_cached[:,n] # get the coordinates of the node
                if coord[3] == Dbound1 # bottom boundary
                    q_upper[n] = 0
                elseif coord[3] == Dbound2 # top boundary
                    q_upper[n] = -1
                end
            end
        elseif ndim_cached == 2
            Dbound1 = 0
            Dbound2 = 1
            iter = 1:size(NodeList_cached,2)
            for n::Int in iter
                coord = NodeList_cached[:,n] # get the coordinates of the node
                if coord[2] == Dbound1 # bottom boundary
                    q_upper[n] = 0
                elseif coord[2] == Dbound2 # top boundary
                    q_upper[n] = -1
                end
            end
        end

        C_uc = C[:,((nNodes_cached)^(ndim_cached-1)+1):((nNodes_cached)^ndim_cached-(nNodes_cached)^(ndim_cached-1))]

    else
        z0Bound = 0
        z1Bound = 50

        rCol = Array{Int}(undef,0)
        iter = 1:size(NodeList_cached,2)
        for nNode::Int in iter
            coord = NodeList_cached[:,nNode]    # get the coordinates of the node
            if coord[3] == z1Bound   # top boundary
                q_upper[3*nNode] = 1     # constraint the z displacement to be -d
                push!(rCol,3*nNode)
            elseif coord[3] == z0Bound   # bottom boundary
                q_lower[3*nNode] = 1     # constraint the z displacement to be -d
                push!(rCol,3*nNode)
            end
        end

        C_uc = C[:,setdiff(1:size(C,2),rCol)]

    end
    return q_upper, q_side, q_lower, C_uc
end

"""
    get_η(t::Float64, F::Float64, R_0::Float64, H_0::Float64, η_0::Float64, n::Float64, K::Float64)
Calculate the shear viscosity of the fluid.
# Arguments:
- `t` : time
- `F::Float64` : force
- `R_0::Float64` : initial radius
- `H_0::Float64` : initial height
- `η_0::Float64` : initial viscosity
- `n::Float64` : power law index
- `K::Float64` : consistency index

# Returns:
- `η::Float64` : shear viscosity
"""
function get_η(t::Number, F::Number, R_0::Number, H_0::Number, η_0::Number, n::Number, K::Number)

    H(t) = H_0*(1+8*H_0^2*F*t/(3*π*η_0*R_0^4))^(-1/4) # height with time
    
    R(t) = R_0*(1+8*H_0^2*F*t/(3*π*η_0*R_0^4))^(1/8) # radius with time

    H_dot(t) = 8/3*(-2*F*H(t)^3/(8*π*η_0*R(t)^4)) # rate of change of height

    γ_dot(t) = H_dot(t)/H(t) # shear rate

    η(t) = K*(abs(γ_dot(t)))^(n-1) # shear viscosity

    return η(t)
end

function def_problem(r::Number, h::Number, ne::Int64, η_0::Float64, ndim::Int64, FunctionClass_u::String, nDof_u::Int64, FunctionClass_p::String, 
                    nDof_p::Int64, β::Float64, F::Float64, control::String, viscosity_type::String, sim_time::Float64, t_steps::Float64)
    n::Float64 = 0.5
    K::Float64 = 2.0
    len_t::Int = round(Int,(sim_time/t_steps)) # number of time steps
    time = collect(Float64, range(start=t_steps, stop=sim_time, step=t_steps))

    if control == "force"
        cParam = -F*ones(Float64, len_t)
        if viscosity_type == "bulk_viscosity"
            η = get_η.(time, F, r, h, η_0, n, K)
        else
            η = [η_0]
        end 
    elseif control == "velocity"
        cParam = -0.02*ones(Float64, len_t)
        η = [η_0]
    else
        throw(ArgumentError("Control type not recognized"))
    end

    # define the model
    stokes = set_model(r, h, ne, η, ndim, FunctionClass_u, nDof_u, FunctionClass_p, nDof_p) # define the model    
    # set the boundary conditions
    q_tp, q_side, q_btm, C_uc = set_boundary_cond(stokes) 
    # define with problem scenario
    squeeze = SqueezeFlow(stokes, [β], [q_tp, q_side, q_btm], C_uc, control, sim_time, t_steps, viscosity_type, cParam)
return stokes, squeeze
end

"""
set_model(ne, NodeList, IEN, ndim, FunctionClass=FunctionClass_p_cached, nDof=1, ID=nothing, Young=1, ν=0.3)

Defines the model for the finite element method.
# Arguments:
- `r::Number`: radius of the cylinder
- `h::Number`: height of the cylinder
- `ne::Int64`: number of elements in each direction
- `η::Number`: dynamic viscosity
- `ndim::Int64`: number of dimensions
- `FunctionClass_u::String`: type of basis functions to be considered for the velocity field (Q1:quadratic or Q2:Lagrange)
- `nDof_u::Int64`: number of degree of freedom per node for the velocity field
- `FunctionClass_p::String`: type of basis functions to be considered for the pressure field (Q1:quadratic or Q2:Lagrange)
- `nDof_p::Int64`: number of degree of freedom per node for the pressure field

# Returns:
- `mdl::Stokes`: model for the finite element method
""" 
function set_model(r::Number, h::Number, ne::Int64, η::Vector{Float64}, ndim::Int64, FunctionClass_u::String, nDof_u::Int64, FunctionClass_p::String, 
                nDof_p::Int64)::Stokes

    mesh_u = meshgrid_cylinder(r, h, ne, FunctionClass=FunctionClass_u)  # generate the mesh grid
    mesh_p = meshgrid_cylinder(r, h, ne, FunctionClass=FunctionClass_p)  # generate the mesh grid

    mdl = Stokes(ndim=ndim, mesh_u=mesh_u, nDof_u=nDof_u, mesh_p=mesh_p, nDof_p=nDof_p, η=η)

    return mdl
end

"""
simulate(x0, x1, y0, y1, z0, z1, ne, η ,ndim, FunctionClass_u, nDof_u, FunctionClass_p, nDof_p, β, camera_matrix, endTime, t_steps, Control, cParam; WRITEDATA=false, filepath=nothing, SIDES::Bool=false)

Simulate the Stokes problem for a given mesh over a given time period. 

# Arguments:
- `x0::Float64` : x coordinate of the lower left corner of the domain
- `x1::Float64` : x coordinate of the upper right corner of the domain
- `y0::Float64` : y coordinate of the lower left corner of the domain
- `y1::Float64` : y coordinate of the upper right corner of the domain
- `z0::Float64` : z coordinate of the lower left corner of the domain
- `z1::Float64` : z coordinate of the upper right corner of the domain
- `ne::Int` : number of elements
- `η::Float64` : viscosity of the fluid
- `ndim::Int` : number of dimensions
- `FunctionClass_u::String` : type of basis function for the velocity field
- `nDof_u::Int` : number of degrees of freedom for the velocity field
- `FunctionClass_p::String` : type of basis function for the pressure field
- `nDof_p::Int` : number of degrees of freedom for the pressure field
- `β::Float64` : parameter for the pressure field
- `camera_matrix::Matrix{Float64}` : camera matrix
- `endTime::Float64` : end time of the simulation
- `t_steps::Float64` : time step of the simulation
- `Control::String` : type of control for the simulation
- `cParam::Float64` : parameter for the control

# Optional arguments:
- `WRITEDATA::Bool` : write the data to a file
- `filepath::String` : path to write the data
- `SIDES::Bool` : if true, only the sides of the cylinder are extracted in the simulation. Simulates the partial observability of the borders.

# Returns:
- `output::Vector{Float64}` : output of the simulation
- `gradList`::Vector{Float64}` : gradient of the solution
- `borderPts2DList::Vector{Float64}` : border points in 2D 
- `splinep::Vector{Float64}` : x coordinates samples of the spline parameters of the border nodes
- `splineq::Vector{Float64}` : y coordinates samples of the spline parameters of the border nodes
- `mdl::Stokes` : model of the simulation
- `pos2D::Vector{Float64}` : position of the mesh in 2D
"""
function simulate(mdl::Stokes, scene::SqueezeFlow, conditions::Conditions)
    @unpack ID, NodeList, top_nodes, bottom_nodes, side_nodes, nNodes = mdl.mesh_u
    @unpack η, nDof_u, nDof_p, ndim = mdl
    @unpack β, viscosity_type, q_d, C_uc, t_steps, sim_time, control = scene
    @unpack camera_matrix, camera_pose, SIDES = conditions

    β_cached::Float64 = β[1]
    viscosity_type_cached::String = viscosity_type
    q_d_cached::Dict{Symbol, Matrix{Float64}} = q_d
    q_d_cached_top::Matrix{Float64} = q_d_cached[:top]
    q_d_cached_btm::Matrix{Float64} = q_d_cached[:bottom]
    q_d_cached_brdr::Matrix{Float64} = q_d_cached[:border]
    C_uc_cached::AbstractMatrix = C_uc
    t_steps_cached::Float64 = t_steps
    sim_time_cached::Float64 = sim_time
    control_cached::String = control

    nNodes_u_cached::Int = nNodes
    nDof_u_cached::Int = nDof_u
    nDof_p_cached::Int = nDof_p
    ndim_cached::Int = ndim
    NodeList_cached::Matrix{Float64} = NodeList
    top_node_list_cached::Vector{Int} = top_nodes # top nodes
    bottom_node_list_cached::Vector{Int} = bottom_nodes # bottom nodes
    side_node_list_cached::Vector{Int} = side_nodes
    ID_cached::Matrix{Int} = ID

    @unpack nNodes = mdl.mesh_p
    nNodes_p_cached::Int = nNodes

    η_cached::Any = η

    camera_matrix_cached::Matrix{Float64} = camera_matrix
    camera_pose_cached::Matrix{Float64} = camera_pose
    SIDES_cached::Bool = SIDES
    
    C_Tu = transpose(C_uc_cached)           # transpose the constraint matrix

    if conditions.filepath != ""
        isnothing(conditions.filepath) || AssertionError("Please provide a filepath to write the data")
        set_file(conditions.filepath)
    end

    time = collect(Float64, range(start=t_steps_cached, stop=sim_time_cached, step=t_steps_cached))
    len_t = length(time)

    μu_btm = 0  
    μu_side = 0
    
    BorderPts2D, SurfacePts2D = extract_borders(NodeList_cached, camera_matrix_cached, camera_pose_cached, side_node_list_cached, nNodes_u_cached, SIDES_cached)
    pi, qi = fit_curve(border=BorderPts2D)
    
    dqdη = zeros(Float64, size(q_d_cached_top))
    dqdβ = zeros(Float64, size(q_d_cached_top))

    displacement = AbstractArray[zeros(Float64,size(NodeList_cached,1),size(NodeList_cached,2))] # store the displacement of the mesh in 3D
    surface_fields = AbstractArray[]
    surface_pts_3D = AbstractArray[vcat(NodeList_cached[:,top_node_list_cached]', NodeList_cached[:,bottom_node_list_cached]', NodeList_cached[:,side_node_list_cached]')'] # store the solution fields of the mesh in 3D
    gradList = AbstractArray[zeros(Float64, size(BorderPts2D,1),size(BorderPts2D,2),2)]                                                 # store the solution fields of the border nodes in 2D 
    pos3D = AbstractArray[NodeList_cached]                                                             # store the solution fields of the mesh in 3D
    pos2D = AbstractArray[SurfacePts2D]                                                                   # store the solution fields of the mesh in 2D
    borderPts2DList = AbstractArray[BorderPts2D]                                                               # store the solution fields of the surfaces in 2D
    splinep = AbstractArray[pi]                                                                            # store the x coordinates samples of the spline parameters of the border nodes
    splineq = AbstractArray[qi]                                                                            # store the y coordinates samples of the spline parameters of the border nodes
    output = Float64[] 
    writeborderList = [vcat(pi', qi')]

    iter::Int = 1
    pr = Progress(len_t; desc= "Simulating with prescribed $(control_cached) ...", showspeed=true)
    if control_cached == "force"
        A_bar = SparseMatrixCSC{Float64,Int}(I, nDof_u_cached*(nNodes_u_cached)^ndim_cached, nDof_u_cached*(nNodes_u_cached)^ndim_cached) # initialize the stiffness matrix
        B = SparseMatrixCSC{Float64,Int}(I, nDof_u_cached*(nNodes_u_cached)^ndim_cached, nDof_p_cached*(nNodes_p_cached)^ndim_cached) # initialize the stiffness matrix
        b = SparseMatrixCSC{Float64,Int}(I, nDof_u_cached*(nNodes_u_cached)^ndim_cached, nDof_u_cached*(nNodes_u_cached)^ndim_cached) # initialize the stiffness matrix
        q_d = zeros(Float64, nDof_u_cached*(nNodes_u_cached)^ndim_cached,1) # initialize the vector of the Dirichlet boundary conditions (for ndof = 1) / Dirichlet boundary conditions upper surface (for ndof > 1)
        A = similar(A_bar)

        A_free = C_Tu*A*C_uc_cached        # extract the free part of the stiffness matrix
        B_free = C_Tu*B             # extract the free part of the stiffness matrix

        dA_freedη = zeros(Float64, size(A_free,1), size(A_free,2))       # extract the free part of the stiffness matrix
        dA_freedβ = zeros(Float64, size(A_free,1), size(A_free,2))            # extract the free part of the stiffness matrix
        dB_free = zeros(Float64, size(B_free))
        zero = zeros(Float64, size(B_free,2),size(B_free,2))

        dAdη = similar(A_bar)
        dAdβ = similar(A_bar)
        dB = zeros(Float64, size(B))
        q = similar(q_d)
        dqfdη = similar(q)
        dqfdβ = similar(q)
        for t in time
            A_bar .= assemble_system_A(mdl)                   # assemble the stiffness matrix
            B .= assemble_system_B(mdl)                   # assemble the stiffness matrix
            b .= apply_boundary_conditions(mdl)           # apply the neumann boundary conditions
        
            q_d .= (μu_btm*q_d_cached_btm + μu_side*q_d_cached_brdr)      # apply the Dirichlet boundary conditions
            
            if viscosity_type_cached == "bulk_viscosity"
                A .= η_cached[iter]*A_bar + β_cached*b
            else
                A .= η_cached[1]*A_bar + β_cached*b
            end
        
            dAdη .= A_bar
            dAdβ .= b
        
            A_free .= C_Tu*A*C_uc_cached        # extract the free part of the stiffness matrix
            # mul!(A_free, C_Tu*A, C_uc_cached) # extract the free part of the stiffness matrix
            B_free .= C_Tu*B             # extract the free part of the stiffness matrix
            # mul!(B_free, C_Tu, B) # extract the free part of the stiffness matrix

            dA_freedη .= C_Tu*dAdη*C_uc_cached        # extract the free part of the stiffness matrix
            # mul!(dA_freedη, C_Tu*dAdη, C_uc_cached) # extract the free part of the stiffness matrix
            dA_freedβ .= C_Tu*dAdβ*C_uc_cached             # extract the free part of the stiffness matrix
            # mul!(dA_freedβ, C_Tu*dAdβ, C_uc_cached) # extract the free part of the stiffness matrix

            M = [A_free B_free C_Tu*A*q_d_cached_top; B_free' zero B'*q_d_cached_top; q_d_cached_top'*A*C_uc_cached q_d_cached_top'*B q_d_cached_top'A*q_d_cached_top] # assemble the system of equations
            dMdη = [dA_freedη dB_free C_Tu*dAdη*q_d_cached_top; dB_free' zero dB'*q_d_cached_top; q_d_cached_top'*dAdη*C_uc_cached q_d_cached_top'*dB q_d_cached_top'dAdη*q_d_cached_top]
            dMdβ = [dA_freedβ dB_free C_Tu*dAdβ*q_d_cached_top; dB_free' zero dB'*q_d_cached_top; q_d_cached_top'*dAdβ*C_uc_cached q_d_cached_top'*dB q_d_cached_top'dAdβ*q_d_cached_top]
            
            r = [-C_Tu*A*q_d; -B'*q_d; scene.cParam[iter].-q_d_cached_top'A*q_d]    # assemble the system of equations
            drdη = -[C_Tu*dAdη*q_d; zeros(Float64, size(B,2),size(q_d,2)); q_d_cached_top'dAdη*q_d] # solve the system of equations
            drdβ = -[C_Tu*dAdβ*q_d; zeros(Float64, size(B,2),size(q_d,2)); q_d_cached_top'dAdβ*q_d] # solve the system of equations

            # invM::Matrix = inv(Matrix(M)) # inverse of the system of equations
            lum = lu(M) # LU decomposition of the system of equations

            # sol = invM*r                 # solve the system of equations
            # dsoldη = invM*(drdη - dMdη*sol) # solve the system of equations
            # dsoldβ = invM*(drdβ - dMdβ*sol) # solve the system of equations

            sol = lum\r # solve the system of equations
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
        
            q .= q_d + C_uc_cached*q_f + μ_tp*q_d_cached_top;              # assemble the solution 
            dqdη .= dqdη + C_uc_cached*dqfdη + dμdη*q_d_cached_top;               # assemble the solution
            dqdβ .= dqdβ + C_uc_cached*dqfdβ + dμdβ*q_d_cached_top;               # assemble the solution

            p = p_f;
            dpdη = dpfdη;                  # assemble the solution
            dpdβ = dpfdβ;                  # assemble the solution

            motion = @views hcat(q[ID_cached[1,:]], q[ID_cached[2,:]], q[ID_cached[3,:]])'* t_steps_cached # extract the motion of the mesh grid
            dmdη_out = @views hcat(dqdη[ID_cached[1,:]], dqdη[ID_cached[2,:]], dqdη[ID_cached[3,:]])'
            dmdβ_out = @views hcat(dqdβ[ID_cached[1,:]], dqdβ[ID_cached[2,:]], dqdβ[ID_cached[3,:]])'

            dmdθ_out = @views cat(dmdη_out,dmdβ_out,dims=3) # concatenate the gradients in to a tensor
            
            NodeList_cached = NodeList_cached + motion # update the mesh grid
            mdl.mesh_u.NodeList = NodeList_cached # update the mesh grid
        
            BorderPts2D, dudθ, SurfacePts2D, ∇SurfacePts2D = extract_borders(NodeList_cached, camera_matrix_cached, camera_pose_cached, side_node_list_cached, GRAD=true, dqdθ=dmdθ_out, SIDES=SIDES_cached)
            pi, qi = fit_curve(border=BorderPts2D)

            mat_nan_inf_check(dudθ[:,:,1])
            mat_nan_inf_check(dudθ[:,:,2])

            # store the solutions in a list
            push!(output, μ_tp*t_steps_cached) # store displacement at the top surface
            push!(displacement, motion)
            push!(surface_fields, motion[:,side_node_list_cached])
            push!(surface_pts_3D, vcat(NodeList_cached[:,top_node_list_cached]', NodeList_cached[:,bottom_node_list_cached]', NodeList_cached[:,side_node_list_cached]')')
            push!(gradList,dudθ)
            push!(pos2D, SurfacePts2D)
            push!(pos3D, NodeList_cached)
            push!(borderPts2DList, BorderPts2D)
            push!(splinep, pi)
            push!(splineq, qi) 
            push!(writeborderList, vcat(pi', qi'))
        
            iter += 1
            next!(pr, showvalues = [(:iterations,iter),(:time,t)])
        end
    elseif control_cached == "velocity"
        for t in time
            A_bar = assemble_system_A(mdl)                   # assemble the stiffness matrix
            B = assemble_system_B(mdl)                   # assemble the stiffness matrix
            b = apply_boundary_conditions(mdl)           # apply the neumann boundary conditions
        
            q_d = (μu_btm*q_d_cached_btm + scene.cParam[iter]*q_d_cached_top + μu_side*q_d_cached_brdr)      # apply the Dirichlet boundary conditions
        
            if viscosity_type_cached == "bulk_viscosity"
                A = η_cached[iter]*A_bar + β_cached*b
            else
                A = η_cached[1]*A_bar + β_cached*b
            end

            dAdη = A_bar
            dAdβ = b
            dB = zeros(Float64, size(B))

            C_Tu = transpose(C_uc_cached)           # transpose the constraint matrix
        
            A_free = C_Tu*A*C_uc_cached        # extract the free part of the stiffness matrix
            B_free = C_Tu*B             # extract the free part of the stiffness matrix

            dA_freedη = C_Tu*dAdη*C_uc_cached        # extract the free part of the stiffness matrix
            dA_freedβ = C_Tu*dAdβ*C_uc_cached             # extract the free part of the stiffness matrix
            dB_free = zeros(Float64, size(B_free))
        
            K_free = [A_free B_free; B_free' zeros(Float64, size(B_free,2),size(B_free,2))]     # assemble the system of equations
            dKdη = [C_Tu*dAdη*C_uc_cached dB_free; dB_free' zeros(Float64, size(B,2),size(B,2))] # assemble the system of equations
            dKdβ = [C_Tu*dAdβ*C_uc_cached dB_free; dB_free' zeros(Float64, size(B,2),size(B,2))] # assemble the system of equations
            
            invK = inv(Matrix(K_free))
        
            r = [C_Tu*A*q_d; B'*q_d]    # assemble the system of equations
            drdη = [C_Tu*dAdη*q_d; zeros(Float64, size(B,2),size(q_d,2))] # solve the system of equations
            drdβ = [C_Tu*dAdβ*q_d; zeros(Float64, size(B,2),size(q_d,2))] # solve the system of equations

            sol = -invK*r                 # solve the system of equations
            dsoldη = -invK*(drdη + dKdη*sol) # solve the system of equations
            dsoldβ = -invK*(drdβ + dKdβ*sol) # solve the system of equations
        
            q_f = sol[1:size(A_free,1)]     # extract the free part of the solution
            dqfdη = dsoldη[1:size(A_free,1)] # extract the free part of the solution
            dqfdβ = dsoldβ[1:size(A_free,1)] # extract the free part of the solution

            p_f = sol[size(A_free,1)+1:end] # extract the free part of the solution
            dpfdη = dsoldη[size(A_free,1)+1:end] # extract the free part of the solution 
            dpfdβ = dsoldβ[size(A_free,1)+1:end] # extract the free part of the solution
        
            q = q_d + C_uc_cached*q_f;              # assemble the solution 
            dqdη = dqdη + C_uc_cached*dqfdη;              # assemble the solution
            dqdβ = dqdβ + C_uc_cached*dqfdβ;              # assemble the solution

            p = p_f;
            dpdη = dpfdη;                  # assemble the solution
            dpdβ = dpfdβ;                  # assemble the solution

            motion = hcat(q[ID_cached[1,:]], q[ID_cached[2,:]], q[ID_cached[3,:]])'*t_steps_cached # get the motion of the mesh
            dmdη_out = hcat(dqdη[ID_cached[1,:]], dqdη[ID_cached[2,:]], dqdη[ID_cached[3,:]])'
            dmdβ_out = hcat(dqdβ[ID_cached[1,:]], dqdβ[ID_cached[2,:]], dqdβ[ID_cached[3,:]])'

            dmdθ_out = cat(dmdη_out,dmdβ_out,dims=3) # concatenate the gradients in to a tensor
            
            NodeList_cached = NodeList_cached + motion # update the mesh grid
            mdl.mesh_u.NodeList = NodeList_cached # update the mesh grid
        
            BorderPts2D, dudθ, SurfacePts2D, ∇SurfacePts2D = extract_borders(NodeList_cached, camera_matrix_cached, camera_pose_cached, side_node_list_cached, GRAD=true, dqdθ=dmdθ_out, SIDES=SIDES_cached)
            pi, qi = fit_curve(border=BorderPts2D)

            # store the solutions in a list
            # push!(output, F_est[1])
            push!(displacement, motion)
            push!(surface_fields, motion[:,side_node_list_cached])
            push!(surface_pts_3D, vcat(NodeList_cached[:,top_node_list_cached]', NodeList_cached[:,bottom_node_list_cached]', NodeList_cached[:,side_node_list_cached]')')
            push!(gradList,dudθ)
            push!(pos2D, SurfacePts2D)
            push!(pos3D, NodeList_cached)
            push!(borderPts2DList, BorderPts2D)
            push!(splinep, pi)
            push!(splineq, qi) 
            push!(writeborderList, vcat(pi', qi'))
        
            iter += 1
            next!(pr, showvalues = [(:iterations,iter),(:time,t)])
        end
    else
            throw(ArgumentError("Control type not unknown"))
    end

    # write the data to a file
    if conditions.ANIMATE
        animate_fields(filepath = string(conditions.filepath,"/Results/images"), fields=surface_pts_3D , IEN=mdl.mesh_u.IEN, BorderNodes2D=borderPts2DList, fields2D=pos2D)
    end
    if conditions.WRITECONTOUR
        write_contour_data(string(conditions.filepath,"/Results"), writeborderList)
    end
    if conditions.WRITEVTK
        write_scene(string(conditions.filepath,"/Results"), pos3D, mdl.mesh_u.IEN, mdl.ne, mdl.ndim, displacement, ID=ID_cached, FunctionClass=mdl.mesh_u.FunctionClass)
    end
    return output, gradList, borderPts2DList, displacement, surface_pts_3D, pos2D
end


