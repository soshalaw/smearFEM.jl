using LinearAlgebra
using SparseArrays
using ProgressMeter
using Parameters
""" 
    assemble_system_A(mdl::Stokes)

Assembles the finite element system. # Returns the global stiffness matrix

# Arguments:
- `mdl::Stokes` : Material model

# Returns:
- `K::SparseMatrixCSC{Float64,Int64}{ndof,ndof}` : sparse stiffness matrix 
"""
function assemble_system_A(mdl::Stokes)::SparseMatrixCSC{Float64,Int64}
    # unpack the model parameters to local variables
    # this is done to avoid the need to pass the model object around
    @unpack ne, ndim, nDof_u = mdl
    @unpack NodeList, IEN, ID, FunctionClass, C_vol, C_top, C_btm, W = mdl.mesh_u    

    IEN_u_cached::Matrix{Int} = IEN
    ID_cached::Matrix{Int} = ID
    ne_cached::Int = ne
    ndim_cached::Int = ndim
    nDof_u_cached::Int = nDof_u
    FunctionClass_u_cached::String = FunctionClass

    @unpack NodeList, IEN, FunctionClass, C_vol, W = mdl.mesh_x

    NodeList_x_cached::Matrix{Float64} = NodeList
    IEN_x_cached::Matrix{Int} = IEN
    FunctionClass_x_cached::String = FunctionClass
    C_vol_x_cached = C_vol
    W_x_cached = W

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
    vol::Float64 = 0.0
    e_iter = 1:ne_cached # iterator for elements loop 
    gpiter = 1:length(wpoints) # iterator for integration loop
    # element loop
    for e::Int in e_iter
        coords::Matrix{Float64} = NodeList_x_cached[:, IEN_x_cached[:, e]]  # Get the coordinates of the nodes of the element
        # integration loop
        for gp::Int in gpiter
            if ndim_cached == 1
                N_u, ΔN = basis_function(x[gp], nothing, nothing, FunctionClass_u_cached) # add function type
            elseif ndim_cached == 2
                if string(FunctionClass_x_cached[1]) == "Q" && string(FunctionClass_u_cached[1]) == "Q"
                    # basis functions for fields
                    N_u, ΔN = basis_function(x[gp], y[gp], FunctionClass_u_cached) 
                    # basis functions for geometry
                    N_x, ΔN_x = N_u, ΔN
                elseif string(FunctionClass_x_cached[1]) == "S" && string(FunctionClass_u_cached[1]) == "S"
                    # basis functions for geometry
                    N_x, ΔN_x = basis_function(x[gp], y[gp], C_vol_x_cached[:,:,e], W_x_cached[IEN_x_cached[:,e]], FunctionClass_x_cached)
                    # basis functions for fields
                    N_u, ΔN_u = basis_function(x[gp], y[gp], C_vol_u_cached[:,:,e], W_u_cached[IEN_u_cached[:,e]], FunctionClass_u_cached)
                elseif string(FunctionClass_x_cached[1]) == "S" && string(FunctionClass_u_cached[1]) == "Q"
                    # basis functions for geometry
                    N_x, ΔN_x = basis_function(x[gp], y[gp], C_vol_x_cached[:,:,e], W_x_cached[IEN_x_cached[:,e]], FunctionClass_x_cached) 
                    # basis functions for fields
                    N_u, ΔN_u = basis_function(x[gp], y[gp], FunctionClass_u_cached) 
                end
            elseif ndim_cached == 3
                if string(FunctionClass_x_cached[1]) == "Q" && string(FunctionClass_u_cached[1]) == "Q"
                    # basis functions for fields
                    N_u, ΔN_u = basis_function(x[gp], y[gp], z[gp], FunctionClass_u_cached) 
                    # basis functions for geometry
                    N_x, ΔN_x = N_u, ΔN_u
                elseif string(FunctionClass_x_cached[1]) == "S" && string(FunctionClass_u_cached[1]) == "S"
                    # basis functions for geometry
                    N_x, ΔN_x = basis_function(x[gp], y[gp], z[gp], C_vol_x_cached[:,:,e], W_x_cached[IEN_x_cached[:,e]], FunctionClass_x_cached)
                    # basis functions for fields
                    N_u, ΔN_u = basis_function(x[gp], y[gp], z[gp], C_vol_u_cached[:,:,e], W_u_cached[IEN_u_cached[:,e]], FunctionClass_u_cached)
                elseif string(FunctionClass_x_cached[1]) == "S" && string(FunctionClass_u_cached[1]) == "Q"
                    # basis functions for geometry
                    N_x, ΔN_x = basis_function(x[gp], y[gp], z[gp], C_vol_x_cached[:,:,e], W_x_cached[IEN_x_cached[:,e]], FunctionClass_x_cached) 
                    # basis functions for fields
                    N_u, ΔN_u = basis_function(x[gp], y[gp], z[gp], FunctionClass_u_cached) 
                end
            end

            dNdX_u = zeros(Float64, size(N_u, 1), ndim_cached) # Gradient of basis functions
            if nDof_u_cached == 2
                B = zeros(Float64, ndim*length(N_u), 3)
            elseif nDof_u_cached == 3
                B = zeros(Float64, ndim*length(N_u), 6)
            end

            szB::Int = size(B, 1)  # Number of basis functions
            Ke = zeros(Float64, szB, szB)  # Element stiffness matrix
            Ke_row::Int, Ke_col::Int = szB, szB
            Ke_len::Int = length(Ke)  

            mul!(Jac, coords, ΔN_x)  # Jacobian matrix [dx/dxi dx/deta; dy/dxi dy/deta]
            w::Float64 = wpoints[gp] * abs(det(Jac))
            vol = vol + w
            invJ .= inv(Jac)  # Inverse of the Jacobian matrix
            mul!(dNdX_u, ΔN_u, invJ)

            if nDof_u_cached == 1
                szN::Int = size(N_u, 1)  # Number of basis functions
                # Loop between basis functions of the element
                for i::Int in 1:szN
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
                
                for iNode::Int in iNodes
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

    @unpack NodeList, IEN, ID, FunctionClass, C_vol, W = mdl.mesh_u

    IEN_u_cached::Matrix{Int} = IEN
    NodeList_cached::Matrix{Float64} = NodeList
    ID_cached::Matrix{Int} = ID
    ne_cached::Int = ne
    ndim_cached::Int = ndim
    nDof_u_cached::Int = nDof_u
    FunctionClass_u_cached::String = FunctionClass
    C_vol_u_cached = C_vol
    W_u_cached = W

    @unpack NodeList, IEN, FunctionClass, C_vol, C_top, C_btm, W = mdl.mesh_x

    NodeList_x_cached::Matrix{Float64} = NodeList
    IEN_x_cached::Matrix{Int} = IEN
    FunctionClass_x_cached::String = FunctionClass
    C_vol_x_cached = C_vol
    W_x_cached = W

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

    e_iter = 1:ne_cached    # iterator for elements loop
    gpiter = 1:length(wpoints)  # iterator for integration loop
    # element loop
    vol = 0.0
    for e::Int in e_iter
        # println("IEN_x_cached[:, $e]: ", IEN_x_cached[:, e])
        coords::Matrix{Float64} = NodeList_x_cached[:, IEN_x_cached[:, e]]  # Get the coordinates of the nodes of the element
        # integration loop
        for gp::Int in gpiter
            if ndim_cached == 1
                N_u, ΔN_u = basis_function(x[gp], nothing, nothing, FunctionClass_u_cached)
                N_p, ΔN_p = basis_function(x[gp], nothing, nothing, FunctionClass_p_cached)
            elseif ndim_cached == 2
                N_u, ΔN_u = basis_function(x[gp], y[gp], nothing, FunctionClass_u_cached)
                N_p, ΔN_p = basis_function(x[gp], y[gp], nothing, FunctionClass_p_cached)
            elseif ndim_cached == 3
                if string(FunctionClass_x_cached[1]) == "Q" && string(FunctionClass_u_cached[1]) == "Q" && string(FunctionClass_p_cached[1]) == "Q"
                    # basis functions for fields
                    N_u, ΔN_u = basis_function(x[gp], y[gp], z[gp], FunctionClass_u_cached)
                    N_p, ΔN_p = basis_function(x[gp], y[gp], z[gp], FunctionClass_p_cached)
                    # basis functions for geometry
                    N_x, ΔN_x = N_u, ΔN_u 
                elseif string(FunctionClass_x_cached[1]) == "S" && string(FunctionClass_u_cached[1]) == "S" && string(FunctionClass_p_cached[1]) == "S"
                    # basis functions for geometry
                    N_x, ΔN_x = basis_function(x[gp], y[gp], z[gp], C_vol_x_cached[:,:,e], W_x_cached[IEN_x_cached[:,e]], FunctionClass_x_cached) 
                    # basis functions for fields
                    N_u, ΔN_u = basis_function(x[gp], y[gp], z[gp], C_vol_u_cached[:,:,e], W_u_cached[IEN_u_cached[:,e]], FunctionClass_u_cached) 
                    N_p, ΔN_p = basis_function(x[gp], y[gp], z[gp], C_vol_p_cached[:,:,e], W_p_cached[IEN_p_cached[:,e]], FunctionClass_p_cached) 
                elseif string(FunctionClass_x_cached[1]) == "S" && string(FunctionClass_u_cached[1]) == "Q" && string(FunctionClass_p_cached[1]) == "Q"
                    # basis functions for geometry
                    N_x, ΔN_x = basis_function(x[gp], y[gp], z[gp], C_vol_x_cached[:,:,e], W_x_cached[IEN_x_cached[:,e]], FunctionClass_x_cached) 
                    # basis functions for fields
                    N_u, ΔN_u = basis_function(x[gp], y[gp], z[gp], FunctionClass_u_cached)
                    N_p, ΔN_p = basis_function(x[gp], y[gp], z[gp], FunctionClass_p_cached)
                end
            end
            dNdX_u = zeros(Float64, size(ΔN_u, 1), ndim_cached) # Gradient of basis functions
            Be = zeros(Float64, 1, 3*size(dNdX_u, 1)) # Gradient of basis functions
            colB::Int = size(Be, 2)  # Number of basis functions
            rowNp::Int = size(N_p, 1)  # Number of basis functions
            Ke = zeros(Float64, rowNp, colB)  # Element stiffness matrix
            Ke_row::Int, Ke_col::Int = rowNp, colB

            mul!(Jac, coords, ΔN_x)  # Jacobian matrix [dx/dxi dx/deta; dy/dxi dy/deta]
            w::Float64 = wpoints[gp] * abs(det(Jac))
            vol = vol + w
            invJ .= inv(Jac)  # Inverse of the Jacobian matrix
            mul!(dNdX_u, ΔN_u, invJ)

            if nDof_u_cached == 1
                szN = size(N_u,1) # number of basis functions
                # loop between basis functions of the element
                for i::Int in 1:szN
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
    # println("volume B:", vol)
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

    e_iter = 1:ne_cached    # iterator for elements loop
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
- `mdl::Stokes` : Material model

# Returns:
- `K: {[ndof,ndof] SparseMatrixCSC{Float64,Int64}}` : sparse stiffness matrix with the boundary conditions applied
- `F: {[ndof] Vector{Float64}}` : force vector
"""
function apply_boundary_conditions(mdl::Stokes)::SparseMatrixCSC{Float64,Int64}

    @unpack ne, ndim, nDof_u = mdl
    @unpack NodeList, IEN_top, IEN_bottom, ID, FunctionClass, C_top, C_btm, W = mdl.mesh_u

    NodeList_u_cached::Matrix{Float64} = NodeList
    IEN_u_top_cached::Matrix{Int} = IEN_top
    IEN_u_btm_cached::Matrix{Int} = IEN_bottom
    ID_cached::Matrix{Int} = ID
    ne_cached::Int = ne
    ndim_cached::Int = ndim
    nDof_u_cached::Int = nDof_u
    FunctionClass_u_cached::String = FunctionClass
    C_top_u_cached = C_top
    C_btm_u_cached = C_btm
    W_u_cached = W

    @unpack NodeList, IEN_top, IEN_bottom, ID, FunctionClass, C_top, C_btm, W = mdl.mesh_x

    NodeList_x_cached::Matrix{Float64} = NodeList
    IEN_x_top_cached::Matrix{Int} = IEN_top
    IEN_x_btm_cached::Matrix{Int} = IEN_bottom
    FunctionClass_x_cached::String = FunctionClass
    C_top_x_cached = C_top
    C_btm_x_cached = C_btm
    W_x_cached = W

    ID_u_rows::Int = size(ID_cached,1)
    IEN_btm_rows::Int = size(IEN_u_btm_cached,1)

    ne_surface = size(IEN_u_top_cached, 2) # number of elements on the surface
    E = zeros(  Int64, ne_surface*(ID_u_rows*IEN_btm_rows)^2*2) # *2 because we have two surfaces
    J = zeros(  Int64, ne_surface*(ID_u_rows*IEN_btm_rows)^2*2) # *2 because we have two surfaces
    V = zeros(Float64, ne_surface*(ID_u_rows*IEN_btm_rows)^2*2) # *2 because we have two surfaces

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
    
    e_iter = 1:size(IEN_u_top_cached, 2)   # iterator for elements loop
    gpiter = 1:length(wpoints) # iterator for integration loop
    # element loop
    A_btm = 0
    A_top = 0
    for e::Int in e_iter
        coords_top::Matrix{Float64} = NodeList_x_cached[:,IEN_x_top_cached[:,e]] # get the coordinates of the nodes of the element
        coords_btm::Matrix{Float64} = NodeList_x_cached[:,IEN_x_btm_cached[:,e]] # get the coordinates of the nodes of the element
        # integration loop
        for gp::Int in gpiter
            if ndim_cached == 2
                N, ΔN = basis_function(x[gp], nothing, nothing, FunctionClass_u_cached)
            elseif ndim_cached == 3
                if string(FunctionClass_x_cached[1]) == "Q" && string(FunctionClass_u_cached[1]) == "Q"
                    # fields
                    N_u_top, ΔN_u_top = basis_function(x[gp], y[gp], nothing, FunctionClass_u_cached) 
                    N_u_btm, ΔN_u_btm = N_u_top, ΔN_u_top
                    # geometry
                    N_x_top, ΔN_x_top = N_u_top, ΔN_u_top
                    N_x_btm, ΔN_x_btm = N_u_top, ΔN_u_top
                elseif string(FunctionClass_x_cached[1]) == "S" && string(FunctionClass_u_cached[1]) == "S"
                    # geometry
                    N_x_top, ΔN_x_top = basis_function(x[gp], y[gp], C_top_x_cached[:,:,e], W_x_cached[IEN_x_top_cached[:,e]], FunctionClass_x_cached)
                    N_x_btm, ΔN_x_btm = basis_function(x[gp], y[gp], C_btm_x_cached[:,:,e], W_x_cached[IEN_x_btm_cached[:,e]], FunctionClass_x_cached)
                    # fields
                    N_u_top, ΔN_u_top = basis_function(x[gp], y[gp], C_top_u_cached[:,:,e], W_u_cached[IEN_u_top_cached[:,e]], FunctionClass_u_cached)
                    N_u_btm, ΔN_u_btm = basis_function(x[gp], y[gp], C_btm_u_cached[:,:,e], W_u_cached[IEN_u_btm_cached[:,e]], FunctionClass_u_cached)
                elseif string(FunctionClass_x_cached[1]) == "S" && string(FunctionClass_u_cached[1]) == "Q"
                    # geometry
                    N_x_top, ΔN_x_top = basis_function(x[gp], y[gp], C_top_x_cached[:,:,e], W_x_cached[IEN_x_top_cached[:,e]], FunctionClass_x_cached) 
                    N_x_btm, ΔN_x_btm = basis_function(x[gp], y[gp], C_btm_x_cached[:,:,e], W_x_cached[IEN_x_btm_cached[:,e]], FunctionClass_x_cached)
                    # fields
                    N_u_top, ΔN_u_top = basis_function(x[gp], y[gp], nothing, FunctionClass_u_cached) 
                    N_u_btm, ΔN_u_btm = N_u_top, ΔN_u_top
                end
            end

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

            for iNode::Int in iNodes
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
            # TODO include in tha assembly function
        end
    end
    # println("Maximum row :",maximum(E))
    # println("Maximum column :",maximum(J))
    # println("A_top = ", A_top)
    # println("A_btm = ", A_btm)
    K = sparse(E,J,V)
    return  K
end

function apply_boundary_conditions_dense(mdl::Stokes)::Matrix{Float64}
    @unpack ne, ndim, nDof_u = mdl
    @unpack NodeList, IEN_top, IEN_bottom, ID, FunctionClass, C_top, C_btm, W = mdl.mesh_u

    NodeList_u_cached::Matrix{Float64} = NodeList
    IEN_u_top_cached::Matrix{Int} = IEN_top
    IEN_u_btm_cached::Matrix{Int} = IEN_bottom
    ID_cached::Matrix{Int} = ID
    ne_cached::Int = ne
    ndim_cached::Int = ndim
    nDof_u_cached::Int = nDof_u
    FunctionClass_u_cached::String = FunctionClass
    C_top_u_cached = C_top
    C_btm_u_cached = C_btm
    W_u_cached = W

    @unpack NodeList, IEN_top, IEN_bottom, ID, FunctionClass, C_top, C_btm, W = mdl.mesh_x

    NodeList_x_cached::Matrix{Float64} = NodeList
    IEN_x_top_cached::Matrix{Int} = IEN_top
    IEN_x_btm_cached::Matrix{Int} = IEN_bottom
    FunctionClass_x_cached::String = FunctionClass
    C_top_x_cached = C_top
    C_btm_x_cached = C_btm
    W_x_cached = W

    ID_u_rows::Int = size(ID_cached,1)
    IEN_btm_rows::Int = size(IEN_u_btm_cached,1)
    ne_surface = size(IEN_u_top_cached, 2) # number of elements on the surface
    nNodes_cached = size(NodeList_x_cached, 2)

    sz = nDof_u_cached*nNodes_cached
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
    
    e_iter = 1:size(IEN_u_top_cached, 2)   # iterator for elements loop
    gpiter = 1:length(wpoints) # iterator for integration loop
    # element loop
    A_btm = 0
    A_top = 0
    for e::Int in e_iter
        coords_top::Matrix{Float64} = NodeList_x_cached[:,IEN_x_top_cached[:,e]] # get the coordinates of the nodes of the element
        coords_btm::Matrix{Float64} = NodeList_x_cached[:,IEN_x_btm_cached[:,e]] # get the coordinates of the nodes of the element
        # integration loop
        for gp::Int in gpiter
            if ndim_cached == 2
                N, ΔN = basis_function(x[gp], nothing, nothing, FunctionClass_u_cached)
            elseif ndim_cached == 3
                if string(FunctionClass_x_cached[1]) == "Q" && string(FunctionClass_u_cached[1]) == "Q"
                    # fields
                    N_u_top, ΔN_u_top = basis_function(x[gp], y[gp], nothing, FunctionClass_u_cached) 
                    N_u_btm, ΔN_u_btm = N_u_top, ΔN_u_top
                    # geometry
                    N_x_top, ΔN_x_top = N_u_top, ΔN_u_top
                    N_x_btm, ΔN_x_btm = N_u_top, ΔN_u_top
                elseif string(FunctionClass_x_cached[1]) == "S" && string(FunctionClass_u_cached[1]) == "S"
                    # geometry
                    N_x_top, ΔN_x_top = basis_function(x[gp], y[gp], C_top_x_cached[:,:,e], W_x_cached[IEN_x_top_cached[:,e]], FunctionClass_x_cached)
                    N_x_btm, ΔN_x_btm = basis_function(x[gp], y[gp], C_btm_x_cached[:,:,e], W_x_cached[IEN_x_btm_cached[:,e]], FunctionClass_x_cached)
                    # fields
                    N_u_top, ΔN_u_top = basis_function(x[gp], y[gp], C_top_u_cached[:,:,e], W_u_cached[IEN_u_top_cached[:,e]], FunctionClass_u_cached)
                    N_u_btm, ΔN_u_btm = basis_function(x[gp], y[gp], C_btm_u_cached[:,:,e], W_u_cached[IEN_u_btm_cached[:,e]], FunctionClass_u_cached)
                elseif string(FunctionClass_x_cached[1]) == "S" && string(FunctionClass_u_cached[1]) == "Q"
                    # geometry
                    N_x_top, ΔN_x_top = basis_function(x[gp], y[gp], C_top_x_cached[:,:,e], W_x_cached[IEN_x_top_cached[:,e]], FunctionClass_x_cached) 
                    N_x_btm, ΔN_x_btm = basis_function(x[gp], y[gp], C_btm_x_cached[:,:,e], W_x_cached[IEN_x_btm_cached[:,e]], FunctionClass_x_cached)
                    # fields
                    N_u_top, ΔN_u_top = basis_function(x[gp], y[gp], nothing, FunctionClass_u_cached) 
                    N_u_btm, ΔN_u_btm = N_u_top, ΔN_u_top
                end
            end

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
            # TODO include in tha assembly function
        end
    end
    return  K
end

"""
set_boundary_cond(mdl::Stokes; DENSE::Bool=false)
Set the Dirichlet boundary conditions for the problem

# Arguments:
- `mdl::Stokes` : Material model
- `DENSE::Bool` : Flag to use dense matrices or sparse matrices

# Returns:
- `q_upper::Vector{Float64}` : vector of the Dirichlet boundary conditions (for ndof = 1) / Dirichlet boundary conditions upper surface (for ndof > 1)
- `q_lower::Vector{Float64}` : vector of the Neumann boundary conditions (for ndof = 1) / Dirichlet boundary conditions lower surface (for ndof > 1)
- `C_uc_cached::SparseMatrixCSC{Float64,Int64}` : onstraint matrix
"""
function set_boundary_cond(mdl::Stokes; DENSE::Bool=false)

    @unpack ndim, nDof_u = mdl
    @unpack NodeList, nNodes, ID = mdl.mesh_u

    NodeList_cached::Matrix{Float64} = NodeList
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
            iter = 1:size(NodeList_cached,2)
            for n::Int in iter
                coord = NodeList_cached[:,n] # get the coordinates of the node
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
            iter = 1:size(NodeList_cached,2)
            for n::Int in iter
                coord = NodeList_cached[:,n] # get the coordinates of the node
                if coord[2] == Dbound1 # bottom boundary
                    q_upper[n] = 0
                    push!(rCol, n)
                elseif coord[2] == Dbound2 # top boundary
                    q_upper[n] = -1
                    push!(rCol, n)
                end
            end
        end

        C_uc = C[:,setdiff(1:size(C,2),rCol)]

    else
        z0Bound = 0
        z1Bound = mdl.mesh_u.h # height of the cylinder

        iter = 1:size(NodeList_cached,2)
        for nNode::Int in iter
            coord = NodeList_cached[:,nNode]    # get the coordinates of the node
            if coord[3] == z1Bound   # top boundary
                q_upper[ID_cached[3,nNode]] = 1     # constraint the z displacement to be -d
                push!(rCol,ID_cached[3,nNode])
            elseif coord[3] == z0Bound   # bottom boundary
                q_lower[ID_cached[3,nNode]] = 1     # constraint the z displacement to be -d
                push!(rCol,ID_cached[3,nNode])
            end
        end

        C_uc = C[:,setdiff(1:size(C,2),rCol)]

    end
    return q_upper, q_side, q_lower, C_uc
end

"""
    get_η_power_law(t::Float64, F::Float64, R_0::Float64, H_0::Float64, η_0::Float64, n::Float64, K::Float64)
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
function get_η_power_law(t::T, F::U, R_0::V, H_0::W, η_0::X, n::Y, K::Z) where {T<:Number,U<:Number,V<:Number,W<:Number,X<:Number,Y<:Number,Z<:Number}

    H(t) = H_0*(1+8*H_0^2*F*t/(3*π*η_0*R_0^4))^(-1/4) # height with time
    
    R(t) = R_0*(1+8*H_0^2*F*t/(3*π*η_0*R_0^4))^(1/8) # radius with time

    H_dot(t) = 8/3*(-2*F*H(t)^3/(8*π*η_0*R(t)^4)) # rate of change of height

    γ_dot(t) = H_dot(t)/H(t) # shear rate

    η(t) = K*(abs(γ_dot(t)))^(n-1) # shear viscosity

    return η(t)
end

function get_η_carraeu(H::T, H_dot::T, F::U, R::V) where {T<:Number,U<:Number,V<:Number}

    η = 8/3*F*H^3/(π*R^4*H_dot)

    return η
end

"""
function def_problem(r::T, h::U, ne::Int64, η_0::V, ndim::Int64, FunctionClass_u::String, nDof_u::Int64, FunctionClass_p::String, 
                    nDof_p::Int64, β::Float64, cParam::Vector{Float64}, control::String, viscosity_type::String, sim_time::W, t_steps::X) where {T<:Number,U<:Number,V<:Number,W<:Number,X<:Number}

Define the conditions and the parameters of the squueze flow problem considered.
"""
function def_problem(r::T, h::U, ne::Int64, η_0::V, ndim::Int64, FunctionClass_u::String, nDof_u::Int64, FunctionClass_p::String, 
                    nDof_p::Int64, FunctionClass_x::String, β::Y, cParam::Vector{Float64}, control::String, viscosity_type::String, sim_time::W, t_steps::X; viscosity_model::String="power_law") where {T<:Number,U<:Number,V<:Number,W<:Number,X<:Number,Y<:Number}
    n::Float64 = 0.9
    K::Float64 = 100.0
    time = collect(Float64, range(start=t_steps, stop=sim_time, step=t_steps))
    len_t::Int = length(time)
    @info "Simulation time: $sim_time, Time step: $t_steps, Number of time steps: $(round(Int, sim_time/t_steps))"
    @info "Length of time array: $(len_t)"
    
    if length(cParam) < len_t   
        @error "Length of the Force vector ($(length(cParam))) is less than length of time array ($(length(time)))"
    end

    η = [η_0]
    if viscosity_type == "bulk_viscosity"
        if viscosity_model == "power_law"
            @info "Using power law viscosity model"
            η = get_η_power_law.(time, -cParam[1:len_t], r, h, η_0, n, K)
        elseif viscosity_model == "carreau"
            @info "Using Carreau viscosity model reading data from file"
            η = [η_0]
        end
    end 

    # define the model
    stokes = set_model(r, h, ne, η, ndim, FunctionClass_u, nDof_u, FunctionClass_p, nDof_p, FunctionClass_x) # define the model    
    # set the boundary conditions
    q_tp, q_side, q_btm, C_uc = set_boundary_cond(stokes) 
    # define with problem scenario
    c = [q_tp, q_side, q_btm]

    squeeze = SqueezeFlow(stokes, [β], [q_tp, q_side, q_btm], C_uc, control, sim_time, t_steps, viscosity_type, cParam)
return stokes, squeeze
end

"""
set_model(ne, NodeList, IEN, ndim, FunctionClass=FunctionClass_p_cached, nDof=1, ID=nothing, Young=1, ν=0.3)

Sets the model for the finite element method.
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
function set_model(r::R, h::H, ne::Int64, η::Vector{Float64}, ndim::Int64, FunctionClass_u::String, nDof_u::Int64, FunctionClass_p::String, 
                nDof_p::Int64, FunctionClass_x::String; GMESH_MESH::Bool=true)::Stokes where {R<:Number,H<:Number}

    filePath = joinpath("home", "soshala", "SMEAR-PhD", "smear-modules", "smearFEM.jl", "cylindergen")

    if FunctionClass_x == "S2"
        mesh_x = meshgrid_cylinder(r, h, ne, FunctionClass=FunctionClass_x, filePath=filePath)  # generate the mesh grid for geometry
        mesh_u = meshgrid_cylinder(r, h, ne, FunctionClass=FunctionClass_u)  # generate the mesh grid
        mesh_p = meshgrid_cylinder(r, h, ne, FunctionClass=FunctionClass_p)  # generate the mesh grid
    else
        if GMESH_MESH == true
            @info "Using Gmsh to generate the mesh"
            filepath_mesh = joinpath("/home", "soshala", "SMEAR-PhD", "smear-modules", "smear-meshes")
            mesh_u = get_gmsh_cylinder(joinpath(filepath_mesh,"cylinder_x_$ne.msh"), nDof_u, r, h, FunctionClass_u)  # generate the mesh grid for velocity field
            mesh_p = get_gmsh_cylinder(joinpath(filepath_mesh,"cylinder_p_$ne.msh"), nDof_p, r, h, FunctionClass_p)  # generate the mesh grid for pressure field
            mesh_x = get_gmsh_cylinder(joinpath(filepath_mesh,"cylinder_x_$ne.msh"), 1, r, h, FunctionClass_x)      # generate the mesh grid for geometry
        else
            @info "Using Julia to generate the mesh"
            mesh_x = meshgrid_cylinder(r, h, ne, FunctionClass=FunctionClass_x)  # generate the mesh grid for geometry
            mesh_u = meshgrid_cylinder(r, h, ne, FunctionClass=FunctionClass_u)  # generate the mesh grid
            mesh_p = meshgrid_cylinder(r, h, ne, FunctionClass=FunctionClass_p)  # generate the mesh grid
        end
    end
    
    mdl = Stokes(ndim=ndim, mesh_x=mesh_x, mesh_u=mesh_u, nDof_u=nDof_u, mesh_p=mesh_p, nDof_p=nDof_p, η=η)

    return mdl
end

"""
simulate(mdl::Stokes, scene::SqueezeFlow, conditions::Conditions)

Simulate the Stokes problem for a given mesh over a given time period. 

# Arguments:
- `mdl::Stokes`: Material model
- `scene::SqueezeFlow`: Parameters for the squeeze flow problem
- `conditions::Conditions` : External environmental conditions

# Returns:
- `output::Vector{Float64}` : output of the simulation
- `gradList`::Vector{Float64}` : gradient of the solution
- `borderPts2DList::Vector{Matrix{Float64}}` : border points in 2D 
- `displacement::Vector{Matrix{Float64}}` : displacement of the nodal points at each timestep
- `surface_pts_3D::Vector{Matrix{Float64}}` : 3D ccordinates of the surface nodes at each timestep
- `pos2D::Vector{Matrix{Float64}}` : 2D projection of the surface nodes at each timestep

"""
function simulate(mdl::Stokes, scene::SqueezeFlow, conditions::Conditions)

    reset_model!(mdl)
    
    @unpack FunctionClass, IEN, IEN_cp, ID, NodeList, C_vol, W = mdl.mesh_x
    FunctionClass_x_cached::String = FunctionClass
    NodeList_x_cached::Matrix{Float64} = NodeList
    IEN_x_cached::Matrix{Int} = IEN
    IEN_x_cp_cached::Matrix{Int} = IEN_cp
    ID_x_cached::Matrix{Int} = ID
    C_vol_x_cached = C_vol
    W_x_cached = W

    @unpack IEN, ID, FunctionClass, top_nodes, bottom_nodes, side_nodes, nNodes = mdl.mesh_u
    FunctionClass_u_cached::String = FunctionClass
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

    @unpack camera_matrix, obj_pose, SIDES = conditions    
    camera_matrix_cached::Matrix{Float64} = camera_matrix
    obj_pose_cached::Matrix{Float64} = obj_pose
    SIDES_cached::Bool = SIDES
    
    NodeList_cached::Matrix{Float64} = NodeList_u_cached
    ID_cached::Matrix{Int} = ID_u_cached
    T = Matrix{Float64}(I, size(NodeList_x_cached,2), size(NodeList_u_cached,2)) # projection matrix from geometry to field mesh

    time = collect(Float64, range(start=t_steps_cached, stop=sim_time_cached, step=t_steps_cached))
    len_t = length(time)

    if FunctionClass_x_cached == "S2" && FunctionClass_u_cached != "S2"
        # ID_cached = ID_x_cached
        T = get_nurbs_2_lagrange_proj(IEN_x_cached, IEN_u_cached, C_vol_x_cached, NodeList_x_cached, W_x_cached)
        NodeList_cached = NodeList_x_cached
    end

    T_ = T'*inv(T*T')
    C_Tu = transpose(C_uc_cached) # transpose the constraint matrix

    if conditions.filepath != ""
        isnothing(conditions.filepath) || AssertionError("Please provide a filepath to write the data")
        set_file(conditions.filepath)
    end

    μu_btm = 0  
    μu_side = 0
    
    NodeList_proj = NodeList_cached*T # project the motion on the geometry mesh grid
            
    BorderPts2D, SurfacePts2D = extract_borders(NodeList_proj, camera_matrix_cached, obj_pose_cached, 13, BorderNodesList=side_node_list_cached)
    pi, qi = fit_curve(border=BorderPts2D)
    
    dqdη = zeros(Float64, size(q_d_cached_top))
    dqdβ = zeros(Float64, size(q_d_cached_top))

    velocity = AbstractArray[zeros(Float64,size(NodeList_proj,1),size(NodeList_proj,2))] # store the velocity of the mesh in 3D
    pressure = AbstractArray[zeros(Float64,size(NodeList_proj,1),1)] # store the pressure of the mesh in 3D
    displacement = AbstractArray[zeros(Float64,size(NodeList_proj,1),size(NodeList_proj,2))] # store the displacement of the mesh in 3D
    surface_fields = AbstractArray[]
    surface_pts_3D = AbstractArray[vcat(NodeList_proj[:,top_node_list_cached]', 
                                        NodeList_proj[:,bottom_node_list_cached]', 
                                        NodeList_proj[:,side_node_list_cached]')'] # store the solution fields of the mesh in 3D
    gradList = AbstractArray[zeros(Float64, size(BorderPts2D,1),size(BorderPts2D,2),2)] # store the solution fields of the border nodes in 2D 
    pos3D = AbstractArray[NodeList_proj]       # store the solution fields of the mesh in 3D
    pos3D_cp = AbstractArray[NodeList_cached]  
    pos2D = AbstractArray[SurfacePts2D]          # store the solution fields of the mesh in 2D
    borderPts2DList = AbstractArray[BorderPts2D] # store the solution fields of the surfaces in 2D
    splinep = AbstractArray[BorderPts2D[1,:]]  # store the x coordinates samples of the spline parameters of the border nodes
    splineq = AbstractArray[BorderPts2D[2,:]]  # store the y coordinates samples of the spline parameters of the border nodes
    output = Float64[] 
    writeborderList = [vcat(pi', qi')]

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
    zero = spzeros(size(B_free,2),size(B_free,2))

    dAdη = similar(_A_bar)
    dAdβ = similar(_A_bar)
    dB = spzeros(size(B))

    q = similar(q_d)
    dqfdη = similar(q)
    dqfdβ = similar(q)

    iter::Int = 1
    pr = progress_guard(len_t; desc= "Simulating with prescribed $(control_cached) ...", showspeed=true)
    if control_cached == "force"

        M = spzeros((size(A_free,1)+size(B_free,2)+1),(size(A_free,2)+size(B_free,2)+1))
        dMdη = spzeros(size(M))
        dMdβ = spzeros(size(M))

        for t in time
            _A_bar .= assemble_system_A(mdl)    # assemble the stiffness matrix
            B .= assemble_system_B(mdl)         # assemble the stiffness matrix
            b .= apply_boundary_conditions(mdl) # apply the neumann boundary conditions
            q_d .= (μu_btm*q_d_cached_btm + μu_side*q_d_cached_brdr) # apply the Dirichlet boundary conditions
   
            if viscosity_type_cached == "bulk_viscosity"
                if length(β_cached) == 1
                    A .= η_cached[iter]*_A_bar + β_cached[1]*b
                    A_bar .= A #η_cached[iter]*_A_bar
                else
                    A .= η_cached[iter]*_A_bar + β_cached[iter]*b
                    A_bar .= A # η_cached[iter]*_A_bar
                end
            else
                A .= η_cached[1]*_A_bar + β_cached[1]*b
                A_bar .= A #η_cached[1]*_A_bar
            end

            dAdη .= _A_bar
            dAdβ .= b
 
            A_free .= C_Tu*A*C_uc_cached # extract the free part of the stiffness matrix
            B_free .= C_Tu*B             # extract the free part of the stiffness matrix

            dA_freedη .= C_Tu*dAdη*C_uc_cached # extract the free part of the stiffness matrix
            dA_freedβ .= C_Tu*dAdβ*C_uc_cached # extract the free part of the stiffness matrix

            M[1:size(A_free,1),1:size(A_free,2)] = A_free
            M[(size(A_free,1)+1):(size(A_free,1)+size(B_free,2)),1:size(A_free,2)] = B_free'
            M[end,1:size(A_free,2)] = q_d_cached_top'*A_bar*C_uc_cached

            M[1:size(A_free,1),(size(A_free,2)+1):(size(A_free,2)+size(B_free,2))] = B_free
            M[(size(A_free,1)+1):(size(A_free,1)+size(B_free,2)),(size(A_free,2)+1):(size(A_free,2)+size(B_free,2))] = zero
            M[end,(size(A_free,2)+1):(size(A_free,2)+size(B_free,2))] = q_d_cached_top'*B

            M[1:size(A_free,1),end] = C_Tu*A*q_d_cached_top
            M[(size(A_free,1)+1):(size(A_free,1)+size(B_free,2)),end] = B'*q_d_cached_top
            M[end,end] = (q_d_cached_top'*A_bar*q_d_cached_top)[end]
            
            @debug begin
                println("Time: ", t)
                println("Iteration: ", iter)
                println("η: ", η_cached[1])
                println("β: ", β_cached[1])
                println("Norm of _A_bar: ", maximum(_A_bar))
                println("Norm of A_bar: ", maximum(A_bar))
                println("Norm of b: ", maximum(b))
                println("Norm of A: ", maximum(A))
                println("Norm of B: ", maximum(B))
                println("Norm of M: ", maximum(M))
                push!(dac_list, norm(q_d_cached_top'*A_bar*C_uc_cached))
                push!(bd_list, norm(q_d_cached_top'*B))
                push!(dad_list, norm((q_d_cached_top'*A_bar*q_d_cached_top)[end]))
                push!(A_list, norm(A_bar))
            end

            dMdη[1:size(A_free,1),1:size(A_free,2)] = dA_freedη
            dMdη[(size(A_free,1)+1):(size(A_free,1)+size(B_free,2)),1:size(A_free,2)] = dB_free'
            dMdη[end,1:size(A_free,2)] = q_d_cached_top'*dAdη*C_uc_cached

            dMdη[1:size(A_free,1),(size(A_free,2)+1):(size(A_free,2)+size(B_free,2))] = dB_free
            dMdη[(size(A_free,1)+1):(size(A_free,1)+size(B_free,2)),(size(A_free,2)+1):(size(A_free,2)+size(B_free,2))] = zero
            dMdη[end,(size(A_free,2)+1):(size(A_free,2)+size(B_free,2))] = q_d_cached_top'*dB

            dMdη[1:size(A_free,1),end] = C_Tu*dAdη*q_d_cached_top
            dMdη[(size(A_free,1)+1):(size(A_free,1)+size(B_free,2)),end] = dB'*q_d_cached_top
            dMdη[end,end] = (q_d_cached_top'dAdη*q_d_cached_top)[end]
            
            dMdβ[1:size(A_free,1),1:size(A_free,2)] = dA_freedβ
            dMdβ[(size(A_free,1)+1):(size(A_free,1)+size(B_free,2)),1:size(A_free,2)] = dB_free'
            dMdβ[end,1:size(A_free,2)] = q_d_cached_top'*dAdβ*C_uc_cached

            dMdβ[1:size(A_free,1),(size(A_free,2)+1):(size(A_free,2)+size(B_free,2))] = dB_free
            dMdβ[(size(A_free,1)+1):(size(A_free,1)+size(B_free,2)),(size(A_free,2)+1):(size(A_free,2)+size(B_free,2))] = zero
            dMdβ[end,(size(A_free,2)+1):(size(A_free,2)+size(B_free,2))] = q_d_cached_top'*dB

            dMdβ[1:size(A_free,1),end] = C_Tu*dAdβ*q_d_cached_top
            dMdβ[(size(A_free,1)+1):(size(A_free,1)+size(B_free,2)),end] = dB'*q_d_cached_top
            dMdβ[end,end] = (q_d_cached_top'dAdβ*q_d_cached_top)[end]

            r = [-C_Tu*A*q_d; -B'*q_d; cParam_cached[iter].-q_d_cached_top'A*q_d]    # assemble the system of equations
            drdη = -[C_Tu*dAdη*q_d; zeros(Float64, size(B,2),size(q_d,2)); q_d_cached_top'dAdη*q_d] # solve the system of equations
            drdβ = -[C_Tu*dAdβ*q_d; zeros(Float64, size(B,2),size(q_d,2)); q_d_cached_top'dAdβ*q_d] # solve the system of equations

            lum = lu(M) # LU decomposition of the system of equations

            sol = lum\Matrix(r)            # solve the system of equations
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
            dqdη .= dqdη + C_uc_cached*dqfdη + dμdη*q_d_cached_top; # assemble the solution
            dqdβ .= dqdβ + C_uc_cached*dqfdβ + dμdβ*q_d_cached_top; # assemble the solution

            p = p_f'       # assemble the solution;
            dpdη = dpfdη'; # assemble the solution
            dpdβ = dpfdβ'; # assemble the solution

            velocity_field = @views hcat(q[ID_cached[1,:]], q[ID_cached[2,:]], q[ID_cached[3,:]])' # reshape the solution to get the velocity field
            dmdη_out_y = @views hcat(dqdη[ID_cached[1,:]], dqdη[ID_cached[2,:]], dqdη[ID_cached[3,:]])'*t_steps_cached
            dmdβ_out_y = @views hcat(dqdβ[ID_cached[1,:]], dqdβ[ID_cached[2,:]], dqdβ[ID_cached[3,:]])'*t_steps_cached
            
            motion_y = velocity_field*t_steps_cached # extract the motion of the mesh grid
            motion =  motion_y # extract the motion of the mesh grid

            NodeList_cached = NodeList_cached + motion # update the mesh grid
            mdl.mesh_x.NodeList = NodeList_cached      # update the mesh grid
            
            NodeList_proj = NodeList_cached # project the motion on the geometry mesh grid
            dmdη_out_proj = dmdη_out_y
            dmdβ_out_proj = dmdβ_out_y
            motion_proj = motion
            
            mat_nan_inf_check(dmdη_out_y)
            mat_nan_inf_check(dmdβ_out_y)
            
            dmdθ_out = @views cat(dmdη_out_proj,dmdβ_out_proj,dims=3) # concatenate the gradients in to a tensor

            BorderPts2D, dudθ, SurfacePts2D, ∇SurfacePts2D = extract_borders(NodeList_proj, camera_matrix_cached, obj_pose_cached, BorderNodesList=side_node_list_cached, GRAD=true, dqdθ=dmdθ_out, SIDES=SIDES_cached)
            pi, qi = fit_curve(border=BorderPts2D)
            
            push!(output, μ_tp*t_steps_cached) # store displacement at the top surface
            push!(velocity, velocity_field) # store the velocity of the mesh in 3D
            push!(pressure, p) # store the pressure of the mesh in 3D
            push!(displacement, motion_proj)
            push!(surface_fields, motion_proj[:,side_node_list_cached])
            push!(surface_pts_3D, vcat(NodeList_proj[:,top_node_list_cached]', NodeList_proj[:,bottom_node_list_cached]', NodeList_proj[:,side_node_list_cached]')')
            push!(gradList,dudθ)
            push!(pos2D, SurfacePts2D)
            push!(pos3D, NodeList_proj)
            push!(pos3D_cp, NodeList_cached)
            push!(borderPts2DList, BorderPts2D)
            push!(splinep, BorderPts2D[1,:])
            push!(splineq, BorderPts2D[2,:])
            push!(writeborderList, vcat(pi', qi'))
        
            iter += 1
            next!(pr, showvalues = [(:iterations,iter),(:time,string(t," seconds"))])
        end
    elseif control_cached == "velocity"
        for t in time
            _A_bar .= assemble_system_A(mdl)    # assemble the stiffness matrix
            B .= assemble_system_B(mdl)         # assemble the stiffness matrix
            b = apply_boundary_conditions(mdl) # apply the neumann boundary conditions
            q_d .= (μu_btm*q_d_cached_btm + cParam_cached[iter]*q_d_cached_top + μu_side*q_d_cached_brdr)      # apply the Dirichlet boundary conditions

            # println("Norm of _A_bar: ", size(_A_bar), " ", norm(_A_bar))
            # println("Norm of b: ", size(b), " ", norm(b))
        
            # display(_A_bar)
            # display(b)

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
            drdη = [C_Tu*dAdη*q_d; zeros(Float64, size(B,2),size(q_d,2))] # solve the system of equations
            drdβ = [C_Tu*dAdβ*q_d; zeros(Float64, size(B,2),size(q_d,2))] # solve the system of equations

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
            dqdη = dqdη + C_uc_cached*dqfdη;   # assemble the solution
            dqdβ = dqdβ + C_uc_cached*dqfdβ;   # assemble the solution

            p = p_f';       # assemble the solution
            dpdη = dpfdη';  # assemble the solution
            dpdβ = dpfdβ';  # assemble the solution

            velocity_field = @views hcat(q[ID_cached[1,:]], q[ID_cached[2,:]], q[ID_cached[3,:]])' # reshape the solution to get the velocity field
            dmdη_out_y = @views hcat(dqdη[ID_cached[1,:]], dqdη[ID_cached[2,:]], dqdη[ID_cached[3,:]])'*t_steps_cached
            dmdβ_out_y = @views hcat(dqdβ[ID_cached[1,:]], dqdβ[ID_cached[2,:]], dqdβ[ID_cached[3,:]])'*t_steps_cached
            
            motion_y = velocity_field*t_steps_cached # extract the motion of the mesh grid
            motion =  motion_y # extract the motion of the mesh grid

            NodeList_cached = NodeList_cached + motion # update the mesh grid
            mdl.mesh_x.NodeList = NodeList_cached      # update the mesh grid
            
            NodeList_proj = NodeList_cached # project the motion on the geometry mesh grid
            dmdη_out_proj = dmdη_out_y
            dmdβ_out_proj = dmdβ_out_y
            motion_proj = motion
            
            mat_nan_inf_check(dmdη_out_y)
            mat_nan_inf_check(dmdβ_out_y)
            
            dmdθ_out = @views cat(dmdη_out_proj,dmdβ_out_proj,dims=3) # concatenate the gradients in to a tensor

            BorderPts2D, dudθ, SurfacePts2D, ∇SurfacePts2D = extract_borders(NodeList_proj, camera_matrix_cached, obj_pose_cached, BorderNodesList=side_node_list_cached, GRAD=true, dqdθ=dmdθ_out, SIDES=SIDES_cached)
            pi, qi = fit_curve(border=BorderPts2D)
            
            # push!(output, μ_tp*t_steps_cached) # store displacement at the top surface
            push!(velocity, velocity_field) # store the velocity of the mesh in 3D
            push!(pressure, p) # store the pressure of the mesh in 3D
            push!(displacement, motion_proj)
            push!(surface_fields, motion_proj[:,side_node_list_cached])
            push!(surface_pts_3D, vcat(NodeList_proj[:,top_node_list_cached]', NodeList_proj[:,bottom_node_list_cached]', NodeList_proj[:,side_node_list_cached]')')
            push!(gradList,dudθ)
            push!(pos2D, SurfacePts2D)
            push!(pos3D, NodeList_proj)
            push!(pos3D_cp, NodeList_cached)
            push!(borderPts2DList, BorderPts2D)
            push!(splinep, BorderPts2D[1,:])
            push!(splineq, BorderPts2D[2,:])
            push!(writeborderList, vcat(pi', qi'))
        
            iter += 1
            next!(pr, showvalues = [(:iterations,iter),(:time,string(t," seconds"))])
        end
    else
            throw(ArgumentError("Control type not unknown"))
    end
    
    if conditions.WRITEVTK
        # write_scene(string(conditions.filepath,"/data"), NodeList_p_cached, mdl.mesh_p.IEN, mdl.ne, mdl.ndim, pressure, ID=ID_cached, FunctionClass=mdl.mesh_p.FunctionClass)
        write_stokes_scene(string(conditions.filepath,"/data"), mdl.mesh_u.NodeList, mdl.mesh_u.IEN, NodeList_p_cached, mdl.mesh_p.IEN, mdl.ne, mdl.ndim, velocity, pressure, pos3D=pos3D)
    end

    # write the data to a file
    if conditions.ANIMATE
        animate_fields(filepath = string(conditions.filepath,"/Results/images/"), Nodes=pos3D , IEN=IEN_u_cached, BorderNodes2D=borderPts2DList, fields2D=pos2D)
        animate_fields(filepath = string(conditions.filepath,"/Results/images/surface"), Nodes=surface_pts_3D)
        if FunctionClass_x_cached == "S2"
            animate_fields(filepath = string(conditions.filepath,"/Results/images/cp"), Nodes=pos3D_cp, IEN=IEN_x_cp_cached)
        end
    end
    if conditions.WRITECONTOUR
        write_data(string(conditions.filepath,"/data/sim_data/contour_data"), writeborderList)
    end
    
    # plt = set_plot(12,legend_column=4)
    # Plots.plot!(plt, time, dac_list, label=L"\delta^{T}_{U}\,A\,C")
    # Plots.plot!(plt, time, bd_list, label=L"\delta^{T}_{U}\,B^{T}")
    # Plots.plot!(plt, time, dad_list, label=L"\delta^{T}_{U}\,A\,\delta_{U}")
    # Plots.plot!(plt, time, A_list, label=L"A")
    # # Plots.ylims!(1e2, 1e8)
    # Plots.savefig(plt, string(conditions.filepath,"/dA_dc.pdf"))
    return output, gradList, borderPts2DList, displacement, surface_pts_3D, pos2D, pos3D, splinep, splineq, velocity, pressure
end