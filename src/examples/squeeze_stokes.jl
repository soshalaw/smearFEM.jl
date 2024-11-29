using LinearAlgebra
using ProgressMeter
using SparseArrays
using LinearAlgebra
using ProgressMeter
# using StaticArrays

""" 
    assemble_system(ne, NodeList, IEN, ndim, FunctionClass="Q1", nDof=1, ID=nothing, Young=1, ν=0.3)

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
function assemble_system_A(mdl::model)

    C = get_cMat("standard",1,0)
    IEN_u_rows = size(mdl.IEN_u,1)

    # (I,J,V) vectors for COO sparse matrix
    if mdl.nDof_u == 1
        E = zeros(  Int64, mdl.ne^mdl.ndim*IEN_u_rows^2)
        J = zeros(  Int64, mdl.ne^mdl.ndim*IEN_u_rows^2)
        V = zeros(Float64, mdl.ne^mdl.ndim*IEN_u_rows^2)
    else
        ID_u_rows = size(mdl.ID_u,1)
        E = zeros(  Int64, mdl.ne^mdl.ndim*((ID_u_rows*IEN_u_rows)^2))
        J = zeros(  Int64, mdl.ne^mdl.ndim*((ID_u_rows*IEN_u_rows)^2))
        V = zeros(Float64, mdl.ne^mdl.ndim*((ID_u_rows*IEN_u_rows)^2))  
    end

    # element loop
    if mdl.ndim == 1
        # gaussian quadrature points for the element [-1,1] 
        ξ, w_ξ = gaussian_quadrature(-1,1)
        
        wpoints =  [w_ξ[1], w_ξ[2]]
        
        x = [ξ[1], ξ[2]]
    elseif mdl.ndim == 2
        # gaussian quadrature points for the element [-1,1]x[-1,1] 
        ξ, w_ξ = gaussian_quadrature(-1,1,nGaussPoints=3)
        η, w_η = gaussian_quadrature(-1,1,nGaussPoints=3)

        x = Float64[]
        y = Float64[]
        wpoints =  Float64[]
        
        n = 1:size(ξ,1)
        m = 1:size(η,1)
        for j in m # loop over η
            for i in n # loop over ξ
                push!(x, ξ[i])
                push!(y, η[j])
                push!(wpoints, w_ξ[i]*w_η[j])
            end
        end

    elseif mdl.ndim == 3
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
        for k in l # loop over ζ
            for j in m # loop over η
                for i in n # loop over ξ
                    push!(x, ξ[i])
                    push!(y, η[j])
                    push!(z, ζ[k])
                    push!(wpoints, w_ξ[i]*w_η[j]*w_ζ[k])
                end
            end
        end
    end

    e_iter = 1:mdl.ne^mdl.ndim
    # integration loop
    gpiter = 1:length(wpoints)
    for gp in gpiter
        if mdl.ndim == 1
            N, ΔN = basis_function(x[gp], nothing, nothing, mdl.FunctionClass_u)
        elseif mdl.ndim == 2
            N, ΔN = basis_function(x[gp], y[gp], nothing, mdl.FunctionClass_u) 
        elseif mdl.ndim == 3
            N, ΔN = basis_function(x[gp], y[gp], z[gp], mdl.FunctionClass_u) 
        end
        # element loop
        for e in e_iter
            coords = mdl.NodeList_u[:,mdl.IEN_u[:,e]] # get the coordinates of the nodes of the element

            Jac  = coords*ΔN # Jacobian matrix [dx/dxi dx/deta; dy/dxi dy/deta]

            w = wpoints[gp]*abs(det(Jac))
            invJ = inv(Jac)
            dNdX = ΔN*invJ

            if mdl.nDof_u == 1
                szN = size(N,1) # number of basis functions
                # loop between basis functions of the element
                for i in 1:szN
                    for j in 1:szN
                        inz = (szN)^2*(e-1) + szN*(i-1) + j # index for the COO sparse matrix
                        E[inz] = mdl.IEN_u[i,e] # row index 
                        J[inz] = mdl.IEN_u[j,e] # column index
                        V[inz] += w*dot(dNdX[i,:],dNdX[j,:])# inner product of the gradient of the basis functions
                    end
                end
            else   
                if mdl.nDof_u == 2
                    B = zeros(Float64, mdl.ndim*length(N), 3)
                    B[1:mdl.nDof_u:end,1] = dNdX[:,1]
                    B[2:mdl.nDof_u:end,2] = dNdX[:,2]
                    B[1:mdl.nDof_u:end,3] = dNdX[:,2]
                    B[2:mdl.nDof_u:end,3] = dNdX[:,1]
                elseif mdl.nDof_u == 3
                    B = zeros(Float64, mdl.ndim*length(N), 6)
                    B[1:mdl.nDof_u:end,1] = dNdX[:,1]
                    B[2:mdl.nDof_u:end,2] = dNdX[:,2]
                    B[3:mdl.nDof_u:end,3] = dNdX[:,3]
                    B[2:mdl.nDof_u:end,4] = dNdX[:,3]
                    B[3:mdl.nDof_u:end,4] = dNdX[:,2]
                    B[1:mdl.nDof_u:end,5] = dNdX[:,3]
                    B[3:mdl.nDof_u:end,5] = dNdX[:,1]
                    B[1:mdl.nDof_u:end,6] = dNdX[:,2]
                    B[2:mdl.nDof_u:end,6] = dNdX[:,1]
                end
                Ke = 2*w*B*C*B' # element stiffness matrix
                
                Ke_row, Ke_col = size(Ke)
                Ke_len = length(Ke)

                # loop between basis functions of the element
                iNodes = 1:Ke_row÷mdl.nDof_u
                jNodes = 1:Ke_col÷mdl.nDof_u
                iDofs = 1:ID_u_rows
                jDofs = 1:ID_u_rows

                for iNode in iNodes
                    for jNode in jNodes
                        for iDof in iDofs
                            for jDof in jDofs
                                i = (iNode-1)*mdl.nDof_u + iDof
                                j = (jNode-1)*mdl.nDof_u + jDof
                                inz = Ke_len*(e-1) + (iNode-1)*mdl.nDof_u*Ke_col + (jNode-1)*mdl.nDof_u^2 + (iDof-1)*mdl.nDof_u + jDof # index for the COO sparse matrix
                                E[inz] = mdl.ID_u[iDof,mdl.IEN_u[iNode,e]] # row index 
                                J[inz] = mdl.ID_u[jDof,mdl.IEN_u[jNode,e]] # column index
                                V[inz] += Ke[i,j] 
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

function assemble_system_A_dense(mdl::model)

    C = get_cMat("standard",1,0)

    IEN_u_rows = size(mdl.IEN_u,1)
    ID_u_rows = size(mdl.ID_u,1)

    if mdl.nDof_u == 1
        sz = mdl.ne^mdl.ndim*IEN_u_rows # no elements x no nodes
        K = zeros(Float64, sz,sz)
    else
        sz = mdl.nDof_u*(2*mdl.ne+1)^mdl.ndim # no elements x no nodes x no dofs
        K = zeros(Float64, sz,sz)  
    end

    # element loop
    if mdl.ndim == 1
        # gaussian quadrature points for the element [-1,1] 
        ξ, w_ξ = gaussian_quadrature(-1,1)
        
        wpoints =  [w_ξ[1], w_ξ[2]]
        
        x = [ξ[1], ξ[2]]
    elseif mdl.ndim == 2
        # gaussian quadrature points for the element [-1,1]x[-1,1] 
        ξ, w_ξ = gaussian_quadrature(-1,1,nGaussPoints=3)
        η, w_η = gaussian_quadrature(-1,1,nGaussPoints=3)

        x = Float64[]
        y = Float64[]
        wpoints =  Float64[]
        
        n = 1:size(ξ,1)
        m = 1:size(η,1)
        for j in m # loop over η
            for i in n # loop over ξ
                push!(x, ξ[i])
                push!(y, η[j])
                push!(wpoints, w_ξ[i]*w_η[j])
            end
        end

    elseif mdl.ndim == 3
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
        for k in l # loop over ζ
            for j in m # loop over η
                for i in n # loop over ξ
                    push!(x, ξ[i])
                    push!(y, η[j])
                    push!(z, ζ[k])
                    push!(wpoints, w_ξ[i]*w_η[j]*w_ζ[k])
                end
            end
        end
    end

    e_iter = 1:mdl.ne^mdl.ndim
    # integration loop
    gpiter = 1:length(wpoints)
    for gp in gpiter
        if mdl.ndim == 1
            N, ΔN = basis_function(x[gp], nothing, nothing, mdl.FunctionClass_u)
        elseif mdl.ndim == 2
            N, ΔN = basis_function(x[gp], y[gp], nothing, mdl.FunctionClass_u) 
        elseif mdl.ndim == 3
            N, ΔN = basis_function(x[gp], y[gp], z[gp], mdl.FunctionClass_u) 
        end
        for e in e_iter
            coords = mdl.NodeList_u[:,mdl.IEN_u[:,e]] # get the coordinates of the nodes of the element

            Jac  = coords*ΔN # Jacobian matrix [dx/dxi dx/deta; dy/dxi dy/deta]

            w = wpoints[gp]*abs(det(Jac))
            invJ = inv(Jac)
            dNdX = ΔN*invJ

            if mdl.nDof_u == 1
                szN = size(N,1) # number of basis functions
                # loop between basis functions of the element
                for i in 1:szN
                    for j in 1:szN
                        inz = (szN)^2*(e-1) + szN*(i-1) + j # index for the COO sparse matrix
                        E[inz] = mdl.IEN_u[i,e] # row index 
                        J[inz] = mdl.IEN_u[j,e] # column index
                        V[inz] += w*dot(dNdX[i,:],dNdX[j,:])# inner product of the gradient of the basis functions
                    end
                end
            else   
                if mdl.nDof_u == 2
                    B = zeros(Float64, mdl.ndim*length(N), 3)
                    B[1:mdl.nDof_u:end,1] = dNdX[:,1]
                    B[2:mdl.nDof_u:end,2] = dNdX[:,2]
                    B[1:mdl.nDof_u:end,3] = dNdX[:,2]
                    B[2:mdl.nDof_u:end,3] = dNdX[:,1]
                elseif mdl.nDof_u == 3
                    B = zeros(Float64, mdl.ndim*length(N), 6)
                    B[1:mdl.nDof_u:end,1] = dNdX[:,1]
                    B[2:mdl.nDof_u:end,2] = dNdX[:,2]
                    B[3:mdl.nDof_u:end,3] = dNdX[:,3]
                    B[2:mdl.nDof_u:end,4] = dNdX[:,3]
                    B[3:mdl.nDof_u:end,4] = dNdX[:,2]
                    B[1:mdl.nDof_u:end,5] = dNdX[:,3]
                    B[3:mdl.nDof_u:end,5] = dNdX[:,1]
                    B[1:mdl.nDof_u:end,6] = dNdX[:,2]
                    B[2:mdl.nDof_u:end,6] = dNdX[:,1]
                end
                Ke = 2*w*B*C*B' # element stiffness matrix
                
                Ke_row, Ke_col = size(Ke)
                # loop between basis functions of the element
                iNodes = 1:Ke_row÷mdl.nDof_u
                jNodes = 1:Ke_col÷mdl.nDof_u
                iDofs = 1:ID_u_rows
                jDofs = 1:ID_u_rows
                for iNode in iNodes
                    for jNode in jNodes
                        for iDof in iDofs
                            for jDof in jDofs
                                i = (iNode-1)*mdl.nDof_u + iDof
                                j = (jNode-1)*mdl.nDof_u + jDof
                                K[mdl.ID[iDof,mdl.IEN[iNode,e]],mdl.ID[jDof,mdl.IEN[jNode,e]]] += Ke[i,j] 
                            end
                        end
                    end
                end
            end
        end
    end
    return K
end

function assemble_system_B(mdl::model)

    # (I,J,V) vectors for COO sparse matrix
    IEN_u_rows = size(mdl.IEN_u,1)
    IEN_p_rows = size(mdl.IEN_p,1)
    ID_u_rows = size(mdl.ID_u,1)

    if mdl.nDof_u == 1
        E = zeros(  Int64, mdl.ne^mdl.ndim*IEN_u_rows*IEN_p_rows)
        J = zeros(  Int64, mdl.ne^mdl.ndim*IEN_u_rows*IEN_p_rows)
        V = zeros(Float64, mdl.ne^mdl.ndim*IEN_u_rows*IEN_p_rows)
    else
        E = zeros(  Int64, mdl.ne^mdl.ndim*(ID_u_rows*IEN_u_rows)*IEN_p_rows)
        J = zeros(  Int64, mdl.ne^mdl.ndim*(ID_u_rows*IEN_u_rows)*IEN_p_rows)
        V = zeros(Float64, mdl.ne^mdl.ndim*(ID_u_rows*IEN_u_rows)*IEN_p_rows)  
    end

    # element loop
    if mdl.ndim == 1
        # gaussian quadrature points for the element [-1,1] 
        ξ, w_ξ = gaussian_quadrature(-1,1)
    
        wpoints = [w_ξ[1], w_ξ[2]]
        
        x = [ξ[1], ξ[2]]
    elseif mdl.ndim == 2
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

    elseif mdl.ndim == 3
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
        for k in l # loop over ζ
            for j in m # loop over η
                for i in n # loop over ξ
                    push!(x, ξ[i])
                    push!(y, η[j])
                    push!(z, ζ[k])
                    push!(wpoints, w_ξ[i]*w_η[j]*w_ζ[k])
                end
            end
        end
    end

    e_iter = 1:mdl.ne^mdl.ndim
    # integration loop
    gpiter = 1:length(wpoints)
    for gp in gpiter
        if mdl.ndim == 1
            N_u, ΔN_u = basis_function(x[gp], nothing, nothing, "Q2")
            N_p, ΔN_p = basis_function(x[gp], nothing, nothing, "Q1")
        elseif mdl.ndim == 2
            N_u, ΔN_u = basis_function(x[gp], y[gp], nothing, "Q2")
            N_p, ΔN_p = basis_function(x[gp], y[gp], nothing, "Q1")
        elseif mdl.ndim == 3
            N_u, ΔN_u = basis_function(x[gp], y[gp], z[gp], "Q2")
            N_p, ΔN_p = basis_function(x[gp], y[gp], z[gp], "Q1")
        end
        for e in e_iter
            coords_u = mdl.NodeList_u[:,mdl.IEN_u[:,e]] # get the coordinates of the nodes of the element

            Jac  = coords_u*ΔN_u # Jacobian matrix [dx/dxi dx/deta; dy/dxi dy/deta]

            w = wpoints[gp]*abs(det(Jac))
            invJ = inv(Jac)
            dNdX = ΔN_u*invJ

            if mdl.nDof_u == 1
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

                Be = zeros(Float64,1,3*size(dNdX,1))
                Be[1:mdl.nDof_u:end] = dNdX[:,1]
                Be[2:mdl.nDof_u:end] = dNdX[:,2]
                Be[3:mdl.nDof_u:end] = dNdX[:,3]
                Ke = -w*N_p*Be # element stiffness matrix
                
                Ke_row, Ke_col = size(Ke)
                # loop between basis functions of the element
                iNodes = 1:Ke_col÷mdl.nDof_u  # column node index
                iDofs = 1:ID_u_rows # column dof index
                jNodes = 1:Ke_row # row index

                for iNode in iNodes 
                    for jNode in jNodes 
                        for iDof in iDofs
                            i = (iNode-1)*mdl.nDof_u + iDof
                            j = jNode
                            inz = length(Ke)*(e-1) + (jNode-1)*Ke_col + (iNode-1)*mdl.nDof_u + iDof # index for the COO sparse matrix
                            E[inz] = mdl.ID_u[iDof,mdl.IEN_u[iNode,e]] # row index 
                            J[inz] = mdl.IEN_p[jNode,e] # column index
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

function assemble_system_B_dense(mdl::model)

    IEN_u_rows = size(mdl.IEN_u,1)
    ID_u_rows = size(mdl.ID_u,1)
    if mdl.nDof_u == 1
        sz = mdl.ne^mdl.ndim*IEN_u_rows # no elements x no nodes
        K = zeros(Float64, sz,sz)
    else
        sz1 = mdl.nDof_u*(2*mdl.ne+1)^mdl.ndim # no elements x no nodes x no dofs
        sz2 = (2*mdl.ne+1)^mdl.ndim
        K = zeros(Float64, sz1, sz2)  
    end

    # element loop
    if mdl.ndim == 1
        # gaussian quadrature points for the element [-1,1] 
        ξ, w_ξ = gaussian_quadrature(-1,1)
    
        wpoints = [w_ξ[1], w_ξ[2]]
        
        x = [ξ[1], ξ[2]]
    elseif mdl.ndim == 2
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

    elseif mdl.ndim == 3
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
        for k in l # loop over ζ
            for j in m # loop over η
                for i in n # loop over ξ
                    push!(x, ξ[i])
                    push!(y, η[j])
                    push!(z, ζ[k])
                    push!(wpoints, w_ξ[i]*w_η[j]*w_ζ[k])
                end
            end
        end
    end

    e_iter = 1:mdl.ne^mdl.ndim
    # integration loop
    gpiter = 1:length(wpoints)
    for gp in gpiter
        if mdl.ndim == 1
            N_u, ΔN_u = basis_function(x[gp], nothing, nothing, "Q2")
            N_p, ΔN_p = basis_function(x[gp], nothing, nothing, "Q1")
        elseif mdl.ndim == 2
            N_u, ΔN_u = basis_function(x[gp], y[gp], nothing, "Q2")
            N_p, ΔN_p = basis_function(x[gp], y[gp], nothing, "Q1")
        elseif mdl.ndim == 3
            N_u, ΔN_u = basis_function(x[gp], y[gp], z[gp], "Q2")
            N_p, ΔN_p = basis_function(x[gp], y[gp], z[gp], "Q1")
        end

        for e in e_iter
            coords_u = mdl.NodeList_u[:,mdl.IEN_u[:,e]] # get the coordinates of the nodes of the element

            Jac  = coords_u*ΔN_u # Jacobian matrix [dx/dxi dx/deta; dy/dxi dy/deta]

            w = wpoints[gp]*abs(det(Jac))
            invJ = inv(Jac)
            dNdX = ΔN_u*invJ

            if mdl.nDof_u == 1
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

                Be = zeros(Float64,1,3*size(dNdX,1))
                Be[1:mdl.nDof_u:end] = dNdX[:,1]
                Be[2:mdl.nDof_u:end] = dNdX[:,2]
                Be[3:mdl.nDof_u:end] = dNdX[:,3]
                Ke = -w*N_p*Be # element stiffness matrix
                
                Ke_row, Ke_col = size(Ke)
                # loop between basis functions of the element
                iNodes = 1:Ke_col÷mdl.nDof_u  # column node index
                iDofs = 1:ID_u_rows # column dof index
                jNodes = 1:Ke_row # row index

                for iNode in iNodes 
                    for jNode in jNodes 
                        for iDof in iDofs
                            i = (iNode-1)*mdl.nDof_u + iDof
                            j = jNode
                            K[mdl.ID_u[iDof,mdl.IEN_u[iNode,e]],mdl.IEN_p[jNode,e]] += Ke[i,j] 
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
function apply_boundary_conditions_stokes(mdl)

    ID_u_rows = size(mdl.ID_u,1)
    IEN_u_btm_rows = size(mdl.IEN_u_btm,1)

    E = zeros(  Int64, mdl.ne^(mdl.ndim-1)*(ID_u_rows*IEN_u_btm_rows)^2*2) # *2 because we have two surfaces
    J = zeros(  Int64, mdl.ne^(mdl.ndim-1)*(ID_u_rows*IEN_u_btm_rows)^2*2) # *2 because we have two surfaces
    V = zeros(Float64, mdl.ne^(mdl.ndim-1)*(ID_u_rows*IEN_u_btm_rows)^2*2) # *2 because we have two surfaces

    if mdl.ndim == 2
        # gaussian quadrature points for the element [-1,1] 
        ξ, w_ξ = gaussian_quadrature(-1,1)

        wpoints = [w_ξ[1], w_ξ[2]]
        
        x = [ξ[1], ξ[2]]
    elseif mdl.ndim == 3
        # gaussian quadrature points for the element [-1,1]x[-1,1] 
        ξ, w_ξ = gaussian_quadrature(-1,1)
        η, w_η = gaussian_quadrature(-1,1)
        
        wpoints = [w_ξ[1]*w_η[1], w_ξ[2]*w_η[1], w_ξ[2]*w_η[2], w_ξ[1]*w_η[2]]
        
        x = [ξ[1], ξ[2], ξ[2], ξ[1]]
        y = [η[1], η[1], η[2], η[2]]
    end

    e_iter = 1:mdl.ne^(mdl.ndim-1)
    # integration loop
    gpiter = 1:length(wpoints)
    for gp in gpiter

        if mdl.ndim == 2
            N, ΔN = basis_function(x[gp], nothing, nothing, mdl.FunctionClass_u)
        elseif mdl.ndim == 3
            N, ΔN = basis_function(x[gp], y[gp], nothing, mdl.FunctionClass_u) 
        end

        # element loop
        for e in e_iter
        
            coords_u_top = mdl.NodeList_u[:,mdl.IEN_u_top[:,e]] # get the coordinates of the nodes of the element
            coords_u_btm = mdl.NodeList_u[:,mdl.IEN_u_btm[:,e]] # get the coordinates of the nodes of the element

            dxdξ_top = coords_u_top*ΔN         # Jacobian matrix [dx/dxi dx/deta; dy/dxi dy/deta; dz/dxi dz/deta]
            dxdξ_btm = coords_u_btm*ΔN         # Jacobian matrix [dx/dxi dx/deta; dy/dxi dy/deta; dz/dxi dz/deta]

            w_top = wpoints[gp]*norm(cross(dxdξ_top[:,1],dxdξ_top[:,2]))     # weight of the quadrature point top surface
            w_btm = wpoints[gp]*norm(cross(dxdξ_btm[:,1],dxdξ_btm[:,2]))     # weight of the quadrature point bottom surface
            
            M = zeros(Float64, 3, mdl.ndim*length(N))
            M[1,1:mdl.nDof_u:end] = N
            M[2,2:mdl.nDof_u:end] = N
            M[3,3:mdl.nDof_u:end] = N

            be = M'*M

            be_row, be_col = size(be)
            len_be = length(be)
            # loop between basis functions of the element
            iNodes = 1:be_row÷mdl.nDof_u
            jNodes = 1:be_col÷mdl.nDof_u
            iDofs = 1:ID_u_rows
            jDofs = 1:ID_u_rows
            for iNode in iNodes
                for jNode in jNodes
                    for iDof in iDofs
                        for jDof in jDofs
                            i = (iNode-1)*mdl.nDof_u + iDof
                            j = (jNode-1)*mdl.nDof_u + jDof
                            
                            inz_btm = len_be*(e-1) + (iNode-1)*mdl.nDof_u*be_col + (jNode-1)*mdl.nDof_u^2 + (iDof-1)*mdl.nDof_u + jDof # index for the COO sparse matrix
                            inz_top = len_be*(mdl.ne)^(mdl.ndim-1)+ len_be*((e-1)) + (iNode-1)*mdl.nDof_u*be_col + (jNode-1)*mdl.nDof_u^2 + (iDof-1)*mdl.nDof_u + jDof # index for the COO sparse matrix
                            
                            E[inz_top] = mdl.ID_u[iDof,mdl.IEN_u_top[iNode,e]]    # row index 
                            J[inz_top] = mdl.ID_u[jDof,mdl.IEN_u_top[jNode,e]]   # column index
                            V[inz_top] += w_top*be[i,j] 

                            E[inz_btm] = mdl.ID_u[iDof,mdl.IEN_u_btm[iNode,e]]    # row index 
                            J[inz_btm] = mdl.ID_u[jDof,mdl.IEN_u_btm[jNode,e]]   # column index
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

function apply_boundary_conditions_stokes_dense(mdl)

    IEN_btm_rows = size(mdl.IEN_btm,1)
    IEN_rows = size(mdl.IEN,1)
    ID_u_rows = size(mdl.ID_u,1)
    sz =    mdl.nDof_u*(2*mdl.ne+1)^mdl.ndim
    K = zeros(Float64, sz,sz)

    if mdl.ndim == 2
        # gaussian quadrature points for the element [-1,1] 
        ξ, w_ξ = gaussian_quadrature(-1,1)

        wpoints = [w_ξ[1], w_ξ[2]]
        
        x = [ξ[1], ξ[2]]
    elseif mdl.ndim == 3
        # gaussian quadrature points for the element [-1,1]x[-1,1] 
        ξ, w_ξ = gaussian_quadrature(-1,1)
        η, w_η = gaussian_quadrature(-1,1)
        
        wpoints = [w_ξ[1]*w_η[1], w_ξ[2]*w_η[1], w_ξ[2]*w_η[2], w_ξ[1]*w_η[2]]
        
        x = [ξ[1], ξ[2], ξ[2], ξ[1]]
        y = [η[1], η[1], η[2], η[2]]
    end 
    e_iter = 1:mdl.ne^(mdl.ndim-1)
    # integration loop
    gpiter = 1:length(wpoints)
    for gp in gpiter

        if mdl.ndim == 2
            N, ΔN = basis_function(x[gp], nothing, nothing, mdl.FunctionClass_u)
        elseif mdl.ndim == 3
            N, ΔN = basis_function(x[gp], y[gp], nothing, mdl.FunctionClass_u) 
        end
        # element loop
        for e in e_iter
        
            coords_u_top = mdl.NodeList_u[:,mdl.IEN_u_top[:,e]] # get the coordinates of the nodes of the element
            coords_u_btm = mdl.NodeList_u[:,mdl.IEN_u_btm[:,e]] # get the coordinates of the nodes of the element

            dxdξ_top = coords_u_top*ΔN         # Jacobian matrix [dx/dxi dx/deta; dy/dxi dy/deta; dz/dxi dz/deta]
            dxdξ_btm = coords_u_btm*ΔN         # Jacobian matrix [dx/dxi dx/deta; dy/dxi dy/deta; dz/dxi dz/deta]

            w_top = wpoints[gp]*norm(cross(dxdξ_top[:,1],dxdξ_top[:,2]))     # weight of the quadrature point top surface
            w_btm = wpoints[gp]*norm(cross(dxdξ_btm[:,1],dxdξ_btm[:,2]))     # weight of the quadrature point bottom surface
            
            M = zeros(Float64, 3, mdl.ndim*length(N))
            M[1,1:mdl.nDof_u:end] = N
            M[2,2:mdl.nDof_u:end] = N
            M[3,3:mdl.nDof_u:end] = N

            be = M'*M

            be_row, be_col = size(be)
            len_be = length(be)

            # loop between basis functions of the element
            iNodes = 1:be_row÷mdl.nDof_u
            jNodes = 1:be_col÷mdl.nDof_u
            iDofs = 1:ID_u_rows
            jDofs = 1:ID_u_rows
            for iNode in iNodes
                for jNode in jNodes
                    for iDof in iDofs
                        for jDof in jDofs
                            i = (iNode-1)*mdl.nDof_u + iDof
                            j = (jNode-1)*mdl.nDof_u + jDof
                            
                            inz_btm = len_be*(e-1) + (iNode-1)*mdl.nDof_u*be_col + (jNode-1)*mdl.nDof_u^2 + (iDof-1)*mdl.nDof_u + jDof # index for the COO sparse matrix
                            inz_top = len_be*(mdl.ne)^(mdl.ndim-1)+ len_be*((e-1)) + (iNode-1)*mdl.nDof_u*be_col + (jNode-1)*mdl.nDof_u^2 + (iDof-1)*mdl.nDof_u + jDof # index for the COO sparse matrix
                            
                            E[inz_top] = mdl.ID_u[iDof,mdl.IEN_u_top[iNode,e]]    # row index 
                            J[inz_top] = mdl.ID_u[jDof,mdl.IEN_u_top[jNode,e]]   # column index
                            V[inz_top] += w_top*be[i,j] 

                            E[inz_btm] = mdl.ID_u[iDof,mdl.IEN_u_btm[iNode,e]]    # row index 
                            J[inz_btm] = mdl.ID_u[jDof,mdl.IEN_u_btm[jNode,e]]   # column index
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
- `C_uc::SparseMatrixCSC{Float64,Int64}` : onstraint matrix
"""
function set_boundary_cond_stokes(NodeList, ne, ndim, FunctionClass, nDof=1)
    if FunctionClass == "Q1"
        q_upper = zeros(nDof*(ne+1)^ndim,1)                  # initialize the vector of the Dirichlet boundary conditions (for ndof = 1) / Dirichlet boundary conditions upper surface (for ndof > 1)
        q_lower = zeros(nDof*(ne+1)^ndim,1)                  # initialize the vector of the Neumann boundary conditions (for ndof = 1) / Dirichlet boundary conditions lower surface (for ndof > 1)
        q_side = zeros(nDof*(ne+1)^ndim,1)                   # initialize the vector of the Neumann boundary conditions (for ndof = 1) / Dirichlet boundary conditions lower surface (for ndof > 1)
        C = sparse(I,nDof*(ne+1)^ndim,nDof*(ne+1)^ndim)      # definition of the constraint matrix
    elseif FunctionClass == "Q2"
        q_upper = zeros(nDof*(2*ne+1)^ndim,1)                # initialize the vector of the Dirichlet boundary conditions (for ndof = 1) / Dirichlet boundary conditions upper surface (for ndof > 1)
        q_lower = zeros(nDof*(2*ne+1)^ndim,1)                # initialize the vector of the Neumann boundary conditions (for ndof = 1) / Dirichlet boundary conditions lower surface (for ndof > 1)
        q_side = zeros(nDof*(2*ne+1)^ndim,1)                  # initialize the vector of the Neumann boundary conditions (for ndof = 1) / Dirichlet boundary conditions lower surface (for ndof > 1)
        C = sparse(I,nDof*(2*ne+1)^ndim,nDof*(2*ne+1)^ndim)  # definition of the constraint matrix
    end
    if nDof == 1
        if ndim == 3
            Dbound1 = 0
            Dbound2 = 1
            iter = 1:size(NodeList,2)
            for n in iter
                coord = NodeList[:,n] # get the coordinates of the node
                if coord[3] == Dbound1 # bottom boundary
                    q_upper[n] = 0
                elseif coord[3] == Dbound2 # top boundary
                    q_upper[n] = -1
                end
            end
        elseif ndim == 2
            Dbound1 = 0
            Dbound2 = 1
            iter = 1:size(NodeList,2)
            for n in iter
                coord = NodeList[:,n] # get the coordinates of the node
                if coord[2] == Dbound1 # bottom boundary
                    q_upper[n] = 0
                elseif coord[2] == Dbound2 # top boundary
                    q_upper[n] = -1
                end
            end
        end

        if FunctionClass == "Q1"
            C_uc = C[:,((ne+1)^(ndim-1)+1):((ne+1)^ndim-(ne+1)^(ndim-1))]
        elseif FunctionClass == "Q2"
            C_uc = C[:,((2*ne+1)^(ndim-1)+1):((2*ne+1)^ndim-(2*ne+1)^(ndim-1))]
        end
    else
        z0Bound = 0
        z1Bound = 1

        rCol = Array{Int}(undef,0)
        iter = 1:size(NodeList,2)
        for nNode in iter
            coord = NodeList[:,nNode]    # get the coordinates of the node
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
Set the Dirichlet boundary conditions for the problem

# Arguments:
- `NodeList::Matrix{Float64}{nNodes,ndim}` : array of nodes
- `ne::Int` : number of elements
- `ndim::Int` : number of dimensions
- `FunctionClass::String` : type of basis function

# Returns:
- `q_upper::Vector{Float64}` : vector of the Dirichlet boundary conditions (for ndof = 1) / Dirichlet boundary conditions upper surface (for ndof > 1)
- `q_lower::Vector{Float64}` : vector of the Neumann boundary conditions (for ndof = 1) / Dirichlet boundary conditions lower surface (for ndof > 1)
- `C_uc::SparseMatrixCSC{Float64,Int64}` : onstraint matrix
"""
function set_boundary_cond_stokes_dense(NodeList, ne, ndim, FunctionClass, nDof=1)
    if FunctionClass == "Q1"
        q_upper = zeros(nDof*(ne+1)^ndim,1)                  # initialize the vector of the Dirichlet boundary conditions (for ndof = 1) / Dirichlet boundary conditions upper surface (for ndof > 1)
        q_lower = zeros(nDof*(ne+1)^ndim,1)                  # initialize the vector of the Neumann boundary conditions (for ndof = 1) / Dirichlet boundary conditions lower surface (for ndof > 1)
        q_side = zeros(nDof*(ne+1)^ndim,1)                   # initialize the vector of the Neumann boundary conditions (for ndof = 1) / Dirichlet boundary conditions lower surface (for ndof > 1)
        C = zeros(I,nDof*(ne+1)^ndim,nDof*(ne+1)^ndim)      # definition of the constraint matrix
    elseif FunctionClass == "Q2"
        q_upper = zeros(nDof*(2*ne+1)^ndim,1)                # initialize the vector of the Dirichlet boundary conditions (for ndof = 1) / Dirichlet boundary conditions upper surface (for ndof > 1)
        q_lower = zeros(nDof*(2*ne+1)^ndim,1)                # initialize the vector of the Neumann boundary conditions (for ndof = 1) / Dirichlet boundary conditions lower surface (for ndof > 1)
        q_side = zeros(nDof*(2*ne+1)^ndim,1)                  # initialize the vector of the Neumann boundary conditions (for ndof = 1) / Dirichlet boundary conditions lower surface (for ndof > 1)
        C = zeros(I,nDof*(2*ne+1)^ndim,nDof*(2*ne+1)^ndim)  # definition of the constraint matrix
    end
    if nDof == 1
        if ndim == 3
            Dbound1 = 0
            Dbound2 = 1
            iter = 1:size(NodeList,2)
            for n in iter
                coord = NodeList[:,n] # get the coordinates of the node
                if coord[3] == Dbound1 # bottom boundary
                    q_upper[n] = 0
                elseif coord[3] == Dbound2 # top boundary
                    q_upper[n] = -1
                end
            end
        elseif ndim == 2
            Dbound1 = 0
            Dbound2 = 1
            iter = 1:size(NodeList,2)
            for n in iter
                coord = NodeList[:,n] # get the coordinates of the node
                if coord[2] == Dbound1 # bottom boundary
                    q_upper[n] = 0
                elseif coord[2] == Dbound2 # top boundary
                    q_upper[n] = -1
                end
            end
        end

        if FunctionClass == "Q1"
            C_uc = C[:,((ne+1)^(ndim-1)+1):((ne+1)^ndim-(ne+1)^(ndim-1))]
        elseif FunctionClass == "Q2"
            C_uc = C[:,((2*ne+1)^(ndim-1)+1):((2*ne+1)^ndim-(2*ne+1)^(ndim-1))]
        end
    else
        z0Bound = 0
        z1Bound = 1

        rCol = Array{Int}(undef,0)
        iter = 1:size(NodeList,2)
        for nNode in iter
            coord = NodeList[:,nNode]    # get the coordinates of the node
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
simulate_stokes(x0, x1, y0, y1, z0, z1, ne, Young, ν, ndim, FunctionClass, nDof, β, CameraMatrix, endTime, tSteps, Control, cParam, cMat; writeData=false, filepath=nothing)

Simulate the Stokes problem

# Arguments:
- `x0::Float64` : x coordinate of the lower left corner of the domain
- `x1::Float64` : x coordinate of the upper right corner of the domain
- `y0::Float64` : y coordinate of the lower left corner of the domain
- `y1::Float64` : y coordinate of the upper right corner of the domain
- `z0::Float64` : z coordinate of the lower left corner of the domain
- `z1::Float64` : z coordinate of the upper right corner of the domain
- `ne::Int` : number of elements
- `Young::Float64` : Young's modulus
- `ν::Float64` : Poisson's ratio
- `ndim::Int` : number of dimensions
- `FunctionClass::String` : type of basis function
- `nDof::Int` : number of degree of freedom per node
- `β::Float64` : penalty parameter
- `CameraMatrix::Matrix{Float64}` : camera matrix
- `endTime::Float64` : end time
- `tSteps::Int` : number of time steps
- `Control::String` : type of control
- `cParam::Float64` : control parameter
- `cMat::Matrix{Float64}` : control matrix

# Optional arguments:
- `writeData::Bool` : write the data to a file
- `filepath::String` : path to write the data

# Returns:
- `fields::Vector{Float64}` : Nodal positions of the mesh
"""
function simulate_stokes(x0, x1, y0, y1, z0, z1, ne, Young, ν, ndim, FunctionClass, nDof, β, CameraMatrix, endTime, tSteps, Control, cParam, cMat; writeData=false, filepath=nothing)

    time = collect(range(start=0,stop=endTime,length=tSteps)) # time vector
    dateTime = Dates.now()
    CameraMatrix = [[8*2048/7.07, 0.0, 2048/2] [0.0, 8*1536/5.3, 1536/2] [0.0, 0.0, 1.0]]
    
    x0 = 0
    x1 = 1
    y0 = 0
    y1 = 1
    z0 = 0
    z1 = 1
    ne = 4
    ndim = 3
    FunctionClass_u = "Q2"
    FunctionClass_p = "Q1"
    nDof_u = ndim  # number of degree of freedom per node
    nDof_p = 1
    β = 100
    ν = 1
    Control = "displacement"
    
    filepath = string("/home/soshala/SMEAR-PhD/SMEAR/Data/sim_experiments/single_simulation/Stokes/fem_runs/",Control,"/",dateTime,"/")
    writeData = true
    
    if writeData
        isnothing(filepath) || AssertionError("Please provide a filepath to write the data")
        set_file(filepath)
    end
    
    μu_btm = 0  
    μu_tp = -0.02
    μu_side = 0
    
    μp_btm = 0
    μp_tp = 0
    μp_side = 0
    
    NodeList_u, IEN_u, ID_u, IEN_u_top, IEN_u_btm, BorderNodesList_u = meshgrid_cube(x0,x1,y0,y1,z0,z1,ne,ndim,FunctionClass=FunctionClass_u)  # generate the mesh grid
    NodeListCylinder = inflate_cylinder(NodeList_u, x0, x1, y0, y1)
    q_tp, q_side, q_btm, C_uc = set_boundary_cond_stokes(NodeList_u, ne, ndim, FunctionClass_u, nDof_u)
    
    NodeList_p, IEN_p, ID_p, IEN_p_top, IEN_p_btm, BorderNodesList_p = meshgrid_cube(x0,x1,y0,y1,z0,z1,ne,ndim,FunctionClass=FunctionClass_p)  # generate the mesh grid
    NodeListCylinderp = inflate_cylinder(NodeList_p, x0, x1, y0, y1)
    
    mdl = def_model("stokes", ne=ne, NodeList=NodeListCylinder, IEN=IEN_u, IEN_top=IEN_u_top, IEN_btm=IEN_u_btm, ndim=ndim, nDof=nDof_u, 
                        FunctionClass=FunctionClass_u, ID=ID_u, IEN_2=IEN_p, IEN_2_top=IEN_p_top, IEN_2_btm=IEN_p_btm, 
                            ndim_2=ndim, nDof_2=nDof_p, FunctionClass_2=FunctionClass_p ) # define the model
    
    state = "init"
    
    BorderPts2D, BorderNodes2D, Nodes2D = extract_borders(NodeListCylinder, CameraMatrix, BorderNodesList_u, state, ne, 2*ne+1)
    pi, qi = fit_curve(border=BorderPts2D)
    
    SideBorders = BorderNodesList_u[1]
    BottomBorders = BorderNodesList_u[2]
    TopBorders = BorderNodesList_u[3]
        
    fields = AbstractArray[]  
    pos3D = AbstractArray[NodeListCylinder]                                                             # store the solution fields of the mesh in 3D
    pos2D = AbstractArray[Nodes2D]                                                                   # store the solution fields of the mesh in 2D
    surfaceNodesList = Float64[NodeList[:,SideBorders] NodeList[:,BottomBorders] NodeList[:,TopBorders]]  # store the solution fields of the surfaces in 3D
    borderPts2DList = AbstractArray[BorderPts2D]                                                               # store the solution fields of the surfaces in 2D
    borderNodeList2D = AbstractArray[BorderNodes2D]                                                       # store the solution fields of the border nodes in 2D
    splinep = AbstractArray[pi]                                                                            # store the x coordinates samples of the spline parameters of the border nodes
    splineq = AbstractArray[qi]                                                                            # store the y coordinates samples of the spline parameters of the border nodes
    output = Float64[] 
    writeborderList = [vcat(pi', qi')]
    
    state = "update"
    iter = 1
    
    pr = Progress(tSteps; desc="Simulating with prescribed $Control ...", showspeed=true)
    for t in time
        A_bar = assemble_system_A(mdl)                   # assemble the stiffness matrix
        B = assemble_system_B(mdl)                   # assemble the stiffness matrix
        b = apply_boundary_conditions_stokes(mdl)           # apply the neumann boundary conditions
    
        q_d = (μu_btm*q_btm + μu_tp*q_tp + μu_side*q_side)      # apply the Dirichlet boundary conditions
    
        A = A_bar + β*b
    
        C_Tu = transpose(C_uc)           # transpose the constraint matrix
    
        A_free = C_Tu*A*C_uc        # extract the free part of the stiffness matrix
        B_free = C_Tu*B             # extract the free part of the stiffness matrix
    
        K_free = [A_free B_free; B_free' zeros(size(B_free,2),size(B_free,2))]     # assemble the system of equations
    
        r = -[C_Tu*A*q_d; B'*q_d]    # assemble the system of equations
        sol = K_free\r                 # solve the system of equations
    
        q_f = sol[1:size(A_free,1)]     # extract the free part of the solution
        p_f = sol[size(A_free,1)+1:end] # extract the free part of the solution 
    
        q = q_d + C_uc*q_f;                 # assemble the solution 
        p = p_f;
    
        motion = [q[ID_u[1,:]] q[ID_u[2,:]] q[ID_u[3,:]]]'
    
        NodeListCylinder = NodeListCylinder + motion*(endTime/tSteps) # update the mesh grid
    
        BorderPts2D, BorderNodes2D, Nodes2D = extract_borders(NodeListCylinder, CameraMatrix, BorderNodesList_u, state)
        surfaceNodesList = [NodeListCylinder[:,SideBorders] NodeListCylinder[:,BottomBorders] NodeListCylinder[:,TopBorders]]
        pi, qi = fit_curve(border=BorderPts2D)
    
        # store the solutions in a list
        # push!(output, F_est[1])
        push!(fields, motion)
        push!(pos2D, Nodes2D)
        push!(pos3D, NodeListCylinder)
        push!(borderPts2DList, BorderPts2D)
        push!(borderNodeList2D, BorderNodes2D)
        push!(splinep, pi)
        push!(splineq, qi) 
        push!(writeborderList, vcat(pi', qi'))
    
        iter += 1
        next!(pr, showvalues = [(:iterations,iter),(:time,t)])
    end
    
    # q_new, IEN_new = rearrange(q, ne, ndim, IEN_u, FunctionClass_u, ID_u) # rearrange the solution
    
    if writeData
        write_scene(string(filepath,"/Results"), NodeListCylinder, IEN_u, ne, ndim, pos3D, ID=ID_u, FunctionClass=FunctionClass_u)
        animate_fields(filepath = string(filepath,"/Results/images"), fields=pos3D , IEN=IEN_u, BorderNodes2D=borderPts2DList, fields2D=pos2D)
        write_contour_data(string(filepath,"/Results"), writeborderList)
    end
    
    return fields
end

