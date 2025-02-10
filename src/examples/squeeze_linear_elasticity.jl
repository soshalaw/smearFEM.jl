using LinearAlgebra
using ProgressMeter
using SparseArrays
using StaticArrays

# set up mesh grid

""" 
    assemble_system(ne, NodeList, IEN, ndim, FunctionClass="Q1", nDof=1, ID=nothing, Young=1, ν=0.3)

Assembles the finite element system. # Returns the global stiffness matrix

# Arguments:
- `ne::Interger`: number of elements in each direction
- `NodeList::Matrix{Float64}{ndim,nNodes}` : coordinates of the nodes
- `IEN::Matrix{Int}{nLocalNodes,nElements}` : connectivity matrix
- `ndim::Interger`: number of dimensions
- `nDof::Interger`: number of degree of freedom per node
- `FunctionClass::String`: type of basis functions to be considered (Q1:quadratic or Q2:Lagrange)
- `ID::Matrix{Int}{nDof,nNodes}` : matrix that maps the global degrees of freedom to the local degrees of freedom
- `Young::Float64`: Young's modulus
- `ν::Float64`: Poisson's ratio

# Returns:
- `K::SparseMatrixCSC{Float64,Int64}{ndof,ndof}` : sparse stiffness matrix 
"""
function assemble_system(mdl::model; GRAD::Bool=false)
    
    # (I,J,V) vectors for COO sparse matrix
    IEN_rows = size(mdl.IEN,1)
    if mdl.nDof == 1
        E = zeros(  Int64, mdl.ne^mdl.ndim*IEN_rows^2)
        J = zeros(  Int64, mdl.ne^mdl.ndim*IEN_rows^2)
        V = zeros(Float64, mdl.ne^mdl.ndim*IEN_rows^2)
    else
        ID_rows = size(mdl.ID,1)
        E = zeros(  Int64, mdl.ne^mdl.ndim*((ID_rows*IEN_rows)^2))
        J = zeros(  Int64, mdl.ne^mdl.ndim*((ID_rows*IEN_rows)^2))
        V = zeros(Float64, mdl.ne^mdl.ndim*((ID_rows*IEN_rows)^2))  

        if GRAD == true
            dE = zeros(  Int64, mdl.ne^mdl.ndim*((ID_rows*IEN_rows)^2))
            dJ = zeros(  Int64, mdl.ne^mdl.ndim*((ID_rows*IEN_rows)^2))
            dV = zeros(Float64, mdl.ne^mdl.ndim*((ID_rows*IEN_rows)^2))  
        end
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
            N, ΔN = basis_function(x[gp], nothing, nothing, mdl.FunctionClass)
        elseif mdl.ndim == 2
            N, ΔN = basis_function(x[gp], y[gp], nothing, mdl.FunctionClass) 
        elseif mdl.ndim == 3
            N, ΔN = basis_function(x[gp], y[gp], z[gp], mdl.FunctionClass) 
        end

        for e in e_iter
            coords = mdl.NodeList[:,mdl.IEN[:,e]] # get the coordinates of the nodes of the element
            Jac  = coords*ΔN # Jacobian matrix [dx/dxi dx/deta; dy/dxi dy/deta]

            w = wpoints[gp]*abs(det(Jac))
            invJ = inv(Jac)
            dNdX = ΔN*invJ

            if mdl.nDof == 1
                szN = size(N,1) # number of basis functions
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
                @argcheck !isnothing(mdl.cMat) "constitutive matrix must be provided for problems with nDof > 1"

                if mdl.nDof == 2
                    B = zeros(Float64, mdl.ndim*length(N), 3)
                    B[1:mdl.nDof:end,1] = dNdX[:,1]
                    B[2:mdl.nDof:end,2] = dNdX[:,2]
                    B[1:mdl.nDof:end,3] = dNdX[:,2]
                    B[2:mdl.nDof:end,4] = dNdX[:,1]

                    @argcheck size(mdl.cMat) == (2,2) "constitutive matrix must be 2x2 for plane stress problems"
                elseif mdl.nDof == 3
                    B = zeros(Float64, mdl.ndim*length(N), 6)
                    B[1:mdl.nDof:end,1] = dNdX[:,1]
                    B[2:mdl.nDof:end,2] = dNdX[:,2]
                    B[3:mdl.nDof:end,3] = dNdX[:,3]
                    B[2:mdl.nDof:end,4] = dNdX[:,3]
                    B[3:mdl.nDof:end,4] = dNdX[:,2]
                    B[1:mdl.nDof:end,5] = dNdX[:,3]
                    B[3:mdl.nDof:end,5] = dNdX[:,1]
                    B[1:mdl.nDof:end,6] = dNdX[:,2]
                    B[2:mdl.nDof:end,6] = dNdX[:,1]

                    @argcheck size(mdl.cMat) == (6,6) "constitutive matrix must be 6x6 for 3D problems"
                end

                Ke = B*mdl.cMat*B'*w # element stiffness matrix
                dKedλ = B*mdl.dcMatdλ*B'*w # element stiffness matrix
        
                # code optimization
                Ke_rows, Ke_cols = size(Ke)
                Ke_len = length(Ke)

                # loop between basis functions of the element
                iNodes = 1:Ke_rows÷mdl.nDof
                jNodes = 1:Ke_cols÷mdl.nDof
                iDofs = 1:ID_rows
                jDofs = 1:ID_rows

                for iNode in iNodes
                    for jNode in jNodes
                        for iDof in iDofs
                            for jDof in jDofs
                                i = (iNode-1)*mdl.nDof + iDof
                                j = (jNode-1)*mdl.nDof + jDof
                                inz = Ke_len*(e-1) + (iNode-1)*mdl.nDof*Ke_cols + (jNode-1)*mdl.nDof^2 + (iDof-1)*mdl.nDof + jDof # index for the COO sparse matrix
                                E[inz] = mdl.ID[iDof,mdl.IEN[iNode,e]] # row index 
                                J[inz] = mdl.ID[jDof,mdl.IEN[jNode,e]] # column index
                                V[inz] += Ke[i,j] 
                                
                                if GRAD == true
                                    dE[inz] = mdl.ID[iDof,mdl.IEN[iNode,e]] # row index 
                                    dJ[inz] = mdl.ID[jDof,mdl.IEN[jNode,e]] # column index
                                    dV[inz] += dKedλ[i,j] 
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    K = sparse(E,J,V)

    if GRAD == true
        dKdλ = sparse(dE,dJ,dV)
        return K, dKdλ
    else
        return K
    end
end

function assemble_system_dense(mdl::model; ∂C = zeros(1,1))

    # (I,J,V) vectors for COO sparse matrix
    IEN_rows = size(mdl.IEN,1)
    ID_rows = size(mdl.ID,1)
    if mdl.nDof == 1
        sz = mdl.ne^mdl.ndim*IEN_rows # no elements x no nodes
        K = zeros(Float64, sz,sz)
    else
        sz = mdl.nDof*(2*mdl.ne+1)^mdl.ndim # no elements x no nodes x no dofs
        K = zeros(Float64, sz,sz)  
    end

    if ∂C != zeros(1,1)
        mdl.cMat = ∂C # replace cMat with ∂C/∂θ_i for getting ∂K/∂θ_i
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
    vol = 0
    
    e_iter = 1:mdl.ne^mdl.ndim
    # integration loop
    gpiter = 1:length(wpoints)
    for gp in gpiter

        if mdl.ndim == 1
            N, ΔN = basis_function(x[gp], nothing, nothing, mdl.FunctionClass)
        elseif mdl.ndim == 2
            N, ΔN = basis_function(x[gp], y[gp], nothing, mdl.FunctionClass) 
        elseif mdl.ndim == 3
            N, ΔN = basis_function(x[gp], y[gp], z[gp], mdl.FunctionClass) 
        end
        # element loop
        for e in e_iter
            coords = mdl.NodeList[:,mdl.IEN[:,e]] # get the coordinates of the nodes of the element

            Jac  = coords*ΔN # Jacobian matrix [dx/dxi dx/deta; dy/dxi dy/deta]

            w = wpoints[gp]*abs(det(Jac))
            invJ = inv(Jac)
            dNdX = ΔN*invJ

            if mdl.nDof == 1
                szN = size(N,1) # number of basis functions
                # loop between basis functions of the element
                for i in 1:szN
                    for j in 1:szN
                        K[mdl.IEN[i,e],mdl.IEN[j,e]] = w*dot(dNdX[i,:],dNdX[j,:])# inner product of the gradient of the basis functions
                    end
                end
            else   
                @argcheck !isnothing(mdl.cMat) "constitutive matrix must be provided for problems with nDof > 1"

                if mdl.nDof == 2
                    B = zeros(Float64, mdl.ndim*length(N), 3)
                    B[1:mdl.nDof:end,1] = dNdX[:,1]
                    B[2:mdl.nDof:end,2] = dNdX[:,2]
                    B[1:mdl.nDof:end,3] = dNdX[:,2]
                    B[2:mdl.nDof:end,4] = dNdX[:,1]

                    @argcheck size(mdl.cMat) == (2,2) "constitutive matrix must be 2x2 for plane stress problems"
                elseif mdl.nDof == 3
                    B = zeros(Float64, mdl.ndim*length(N), 6)
                    B[1:mdl.nDof:end,1] = dNdX[:,1]
                    B[2:mdl.nDof:end,2] = dNdX[:,2]
                    B[3:mdl.nDof:end,3] = dNdX[:,3]
                    B[2:mdl.nDof:end,4] = dNdX[:,3]
                    B[3:mdl.nDof:end,4] = dNdX[:,2]
                    B[1:mdl.nDof:end,5] = dNdX[:,3]
                    B[3:mdl.nDof:end,5] = dNdX[:,1]
                    B[1:mdl.nDof:end,6] = dNdX[:,2]
                    B[2:mdl.nDof:end,6] = dNdX[:,1]

                    @argcheck size(mdl.cMat) == (6,6) "constitutive matrix must be 6x6 for 3D problems"
                end

                Ke = B*mdl.cMat*B'*w # element stiffness matrix
        
                # code optimization
                Ke_rows, Ke_cols = size(Ke)
                # Ke_len = length(Ke)

                # loop between basis functions of the element
                iNodes = 1:Ke_rows÷mdl.nDof
                jNodes = 1:Ke_cols÷mdl.nDof
                iDofs = 1:ID_rows
                jDofs = 1:ID_rows

                for iNode in iNodes
                    for jNode in jNodes
                        for iDof in iDofs
                            for jDof in jDofs
                                i = (iNode-1)*mdl.nDof + iDof
                                j = (jNode-1)*mdl.nDof + jDof
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

""" Apply the Neumann slip boundary conditions to the global stiffness matrix
# Arguments:
K: {[ndof,ndof] SparseMatrixCSC{Float64,Int64}} : sparse stiffness matrix 
ID: {[nNodes,nDof] Matrix{Int}} : matrix that maps the global degrees of freedom to the local degrees of freedom
q_d: {[ndof] Vector{Float64}} : Dirichlet boundary conditions
q_n: {[ndof] Vector{Float64}} : Neumann boundary conditions

# Returns:
K: {[ndof,ndof] SparseMatrixCSC{Float64,Int64}} : sparse stiffness matrix with the boundary conditions applied
F: {[ndof] Vector{Float64}} : force vector
"""
function set_slip_conditions(mdl::model)

    IEN_btm_rows = size(mdl.IEN_btm,1)
    ID_rows = size(mdl.ID,1)
    E = zeros(  Int64, mdl.ne^(mdl.ndim-1)*(ID_rows*IEN_btm_rows)^2*2) # *2 because we have two surfaces
    J = zeros(  Int64, mdl.ne^(mdl.ndim-1)*(ID_rows*IEN_btm_rows)^2*2) # *2 because we have two surfaces
    V = zeros(Float64, mdl.ne^(mdl.ndim-1)*(ID_rows*IEN_btm_rows)^2*2) # *2 because we have two surfaces

    if mdl.ndim == 2
        # gaussian quadrature points for the element [-1,1] 
        ξ, w_ξ = gaussian_quadrature(-1,1)

        wpoints =  [w_ξ[1], w_ξ[2]]
        
        x = [ξ[1], ξ[2]]
    elseif mdl.ndim == 3            
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
    end 

    e_iter = 1:mdl.ne^(mdl.ndim-1)
    # integration loop
    gpiter = 1:length(wpoints)
    for gp in gpiter
        # element loop
        for e in e_iter
            if mdl.ndim == 2
                N, ΔN = basis_function(x[gp], nothing, nothing, mdl.FunctionClass) 
            elseif mdl.ndim == 3
                N, ΔN = basis_function(x[gp], y[gp], nothing, mdl.FunctionClass) 
            end
        
            coords_top = mdl.NodeList[:,mdl.IEN_top[:,e]] # get the coordinates of the nodes of the element
            coords_btm = mdl.NodeList[:,mdl.IEN_btm[:,e]] # get the coordinates of the nodes of the element

            dxdξ_top = coords_top*ΔN         # Jacobian matrix [dx/dxi dx/deta; dy/dxi dy/deta; dz/dxi dz/deta]
            dxdξ_btm = coords_btm*ΔN         # Jacobian matrix [dx/dxi dx/deta; dy/dxi dy/deta; dz/dxi dz/deta]

            w_top = wpoints[gp]*norm(cross(dxdξ_top[:,1],dxdξ_top[:,2]))     # weight of the quadrature point top surface
            w_btm = wpoints[gp]*norm(cross(dxdξ_btm[:,1],dxdξ_btm[:,2]))     # weight of the quadrature point bottom surface
            
            M = zeros(3, mdl.ndim*length(N))
            M[1,1:mdl.nDof:end] = N
            M[2,2:mdl.nDof:end] = N
            M[3,3:mdl.nDof:end] = N

            be = M'*M
            be_cols, be_rows = size(be)
            be_len = length(be)

            # loop between basis functions of the element
            iNodes = 1:be_cols÷mdl.nDof
            jNodes = 1:be_rows÷mdl.nDof
            iDofs = 1:ID_rows
            jDofs = 1:ID_rows
            for iNode in iNodes
                for jNode in jNodes
                    for iDof in iDofs
                        for jDof in jDofs
                            i = (iNode-1)*mdl.nDof + iDof
                            j = (jNode-1)*mdl.nDof + jDof
                            
                            inz_btm = be_len*(e-1) + (iNode-1)*mdl.nDof*size(be,2) + (jNode-1)*mdl.nDof^2 + (iDof-1)*mdl.nDof + jDof # index for the COO sparse matrix
                            inz_top = be_len*(mdl.ne)^(mdl.ndim-1)+ be_len*((e-1)) + (iNode-1)*mdl.nDof*size(be,2) + (jNode-1)*mdl.nDof^2 + (iDof-1)*mdl.nDof + jDof # index for the COO sparse matrix
                            
                            E[inz_top] = mdl.ID[iDof, mdl.IEN_top[iNode,e]]    # row index 
                            J[inz_top] = mdl.ID[jDof, mdl.IEN_top[jNode,e]]   # column index
                            V[inz_top] += w_top*be[i,j] 

                            E[inz_btm] = mdl.ID[iDof, mdl.IEN_btm[iNode,e]]    # row index 
                            J[inz_btm] = mdl.ID[jDof, mdl.IEN_btm[jNode,e]]   # column index
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

function set_slip_conditions_dense(mdl::model)

    IEN_btm_rows = size(mdl.IEN_btm,1)
    IEN_rows = size(mdl.IEN,1)
    ID_rows = size(mdl.ID,1)
    sz =    mdl.nDof*(2*mdl.ne+1)^mdl.ndim
    K = zeros(Float64, sz,sz)

    if mdl.ndim == 2
        # gaussian quadrature points for the element [-1,1] 
        ξ, w_ξ = gaussian_quadrature(-1,1)

        wpoints =  [w_ξ[1], w_ξ[2]]
        
        x = [ξ[1], ξ[2]]
    elseif mdl.ndim == 3            
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
    end 

    e_iter = 1:mdl.ne^(mdl.ndim-1)
    # integration loop
    gpiter = 1:length(wpoints)
    for gp in gpiter

        if mdl.ndim == 2
            N, ΔN = basis_function(x[gp], nothing, nothing, mdl.FunctionClass)
        elseif mdl.ndim == 3
            N, ΔN = basis_function(x[gp], y[gp], nothing, mdl.FunctionClass) 
        end

        for e in e_iter
        
            coords_top = mdl.NodeList[:,mdl.IEN_top[:,e]] # get the coordinates of the nodes of the element
            coords_btm = mdl.NodeList[:,mdl.IEN_btm[:,e]] # get the coordinates of the nodes of the element

            dxdξ_top = coords_top*ΔN         # Jacobian matrix [dx/dxi dx/deta; dy/dxi dy/deta; dz/dxi dz/deta]
            dxdξ_btm = coords_btm*ΔN         # Jacobian matrix [dx/dxi dx/deta; dy/dxi dy/deta; dz/dxi dz/deta]

            w_top = wpoints[gp]*norm(cross(dxdξ_top[:,1],dxdξ_top[:,2]))     # weight of the quadrature point top surface
            w_btm = wpoints[gp]*norm(cross(dxdξ_btm[:,1],dxdξ_btm[:,2]))     # weight of the quadrature point bottom surface
            
            M = zeros(3, mdl.ndim*length(N))
            M[1,1:mdl.nDof:end] = N
            M[2,2:mdl.nDof:end] = N
            M[3,3:mdl.nDof:end] = N

            be = M'*M
            be_cols, be_rows = size(be)
            be_len = length(be)

            # loop between basis functions of the element
            iNodes = 1:be_cols÷mdl.nDof
            jNodes = 1:be_rows÷mdl.nDof
            iDofs = 1:ID_rows
            jDofs = 1:ID_rows
            for iNode in iNodes
                for jNode in jNodes
                    for iDof in iDofs
                        for jDof in jDofs
                            i = (iNode-1)*mdl.nDof + iDof
                            j = (jNode-1)*mdl.nDof + jDof
                            
                            K[mdl.ID[iDof, mdl.IEN_top[iNode,e]],mdl.ID[jDof, mdl.IEN_top[jNode,e]]] += w_top*be[i,j] # top surface
                            K[mdl.ID[iDof, mdl.IEN_btm[iNode,e]],mdl.ID[jDof, mdl.IEN_btm[jNode,e]]] += w_btm*be[i,j] # bottom surface
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
function set_boundary_conditions(mdl::model)
    tDof = mdl.nDof*size(mdl.NodeList,2)     # Total number of DOFs
    q_upper = zeros(tDof,1)                  # initialize the vector of the Dirichlet boundary conditions (for mdl.nDof = 1) / Dirichlet boundary conditions upper surface (for mdl.nDof > 1)
    q_lower = zeros(tDof,1)                  # initialize the vector of the mdl.neumann boundary conditions (for mdl.nDof = 1) / Dirichlet boundary conditions lower surface (for mdl.nDof > 1)
    C = sparse(I,tDof,tDof)      # definition of the constraint matrix

    if mdl.nDof == 1
        if mdl.ndim == 3
            Dbound1 = 0
            Dbound2 = 1
            iter = 1:size(mdl.NodeList,2)
            for n in iter
                coord = mdl.NodeList[:,n] # get the coordinates of the node
                if coord[3] == Dbound1 # bottom boundary
                    q_upper[n] = 0
                elseif coord[3] == Dbound2 # top boundary
                    q_upper[n] = -1
                end
            end
        elseif mdl.ndim == 2
            Dbound1 = 0
            Dbound2 = 1
            iter = 1:size(mdl.NodeList,2)
            for n in iter
                coord = mdl.NodeList[:,n] # get the coordinates of the node
                if coord[2] == Dbound1 # bottom boundary
                    q_upper[n] = 0
                elseif coord[2] == Dbound2 # top boundary
                    q_upper[n] = -1
                end
            end
        end

        if mdl.FunctionClass == "Q1"
            C_uc = C[:,((mdl.ne+1)^(mdl.ndim-1)+1):((mdl.ne+1)^mdl.ndim-(mdl.ne+1)^(mdl.ndim-1))]
        elseif mdl.FunctionClass == "Q2"
            C_uc = C[:,((2*mdl.ne+1)^(mdl.ndim-1)+1):((2*mdl.ne+1)^mdl.ndim-(2*mdl.ne+1)^(mdl.ndim-1))]
        end

    elseif mdl.nDof == 3
        z0Bound = 0
        z1Bound = 1

        rCol = Array{Int}(undef,0)
        iter = 1:size(mdl.NodeList,2)
        for nNode in iter
            coord = mdl.NodeList[:,nNode]    # get the coordinates of the node
            if coord[3] == z0Bound       # bottom boundary
                q_lower[3*nNode] = 1     # constraint the z displacement to be zero at the bottom boundary
                push!(rCol,3*nNode)
            elseif coord[3] == z1Bound   # top boundary
                q_upper[3*nNode] = 1     # constraint the z displacement to be -d
                push!(rCol,3*nNode)
            end
        end

        C_uc = C[:,setdiff(1:size(C,2),rCol)]
    end
    return q_upper, q_lower, C_uc
end

function set_boundary_conditions_dense(mdl::model)
    tDof = mdl.nDof*size(mdl.NodeList,2)    # Total number od Degrees od freedom
    q_upper = zeros(tDof,1)                 # initialize the vector of the Dirichlet boundary conditions (for mdl.nDof = 1) / Dirichlet boundary conditions upper surface (for mdl.nDof > 1)
    q_lower = zeros(tDof,1)                 # initialize the vector of the mdl.neumann boundary conditions (for mdl.nDof = 1) / Dirichlet boundary conditions lower surface (for mdl.nDof > 1)
    C = Matrix(I,tDof,tDof)                 # Definition of the constraint matrix
    
    if mdl.nDof == 1
        if mdl.ndim == 3
            Dbound1 = 0
            Dbound2 = 1
            iter = 1:size(mdl.NodeList,2)
            for n in iter
                coord = mdl.NodeList[:,n] # get the coordinates of the node
                if coord[3] == Dbound1 # bottom boundary
                    q_upper[n] = 0
                elseif coord[3] == Dbound2 # top boundary
                    q_upper[n] = -1
                end
            end
        elseif mdl.ndim == 2
            Dbound1 = 0
            Dbound2 = 1
            iter = 1:size(mdl.NodeList,2)
            for n in iter
                coord = mdl.NodeList[:,n] # get the coordinates of the node
                if coord[2] == Dbound1 # bottom boundary
                    q_upper[n] = 0
                elseif coord[2] == Dbound2 # top boundary
                    q_upper[n] = -1
                end
            end
        end

        if mdl.FunctionClass == "Q1"
            C_uc = C[:,((mdl.ne+1)^(mdl.ndim-1)+1):((mdl.ne+1)^mdl.ndim-(mdl.ne+1)^(mdl.ndim-1))]
        elseif mdl.FunctionClass == "Q2"
            C_uc = C[:,((2*mdl.ne+1)^(mdl.ndim-1)+1):((2*mdl.ne+1)^mdl.ndim-(2*mdl.ne+1)^(mdl.ndim-1))]
        end

    else
        z0Bound = 0
        z1Bound = 1

        rCol = Array{Int}(undef,0)
        iter = 1:size(mdl.NodeList,2)
        for nNode in iter
            coord = mdl.NodeList[:,nNode]    # get the coordinates of the node
            if coord[3] == z0Bound       # bottom boundary
                q_lower[3*nNode] = 1     # constraint the z displacement to be zero at the bottom boundary
                push!(rCol,3*nNode)
            elseif coord[3] == z1Bound   # top boundary
                q_upper[3*nNode] = 1     # constraint the z displacement to be -d
                push!(rCol,3*nNode)
            end
        end

        C_uc = C[:,setdiff(1:size(C,2),rCol)]
    end
    return q_upper, q_lower, C_uc
end

""" 
    get_cMat(type; λ=nothing, μ=nothing, Young=nothing, ν=nothing)

Returns the stiffness matrix for a given type of material

# Arguments:
- `type::String`: type of constitutive matrix to be considered (lame or standard)
- `λ::Float64`: Lame's first parameter
- `μ::Float64`: Sheer modulus
- `Young::Float64`: Young's modulus
- `ν::Float64`: Poisson's ratio

# Returns:
- `cMat::Matrix{Float64}`: constitutive matrix
"""
function get_cMat(type="standard", c1=nothing, c2=nothing)

    if type == "lame"
        cMat =  [[ 2*c2+c1  c1    c1    0  0  0]; 
                        [  c1   2*c2+c1  c1    0  0  0]; 
                        [  c1     c1   2*c2+c1 0  0  0]; 
                        [  0      0      0     c2 0  0]; 
                        [  0      0      0     0 c2  0]
                    [  0      0      0     0  0 c2]]  # constitutive matrix
    elseif type == "standard"
        cMat =  (c1/((1+c2)*(1-2*c2)))*[[1-c2 c2 c2   0      0      0]; 
                                                [c2 1-c2 c2   0      0      0]; 
                                                [c2 c2 1-c2   0      0      0]; 
                                                [0 0  0 (1-2*c2)/2 0      0]; 
                                                [0 0  0    0   (1-2*c2)/2 0]; 
                                                [0 0  0    0      0   (1-2*c2)/2]]  # constitutive matrix   
    end
    return cMat
end

function get_volume(mdl::model)

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

    vol = 0
    e_iter = 1:mdl.ne^mdl.ndim
    # integration loop
    gpiter = 1:length(wpoints)
    for gp in gpiter
        if mdl.ndim == 1
            N, ΔN = basis_function(x[gp], nothing, nothing, mdl.FunctionClass)
        elseif mdl.ndim == 2
            N, ΔN = basis_function(x[gp], y[gp], nothing, mdl.FunctionClass) 
        elseif mdl.ndim == 3
            N, ΔN = basis_function(x[gp], y[gp], z[gp], mdl.FunctionClass) 
        end
        # element loop
        for e in e_iter
            coords = mdl.NodeList[:,mdl.IEN[:,e]] # get the coordinates of the nodes of the element

            Jac  = coords*ΔN # Jacobian matrix [dx/dxi dx/deta; dy/dxi dy/deta]

            w = wpoints[gp]*abs(det(Jac))
            vol += w
        end
    end
end

"""
    simulate(x0, x1, y0, y1, z0, z1, ne, Young, ν, ndim, FunctionClass, nDof, β, CameraMatrix, endTime, tSteps, Control, cParam, cMat; writeData=false, filepath=nothing)

Simulate the deformation of a cylindrical under compression

# Arguments:
- `x0::Float64` : x-coordinate of the lower left corner of the domain
- `x1::Float64` : x-coordinate of the upper right corner of the domain
- `y0::Float64` : y-coordinate of the lower left corner of the domain
- `y1::Float64` : y-coordinate of the upper right corner of the domain
- `z0::Float64` : z-coordinate of the lower left corner of the domain
- `z1::Float64` : z-coordinate of the upper right corner of the domain
- `ne::Int` : number of elements
- `Young::Float64` : Young's modulus
- `ν::Float64` : Poisson's ratio
- `ndim::Int` : number of dimensions
- `FunctionClass::String` : type of basis function
- `nDof::Int` : number of degree of freedom per node
- `β::Float64` : friction parameter
- `CameraMatrix::Matrix{Float64}` : camera matrix
- `endTime::Float64` : end time
- `tSteps::Int` : number of time steps to be taken
- `Control::String` : type of control (force or displacement)
- `cParam::Vector{Float64}` : control parameter (force or displacement prescribed at the top surface per time step)
- `cMat::Matrix{Float64}` : control matrix
- `writeData::Bool` : write the data to a file
- `filepath::String` : path to the file
"""
function simulate(x0, x1, y0, y1, z0, z1, ne, Young, ν, ndim, FunctionClass, nDof, β, CameraMatrix, endTime, tSteps, Control, cParam, cMat; writeData=false, 
                    filepath=nothing, SIDES::Bool=false, dcdλ::Matrix{Float64}=zeros(Float64,1,1))

    time = collect(range(start=0,stop=endTime,length=tSteps)) # time vector

    NodeList, IEN, ID, IEN_top, IEN_btm, BorderNodesList = meshgrid_cube(x0,x1,y0,y1,z0,z1,ne,ndim,FunctionClass=FunctionClass)  # generate the mesh grid
    NodeListCylinder = inflate_cylinder(NodeList, x0, x1, y0, y1)                                 # inflate the sphere to a unit sphere


    mdl = def_model("linear_elasticity", ne=ne, NodeList=NodeListCylinder, IEN=IEN, IEN_top=IEN_top, IEN_btm=IEN_btm, ndim=ndim, nDof=nDof, ID = ID,
                     FunctionClass=FunctionClass, Young=Float64(Young), ν=ν, cMat=cMat,dcMatdλ=dcdλ)
    
    q_tp, q_btm, C_uc = set_boundary_conditions(mdl)

    BorderPts2D, BorderNodes2D, Nodes2D = extract_borders(NodeListCylinder, CameraMatrix, BorderNodesList, ne, 2*ne+1, SIDES)
    pi, qi = fit_curve(border=BorderPts2D)

    SideBorders = BorderNodesList[1]
    BottomBorders = BorderNodesList[2]
    TopBorders = BorderNodesList[3]
        
    fields = AbstractArray[]  
    dfields = AbstractArray[]  
    pos3D = AbstractArray[NodeListCylinder]                                                             # store the solution fields of the mesh in 3D
    pos2D = AbstractArray[Nodes2D]                                                                   # store the solution fields of the mesh in 2D
    surfaceNodesList = Float64[NodeList[:,SideBorders] NodeList[:,BottomBorders] NodeList[:,TopBorders]]  # store the solution fields of the surfaces in 3D
    borderPts2DList = AbstractArray[BorderPts2D]                                                               # store the solution fields of the surfaces in 2D
    borderNodeList2D = AbstractArray[BorderNodes2D]                                                       # store the solution fields of the border nodes in 2D
    splinep = AbstractArray[pi]                                                                            # store the x coordinates samples of the spline parameters of the border nodes
    splineq = AbstractArray[qi]                                                                            # store the y coordinates samples of the spline parameters of the border nodes
    output = Float64[] 
    writeborderList = [vcat(pi', qi')]

    # K, dKdλ = assemble_system(mdl,GRAD=true)             # assemble the stiffness matrix
    K = assemble_system(mdl)  
    b = set_slip_conditions(mdl)   # apply the neumann boundary conditions

    μ_btm = 0      
    iter = 1
    pr = Progress(tSteps; desc= "Simulating with prescribed $Control ...", showspeed=true)
    if Control == "force"
        for t in time

            q_btm = μ_btm*q_btm
            
            K_bar = K + β*b

            C_T = transpose(C_uc)              # transpose the constraint matrix

            M = [C_T*K_bar*C_uc C_T*K_bar*q_tp; q_tp'*K_bar*C_uc q_tp'*K_bar*q_tp] # assemble the system of equations]

            sol = M\([-C_T*K_bar*q_btm; cParam[iter].-q_tp'*K_bar*q_btm]) # solve the system of equations

            q_f = sol[1:end-1]        # solve the system of equations
            μ_tp = sol[end]           # solve the system of equations
            q = q_btm + C_uc*q_f + μ_tp*q_tp;                # assemble the solution
            # dqdλ = -K_free\(dKdλ*q)

            # post process the solution
            motion = [q[ID[1,:]] q[ID[2,:]] q[ID[3,:]]]'    # update the nodal positions
            # dmotion = [dqdλ[ID[1,:]] dqdλ[ID[2,:]] dqdλ[ID[3,:]]]'
            NodeListCylinder = NodeListCylinder + motion    # update the node coordinates

            BorderPts2D, BorderNodes2D, Nodes2D = extract_borders(NodeListCylinder, CameraMatrix, BorderNodesList, SIDES)
            surfaceNodesList = [NodeListCylinder[:,SideBorders] NodeListCylinder[:,BottomBorders] NodeListCylinder[:,TopBorders]]
            pi, qi = fit_curve(border=BorderPts2D)

            push!(output, μ_tp)
            push!(fields, motion)
            # push!(dfields,dmotion)
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
    elseif Control == "displacement"
        pr = Progress(tSteps; desc="Simulating with prescribed $Control ...", showspeed=true)
        for t in time
        
            q_d = (μ_btm*q_btm + cParam[iter]*q_tp)                  # apply the Dirichlet boundary conditions

            K_bar = K + β*b

            C_T = transpose(C_uc)           # transpose the constraint matrix
            K_free = C_T*K_bar*C_uc         # extract the free part of the stiffness matrix

            q_f = K_free\(C_T*(-K_bar*q_d))         # solve the system of equations
            q = q_d + C_uc*q_f;                 # assemble the solution 

            # dqdλ = -K_bar\(dKdλ*q)

            # post process the solution
            f_R = K_bar*q
            motion = [q[ID[1,:]] q[ID[2,:]] q[ID[3,:]]]'
            # dmotion = [dqdλ[ID[1,:]] dqdλ[ID[2,:]] dqdλ[ID[3,:]]]'
            F_est = q_tp'*f_R                                         # calculate the reaction force at the top surface F = Σf^{tp}_{iR} = q_tp'*f_R
            NodeListCylinder = NodeListCylinder + motion            # update the node coordinates

            BorderPts2D, BorderNodes2D, Nodes2D = extract_borders(NodeListCylinder, CameraMatrix, BorderNodesList, SIDES)
            surfaceNodesList = [NodeListCylinder[:,SideBorders] NodeListCylinder[:,BottomBorders] NodeListCylinder[:,TopBorders]]
            pi, qi = fit_curve(border=BorderPts2D)

            # store the solutions in a list
            push!(output, F_est[1])
            push!(fields, motion)
            # push!(dfields,dmotion)
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
    end

    if writeData
        write_scene(string(filepath,"/Results"), NodeListCylinder, IEN, ne, ndim, pos3D, ID=ID, FunctionClass=FunctionClass)
        animate_fields(filepath = string(filepath,"/Results/images"), fields=pos3D , IEN=IEN, BorderNodes2D=borderPts2DList, fields2D=pos2D)
        write_contour_data(string(filepath,"/Results"), writeborderList)
    end

    return output, borderPts2DList, borderNodeList2D, splinep, splineq, mdl
end