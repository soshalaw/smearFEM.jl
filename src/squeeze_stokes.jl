using LinearAlgebra
using ProgressMeter
using SparseArrays

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
- `ID::Matrix{Int}{nNodes,nDof}` : matrix that maps the global degrees of freedom to the local degrees of freedom
- `Young::Float64`: Young's modulus
- `ν::Float64`: Poisson's ratio

# Returns:
- `K::SparseMatrixCSC{Float64,Int64}{ndof,ndof}` : sparse stiffness matrix 
"""
function assemble_system_A(mdl::model)

    # (I,J,V) vectors for COO sparse matrix
    if mdl.nDof == 1
        E = zeros(  Int64, mdl.ne^mdl.ndim*size(mdl.IEN_u,2)^2)
        J = zeros(  Int64, mdl.ne^mdl.ndim*size(mdl.IEN_u,2)^2)
        V = zeros(Float64, mdl.ne^mdl.ndim*size(mdl.IEN_u,2)^2)
    else
        E = zeros(  Int64, mdl.ne^mdl.ndim*((size(mdl.ID,2)*size(mdl.IEN_u,2))^2))
        J = zeros(  Int64, mdl.ne^mdl.ndim*((size(mdl.ID,2)*size(mdl.IEN_u,2))^2))
        V = zeros(Float64, mdl.ne^mdl.ndim*((size(mdl.ID,2)*size(mdl.IEN_u,2))^2))  
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

        x = []
        y = []
        wpoints = []
        
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

        x = []
        y = []
        z = []
        wpoints = []
        
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

    for e in 1:mdl.ne^mdl.ndim
        coords = mdl.NodeList[:,mdl.IEN_u[e,:]] # get the coordinates of the nodes of the element
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
                        E[inz] = mdl.IEN_u[e,i] # row index 
                        J[inz] = mdl.IEN_u[e,j] # column index
                        V[inz] += w*dot(dNdX[i,:],dNdX[j,:])# inner product of the gradient of the basis functions
                    end
                end
            else   

                divN = sum(ΔN, dims=2)

                B = zeros(3,3*length(divN))
                B[1,1:mdl.nDof:end] = divN
                B[2,2:mdl.nDof:end] = divN
                B[3,3:mdl.nDof:end] = divN
                Ke = w*B'*B # element stiffness matrix
        
                # loop between basis functions of the element
                iNodes = 1:size(Ke,1)÷mdl.nDof
                jNodes = 1:size(Ke,2)÷mdl.nDof
                iDofs = 1:size(mdl.ID,2)
                jDofs = 1:size(mdl.ID,2)
                for iNode in iNodes
                    for jNode in jNodes
                        for iDof in iDofs
                            for jDof in jDofs
                                i = (iNode-1)*mdl.nDof + iDof
                                j = (jNode-1)*mdl.nDof + jDof
                                inz = length(Ke)*(e-1) + (iNode-1)*mdl.nDof*size(Ke,2) + (jNode-1)*mdl.nDof^2 + (iDof-1)*mdl.nDof + jDof # index for the COO sparse matrix
                                E[inz] = mdl.ID[mdl.IEN_u[e,iNode],iDof] # row index 
                                J[inz] = mdl.ID[mdl.IEN_u[e,jNode], jDof] # column index
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

function assemble_system_B(mdl::model)

    # (I,J,V) vectors for COO sparse matrix
    if mdl.nDof == 1
        E = zeros(  Int64, mdl.ne^mdl.ndim*size(mdl.IEN_u,2)*size(mdl.IEN_p,2))
        J = zeros(  Int64, mdl.ne^mdl.ndim*size(mdl.IEN_u,2)*size(mdl.IEN_p,2))
        V = zeros(Float64, mdl.ne^mdl.ndim*size(mdl.IEN_u,2)*size(mdl.IEN_p,2))
    else
        E = zeros(  Int64, mdl.ne^mdl.ndim*(size(mdl.ID,2)*size(mdl.IEN_u,2))*size(mdl.IEN_p,2))
        J = zeros(  Int64, mdl.ne^mdl.ndim*(size(mdl.ID,2)*size(mdl.IEN_u,2))*size(mdl.IEN_p,2))
        V = zeros(Float64, mdl.ne^mdl.ndim*(size(mdl.ID,2)*size(mdl.IEN_u,2))*size(mdl.IEN_p,2))  
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

        x = []
        y = []
        wpoints = []
        
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

        x = []
        y = []
        z = []
        wpoints = []
        
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

    for e in 1:mdl.ne^mdl.ndim
        coords_u = mdl.NodeList[:,mdl.IEN_u[e,:]] # get the coordinates of the nodes of the element
        coords_p = mdl.NodeList[:,mdl.IEN_p[e,:]] # get the coordinates of the nodes of the element
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

            Jac  = coords_u*ΔN_u # Jacobian matrix [dx/dxi dx/deta; dy/dxi dy/deta]

            w = wpoints[gp]*abs(det(Jac))
            invJ = inv(Jac)
            dNdX = ΔN_u*invJ

            if mdl.nDof == 1
                szN = size(N_u,1) # number of basis functions
                # loop between basis functions of the element
                for i in 1:szN
                    for j in 1:szN
                        inz = (szN)^2*(e-1) + szN*(i-1) + j # index for the COO sparse matrix
                        E[inz] = mdl.IEN[e,i] # row index 
                        J[inz] = mdl.IEN[e,j] # column index
                        V[inz] += w*dot(dNdX[i,:],dNdX[j,:])# inner product of the gradient of the basis functions
                    end
                end
            else   
                divN = sum(dNdX, dims=2)

                B = zeros(3,3*length(divN))
                B[1,1:mdl.nDof:end] = divN
                B[2,2:mdl.nDof:end] = divN
                B[3,3:mdl.nDof:end] = divN

                Ke = [N_p N_p N_p]*B # element stiffness matrix

                # loop between basis functions of the element
                jNodes = 1:size(Ke,1) # row index
                iNodes = 1:size(Ke,2)÷mdl.nDof  # column node index
                iDofs = 1:size(mdl.ID,2) # column dof index

                for iNode in iNodes 
                    for jNode in jNodes 
                        for iDof in iDofs
                            i = (iNode-1)*mdl.nDof + iDof
                            j = jNode
                            inz = length(Ke)*(e-1) + (jNode-1)*size(Ke,2) + (iNode-1)*mdl.nDof + iDof # index for the COO sparse matrix
                            E[inz] = mdl.ID[mdl.IEN_u[e,iNode],iDof] # row index 
                            J[inz] = mdl.IEN_p[e,jNode] # column index
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
function apply_boundary_conditions_stokes(mdl)

    E = zeros(  Int64, mdl.ne^(mdl.ndim-1)*(size(mdl.ID,2)*size(mdl.IEN_u_btm,2))^2*2) # *2 because we have two surfaces
    J = zeros(  Int64, mdl.ne^(mdl.ndim-1)*(size(mdl.ID,2)*size(mdl.IEN_u_btm,2))^2*2) # *2 because we have two surfaces
    V = zeros(Float64, mdl.ne^(mdl.ndim-1)*(size(mdl.ID,2)*size(mdl.IEN_u_btm,2))^2*2) # *2 because we have two surfaces

    for e in 1:mdl.ne^(mdl.ndim-1)
    
        coords_u_top = mdl.NodeList[:,mdl.IEN_u_top[e,:]] # get the coordinates of the nodes of the element
        coords_u_btm = mdl.NodeList[:,mdl.IEN_u_btm[e,:]] # get the coordinates of the nodes of the element

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

        for gp in 1:2^(mdl.ndim-1)

            if mdl.ndim == 2
                N, ΔN = basis_function(x[gp], nothing, nothing, mdl.FunctionClass)
            elseif mdl.ndim == 3
                N, ΔN = basis_function(x[gp], y[gp], nothing, mdl.FunctionClass) 
            end

            dxdξ_top = coords_u_top*ΔN         # Jacobian matrix [dx/dxi dx/deta; dy/dxi dy/deta; dz/dxi dz/deta]
            dxdξ_btm = coords_u_btm*ΔN         # Jacobian matrix [dx/dxi dx/deta; dy/dxi dy/deta; dz/dxi dz/deta]

            w_top = wpoints[gp]*norm(cross(dxdξ_top[:,1],dxdξ_top[:,2]))     # weight of the quadrature point top surface
            w_btm = wpoints[gp]*norm(cross(dxdξ_btm[:,1],dxdξ_btm[:,2]))     # weight of the quadrature point bottom surface
            
            M = zeros(3, mdl.ndim*length(N))
            M[1,1:mdl.nDof:end] = N
            M[2,2:mdl.nDof:end] = N
            M[3,3:mdl.nDof:end] = N

            be = M'*M

            # loop between basis functions of the element
            iNodes = 1:size(be,1)÷mdl.nDof
            jNodes = 1:size(be,2)÷mdl.nDof
            iDofs = 1:size(mdl.ID,2)
            jDofs = 1:size(mdl.ID,2)
            for iNode in iNodes
                for jNode in jNodes
                    for iDof in iDofs
                        for jDof in jDofs
                            i = (iNode-1)*mdl.nDof + iDof
                            j = (jNode-1)*mdl.nDof + jDof
                            
                            inz_btm = length(be)*(e-1) + (iNode-1)*mdl.nDof*size(be,2) + (jNode-1)*mdl.nDof^2 + (iDof-1)*mdl.nDof + jDof # index for the COO sparse matrix
                            inz_top = length(be)*(mdl.ne)^(mdl.ndim-1)+ length(be)*((e-1)) + (iNode-1)*mdl.nDof*size(be,2) + (jNode-1)*mdl.nDof^2 + (iDof-1)*mdl.nDof + jDof # index for the COO sparse matrix
                            
                            E[inz_top] = mdl.ID[mdl.IEN_u_top[e,iNode],iDof]    # row index 
                            J[inz_top] = mdl.ID[mdl.IEN_u_top[e,jNode], jDof]   # column index
                            V[inz_top] += w_top*be[i,j] 

                            E[inz_btm] = mdl.ID[mdl.IEN_u_btm[e,iNode],iDof]    # row index 
                            J[inz_btm] = mdl.ID[mdl.IEN_u_btm[e,jNode], jDof]   # column index
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
        C = sparse(I,nDof*(ne+1)^ndim,nDof*(ne+1)^ndim)      # definition of the constraint matrix
    elseif FunctionClass == "Q2"
        q_upper = zeros(nDof*(2*ne+1)^ndim,1)                # initialize the vector of the Dirichlet boundary conditions (for ndof = 1) / Dirichlet boundary conditions upper surface (for ndof > 1)
        q_lower = zeros(nDof*(2*ne+1)^ndim,1)                # initialize the vector of the Neumann boundary conditions (for ndof = 1) / Dirichlet boundary conditions lower surface (for ndof > 1)
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