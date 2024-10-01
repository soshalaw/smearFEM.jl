using LinearAlgebra
using SparseArrays
using ArgCheck

function greet_fem()
    println("Hello, I am the FEM module")
end

mutable struct model
    ne::Int64
    NodeList::Matrix{Float64}
    IEN::Matrix{Int}
    IEN_top::Matrix{Int}
    IEN_btm::Matrix{Int}
    IEN_border::Matrix{Int}
    ndim::Int64
    nDof::Int64
    FunctionClass::String
    ID::Matrix{Int}
    Young::Float64
    ν::Float64
    cMat::Matrix{Float64}
    dcMatdλ::Matrix{Float64}
    dcMatdμ::Matrix{Float64}
end


function def_model(; ne::Int64 = 1, 
    NodeList::Matrix{Float64} = [0.0 1.0 1.1 0.1], 
    IEN::Matrix{Int} = [1 2 3 4], 
    ndim::Int64 = 2, 
    nDof::Int64 = 1, 
    FunctionClass::String = "Q1", 
    Young = 1.0, 
    ν = 0.3, 
    ID::Matrix{Int} = [1 2 3 4],
    cMat::Matrix{Float64} = [0.0 1.0 1.1 0.1],
    IEN_top::Matrix{Int} = [1 2 3 4],
    IEN_btm::Matrix{Int} = [1 2 3 4],
    IEN_border::Matrix{Int} = [1 2 3 4],
    dcMatdλ::Matrix{Float64} = [0.0 1.0 1.1 0.1],
    dcMatdμ::Matrix{Float64} = [0.0 1.0 1.1 0.1]
    )

    model(ne, NodeList, IEN, IEN_top, IEN_btm, IEN_border, ndim, nDof, FunctionClass, ID, Young, ν, cMat, dcMatdλ, dcMatdμ)
end

params(mdl::model) = "ne = $(mdl.ne), ndim = $(mdl.ndim), nDof = $(mdl.nDof), FunctionClass = $(mdl.FunctionClass), Young = $(mdl.Young), ν = $(mdl.ν)"

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
        cMat = [[ 2*c2+c1  c1    c1    0 0 0]; 
                [  c1   2*c2+c1  c1    0 0 0]; 
                [  c1     c1   2*c2+c1 0 0 0]; 
                [  0     0    0    c2 0 0]; 
                [  0     0    0    0 c2 0]
                [  0     0    0    0 0 c2]]  # constitutive matrix
    elseif type == "standard"
        cMat = (c1/((1+c2)*(1-2*c2)))*[[1-c2 c2 c2   0      0      0]; 
                                        [c2 1-c2 c2   0      0      0]; 
                                        [c2 c2 1-c2   0      0      0]; 
                                        [0 0  0 (1-2*c2)/2 0      0]; 
                                        [0 0  0    0   (1-2*c2)/2 0]; 
                                        [0 0  0    0      0   (1-2*c2)/2]]  # constitutive matrix   
    end
    return cMat
end

""" 
    smearFEM.gaussian_quadrature(a,b,nGaussPoints)

Compute the nodes and weights for the Gaussian quadrature of order 2
    
# Arguments:    
- `a,b::Integer` : the limits of the integration interval
- `nGaussPoints::Integer` : number of Gauss points to be considered (2 or 3)

# Returns:    
- `ξ::Vector{Float64}{,nGaussPoints}`: nodes.
- `w::Vector{Float64}{,nGaussPoints}`: weights 
"""
function gaussian_quadrature(a,b,nGaussPoints=2)
  
    if nGaussPoints == 2
        ξ = [-(b-a)/(2*sqrt(3))+(b+a)/2, (b-a)/(2*sqrt(3))+(b+a)/2]
        w = [(b-a)/2, (b-a)/2]
    elseif nGaussPoints == 3
        ξ = [-(b-a)/(2*sqrt(5/3))+(b+a)/2, 0, (b-a)/(2*sqrt(5/3))+(b+a)/2]
        w = [(b-a)/2*5/9, (b-a)/2*8/9, (b-a)/2*5/9]
    end
    return ξ, w
end

""" 
    basis_function(ξ,η=nothing,ζ=nothing,FunctionClass = "Q1")

Define the basis functions and the gradients for a master element

# Arguments:
- `ξ::Float64`: ξ coordinate of the point where the basis function is evaluated
- `η::Float64`: η coordinate of the point where the basis function is evaluated
- `ζ::Float64`: ζ coordinate of the point where the basis function is evaluated
- `FunctionClass::String`: type of basis functions to be considered (Q1:quadratic or Q2:Lagrange)

# Returns:
- `N::Vector{Float64}{,ndof}`: basis functions
- `Delta_N::Matrix{Float64}{ndof,ndim}`: gradient of the basis functions 
"""
function basis_function(ξ,η=nothing,ζ=nothing,FunctionClass = "Q1")
    
    if FunctionClass == "Q1"
        if !isnothing(ζ) # Considering a 3D master element
            # basis functions
            N = [(1-ξ)*(1-η)*(1-ζ)/8, 
                (1+ξ)*(1-η)*(1-ζ)/8, 
                (1+ξ)*(1+η)*(1-ζ)/8, 
                (1-ξ)*(1+η)*(1-ζ)/8, 
                (1-ξ)*(1-η)*(1+ζ)/8, 
                (1+ξ)*(1-η)*(1+ζ)/8, 
                (1+ξ)*(1+η)*(1+ζ)/8, 
                (1-ξ)*(1+η)*(1+ζ)/8]

            # gradient of the basis functions
            ΔN = [[-(1-η)*(1-ζ)/8, (1-η)*(1-ζ)/8, (1+η)*(1-ζ)/8, -(1+η)*(1-ζ)/8, -(1-η)*(1+ζ)/8, (1-η)*(1+ζ)/8, (1+η)*(1+ζ)/8, -(1+η)*(1+ζ)/8] [-(1-ξ)*(1-ζ)/8, -(1+ξ)*(1-ζ)/8, (1+ξ)*(1-ζ)/8, (1-ξ)*(1-ζ)/8, -(1-ξ)*(1+ζ)/8, -(1+ξ)*(1+ζ)/8, (1+ξ)*(1+ζ)/8, (1-ξ)*(1+ζ)/8] [-(1-ξ)*(1-η)/8, -(1+ξ)*(1-η)/8, -(1+ξ)*(1+η)/8, -(1-ξ)*(1+η)/8, (1-ξ)*(1-η)/8, (1+ξ)*(1-η)/8, (1+ξ)*(1+η)/8, (1-ξ)*(1+η)/8]]  # [dN/dxi dN/deta dN/dzeta]
        elseif !isnothing(η) # Considering a 2D master element
            # basis functions
            N = [(1-ξ)*(1-η)/4, (ξ+1)*(1-η)/4, (1+ξ)*(η+1)/4, (1-ξ)*(1+η)/4]

            # gradient of the basis functions
            ΔN = [[-(1-η)/4 , (1-η)/4, (η+1)/4, -(1+η)/4] [-(1-ξ)/4, -(ξ+1)/4, (1+ξ)/4, (1-ξ)/4] ]
        else # Considering a 1D master element
            # basis functions
            N  = [0.5-0.5*ξ, 0.5+0.5*ξ]

            # gradient of the basis functions
            ΔN = [-0.5 0.5]
        end
    elseif FunctionClass == "Q2"
        if !isnothing(η) # Considering a 2D master element
            # basis functions
            N = [(1-ξ)*ξ*(1-η)*η/4, 
                -ξ*(1+ξ)*(1-η)*η/4, 
                ξ*(1+ξ)*η*(1+η)/4, 
                -(1-ξ)*ξ*η*(1+η)/4,
                -(1-ξ)*(1+ξ)*(1-η)*η/2,
                ξ*(1+ξ)*(1-η)*(1+η)/2,
                (1-ξ)*(1+ξ)*η*(1+η)/2,
                -(1-ξ)*ξ*(1-η)*(1+η)/2,
                (1-ξ)*(1+ξ)*(1-η)*(1+η)] # gradient of the basis functions

            Delta_Nx = [(1-2*ξ)*(1-η)*η/4, 
                -(1+2*ξ)*(1-η)*η/4, 
                (1+2*ξ)*η*(1+η)/4, 
                -(1-2*ξ)*η*(1+η)/4,
                ξ*(1-η)*η,
                (1+2*ξ)*(1-η)*(1+η)/2,
                -ξ*η*(1+η),
                -(1-2*ξ)*(1-η)*(1+η)/2,
                -2*ξ*(1-η)*(1+η)] # dN/dξ gradient of the basis functions
            
            Delta_Ny = [(1-ξ)*ξ*(1-2*η)/4, 
                -ξ*(1+ξ)*(1-2*η)/4, 
                ξ*(1+ξ)*(1+2*η)/4, 
                -(1-ξ)*ξ*(1+2*η)/4,
                -(1-ξ)*(1+ξ)*(1-2*η)/2,
                -ξ*(1+ξ)*η,
                (1-ξ)*(1+ξ)*(1+2*η)/2,
                (1-ξ)*ξ*η,
                -(1-ξ)*(1+ξ)*2*η] # dN/dη gradient of the basis functions
        
            ΔN = [Delta_Nx Delta_Ny] # [dN/dξ dN/dη]
        end
    end
    return N, ΔN
end

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
function assemble_system(mdl::model)

    # (I,J,V) vectors for COO sparse matrix
    if mdl.nDof == 1
        E = zeros(  Int64, mdl.ne^mdl.ndim*size(mdl.IEN,2)^2)
        J = zeros(  Int64, mdl.ne^mdl.ndim*size(mdl.IEN,2)^2)
        V = zeros(Float64, mdl.ne^mdl.ndim*size(mdl.IEN,2)^2)
    else
        E = zeros(  Int64, mdl.ne^mdl.ndim*((size(mdl.ID,2)*2^mdl.ndim)^2))
        J = zeros(  Int64, mdl.ne^mdl.ndim*((size(mdl.ID,2)*2^mdl.ndim)^2))
        V = zeros(Float64, mdl.ne^mdl.ndim*((size(mdl.ID,2)*2^mdl.ndim)^2))  
    end

    # element loop
    if mdl.ndim == 1
        # gaussian quadrature points for the element [-1,1] 
        ξ, w_ξ = gaussian_quadrature(-1,1)
        
        wpoints = [w_ξ[1], w_ξ[2]]
        
        x = [ξ[1], ξ[2]]
    elseif mdl.ndim == 2
        # gaussian quadrature points for the element [-1,1]x[-1,1] 
        ξ, w_ξ = gaussian_quadrature(-1,1)
        η, w_η = gaussian_quadrature(-1,1)
        
        wpoints = [w_ξ[1]*w_η[1], w_ξ[2]*w_η[1], w_ξ[2]*w_η[2], w_ξ[1]*w_η[2]]
        
        x = [ξ[1], ξ[2], ξ[2], ξ[1]]
        y = [η[1], η[1], η[2], η[2]]

    elseif mdl.ndim == 3
        # gaussian quadrature points for the element [-1,1]x[-1,1]x[-1,1] 
        ξ, w_ξ = gaussian_quadrature(-1,1)
        η, w_η = gaussian_quadrature(-1,1)
        ζ, w_ζ = gaussian_quadrature(-1,1)
        
        wpoints = [w_ξ[1]*w_η[1]*w_ζ[1], w_ξ[2]*w_η[1]*w_ζ[1], w_ξ[2]*w_η[2]*w_ζ[1], w_ξ[1]*w_η[2]*w_ζ[1], w_ξ[1]*w_η[1]*w_ζ[2], w_ξ[2]*w_η[1]*w_ζ[2], w_ξ[2]*w_η[2]*w_ζ[2], w_ξ[1]*w_η[2]*w_ζ[2]]
        
        x = [ξ[1], ξ[2], ξ[2], ξ[1], ξ[1], ξ[2], ξ[2], ξ[1]]
        y = [η[1], η[1], η[2], η[2], η[1], η[1], η[2], η[2]]
        z = [ζ[1], ζ[1], ζ[1], ζ[1], ζ[2], ζ[2], ζ[2], ζ[2]]
    end

    for e in 1:mdl.ne^mdl.ndim
        coords = mdl.NodeList[:,mdl.IEN[e,:]] # get the coordinates of the nodes of the element
        
        # integration loop
        for gp in 1:2^mdl.ndim
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
                        E[inz] = mdl.IEN[e,i] # row index 
                        J[inz] = mdl.IEN[e,j] # column index
                        V[inz] += w*dot(dNdX[i,:],dNdX[j,:])# inner product of the gradient of the basis functions
                    end
                end
            else   
                @argcheck !isnothing(mdl.cMat) "constitutive matrix must be provided for problems with nDof > 1"

                if mdl.nDof == 2
                    B = zeros(3, mdl.ndim*length(N))
                    B[1,1:mdl.nDof:end] = dNdX[:,1]
                    B[2,2:mdl.nDof:end] = dNdX[:,2]
                    B[3,1:mdl.nDof:end] = dNdX[:,2]
                    B[3,2:mdl.nDof:end] = dNdX[:,1]

                    @argcheck size(mdl.cMat) == (2,2) "constitutive matrix must be 2x2 for plane stress problems"
                elseif mdl.nDof == 3
                    B = zeros(6, mdl.ndim*length(N))
                    B[1,1:mdl.nDof:end] = dNdX[:,1]
                    B[2,2:mdl.nDof:end] = dNdX[:,2]
                    B[3,3:mdl.nDof:end] = dNdX[:,3]
                    B[4,2:mdl.nDof:end] = dNdX[:,3]
                    B[4,3:mdl.nDof:end] = dNdX[:,2]
                    B[5,1:mdl.nDof:end] = dNdX[:,3]
                    B[5,3:mdl.nDof:end] = dNdX[:,1]
                    B[6,1:mdl.nDof:end] = dNdX[:,2]
                    B[6,2:mdl.nDof:end] = dNdX[:,1]

                    @argcheck size(mdl.cMat) == (6,6) "constitutive matrix must be 6x6 for 3D problems"
                end

                Ke = B'*mdl.cMat*B*w # element stiffness matrix
        
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
                                E[inz] = mdl.ID[mdl.IEN[e,iNode],iDof] # row index 
                                J[inz] = mdl.ID[mdl.IEN[e,jNode], jDof] # column index
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

function apply_boundary_conditions(mdl)
    """ Apply the Neumann slip boundary conditions to the global stiffness matrix

        Parameters:
        K: {[ndof,ndof] SparseMatrixCSC{Float64,Int64}} : sparse stiffness matrix 
        ID: {[nNodes,nDof] Matrix{Int}} : matrix that maps the global degrees of freedom to the local degrees of freedom
        q_d: {[ndof] Vector{Float64}} : Dirichlet boundary conditions
        q_n: {[ndof] Vector{Float64}} : Neumann boundary conditions

        Returns:
        K: {[ndof,ndof] SparseMatrixCSC{Float64,Int64}} : sparse stiffness matrix with the boundary conditions applied
        F: {[ndof] Vector{Float64}} : force vector
    """

    E = zeros(  Int64, mdl.ne^(mdl.ndim-1)*(size(mdl.ID,2)*size(mdl.IEN_btm,2))^2*2) # *2 because we have two surfaces
    J = zeros(  Int64, mdl.ne^(mdl.ndim-1)*(size(mdl.ID,2)*size(mdl.IEN_btm,2))^2*2) # *2 because we have two surfaces
    V = zeros(Float64, mdl.ne^(mdl.ndim-1)*(size(mdl.ID,2)*size(mdl.IEN_btm,2))^2*2) # *2 because we have two surfaces

    for e in 1:mdl.ne^(mdl.ndim-1)
    
        coords_top = mdl.NodeList[:,mdl.IEN_top[e,:]] # get the coordinates of the nodes of the element
        coords_btm = mdl.NodeList[:,mdl.IEN_btm[e,:]] # get the coordinates of the nodes of the element

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

            dxdξ_top = coords_top*ΔN         # Jacobian matrix [dx/dxi dx/deta; dy/dxi dy/deta; dz/dxi dz/deta]
            dxdξ_btm = coords_btm*ΔN         # Jacobian matrix [dx/dxi dx/deta; dy/dxi dy/deta; dz/dxi dz/deta]

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
                            
                            E[inz_top] = mdl.ID[mdl.IEN_top[e,iNode],iDof]    # row index 
                            J[inz_top] = mdl.ID[mdl.IEN_top[e,jNode], jDof]   # column index
                            V[inz_top] += w_top*be[i,j] 

                            E[inz_btm] = mdl.ID[mdl.IEN_btm[e,iNode],iDof]    # row index 
                            J[inz_btm] = mdl.ID[mdl.IEN_btm[e,jNode], jDof]   # column index
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