using LinearAlgebra
using SparseArrays

function greet_fem()
    println("Hello, I am the FEM module")
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
function assemble_system(ne, NodeList, IEN, ndim, FunctionClass="Q1", nDof=1, ID=nothing, Young=1, ν=0.3)

    # (I,J,V) vectors for COO sparse matrix
    if nDof == 1
        E = zeros(  Int64, ne^ndim*size(IEN,2)^2)
        J = zeros(  Int64, ne^ndim*size(IEN,2)^2)
        V = zeros(Float64, ne^ndim*size(IEN,2)^2)
    else
        E = zeros(  Int64, ne^ndim*((size(ID,2)*2^ndim)^2))
        J = zeros(  Int64, ne^ndim*((size(ID,2)*2^ndim)^2))
        V = zeros(Float64, ne^ndim*((size(ID,2)*2^ndim)^2))  
    end

    # element loop
    if ndim == 1
        # gaussian quadrature points for the element [-1,1] 
        ξ, w_ξ = gaussian_quadrature(-1,1)
        
        wpoints = [w_ξ[1], w_ξ[2]]
        
        x = [ξ[1], ξ[2]]
    elseif ndim == 2
        # gaussian quadrature points for the element [-1,1]x[-1,1] 
        ξ, w_ξ = gaussian_quadrature(-1,1)
        η, w_η = gaussian_quadrature(-1,1)
        
        wpoints = [w_ξ[1]*w_η[1], w_ξ[2]*w_η[1], w_ξ[2]*w_η[2], w_ξ[1]*w_η[2]]
        
        x = [ξ[1], ξ[2], ξ[2], ξ[1]]
        y = [η[1], η[1], η[2], η[2]]

    elseif ndim == 3
        # gaussian quadrature points for the element [-1,1]x[-1,1]x[-1,1] 
        ξ, w_ξ = gaussian_quadrature(-1,1)
        η, w_η = gaussian_quadrature(-1,1)
        ζ, w_ζ = gaussian_quadrature(-1,1)
        
        wpoints = [w_ξ[1]*w_η[1]*w_ζ[1], w_ξ[2]*w_η[1]*w_ζ[1], w_ξ[2]*w_η[2]*w_ζ[1], w_ξ[1]*w_η[2]*w_ζ[1], w_ξ[1]*w_η[1]*w_ζ[2], w_ξ[2]*w_η[1]*w_ζ[2], w_ξ[2]*w_η[2]*w_ζ[2], w_ξ[1]*w_η[2]*w_ζ[2]]
        
        x = [ξ[1], ξ[2], ξ[2], ξ[1], ξ[1], ξ[2], ξ[2], ξ[1]]
        y = [η[1], η[1], η[2], η[2], η[1], η[1], η[2], η[2]]
        z = [ζ[1], ζ[1], ζ[1], ζ[1], ζ[2], ζ[2], ζ[2], ζ[2]]
    end

    for e in 1:ne^ndim
        coords = NodeList[:,IEN[e,:]] # get the coordinates of the nodes of the element
        
        # integration loop
        for gp in 1:2^ndim
            if ndim == 1
                N, ΔN = basis_function(x[gp], nothing, nothing, FunctionClass)
            elseif ndim == 2
                N, ΔN = basis_function(x[gp], y[gp], nothing, FunctionClass) 
            elseif ndim == 3
                N, ΔN = basis_function(x[gp], y[gp], z[gp], FunctionClass) 
            end

            Jac  = coords*ΔN # Jacobian matrix [dx/dxi dx/deta; dy/dxi dy/deta]

            w = wpoints[gp]*abs(det(Jac))
            invJ = inv(Jac)
            dNdX = ΔN*invJ
            
            if nDof == 1
                szN = size(N,1) # number of basis functions
                # loop between basis functions of the element
                for i in 1:szN
                    for j in 1:szN
                        inz = (szN)^2*(e-1) + szN*(i-1) + j # index for the COO sparse matrix
                        E[inz] = IEN[e,i] # row index 
                        J[inz] = IEN[e,j] # column index
                        V[inz] += w*dot(dNdX[i,:],dNdX[j,:])# inner product of the gradient of the basis functions
                    end
                end
            else   
                if nDof == 2
                    B = zeros(3, ndim*length(N))
                    B[1,1:nDof:end] = dNdX[:,1]
                    B[2,2:nDof:end] = dNdX[:,2]
                    B[3,1:nDof:end] = dNdX[:,2]
                    B[3,2:nDof:end] = dNdX[:,1]

                    cMat = [[Young/(1-ν^2) ν*Young/(1-ν^2) 0]; [ν*Young/(1-ν^2) Young/(1-ν^2) 0]; [0 0 Young/(2*(1+ν))]] # constitutive matrix for plane stress
                elseif nDof == 3
                    B = zeros(6, ndim*length(N))
                    B[1,1:nDof:end] = dNdX[:,1]
                    B[2,2:nDof:end] = dNdX[:,2]
                    B[3,3:nDof:end] = dNdX[:,3]
                    B[4,2:nDof:end] = dNdX[:,3]
                    B[4,3:nDof:end] = dNdX[:,2]
                    B[5,1:nDof:end] = dNdX[:,3]
                    B[5,3:nDof:end] = dNdX[:,1]
                    B[6,1:nDof:end] = dNdX[:,2]
                    B[6,2:nDof:end] = dNdX[:,1]

                    cMat = [[1-ν ν ν 0 0 0]; [ν 1-ν ν 0 0 0]; [ν ν 1-ν 0 0 0]; [0 0 0 (1-2*ν)/2 0 0]; [0 0 0 0 (1-2*ν)/2 0]; [0 0 0 0 0 (1-2*ν)/2]]*(Young/((1+ν)*(1-2*ν))) # constitutive matrix
                end

                Ke = B'*cMat*B*w # element stiffness matrix
        
                # loop between basis functions of the element
                for iNode in 1:size(Ke,1)÷nDof
                    for jNode in 1:size(Ke,2)÷nDof
                        for iDof in 1:size(ID,2)
                            for jDof in 1:size(ID,2)
                                i = (iNode-1)*nDof + iDof
                                j = (jNode-1)*nDof + jDof
                                inz = length(Ke)*(e-1) + (iNode-1)*nDof*size(Ke,2) + (jNode-1)*nDof^2 + (iDof-1)*nDof + jDof # index for the COO sparse matrix
                                E[inz] = ID[IEN[e,iNode],iDof] # row index 
                                J[inz] = ID[IEN[e,jNode], jDof] # column index
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