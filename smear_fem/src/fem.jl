module fem
    using LinearAlgebra
    using SparseArrays

    function greet()
        println("Hello, I am the FEM module")
    end

    function gaussian_quadrature(a,b)
        """ Compute the nodes and weights for the Gaussian quadrature of order 2
            a, b are the limits of the integration interval
            
            Returns:    
            ξ: nodes [nGaussPoints]-element Vector{Float64}
            w: weights [nGaussPoints]-element Vector{Float64}
        """

        ξ = [-(b-a)/(2*sqrt(3))+(b+a)/2, (b-a)/(2*sqrt(3))+(b+a)/2]
        w = [(b-a)/2, (b-a)/2]
        return ξ, w
    end

    function basis_function(ξ,η=nothing,ζ=nothing)
        """ Define the basis functions and the gradients for a master element
            xi, eta, zeta are the coordinates of the point where the basis function is evaluated
            FunctionClass is the type of basis functions to be considered
            
            Returns:
            N: basis functions {[ndof] Vector{Float64}}
            Delta_N: gradient of the basis functions {[ndof,ndim] Matrix{Float64}}
        """
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
        return N, ΔN
    end

    function assemble_system(ne, NodeList, IEN, ndim)
        """ Assembles the finite element system. Returns the global stiffness matrix
    
            Returns:
            K: sparse stiffness matrix {[ndof,ndof] SparseMatrixCSC{Float64,Int64}}
        """
        # (I,J,V) vectors for COO sparse matrix

        E = zeros(Int64, (4*ne)^ndim)
        J = zeros(Int64, (4*ne)^ndim)
        V = zeros(Float64, (4*ne)^ndim)

        # element loop
        if ndim == 2
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

                if ndim == 2
                    N, ΔN = basis_function(x[gp], y[gp])
                elseif ndim == 3    
                    N, ΔN = basis_function(x[gp], y[gp], z[gp])
                end
    
                Jac  = coords*ΔN # Jacobian matrix [dx/dxi dx/deta; dy/dxi dy/deta]
    
                w = wpoints[gp]*abs(det(Jac))
                invJ = inv(Jac)
                dNdX = ΔN*invJ
                
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
            end
        end
    
        K = sparse(E,J,V)
    
        return K
    end

end # module fem
