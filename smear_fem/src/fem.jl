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

    function basis_function_glob(x1, x2, y1, y2, x, y, h)
        """ Define local basis functions and the gradients for a square element [x1,x2]x[y1,y2]
            x1, x2, y1, y2 are the coordinates of the lower left corner of the element
            x, y are the coordinates of the point where the basis function is evaluated
            h is the size of the element
            
            Returns:
            N: basis functions {[ndof] Vector{Float64}}
            Delta_N: gradient of the basis functions {[ndof,ndim] Matrix{Float64}}
        """

        N = [(x-x1)*(y-y1)/(h^2), (x-x1)*(y2-y)/(h^2), (x2-x)*(y-y1)/(h^2), (x2-x)*(y2-y)/(h^2)]

        Delta_N = [[(y-y1)/(h^2) , (y2-y)/(h^2), -(y-y1)/(h^2), -(y2-y)/(h^2)] [(x-x1)/(h^2), -(x-x1)/(h^2), (x2-x)/(h^2), -(x2-x)/(h^2)]]

        return N, Delta_N
    end 

    function basis_function(x, y)
        
        """ Define local basis functions and the gradients for a square element [-1,1]x[-1,1]
            x, y are the coordinates of the point where the basis function is evaluated
            
            Returns:
            N: basis functions {[ndof] Vector{Float64}} 
            Delta_N: gradient of the basis functions {[ndof,ndim] Matrix{Float64}}
        """
        
            # basis functions
            N = [(1-x)*(1-y)/4, (x+1)*(1-y)/4, (1+x)*(y+1)/4, (1-x)*(1+y)/4]
        
            # gradient of the basis functions
            Delta_N = [[-(1-y)/4 , (1-y)/4, (y+1)/4, -(1+y)/4] [-(1-x)/4, -(x+1)/4, (1+x)/4, (1-x)/4] ]
        
            return N, Delta_N
    end 

    function assemble_system(ne, NodeList, IEN)
        """ Assembles the finite element system. Returns the gloval stiffness matrix

            Returns:
            K: sparse stiffness matrix [ndof,ndof] SparseMatrixCSC{Float64,Int64}
        """
        # (I,J,V) vectors for COO sparse matrix
        E = zeros(Int64, 16*ne*ne)
        J = zeros(Int64, 16*ne*ne)
        V = zeros(Float64, 16*ne*ne)

        v = 0
        # element loop
        b = zeros((ne+1)*(ne+1),1)

        for e in 1:ne*ne

            # gaussian quadrature points for the element [-1,1]x[-1,1] 
            xi, w_xi = gaussian_quadrature(-1,1)
            eta, w_eta = gaussian_quadrature(-1,1)

            wpoints = [w_xi[1]*w_eta[1], w_xi[2]*w_eta[1], w_xi[2]*w_eta[2], w_xi[1]*w_eta[2]]
            
            x = [xi[1], xi[2], xi[2], xi[1]]
            y = [eta[1], eta[1], eta[2], eta[2]]

            coords = NodeList[:,IEN[e,:]] # get the coordinates of the nodes of the element
            
            # integration loop
            for gp in 1:4
                N, Delta_N = basis_function(x[gp],y[gp]) 

                Jac  = coords*Delta_N # Jacobian matrix [dx/dxi dx/deta; dy/dxi dy/deta]

                w = wpoints[gp]*abs(det(Jac))
                invJ = inv(Jac)
                dNdX = Delta_N*invJ
                
                # loop between basis functions of the element
                for i in 1:4
                    for j in 1:4
                        inz = 16*(e-1) + 4*(i-1) + j # index for COO sparse matrix
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