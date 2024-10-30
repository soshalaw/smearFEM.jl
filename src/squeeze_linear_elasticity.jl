using LinearAlgebra
using ProgressMeter
using SparseArrays

# set up mesh grid

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
        E = zeros(  Int64, mdl.ne^mdl.ndim*((size(mdl.ID,2)*size(mdl.IEN,2))^2))
        J = zeros(  Int64, mdl.ne^mdl.ndim*((size(mdl.ID,2)*size(mdl.IEN,2))^2))
        V = zeros(Float64, mdl.ne^mdl.ndim*((size(mdl.ID,2)*size(mdl.IEN,2))^2))  
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
    vol = 0
    for e in 1:mdl.ne^mdl.ndim
        coords = mdl.NodeList[:,mdl.IEN[e,:]] # get the coordinates of the nodes of the element
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
function apply_boundary_conditions(mdl)

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
function setboundaryCond(NodeList, ne, ndim, FunctionClass, nDof=1)
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
        cMat = [[ 2*c2+c1  c1    c1    0  0  0]; 
                [  c1   2*c2+c1  c1    0  0  0]; 
                [  c1     c1   2*c2+c1 0  0  0]; 
                [  0      0      0     c2 0  0]; 
                [  0      0      0     0 c2  0]
                [  0      0      0     0  0 c2]]  # constitutive matrix
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
function simulate(x0, x1, y0, y1, z0, z1, ne, Young, ν, ndim, FunctionClass, nDof, β, CameraMatrix, endTime, tSteps, Control, cParam, cMat; writeData=false, filepath=nothing)

    time = collect(range(start=0,stop=endTime,length=tSteps)) # time vector

    NodeList, IEN, ID, IEN_top, IEN_btm, BorderNodesList = meshgrid_cube(x0,x1,y0,y1,z0,z1,ne,ndim,FunctionClass=FunctionClass)  # generate the mesh grid
    NodeListCylinder = inflate_cylinder(NodeList, x0, x1, y0, y1)                                 # inflate the sphere to a unit sphere
    q_tp, q_btm, C_uc = setboundaryCond(NodeList, ne, ndim, FunctionClass, nDof)

    mdl = def_model("linear_elasticity", ne=ne, NodeList=NodeList, IEN=IEN, IEN_top=IEN_top, IEN_btm=IEN_btm, ndim=ndim, nDof=nDof, FunctionClass=FunctionClass, ID=ID, Young=Float64(Young), ν=ν, cMat=cMat)
                
    state = "init"

    BorderPts2D, BorderNodes2D, Nodes2D = extract_borders(NodeListCylinder, CameraMatrix, BorderNodesList, state, ne, 2*ne+1)
    pi, qi = fit_curve(border=BorderPts2D)

    SideBorders = BorderNodesList[1]
    BottomBorders = BorderNodesList[2]
    TopBorders = BorderNodesList[3]
        
    fields = []  
    pos3D = [NodeListCylinder]                                                             # store the solution fields of the mesh in 3D
    pos2D = [Nodes2D]                                                                   # store the solution fields of the mesh in 2D
    surfaceNodesList = [NodeList[:,SideBorders] NodeList[:,BottomBorders] NodeList[:,TopBorders]]  # store the solution fields of the surfaces in 3D
    borderPts2DList = [BorderPts2D]                                                               # store the solution fields of the surfaces in 2D
    borderNodeList2D = [BorderNodes2D]                                                       # store the solution fields of the border nodes in 2D
    splinep = [pi]                                                                            # store the x coordinates samples of the spline parameters of the border nodes
    splineq = [qi]                                                                            # store the y coordinates samples of the spline parameters of the border nodes
    output = []
    writeborderList = [vcat(pi', qi')]

    state = "update"
    μ_btm = 0      
    iter = 1

    pr = Progress(tSteps; desc= "Simulating with prescribed $Control ...", showspeed=true)
    if Control == "force"
        for t in time

            K = assemble_system(mdl)                   # assemble the stiffness matrix
            b = apply_boundary_conditions(mdl)   # apply the neumann boundary conditions
            q_btm = μ_btm*q_btm
            
            K_bar = K + β*b

            C_T = transpose(C_uc)              # transpose the constraint matrix

            M = [C_T*K_bar*C_uc C_T*K_bar*q_tp; q_tp'*K_bar*C_uc q_tp'*K_bar*q_tp] # assemble the system of equations]

            invM = inv(Matrix(M))                      # invert the system of equations

            sol = invM*[-C_T*K_bar*q_btm; cParam[iter].-q_tp'*K_bar*q_btm] # solve the system of equations

            q_f = sol[1:end-1]        # solve the system of equations
            μ_tp = sol[end]           # solve the system of equations
            q = q_btm + C_uc*q_f + μ_tp*q_tp;                # assemble the solution

            # post process the solution
            motion = [q[ID[:,1]] q[ID[:,2]] q[ID[:,3]]]'    # update the nodal positions
            NodeListCylinder = NodeListCylinder + motion    # update the node coordinates

            BorderPts2D, BorderNodes2D, Nodes2D = extract_borders(NodeListCylinder, CameraMatrix, BorderNodesList, state)
            surfaceNodesList = [NodeListCylinder[:,SideBorders] NodeListCylinder[:,BottomBorders] NodeListCylinder[:,TopBorders]]
            pi, qi = fit_curve(border=BorderPts2D)

            push!(output, μ_tp)
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
    elseif Control == "displacement"
        pr = Progress(tSteps; desc="Simulating with prescribed $Control ...", showspeed=true)
        for t in time
        
            K = assemble_system(mdl)                   # assemble the stiffness matrix
            b = apply_boundary_conditions(mdl)       # apply the neumann boundary conditions
            q_d = (μ_btm*q_btm + cParam[iter]*q_tp)                  # apply the Dirichlet boundary conditions

            K_bar = K + β*b

            C_T = transpose(C_uc)           # transpose the constraint matrix
            K_free = C_T*K_bar*C_uc         # extract the free part of the stiffness matrix

            invK = inv(Matrix(K_free))

            q_f = invK*C_T*(-K_bar*q_d)         # solve the system of equations
            q = q_d + C_uc*q_f;                 # assemble the solution 

            # post process the solution
            f_R = K_bar*q
            motion = [q[ID[:,1]] q[ID[:,2]] q[ID[:,3]]]'
            F_est = q_tp'*f_R                                         # calculate the reaction force at the top surface F = Σf^{tp}_{iR} = q_tp'*f_R
            NodeListCylinder = NodeListCylinder + motion              # update the node coordinates
            # q_new, IEN_new = rearrange(q, ne, ndim, IEN, FunctionClass)
            BorderPts2D, BorderNodes2D, Nodes2D = extract_borders(NodeListCylinder, CameraMatrix, BorderNodesList, state)
            surfaceNodesList = [NodeListCylinder[:,SideBorders] NodeListCylinder[:,BottomBorders] NodeListCylinder[:,TopBorders]]
            pi, qi = fit_curve(border=BorderPts2D)

            # store the solutions in a list
            push!(output, F_est[1])
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
    end

    if writeData
        write_scene(string(filepath,"/Results"), NodeList, IEN, ne, ndim, pos3D)
        animate_fields(filepath = string(filepath,"/Results/images"), fields=pos3D , IEN=IEN, BorderNodes2D=borderPts2DList, fields2D=pos2D)
        writeCSV(string(filepath,"/Results"), writeborderList)
    end
    return output, borderPts2DList, borderNodeList2D, splinep, splineq, mdl
end

"""
    set_file(filepath)

Create the directories to store the data

# Arguments:
- `filepath::String` : path to the file
"""
function set_file(filepath)
    if !isdir(filepath)
        mkdir(filepath)
        mkdir(string(filepath,"/Results"))
        mkdir(string(filepath,"/Results/contour_data"))
        mkdir(string(filepath,"/Results/vtkFiles"))
        mkdir(string(filepath,"/Results/images"))
        mkdir(string(filepath,"/Results/cost"))
        mkdir(string(filepath,"/Results/cost_lame"))
    end
end

"""
    test(x0, x1, y0, y1, z0, z1, ne, Young, ν, ndim, FunctionClass, nDof, β, CameraMatrix, endTime, tSteps, Control; writeData=false, filepath=nothing, mode = "standard")

Test the simulation

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
- `writeData::Bool` : write the data to a file
- `filepath::String` : path to the file
- `mode::String` : type of constitutive matrix
"""
function test(x0, x1, y0, y1, z0, z1, ne, Young, ν, ndim, FunctionClass, nDof, β, CameraMatrix, endTime, tSteps, Control; writeData=false, filepath=nothing, mode = "standard")
    
    if writeData
        isnothing(filepath) || AssertionError("Please provide a filepath to write the data")
        set_file(filepath)
    end

    μ_tp = 0.3

    if Control == "force"
        cParam = -0.25*ones(tSteps)
    elseif Control == "displacement"
        cParam = -(μ_tp/tSteps)*ones(tSteps)
    end

    cMat = get_cMat(mode, Young, ν)

    μ_list, simBorderPts, simBorderNodes, splinex, spliney, mdl = simulate(x0, x1, y0, y1, z0, z1, ne, Young, ν, ndim, FunctionClass, nDof, β, CameraMatrix, endTime, tSteps, Control, cParam, cMat, writeData=writeData, filepath=filepath)

    if !writeData
        println("comparing the simulated and observed data")
        isnothing(filepath) || AssertionError("Please provide a filepath to read the data")
        # animate_fields(filepath = string(filepath,"/Results"), p=splinex, q=spliney, pObs=splinexObs, qObs=splineyObs)
        obsBorderPts, splinexObs, splineyObs = readCSV(string(filepath,"/Results/contour_data"))                         # read the observation data
        
        # animate_fields(filepath = string(filepath,"/Results/cost"), BorderNodes2D=simBorderPts, p=splinex, q=spliney, pObs=splinexObs, qObs=splineyObs) # animate the fields

        pairs = match_points(simBorderPts[1], [splinexObs[1],splineyObs[1]]) # match the points using the first border

        # plot_matches(simBorderPts, splinex, spliney, splinexObs, splineyObs, pairs, string(filepath,"/Results/cost"))

        # test the closest point function
        d_h, xObsintlst, xSimintlst, ySimintlst = height_sample(simBorderPts, obsBorderPts)
        # plot_matches_h(xObsintlst, ySimintlst, xSimintlst, splinex, spliney, splinexObs, splineyObs, string(filepath,"/Results/cost"))

        d_cp = closest_point(simBorderPts, obsBorderPts, pairs)
        # d_cp = zeros(length(d_h))

        return d_h, d_cp
    else
        return 0, 0
    end
end

"""
    write_sim_data(x0, x1, y0, y1, z0, z1, ne, Young, ν, ndim, FunctionClass, nDof, β, CameraMatrix, endTime, tSteps, Control, filename; mode = "standard")

Initialize the simulation and write the data to a file

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
- `filename::String` : path to the file
"""
function write_sim_data(x0, x1, y0, y1, z0, z1, ne, Youngtst, νtst, ndim, FunctionClass, nDof, β, CameraMatrix, endTime, tSteps, Control, filename; mode = "standard")
    writeData = true
    hcost, cpCost = test(x0, x1, y0, y1, z0, z1, ne, Youngtst, νtst, ndim, FunctionClass, nDof, β, CameraMatrix, endTime, tSteps, Control, writeData=writeData, filepath=filename, mode = mode)
end

"""
    function initialize_mesh_test(x0, x1, y0, y1, z0, z1, ne, ndim, FunctionClass, CameraMatrix, filepath=nothing)
        
Initialize the mesh and write the data to a file

# Arguments:
- `x0::Float64` : x-coordinate of the lower left corner of the domain
- `x1::Float64` : x-coordinate of the upper right corner of the domain
- `y0::Float64` : y-coordinate of the lower left corner of the domain
- `y1::Float64` : y-coordinate of the upper right corner of the domain
- `z0::Float64` : z-coordinate of the lower left corner of the domain
- `z1::Float64` : z-coordinate of the upper right corner of the domain
- `ne::Int` : number of elements
- `ndim::Int` : number of dimensions
- `FunctionClass::String` : type of basis function
- `CameraMatrix::Matrix{Float64}` : camera matrix
- `filepath::String` : path to the file
"""
function initialize_mesh_test(x0, x1, y0, y1, z0, z1, ne, ndim, FunctionClass, CameraMatrix, filepath=nothing)

    isnothing(filepath) || AssertionError("Please provide a filepath to write the data")
    set_file(filepath)
    
    NodeList, IEN, ID, IEN_top, IEN_btm, BorderNodesList = meshgrid_cube(x0,x1,y0,y1,z0,z1,ne,ndim,FunctionClass=FunctionClass)  # generate the mesh grid
    NodeListCylinder = inflate_cylinder(NodeList, x0, x1, y0, y1)                                 # inflate the sphere to a unit sphere

    state = "init"

    BorderPts2D, BorderNodes2D, Nodes2D = extract_borders(NodeListCylinder, CameraMatrix, BorderNodesList, state, ne, (2*ne+1))
    pi, qi = fit_curve(border=BorderPts2D)

    SideBorders = BorderNodesList[1]
    BottomBorders = BorderNodesList[2]
    TopBorders = BorderNodesList[3]
        
    fields = [NodeListCylinder]                                                               # store the solution fields of the mesh in 3D
    pos2D = [Nodes2D]                                                                   # store the solution fields of the mesh in 2D
    surfaceNodesList = [NodeList[:,SideBorders] NodeList[:,BottomBorders] NodeList[:,TopBorders]]  # store the solution fields of the surfaces in 3D
    borderPts2DList = [BorderPts2D]                                                               # store the solution fields of the surfaces in 2D
    borderNodeList2D = [BorderNodes2D]                                                       # store the solution fields of the border nodes in 2D
    splinep = [pi]                                                                            # store the x coordinates samples of the spline parameters of the border nodes
    splineq = [qi]                                                                            # store the y coordinates samples of the spline parameters of the border nodes
    writeborderList = [vcat(pi', qi')]

    animate_fields(filepath = string(filepath,"/Results/images"), fields=pos3D , IEN=IEN, BorderNodes2D=borderPts2DList, fields2D=pos2D)

end

function simulate_single_tstep(x0, x1, y0, y1, z0, z1, ne, Young, ν, ndim, FunctionClass, nDof, β, μ_tp, μ_btm)

    mode = "standard"

    cMat = get_cMat(mode, Young, ν)

    NodeList, IEN, ID, IEN_top, IEN_btm, BorderNodesList = meshgrid_cube(x0,x1,y0,y1,z0,z1,ne,ndim,FunctionClass=FunctionClass)  # generate the mesh grid
    NodeListCylinder = inflate_cylinder(NodeList, x0, x1, y0, y1)                                 # inflate the sphere to a unit sphere
    q_tp, q_btm, C_uc = setboundaryCond(NodeList, ne, ndim, FunctionClass, nDof)

    mdl = def_model("linear_elasticity", ne=ne, NodeList=NodeList, IEN=IEN, IEN_top=IEN_top, IEN_btm=IEN_btm, ndim=ndim, nDof=nDof, 
                    FunctionClass=FunctionClass, ID=ID, Young=Float64(Young), ν=ν, cMat=cMat)
        
    K = assemble_system(mdl)                   # assemble the stiffness matrix
    b = apply_boundary_conditions(mdl)       # apply the neumann boundary conditions
    q_d = (μ_btm*q_btm + μ_tp*q_tp)                  # apply the Dirichlet boundary conditions

    K_bar = K + β*b

    C_T = transpose(C_uc)           # transpose the constraint matrix
    K_free = C_T*K_bar*C_uc         # extract the free part of the stiffness matrix

    invK = inv(Matrix(K_free))

    q_f = invK*C_T*(-K_bar*q_d)         # solve the system of equations
    q = q_d + C_uc*q_f;                 # assemble the solution 

    q_out = [q[ID[:,1]] q[ID[:,2]] q[ID[:,3]]]'

    return q_out, mdl
end
