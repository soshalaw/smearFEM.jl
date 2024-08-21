using LinearAlgebra
using ProgressMeter
using SparseArrays
using Plots
using DelimitedFiles

using smearFEM

# set up mesh grid
function meshgrid(x0,x1,y0,y1,z0,z1,ne,ndim)

    NodeList = zeros(ndim,(ne+1)^ndim)
    IEN = zeros(Int64,ne^ndim,2^ndim) # IEN for the 3D mesh
    IEN_top = zeros(Int64,ne^(ndim-1),2^(ndim-1)) # IEN for the top surface
    IEN_btm = zeros(Int64,ne^(ndim-1),2^(ndim-1)) # IEN for the bottom surface
    IEN_side = zeros(Int64,ne^(ndim-1),2^(ndim-1)) # IEN for the side surfaces
    ID = zeros(Int64,(ne+1)^ndim,ndim)
    
    BorderNodes = []
    BottomBorderNodes = []
    TopBorderNodes = []

    if ndim == 2
        x = collect(range(x0, x1, length=ne+1))
        y = collect(range(y0, y1, length=ne+1))
    
        m = 1
        for j in 1:ne+1 # y direction
            for i in 1:ne+1 # x direction
                NodeList[1,m] = x[i]
                NodeList[2,m] = y[j]
                for l in 1:ndim
                    ID[m,l] = ndim*(m-1) + l
                end
                if i == 1 || i == ne+1 # populate the BorderNodes with the nodes on the left and right boundaries
                    push!(BorderNodes,m)
                end
                m = m + 1
            end
        end 
        
        n = 1
        for j in 1:ne # y direction
            for i in 1:ne # x direction
                IEN[n,1] = (j-1)*(ne+1) + i
                IEN[n,2] = (j-1)*(ne+1) + i + 1
                IEN[n,3] = j*(ne+1) + i + 1
                IEN[n,4] = j*(ne+1) + i
                if j == 1 # populate the IEN for the bottom surface
                    IEN_btm[i,1] = IEN[n,1]
                    IEN_btm[i,2] = IEN[n,2]
                elseif j == ne # populate the IEN for the top surface
                    IEN_top[i,1] = IEN[n,4]
                    IEN_top[i,2] = IEN[n,3]
                end
                n = n + 1
            end
        end

    elseif ndim == 3

        x = collect(range(x0, x1, length=ne+1))
        y = collect(range(y0, y1, length=ne+1))
        z = collect(range(z0, z1, length=ne+1))
        
        m = 1
        for k in 1:ne+1 # z direction
            for j in 1:ne+1 # y direction
                for i in 1:ne+1 # x direction
                    NodeList[1,m] = x[i]
                    NodeList[2,m] = y[j]
                    NodeList[3,m] = z[k]
                    for l in 1:ndim
                        ID[m,l] = ndim*(m-1) + l
                    end
                    if (i == 1 || i == ne+1 || j == 1 || j == ne+1)
                        push!(BorderNodes,m) # populate the BorderNodes with the nodes on the boundaries (excluding the top and bottom surfaces)
                    elseif k == 1
                        push!(BottomBorderNodes,m) # populate the BorderNodes with the nodes on the top and bottom surfaces
                    elseif k == ne + 1
                        push!(TopBorderNodes,m) # populate the BorderNodes with the nodes on the top and bottom surfaces
                    end
                    m = m + 1
                end
            end
        end
        
        n = 1       # element number on 3D mesh
        nt = 1      # element number on 2D mesh of top surface
        nb = 1      # element number on 2D mesh of bottom surface
        for k in 1:ne            # z direction
            for j in 1:ne        # y direction
                for i in 1:ne    # x direction
                    IEN[n,1] = (k-1)*(ne+1)^2 + (j-1)*(ne+1) + i
                    IEN[n,2] = (k-1)*(ne+1)^2 + (j-1)*(ne+1) + i + 1
                    IEN[n,3] = (k-1)*(ne+1)^2 + j*(ne+1) + i + 1
                    IEN[n,4] = (k-1)*(ne+1)^2 + j*(ne+1) + i
                    IEN[n,5] = k*(ne+1)^2 + (j-1)*(ne+1) + i
                    IEN[n,6] = k*(ne+1)^2 + (j-1)*(ne+1) + i + 1
                    IEN[n,7] = k*(ne+1)^2 + j*(ne+1) + i + 1
                    IEN[n,8] = k*(ne+1)^2 + j*(ne+1) + i
                    if k == 1 # populate the IEN for the bottom surface
                        IEN_btm[nb,1] = IEN[n,1]
                        IEN_btm[nb,2] = IEN[n,2]
                        IEN_btm[nb,3] = IEN[n,3]
                        IEN_btm[nb,4] = IEN[n,4]
                        nb = nb + 1
                    elseif k == ne # populate the IEN for the top surface
                        IEN_top[nt,1] = IEN[n,5]
                        IEN_top[nt,2] = IEN[n,6]
                        IEN_top[nt,3] = IEN[n,7]
                        IEN_top[nt,4] = IEN[n,8]
                        nt = nt + 1
                    # elseif j == 1 || j == ne || i == 1 || i == ne
                    #     IEN_side[n,1] = IEN[n,1]
                    #     IEN_side[n,2] = IEN[n,2]
                    #     IEN_side[n,3] = IEN[n,3]
                    #     IEN_side[n,4] = IEN[n,4]
                    #     IEN_side[n,5] = IEN[n,5]
                    #     IEN_side[n,6] = IEN[n,6]
                    #     IEN_side[n,7] = IEN[n,7]
                    #     IEN_side[n,8] = IEN[n,8]
                    end
                    n = n + 1
                end
            end
        end
    end
    return NodeList, IEN, ID, IEN_top, IEN_btm, [BorderNodes, BottomBorderNodes, TopBorderNodes]
end


function setboundaryCond(NodeList, ne, ndim, FunctionClass, nDof=1)
    """
        Set the Dirichlet boundary conditions for the problem

        Parameters:
        NodeList: {[ndim,nNodes] Matrix{Float64}} : matrix of the coordinates of the nodes
        ne: {Int} : number of elements
        ndim: {Int} : number of dimensions
        FunctionClass: {String} : type of basis function
        d: {Float64} : displacement increment
        nDOf: {Int} : number of degree of freedom per node

        Returns:
        q_d: {[ndof] Vector{Float64}} : Dirichlet boundary conditions
    """

    if FunctionClass == "Q1"
        q_upper = zeros(nDof*(ne+1)^ndim,1)              # initialize the vector of the Dirichlet boundary conditions upper surface
        q_lower = zeros(nDof*(ne+1)^ndim,1)              # initialize the vector of the Dirichlet boundary conditions lower surface
        C = sparse(I,ndim*(ne+1)^ndim,ndim*(ne+1)^ndim)  # definition of the constraint matrix
    end

    z0Bound = 0
    z1Bound = 1

    rCol = Array{Int}(undef,0)

    for nNode in 1:size(NodeList,2)
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

    return q_upper, q_lower, C_uc
end

function apply_boundary_conditions(ne, NodeList, IEN, IEN_top, IEN_btm, ndim, FunctionClass, ID=None, nDof=3)
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

    E = zeros(  Int64, ne^(ndim-1)*(size(ID,2)*size(IEN_btm,2))^2*2) # *2 because we have two surfaces
    J = zeros(  Int64, ne^(ndim-1)*(size(ID,2)*size(IEN_btm,2))^2*2) # *2 because we have two surfaces
    V = zeros(Float64, ne^(ndim-1)*(size(ID,2)*size(IEN_btm,2))^2*2) # *2 because we have two surfaces

    for e in 1:ne^(ndim-1)
    
        coords_top = NodeList[:,IEN_top[e,:]] # get the coordinates of the nodes of the element
        coords_btm = NodeList[:,IEN_btm[e,:]] # get the coordinates of the nodes of the element

        if ndim == 2
            # gaussian quadrature points for the element [-1,1] 
            ξ, w_ξ = gaussian_quadrature(-1,1)

            wpoints = [w_ξ[1], w_ξ[2]]
            
            x = [ξ[1], ξ[2]]
        elseif ndim == 3
            # gaussian quadrature points for the element [-1,1]x[-1,1] 
            ξ, w_ξ = gaussian_quadrature(-1,1)
            η, w_η = gaussian_quadrature(-1,1)
            
            wpoints = [w_ξ[1]*w_η[1], w_ξ[2]*w_η[1], w_ξ[2]*w_η[2], w_ξ[1]*w_η[2]]
            
            x = [ξ[1], ξ[2], ξ[2], ξ[1]]
            y = [η[1], η[1], η[2], η[2]]
        end 

        for gp in 1:2^(ndim-1)

            if ndim == 2
                N, ΔN = basis_function(x[gp], nothing, nothing, FunctionClass)
            elseif ndim == 3
                N, ΔN = basis_function(x[gp], y[gp], nothing, FunctionClass) 
            end

            dxdξ_top = coords_top*ΔN         # Jacobian matrix [dx/dxi dx/deta; dy/dxi dy/deta; dz/dxi dz/deta]
            dxdξ_btm = coords_btm*ΔN         # Jacobian matrix [dx/dxi dx/deta; dy/dxi dy/deta; dz/dxi dz/deta]

            w_top = wpoints[gp]*norm(cross(dxdξ_top[:,1],dxdξ_top[:,2]))     # weight of the quadrature point top surface
            w_btm = wpoints[gp]*norm(cross(dxdξ_btm[:,1],dxdξ_btm[:,2]))     # weight of the quadrature point bottom surface
            
            M = zeros(3, ndim*length(N))
            M[1,1:nDof:end] = N
            M[2,2:nDof:end] = N
            M[3,3:nDof:end] = N

            be = M'*M

            # loop between basis functions of the element
            for iNode in 1:size(be,1)÷nDof
                for jNode in 1:size(be,2)÷nDof
                    for iDof in 1:size(ID,2)
                        for jDof in 1:size(ID,2)
                            i = (iNode-1)*nDof + iDof
                            j = (jNode-1)*nDof + jDof
                            
                            inz_btm = length(be)*(e-1) + (iNode-1)*nDof*size(be,2) + (jNode-1)*nDof^2 + (iDof-1)*nDof + jDof # index for the COO sparse matrix
                            inz_top = length(be)*(ne)^(ndim-1)+ length(be)*((e-1)) + (iNode-1)*nDof*size(be,2) + (jNode-1)*nDof^2 + (iDof-1)*nDof + jDof # index for the COO sparse matrix
                            
                            E[inz_top] = ID[IEN_top[e,iNode],iDof]    # row index 
                            J[inz_top] = ID[IEN_top[e,jNode], jDof]   # column index
                            V[inz_top] += w_top*be[i,j] 

                            E[inz_btm] = ID[IEN_btm[e,iNode],iDof]    # row index 
                            J[inz_btm] = ID[IEN_btm[e,jNode], jDof]   # column index
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

function simulate(x0, x1, y0, y1, z0, z1, ne, Young, ν, ndim, FunctionClass, nDof, β, endTime, tSteps, Control, cParam, writeData=false)

    time = collect(range(0,endTime,tSteps)) # time vector

    NodeList, IEN, ID, IEN_top, IEN_btm, BorderNodesList = meshgrid(x0,x1,y0,y1,z0,z1,ne,ndim)  # generate the mesh grid
    NodeListCylinder = inflate_sphere(NodeList, x0, x1, y0, y1)                                 # inflate the sphere to a unit sphere
    q_tp, q_btm, C_uc = setboundaryCond(NodeList, ne, ndim, FunctionClass, nDof)

    state = "init"

    BorderNodes2D, NodeList2D = extract_borders(NodeListCylinder, CameraMatrix, BorderNodesList, state, ne)
    pi, qi = fit_curve(BorderNodes2D)

    SideBorders = BorderNodesList[1]
    BottomBorders = BorderNodesList[2]
    TopBorders = BorderNodesList[3]
        
    fields = [NodeListCylinder]                                                               # store the solution fields of the mesh in 3D
    fields2D = [NodeList2D]                                                                   # store the solution fields of the mesh in 2D
    BorderNodes = [NodeList[:,SideBorders] NodeList[:,BottomBorders] NodeList[:,TopBorders]]  # store the solution fields of the surfaces in 3D
    borderfields2D = [BorderNodes2D]                                                          # store the solution fields of the surfaces in 2D
    splinep = [pi]                                                                            # store the x coordinates samples of the spline parameters of the border nodes
    splineq = [qi]                                                                            # store the y coordinates samples of the spline parameters of the border nodes
    output = []

    state = "update"
    μ_btm = 0      
    iter = 1

    if Control == "force"
        @showprogress "Simulating with prescribed $Control ..." for t in time

            K = assemble_system(ne, NodeListCylinder, IEN, ndim, FunctionClass, nDof, ID, Young, ν)                   # assemble the stiffness matrix
            b = apply_boundary_conditions(ne, NodeListCylinder, IEN, IEN_top, IEN_btm, ndim, FunctionClass, ID)   # apply the neumann boundary conditions
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

            BorderNodes2D, NodeList2D = extract_borders(NodeListCylinder, CameraMatrix, BorderNodesList, state)
            BorderNodes = [NodeListCylinder[:,SideBorders] NodeListCylinder[:,BottomBorders] NodeListCylinder[:,TopBorders]]
            pi, qi = fit_curve(BorderNodes2D)

            push!(output, μ_tp)
            push!(fields, NodeListCylinder)
            push!(fields2D, NodeList2D)
            push!(borderfields2D, BorderNodes2D)
            push!(splinep, pi)
            push!(splineq, qi) 

            iter += 1
        end

    elseif Control == "displacement"

        @showprogress "Simulating with prescribed $Control ..." for t in time
        
            K = assemble_system(ne, NodeListCylinder, IEN, ndim, FunctionClass, nDof, ID, Young, ν)                   # assemble the stiffness matrix
            b = apply_boundary_conditions(ne, NodeListCylinder, IEN, IEN_top, IEN_btm, ndim, FunctionClass, ID)       # apply the neumann boundary conditions
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

            BorderNodes2D, NodeList2D = extract_borders(NodeListCylinder, CameraMatrix, BorderNodesList, state)
            BorderNodes = [NodeListCylinder[:,SideBorders] NodeListCylinder[:,BottomBorders] NodeListCylinder[:,TopBorders]]
            pi, qi = fit_curve(BorderNodes2D)

            # store the solutions in a list
            push!(output, F_est[1])
            push!(fields, NodeListCylinder)
            push!(fields2D, NodeList2D)
            push!(borderfields2D, BorderNodes2D)
            push!(splinep, pi)
            push!(splineq, qi) 

            iter += 1
        end
    end

    if writeData
        fileName = "vtkFiles/squeeze_flow"
        write_scene(fileName, NodeList, IEN, ne, ndim, fields)
        animate_fields(fields,fields2D,borderfields2D,IEN, splinep, splineq)
    end

    return output
end

function test()
    # test case 
    x0 = 0
    x1 = 1
    y0 = 0
    y1 = 1
    z0 = 0
    z1 = 1
    ne = 8
    Young = 40
    ν = 0.4
    ndim = 3
    FunctionClass = "Q1"
    nDof = ndim  # number of degree of freedom per node
    β = 100

    Control = "force"
    cParam = F*ones(tSteps)
    μ_list = simulate(x0, x1, y0, y1, z0, z1, ne, Young, ν, ndim, FunctionClass, nDof, β, endTime, tSteps, Control, cParam)

    Control = "displacement"
    F_list = simulate(x0, x1, y0, y1, z0, z1, ne, Young, ν, ndim, FunctionClass, nDof, β, endTime, tSteps, Control, μ_list)

    isapprox(F_list, cParam, atol=1e-3)

end

test()
