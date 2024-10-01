using LinearAlgebra
using ProgressMeter
using SparseArrays
using Plots

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

function simulate(x0, x1, y0, y1, z0, z1, ne, Young, ν, ndim, FunctionClass, nDof, β, CameraMatrix, endTime, tSteps, Control, cParam; cMat, writeData=true, csv_path)

    time = collect(range(start=0,stop=endTime,length=tSteps)) # time vector

    NodeList, IEN, ID, IEN_top, IEN_btm, BorderNodesList = meshgrid(x0,x1,y0,y1,z0,z1,ne,ndim)  # generate the mesh grid
    NodeListCylinder = inflate_sphere(NodeList, x0, x1, y0, y1)                                 # inflate the sphere to a unit sphere
    q_tp, q_btm, C_uc = setboundaryCond(NodeList, ne, ndim, FunctionClass, nDof)

    mdl = def_model(ne=ne, NodeList=NodeList, IEN=IEN, IEN_top=IEN_top, IEN_btm=IEN_btm, ndim=ndim, nDof=nDof, FunctionClass=FunctionClass, ID=ID, Young=Float64(Young), ν=ν, cMat=cMat)
                
    state = "init"

    BorderNodes2D, NodeList2D = extract_borders(NodeListCylinder, CameraMatrix, BorderNodesList, state, ne)
    pi, qi = fit_curve(border=BorderNodes2D)

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

            BorderNodes2D, NodeList2D = extract_borders(NodeListCylinder, CameraMatrix, BorderNodesList, state)
            BorderNodes = [NodeListCylinder[:,SideBorders] NodeListCylinder[:,BottomBorders] NodeListCylinder[:,TopBorders]]
            # border = average_pts(BorderNodes2D, 2048/2)
            pi, qi = fit_curve(border=BorderNodes2D)

            push!(output, μ_tp)
            push!(fields, NodeListCylinder)
            push!(fields2D, NodeList2D)
            push!(borderfields2D, border)
            push!(splinep, pi)
            push!(splineq, qi) 

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

            BorderNodes2D, NodeList2D = extract_borders(NodeListCylinder, CameraMatrix, BorderNodesList, state)
            BorderNodes = [NodeListCylinder[:,SideBorders] NodeListCylinder[:,BottomBorders] NodeListCylinder[:,TopBorders]]
            pi, qi = fit_curve(border=BorderNodes2D)

            # store the solutions in a list
            push!(output, F_est[1])
            push!(fields, NodeListCylinder)
            push!(fields2D, NodeList2D)
            push!(borderfields2D, BorderNodes2D)
            push!(splinep, pi)
            push!(splineq, qi) 

            iter += 1
            next!(pr, showvalues = [(:iterations,iter),(:time,t)])
        end
    end

    if writeData
        fileName = "vtkFiles/squeeze_flow"
        write_scene(fileName, NodeList, IEN, ne, ndim, fields)
        animate_fields(fields=fields , IEN=IEN, BorderNodes2D=borderfields2D, fields2D=fields2D, p=splinep, q=splineq)
        writeCSV(csv_path, borderfields2D)
    end;

    return output, borderfields2D, splinep, splineq, mdl
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
    ndim = 3
    FunctionClass = "Q1"
    nDof = ndim  # number of degree of freedom per node
    β = 100
    CameraMatrix = [[8*2048/7.07, 0.0, 2048/2] [0.0, 8*1536/5.3, 1536/2] [0.0, 0.0, 1.0]]
    endTime = 20
    tSteps = 45
    Control = "displacement" # "force" or "displacement"
    writeData = false
    Young = 40
    ν = 0.4

    csv_path = "/home/soshala/SMEAR-PhD/SMEAR/Data/experiment_01/Results/contour_data"
   
    if Control == "force"
        cParam = -0.25*ones(tSteps)
    elseif Control == "displacement"
        μ_tp = 0.25
        cParam = -(μ_tp/tSteps)*ones(tSteps)
    end

    cMat = get_cMat("standard", Young = Young, ν = ν)

    μ_list, simborderfields, splinex, spliney = simulate(x0, x1, y0, y1, z0, z1, ne, Young, ν, ndim, FunctionClass, nDof, β, CameraMatrix, endTime, tSteps, Control, cParam; cMat, writeData=writeData, csv_path=csv_path)

end

test()