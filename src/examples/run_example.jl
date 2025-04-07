using LinearSolve
using SparseArrays
using IterativeSolvers

""" 
simulate_single_tstep(x0, x1, y0, y1, z0, z1, ne, Young, ν, ndim, FunctionClass, nDof, β, μ_tp, μ_btm; DENSE=false, CG=false)

Simulate the deformation of the mesh for a single time step using the linear elasticity model

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
- `μ_tp::Float64` : top boundary condition
- `μ_btm::Float64` : bottom boundary condition
- `DENSE::Bool` : use dense matrices
- `CG::Bool` : use conjugate gradient method to solve the system of equations
"""
function simulate_single_tstep(x0, x1, y0, y1, z0, z1, ne, c1, c2, ndim, FunctionClass, nDof, β, μ_tp, μ_btm; mode::String="lame", GRAD::Bool=false, DENSE::Bool=false, CG::Bool=false)

    cMat = get_cMat(mode, c1, c2)
    dcdλ = get_cMat(mode, 1.0 , 0.0)

    NodeList, IEN, ID, IEN_top, IEN_btm, BorderNodesList = meshgrid_cube(x0,x1,y0,y1,z0,z1,ne,ndim,FunctionClass=FunctionClass)  # generate the mesh grid
    NodeListCylinder = inflate_cylinder(NodeList, x0, x1, y0, y1) # inflate the sphere to a unit sphere
    
    mdl = def_model("linear_elasticity", ne=ne, NodeList=NodeListCylinder, IEN=IEN, IEN_top=IEN_top, IEN_btm=IEN_btm, ndim=ndim, nDof=nDof, ID = ID,
                        FunctionClass=FunctionClass, Young=Float64(c1), ν=c2, cMat=cMat, dcMatdλ=dcdλ)

    if DENSE == true
        q_tp, q_btm, C_uc = set_boundary_conditions_dense(mdl)
        K = assemble_system_dense(mdl)                   # assemble the stiffness matrix
        b = set_slip_conditions_dense(mdl)         # apply the neumann boundary conditions
    else
        q_tp, q_btm, C_uc = set_boundary_conditions(mdl)
        b = set_slip_conditions(mdl)         # apply the neumann boundary conditions
        if GRAD
            K, dKdλ = assemble_system(mdl, GRAD=true) 
            dKdβ = b             
        else
            K = assemble_system(mdl)                   # assemble the stiffness matrix
        end
        
    end
    
    q_d = (μ_btm*q_btm + μ_tp*q_tp)            # apply the Dirichlet boundary conditions
    K_bar = K + β*b
    C_T = transpose(C_uc)                      # transpose the constraint matrix
    K_free = C_T*K_bar*C_uc                    # extract the free part of the stiffness matrix

    if CG
        q_f = cg(K_free, C_T*(-K_bar*q_d))
    else
        q_f = K_free\(C_T*(-K_bar*q_d))         # solve the system of equations
    end

    q = q_d + C_uc*q_f                 # assemble the solution 
    q_out = hcat([q[ID[1,:]] q[ID[2,:]] q[ID[3,:]]])'    # update the nodal positions

    if GRAD && CG==false

        dqfdλ = -K_free\(C_T*dKdλ*C_uc*q_f + C_T*dKdλ*q_d)
        dqfdβ = -K_free\(C_T*dKdβ*C_uc*q_f + C_T*dKdβ*q_d)

        dqdλ = zeros(size(q_d)) + C_uc*dqfdλ
        dqdβ = zeros(size(q_d)) + C_uc*dqfdβ

        dqdλ_out = hcat(dqdλ[ID[1,:]], dqdλ[ID[2,:]], dqdλ[ID[3,:]])'
        dqdβ_out = hcat(dqdβ[ID[1,:]], dqdβ[ID[2,:]], dqdβ[ID[3,:]])'

        dqdθ_out = cat(dqdλ_out,dqdβ_out,dims=(3,3)) # concatenate the gradients in to a tensor
        return q_out, dqdθ_out, mdl
    else
        return q_out, mdl
    end
end

"""
    simulate_single_tstep_stokes(x0, x1, y0, y1, z0, z1, ne, ν, ndim, FunctionClass_u, FunctionClass_p, nDof_u, nDof_p, β, μu_tp, μu_btm, μu_side; DENSE=false)
    
Simulate the deformation of the mesh for a single time step using the Stokes model

# Arguments:
- `x0::Float64` : x-coordinate of the lower left corner of the domain
- `x1::Float64` : x-coordinate of the upper right corner of the domain
- `y0::Float64` : y-coordinate of the lower left corner of the domain
- `y1::Float64` : y-coordinate of the upper right corner of the domain
- `z0::Float64` : z-coordinate of the lower left corner of the domain
- `z1::Float64` : z-coordinate of the upper right corner of the domain
- `ne::Int` : number of elements
- `ν::Float64` : Poisson's ratio
- `ndim::Int` : number of dimensions
- `FunctionClass_u::String` : type of basis function for the velocity field
- `FunctionClass_p::String` : type of basis function for the pressure field
- `nDof_u::Int` : number of degree of freedom per node for the velocity field
- `nDof_p::Int` : number of degree of freedom per node for the pressure field
- `β::Float64` : friction parameter
- `μu_tp::Float64` : top boundary condition
- `μu_btm::Float64` : bottom boundary condition
- `μu_side::Float64` : side boundary condition
- `DENSE::Bool` : use dense matrix
"""
function simulate_single_tstep_stokes(x0, x1, y0, y1, z0, z1, ne, η, ndim, FunctionClass_u, FunctionClass_p, nDof_u, nDof_p, β, μu_tp, μu_btm, μu_side; GRAD::Bool=false, DENSE=false)
    
    NodeList_u, IEN_u, ID_u, IEN_u_top, IEN_u_btm, BorderNodesList_u = meshgrid_cube(x0,x1,y0,y1,z0,z1,ne,ndim,FunctionClass=FunctionClass_u)  # generate the mesh grid
    NodeListCylinder = inflate_cylinder(NodeList_u, x0, x1, y0, y1)
    q_tp, q_side, q_btm, C_uc = set_boundary_cond_stokes(NodeList_u, ne, ndim, FunctionClass_u, nDof_u)
    
    NodeList_p, IEN_p, ID_p, IEN_p_top, IEN_p_btm, BorderNodesList_p = meshgrid_cube(x0,x1,y0,y1,z0,z1,ne,ndim,FunctionClass=FunctionClass_p)  # generate the mesh grid
    
    mdl = def_model("stokes", ne=ne, NodeList=NodeListCylinder, IEN=IEN_u, IEN_top=IEN_u_top, IEN_btm=IEN_u_btm, ndim=ndim, nDof=nDof_u, FunctionClass=FunctionClass_u, 
                    ID=ID_u, IEN_2=IEN_p, IEN_2_top=IEN_p_top, IEN_2_btm=IEN_p_btm, nDof_2=nDof_p, FunctionClass_2=FunctionClass_p ) # define the model
                            
    if DENSE == true
        A_bar = assemble_system_A(mdl)                   # assemble the stiffness matrix
        B = assemble_system_B(mdl)                   # assemble the stiffness matrix
        b = apply_boundary_conditions_stokes(mdl)           # apply the neumann boundary conditions
    else
        A_bar = assemble_system_A(mdl)                   # assemble the stiffness matrix
        B = assemble_system_B(mdl)                   # assemble the stiffness matrix
        b = apply_boundary_conditions_stokes(mdl)           # apply the neumann boundary conditions
    end

    q_d = (μu_btm*q_btm + μu_tp*q_tp + μu_side*q_side)      # apply the Dirichlet boundary conditions

    if GRAD
        dAdη = A_bar
        dAdβ = b        
    end

    A = η*A_bar + β*b

    C_Tu = transpose(C_uc)           # transpose the constraint matrix

    A_free = C_Tu*A*C_uc        # extract the free part of the stiffness matrix
    B_free = C_Tu*B             # extract the free part of the stiffness matrix

    K_free = [A_free B_free; B_free' zeros(size(B_free,2),size(B_free,2))]     # assemble the system of equations

    r = [C_Tu*A*q_d; B'*q_d]    # assemble the system of equations
    sol = -K_free\r                 # solve the system of equations

    q_f = sol[1:size(A_free,1)]     # extract the free part of the solution
    p_f = sol[size(A_free,1)+1:end] # extract the free part of the solution 

    q = q_d + C_uc*q_f;                 # assemble the solution 
    p = p_f;

    q_out = [q[ID_u[1,:]] q[ID_u[2,:]] q[ID_u[3,:]]]'

    if GRAD
        dKdη = [C_Tu*dAdη*C_uc zeros(size(B_free)); zeros(size(B_free')) zeros(size(B,2),size(B,2))] # assemble the system of equations
        dKdβ = [C_Tu*dAdβ*C_uc zeros(size(B_free)); zeros(size(B_free')) zeros(size(B,2),size(B,2))] # assemble the system of equations
    
        drdη = [C_Tu*dAdη*q_d; zeros(size(B,2),size(q_d,2))] # solve the system of equations
        drdβ = [C_Tu*dAdβ*q_d; zeros(size(B,2),size(q_d,2))] # solve the system of equations

        dsoldη = -K_free\(drdη + dKdη*sol) # solve the system of equations
        dsoldβ = -K_free\(drdβ + dKdβ*sol) # solve the system of equations

        dqfdη = dsoldη[1:size(A_free,1)] # extract the free part of the solution
        dqfdβ = dsoldβ[1:size(A_free,1)] # extract the free part of the solution

        dpfdη = dsoldη[size(A_free,1)+1:end] # extract the free part of the solution 
        dpfdβ = dsoldβ[size(A_free,1)+1:end] # extract the free part of the solution

        dqdη = C_uc*dqfdη;              # assemble the solution
        dqdβ = C_uc*dqfdβ;              # assemble the solution

        dpdη = dpfdη;                  # assemble the solution
        dpdβ = dpfdβ;                  # assemble the solution

        dqdη_out = hcat(dqdη[ID_u[1,:]], dqdη[ID_u[2,:]], dqdη[ID_u[3,:]])'
        dqdβ_out = hcat(dqdβ[ID_u[1,:]], dqdβ[ID_u[2,:]], dqdβ[ID_u[3,:]])'

        dqdθ_out = cat(dqdη_out,dqdβ_out,dims=(3,3)) # concatenate the gradients in to a tensor
        return q_out, dqdθ_out, mdl
    else
        return q_out, mdl
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
function write_sim_data(x0, x1, y0, y1, z0, z1, ne, Youngtst, νtst, ndim, FunctionClass, nDof, β, CameraMatrix, endTime, tSteps, Control::String, filename::String; mode::String = "lame")
    writeData = true
    μ_list, simBorderPts, simBorderNodes, splinex, spliney, mdl, SurfacePt2D = test(x0, x1, y0, y1, z0, z1, ne, Youngtst, νtst, ndim, FunctionClass, nDof, β, CameraMatrix, endTime, tSteps, Control, writeData=writeData, filepath=filename, mode = mode)

    return simBorderPts, simBorderNodes, splinex, spliney
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
function test(x0, x1, y0, y1, z0, z1, ne, c1, c2, ndim::Int64, FunctionClass::String, nDof::Int64, β, CameraMatrix, endTime, tSteps, Control::String; 
              writeData::Bool=false, filepath=nothing, mode::String = "lame", SIDES::Bool=false)
    
    if writeData
        isnothing(filepath) || AssertionError("Please provide a filepath to write the data")
        set_file(filepath)
    end

    if Control == "force"
        cParam = -3*ones(round(Int,endTime/tSteps))
    elseif Control == "displacement"
        μ_tp = -0.03 
        # cParam = collect(range(start=tSteps,stop=endTime,step=tSteps))*μ_tp 
        cParam = μ_tp*ones(round(Int,endTime/tSteps))
    end

    cMat = get_cMat(mode, c1, c2) # E, ν or λ, μ
    dcdλ = get_cMat(mode, 1.0 , 0.0)
    # dcdβ = get_cMat(mode, 0.0 , 1.0)

    μ_list, gradList, borderPts2DList, splinep, splineq, mdl, SurfacePt2D = simulate(x0, x1, y0, y1, z0, z1, ne, c1, c2, ndim, FunctionClass, nDof, β, 
                                                                        CameraMatrix, endTime, tSteps, Control, cParam, cMat, writeData=writeData, 
                                                                        filepath=filepath, SIDES=SIDES, dcdλ=dcdλ)

    return μ_list, gradList, borderPts2DList, splinep, splineq, mdl, SurfacePt2D
end

"""
    compare(x0, x1, y0, y1, z0, z1, ne, Young, ν, ndim::Int64, FunctionClass::String, nDof::Int64, β, CameraMatrix, endTime, tSteps, Control::String, 
            mode::String, ObsData, PLOT::Bool=false, filepath=nothing)

Compare the simulation with the observation data

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
- `mode::String` : type of constitutive matrix
- `ObsData::Tuple{Array{Float64,2},Array{Float64,1},Array{Float64,1}}` : observation data
- `PLOT::Bool` : plot the data
- `filepath::String` : path to the file
"""
function compare(x0, x1, y0, y1, z0, z1, ne, Young, ν, ndim::Int64, FunctionClass::String, nDof::Int64, β, CameraMatrix, endTime, tSteps, Control::String, 
                mode::String, ObsData, SIDES::Bool=false, PLOT::Bool=false, filepath=nothing)

    μ_list, gradList, simBorderPts, splinex, spliney, mdl, SurfacePt2D = test(x0, x1, y0, y1, z0, z1, ne, Young, ν, ndim, FunctionClass, nDof, β, CameraMatrix, 
                                                                       endTime, tSteps, Control, mode=mode, SIDES=SIDES)    

    obsBorderPts = ObsData[1]
    splinexObs = ObsData[2]
    splineyObs = ObsData[3]
    
    # test the closest point function
    
    d_cp, pairs = closest_point(simBorderPts, obsBorderPts)
    
    if PLOT
        @assert !isnothing(filepath) "File path not provided"
        plot_matches(simBorderPts, splinex, spliney, splinexObs, splineyObs, pairs, string(filepath,"/Results/images/"))
    end
    d_h = 0
    return d_h, d_cp
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
function initialize_mesh_test(x0, x1, y0, y1, z0, z1, ne, ndim, FunctionClass, CameraMatrix, filepath=nothing, SIDES::Bool=false)

    isnothing(filepath) || AssertionError("Please provide a filepath to write the data")
    set_file(filepath)
    
    NodeList, IEN, ID, IEN_top, IEN_btm, BorderNodesList = meshgrid_cube(x0,x1,y0,y1,z0,z1,ne,ndim,FunctionClass=FunctionClass)  # generate the mesh grid
    NodeListCylinder = inflate_cylinder(NodeList, x0, x1, y0, y1)                                 # inflate the sphere to a unit sphere

    BorderPts2D, Nodes2D = extract_borders(NodeListCylinder, CameraMatrix, BorderNodesList, 2*ne+1, SIDES)
    pi, qi = fit_curve(border=BorderPts2D)

    SideBorders = BorderNodesList[1]
    BottomBorders = BorderNodesList[2]
    TopBorders = BorderNodesList[3]
        
    pos3D = AbstractArray[NodeListCylinder]                                                             # store the solution fields of the mesh in 3D
    pos2D = AbstractArray[Nodes2D]                                                                   # store the solution fields of the mesh in 2D
    surfaceNodesList = Float64[NodeList[:,SideBorders] NodeList[:,BottomBorders] NodeList[:,TopBorders]]  # store the solution fields of the surfaces in 3D
    borderPts2DList = AbstractArray[BorderPts2D]                                                               # store the solution fields of the surfaces in 2D
    splinep = AbstractArray[pi]                                                                            # store the x coordinates samples of the spline parameters of the border nodes
    splineq = AbstractArray[qi]                                                                            # store the y coordinates samples of the spline parameters of the border nodes 
    writeborderList = [vcat(pi', qi')]

    animate_fields(filepath = string(filepath,"/Results/images"), fields=pos3D , IEN=IEN, BorderNodes2D=borderPts2DList, fields2D=pos2D)
    write_contour_data(string(filepath,"/Results"), writeborderList)

    return borderPts2DList, surfaceNodesList, splinep, splineq
end

"""
    readData(filepath::String; NOISE::Bool=false, nProfile::Int64=1, nFactor=0)

Read the data from a file

# Arguments:
- `filepath::String` : path to the file
- `NOISE::Bool` : add noise to the data
- `nProfile::Int` : noise profile to be added to the data (1: uniform, 2: normal)
- `nFactor::Int` : noise level
"""
function readData(filepath::String)
    obsBorderPts, splinexObs, splineyObs, pd = read_csv(string(filepath,"/Results/contour_data"))   
    animate_fields(filepath = string(filepath,"/Results/cost"),pObs=splinexObs, qObs=splineyObs) # animate the fields
    plot(x->pdf(pd, x))
    savefig(string(filepath,"/Results"))
    return obsBorderPts, splinexObs, splineyObs
end