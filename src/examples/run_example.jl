using LinearSolve
using SparseArrays
using IterativeSolvers

"""
    simulate_single_tstep(r, h, ne, c1, c2, ndim, FunctionClass, nDof, β, μ_tp, μ_btm; mode="lame", GRAD=false, DENSE=false, CG=false)

Simulates the deformation of the mesh for a single time step using the linear elasticity model.

# Arguments
- `r::Number`: Radius of the cylinder.
- `h::Number`: Height of the cylinder.
- `ne::Int64`: Number of elements.
- `c1::Number`: First material parameter (e.g., Young's modulus or Lame's first parameter).
- `c2::Number`: Second material parameter (e.g., Poisson's ratio or Lame's second parameter).
- `ndim::Int64`: Number of dimensions.
- `FunctionClass::String`: Type of basis function.
- `nDof::Int64`: Number of degrees of freedom per node.
- `β::Number`: Friction parameter.
- `μ_tp::Number`: Top boundary condition.
- `μ_btm::Number`: Bottom boundary condition.

# Keyword Arguments
- `mode::String`: Type of constitutive matrix ("lame" by default).
- `GRAD::Bool`: Whether to compute gradients (default: `false`).
- `DENSE::Bool`: Whether to use dense matrices (default: `false`).
- `CG::Bool`: Whether to use the conjugate gradient method (default: `false`).

# Returns
- `q_out::Matrix{Float64}`: Displacement fields.
- `dqdθ_out::Matrix{Float64}` (if `GRAD=true`): Gradients wrt to model parameters at the nodal points.
- `mdl::model`: Model object containing mesh and material properties.
"""
function simulate_single_tstep(r::Number, h::Number, ne::Int64, c1::Number, c2::Number, ndim::Int64, FunctionClass::String, nDof::Int64, β::Number, μ_tp::Number, 
                            μ_btm::Number; mode::String="lame", GRAD::Bool=false, DENSE::Bool=false, CG::Bool=false)

    cMat = get_cMat(mode, c1, c2)
    dcdλ = get_cMat(mode, 1.0 , 0.0)

    NodeList, IEN, ID, IEN_top, IEN_btm, BorderNodesList = meshgrid_cylinder(r, h, ne, FunctionClass=FunctionClass)  # generate the mesh grid
    
    mdl = def_model("linear_elasticity", ne=ne, NodeList=NodeList, IEN=IEN, IEN_top=IEN_top, IEN_btm=IEN_btm, ndim=ndim, nDof=nDof, ID = ID,
                        FunctionClass=FunctionClass, θ1=c1, θ2=c2, cMat=cMat, dcMatdθ1=dcdλ)

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
    simulate_single_tstep_stokes(r, h, ne, η, ndim, FunctionClass_u, FunctionClass_p, nDof_u, nDof_p, β, μu_tp, μu_btm, μu_side; GRAD=false, DENSE=false)

Simulates the deformation of the mesh for a single time step using the Stokes model.

# Arguments
- `r::Number`: Radius of the cylinder.
- `h::Number`: Height of the cylinder.
- `ne::Int64`: Number of elements.
- `η::Number`: Viscosity parameter.
- `ndim::Int64`: Number of dimensions.
- `FunctionClass_u::String`: Basis function for the velocity field.
- `FunctionClass_p::String`: Basis function for the pressure field.
- `nDof_u::Int64`: Degrees of freedom per node for the velocity field.
- `nDof_p::Int64`: Degrees of freedom per node for the pressure field.
- `β::Number`: Friction parameter.
- `μu_tp::Number`: Top boundary condition.
- `μu_btm::Number`: Bottom boundary condition.
- `μu_side::Number`: Side boundary condition.

# Keyword Arguments
- `GRAD::Bool`: Whether to compute gradients (default: `false`).
- `DENSE::Bool`: Whether to use dense matrices (default: `false`).

# Returns
- `q_out::Matrix{Float64}`: Displacement fields.
- `dqdθ_out::Matrix{Float64}` (if `GRAD=true`): Gradients wrt to model parameters at the nodal points.
- `mdl::model`: Model object containing mesh and material properties.
"""
function simulate_single_tstep_stokes(r::Number, h::Number, ne::Int64, η::Number, ndim::Int64, FunctionClass_u::String, FunctionClass_p::String, nDof_u::Int64,
                                    nDof_p::Int64, β::Number, μu_tp::Number, μu_btm::Number, μu_side::Number; FunctionClass_x::String=FunctionClass_u, GRAD::Bool=false, DENSE::Bool=false)
    
    
    filePath = "/home/soshala/SMEAR-PhD/smear-modules/smearFEM.jl/cylindergen"

    mesh_x = meshgrid_cylinder(r, h, ne, FunctionClass=FunctionClass_x, filePath=filePath)  # generate the mesh grid for geometry
    mesh_u = meshgrid_cylinder(r, h, ne, FunctionClass=FunctionClass_u)  # generate the mesh grid
    mesh_p = meshgrid_cylinder(r, h, ne, FunctionClass=FunctionClass_p)  # generate the mesh grid

    mdl = Stokes(ndim=ndim, mesh_x=mesh_x, mesh_u=mesh_u, nDof_u=nDof_u, mesh_p=mesh_p, nDof_p=nDof_p, η=[η])

    ID_u = mdl.mesh_u.ID
    
    q_tp, q_side, q_btm, C_uc = set_boundary_cond(mdl)

    if DENSE == true
        A_bar = assemble_system_A(mdl)               # assemble the stiffness matrix
        B = assemble_system_B(mdl)                   # assemble the stiffness matrix
        b = apply_boundary_conditions(mdl)           # apply the neumann boundary conditions
    else
        A_bar = assemble_system_A(mdl)               # assemble the stiffness matrix
        B = assemble_system_B(mdl)                   # assemble the stiffness matrix
        b = apply_boundary_conditions(mdl)           # apply the neumann boundary conditions
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
    sol = -K_free\Matrix(r)                 # solve the system of equations

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
    write_sim_data(r, h, ne, Young, ν, ndim, FunctionClass, nDof, β, CameraMatrix, endTime, tSteps, Control, filename; mode="lame")

Initializes the simulation and writes the data to a file.

# Arguments
- `r::Number`: Radius of the cylinder.
- `h::Number`: Height of the cylinder.
- `ne::Int64`: Number of elements.
- `Young::Number`: Young's modulus.
- `ν::Number`: Poisson's ratio.
- `ndim::Int64`: Number of dimensions.
- `FunctionClass::String`: Type of basis function.
- `nDof::Int64`: Number of degrees of freedom per node.
- `β::Number`: Friction parameter.
- `CameraMatrix::AbstractMatrix{Float64}`: Camera matrix.
- `endTime::Number`: End time of the simulation.
- `tSteps::Number`: Number of time steps.
- `Control::String`: Type of control ("force" or "displacement").
- `filename::String`: Path to the output file.

# Keyword Arguments
- `mode::String`: Type of constitutive matrix ("lame" by default).

# Returns
None.
"""
function write_sim_data(_model::AbstractModel, _scene::AbstractScenario, camera_matrix::AbstractMatrix{Float64}, camera_pose::AbstractMatrix{Float64}, 
                        filepath::String="nothing"; ANIMATE=true, WRITECONTOUR=true, RENDER=true, WRITEVTK=true)

    # copy the model and scene to avoid modifying the original objects
    model::AbstractModel = deepcopy(_model)
    scene::AbstractScenario = deepcopy(_scene)
    
    # set the model parameters
    conditions = Conditions(ANIMATE=ANIMATE, WRITECONTOUR=WRITECONTOUR, RENDER=RENDER, WRITEVTK=WRITEVTK, camera_matrix=camera_matrix, camera_pose=camera_pose, filepath=filepath)
    
    # remove the previous results folder if it exists
    if isdir(string(filepath))
        @info "Removing previous results folder: $filepath"
        rm(string(filepath), recursive=true, force=true) # remove the previous results folder if it exists
    end

    
    h_, gradList, borderPts2DList, fields, pos3D, pos2D, _ = simulate(model, scene, conditions) # run the simulation
    
    h = get_height(h_, model.mesh_u.h) # get the mesh height with time
    
    params = Dict("r"=>model.mesh_u.r, "h"=>model.mesh_u.h, "η" => model.η, "β" => scene.β, "camera_matrix" => conditions.camera_matrix, "camera_pose" => conditions.camera_pose, 
                    "control_type"=>scene.control, "cParam"=>scene.cParam, "simulation_time" => scene.sim_time, "time_steps" => scene.t_steps, 
                    "viscosity_type"=>scene.viscosity_type)

    # write results to files
    write_json(string(filepath,"/data/sim_params"), params)
    write_csv(string(filepath,"/data/h"), h)
    write_data(string(conditions.filepath,"/data/sim_data/2D_surface_points"), pos2D)
    write_data(string(filepath,"/data/sim_data/node_points"), pos3D)
    write_data(string(filepath,"/data/sim_data/displacement_fields"), fields)
    write_data(string(filepath,"/data/sim_data/2D_border_points"), borderPts2DList)

    @info "Data written to $filepath"

end

"""
    test(r, h, ne, c1, c2, ndim, FunctionClass, nDof, β, CameraMatrix, endTime, tSteps, Control; writeData=false, filepath="nothing", mode="lame", SIDES=false)

Runs a test simulation.

# Arguments
- `r::Number`: Radius of the cylinder.
- `h::Number`: Height of the cylinder.
- `ne::Int64`: Number of elements.
- `c1::Number`: First material parameter.
- `c2::Number`: Second material parameter.
- `ndim::Int64`: Number of dimensions.
- `FunctionClass::String`: Type of basis function.
- `nDof::Int64`: Number of degrees of freedom per node.
- `β::Number`: Friction parameter.
- `CameraMatrix::AbstractMatrix{Float64}`: Camera matrix.
- `endTime::Number`: End time of the simulation.
- `tSteps::Number`: Number of time steps.
- `Control::String`: Type of control ("force" or "displacement").

# Keyword Arguments
- `writeData::Bool`: Whether to write data to a file (default: `false`).
- `filepath::String`: Path to the output file (default: `"nothing"`).
- `mode::String`: Type of constitutive matrix ("lame" by default).
- `SIDES::Bool`: Whether to include side boundaries (default: `false`).

# Returns
- `μ_list::Vector{Float64}`: Force applied on the top boundary in the force constrolled scenario / Force applied on the top boundary in the displacement controlled scenario.
- `gradList::Vector{Matrix{Float64}}`: List of gradients of the border points at each timestep.
- `borderPts2DList::Vector{Matrix{Float64}}`: List of 2D border points at each timestep.
- `splinep::Vector{Vector{Float64}}`: x coordinates of the border observation at each timestep interpolated.
- `splineq::Vector{Vector{Float64}}`: y coordinates of the border observation at each timestep interpolated.
- `mdl::model`: Model object containing mesh and material properties.
- `SurfacePt2D::Vector{Matrix{Float64}}`: List of surface nodes of the mesh at each timestep.
"""
function test(r::Number, h::Number, ne::Int64, c1::Float64, c2::Float64, ndim::Int64, FunctionClass::String, nDof::Int64, β::Number, CameraMatrix::AbstractMatrix{Float64}, 
            endTime::Number, tSteps::Number, Control::String; writeData::Bool=false, filepath::String="nothing", mode::String = "lame", SIDES::Bool=false)
    
    if writeData
        isnothing(filepath) || AssertionError("Please provide a filepath to write the data")
        set_file(filepath)
    end

    time = collect(Float64, range(0, stop=endTime, step=tSteps))
    len_t = length(time)

    if Control == "force"
        cParam = -3*ones(len_t)
    elseif Control == "displacement"
        μ_tp = -0.1
        cParam = μ_tp*ones(len_t)
    end

    cMat = get_cMat(c1, c2, type=mode) # E, ν or λ, μ
    dcdλ = get_cMat(1.0 , 0.0, type=mode)
    dcdβ = get_cMat(0.0 , 1.0, type=mode)

    dcdθl = cat(dcdλ, dcdβ, dims=3)

    μ_list, gradList, borderPts2DList, splinep, splineq, mdl, SurfacePt2D = simulate(r, h, ne, c1, c2, ndim, FunctionClass, nDof, β, CameraMatrix, time, 
                                                                                Control, cParam, cMat, writeData=writeData, filepath=filepath, SIDES=SIDES, 
                                                                                dcdλ=dcdλ)

    return μ_list, gradList, borderPts2DList, splinep, splineq, mdl, SurfacePt2D
end

"""
    compare(r, h, ne, Young, ν, ndim, FunctionClass, nDof, β, CameraMatrix, endTime, tSteps, Control, mode, ObsData; SIDES=false, PLOT=false, filepath="nothing")

Compares the simulation results with observation data.

# Arguments
- `r::Number`: Radius of the cylinder.
- `h::Number`: Height of the cylinder.
- `ne::Int64`: Number of elements.
- `Young::Number`: Young's modulus.
- `ν::Number`: Poisson's ratio.
- `ndim::Int64`: Number of dimensions.
- `FunctionClass::String`: Type of basis function.
- `nDof::Int64`: Number of degrees of freedom per node.
- `β::Number`: Friction parameter.
- `CameraMatrix::AbstractMatrix{Float64}`: Camera matrix.
- `endTime::Number`: End time of the simulation.
- `tSteps::Number`: Number of time steps.
- `Control::String`: Type of control ("force" or "displacement").
- `mode::String`: Type of constitutive matrix.
- `ObsData::AbstractArray`: Observation data.

# Keyword Arguments
- `SIDES::Bool`: Whether to include side boundaries (default: `false`).
- `PLOT::Bool`: Whether to generate plots (default: `false`).
- `filepath::String`: Path to the output file (default: `"nothing"`).

# Returns
- `d_cp::Vector{Float64}`: Closest point distances between borders at each timestep.
"""
function compare(r::Number, h::Number, ne::Int64, Young::Number, ν::Number, ndim::Int64, FunctionClass::String, nDof::Int64, β::Number, CameraMatrix::AbstractMatrix{Float64}, 
                endTime::Number, tSteps::Number, Control::String, mode::String, ObsData::AbstractArray, SIDES::Bool=false, PLOT::Bool=false, 
                filepath::String="nothing")

    μ_list, gradList, simBorderPts, splinex, spliney, mdl, SurfacePt2D = test(r, h, ne, Young, ν, ndim, FunctionClass, nDof, β, CameraMatrix, 
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
    initialize_mesh_test(r, h, ne, FunctionClass, CameraMatrix; filepath="nothing", SIDES=false)

Initializes the mesh and writes the data to a file.

# Arguments
- `r::Number`: Radius of the cylinder.
- `h::Number`: Height of the cylinder.
- `ne::Int64`: Number of elements.
- `FunctionClass::String`: Type of basis function.
- `CameraMatrix::AbstractMatrix{Float64}`: Camera matrix.

# Keyword Arguments
- `filepath::String`: Path to the output file (default: `"nothing"`).
- `SIDES::Bool`: Whether to include side boundaries (default: `false`).

# Returns
- `borderPts2DList::Vector{Matrix{Float64}}`: List of 2D border points at each timestep.
- `surfaceNodesList::Vector{Matrix{Float64}}`: List of surface nodes of the mesh at each timestep.
- `splinep::Vector{Vector{Float64}}`: x coordinates of the border observation at each timestep interpolated.
- `splineq::Vector{Vector{Float64}}`: y coordinates of the border observation at each timestep interpolated.
"""
function initialize_mesh(r::Number, h::Number, ne::Int64, FunctionClass::String, camera_matrix::AbstractMatrix{Float64}, camera_pose::AbstractMatrix{Float64},filepath::String="nothing", 
                            SIDES::Bool=false)

    isnothing(filepath) || AssertionError("Please provide a filepath to write the data")
    set_file(filepath)
    
    mesh = meshgrid_cylinder(r, h, ne, FunctionClass=FunctionClass)

    BorderPts2D, SurfacePts2D = extract_borders(mesh.NodeList, camera_matrix, camera_pose, mesh.nNodes, BorderNodesList= mesh.side_nodes)
    pi, qi = fit_curve(border=BorderPts2D)
        
                                               # store the solution fields of the border nodes in 2D 
    pos3D = AbstractArray[mesh.NodeList]                                                             # store the solution fields of the mesh in 3D
    surface_pts_3D = [vcat(mesh.NodeList[:,mesh.top_nodes]', mesh.NodeList[:,mesh.bottom_nodes]', mesh.NodeList[:,mesh.side_nodes]')'] # store the solution fields of the mesh in 3D
    pos2D = AbstractArray[SurfacePts2D]                                                                   # store the solution fields of the mesh in 2D
    borderPts2DList = AbstractArray[BorderPts2D]                                                          # store the solution fields of the surfaces in 2D
    splinep = AbstractArray[pi]                                                                           # store the x coordinates samples of the spline parameters of the border nodes
    splineq = AbstractArray[qi]                                                                           # store the y coordinates samples of the spline parameters of the border nodes 
    writeborderList = [vcat(pi', qi')]

    animate_fields(filepath = string(filepath,"/Results/images"), SurfaceNodes3D=surface_pts_3D , IEN=mesh.IEN, BorderNodes2D=borderPts2DList, fields2D=pos2D)
    write_data(string(filepath,"/Results"), writeborderList)

    return borderPts2DList, pos2D, splinep, splineq
end

"""
    readData(filepath)

Reads the data from a file.

# Arguments
- `filepath::String`: Path to the file.

# Returns
- `obsBorderPts::Vector{Matrix{Float64}}`: list of observation data.
- `splinexObs::Vector{Vector{Float64}}`: x-coordinates of the border observation at each time step.
- `splineyObs::Vector{Vector{Float64}}`: y-coordinates of the border observation at each time step.
"""
function readData(filepath::String)
    obsBorderPts, splinexObs, splineyObs, pd = read_csv(string(filepath,"/Results/contour_data"))   
    animate_fields(filepath = string(filepath,"/Results/cost"),pObs=splinexObs, qObs=splineyObs) # animate the fields
    plot(x->pdf(pd, x))
    savefig(string(filepath,"/Results"))
    return obsBorderPts, xObs, yObs
end