using Plots
using smearFEM
using FileIO, Images

function initialize_meshes()
    # Initialization code here
    gt_path = "/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/ground_truth/physical_data/5"
    t_obs = read_perception_data(string(gt_path, "/Results/sequence.hdf5"))

    _ObsDataList, _splinexObs, _splineyObs = read_csv(string(gt_path, "/data/sim_data/contour_data"))
    meta_data = read_json(string(gt_path, "/data/video_metadata.json"))

    frame_rate = round(Int, meta_data["frame_rate"])
    frame_width = meta_data["frame_width"]
    frame_height = meta_data["frame_height"]
    compression_frames = meta_data["compressed_frames"]

    sim_time_exp = 10.0  # seconds
    steps_exp = sim_time_exp * frame_rate
    println("steps exp: $steps_exp")
    t_steps_exp = 1/frame_rate
    println("time steps exp: $t_steps_exp")

    ndim::Int = 3
    nDof_p::Int = 1  # number of degree of freedom per node
    nDof_u::Int = ndim  # number of degree of freedom per node

    # test case 
    r::Float64 = 25  # radius of the cylinder in mm
    h::Float64 = 40  # height of the cylinder in mm
    F_ext::Float64 = 1*9.812*1e3 # force applied to the cylinder in N
    F = -F_ext*ones(Float64, round(Int, steps_exp)) # force applied to the cylinder in N
    ne_exp = 4
    FunctionClass_u = "Q2"
    FunctionClass_p = "Q1"
    FunctionClass_x = "Q2"
    control = "force"            # "force" or "velocity"
    gt_viscosity_type = "constant"  # "constant" or "bulk_viscosity"
    β_start::Float64 = 2.0       # initial guess for the inverse of permeability
    η_start::Float64 = 1.0     # initial guess for the viscosity

    camera_matrix = [[2.39642674e+03, 0.0, 1.00429248e+03] [0.0, 2.40565353e+03, 7.57028161e+02] [0.0, 0.0, 1.0]]'
    obj_pose = Float64.(t_obs)

    filepath = "/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/sim_experiments/visualization_test/"
    set_file(filepath)
    SIDES = false

    obj_pose = get_pose(t_obs)  # Example usage of get_pose function
    obj_pose_ = zeros(Float64, 4,4)
    obj_pose_[1,1] = -1.0
    obj_pose_[2,3] = -1.0
    obj_pose_[3,2] = -1.0
    obj_pose_[1:3,4] = obj_pose[1:3,4]

    ObsDataList = _ObsDataList[compression_frames]
    splinexObs = _splinexObs[compression_frames]
    splineyObs = _splineyObs[compression_frames]
    valid_frames, outlier_frames = detect_outlier_observations(ObsDataList)

    println("Valid frames: ", length(valid_frames), " out of ", length(ObsDataList))
    println("Outlier frames: ", length(outlier_frames))

    # borderPts2DList, pos2D, splinep, splineq = initialize_mesh(r, h, ne_exp, FunctionClass_u, camera_matrix, obj_pose, filepath, SIDES)
    model, scene = def_problem(r, h, ne_exp, η_start, ndim, FunctionClass_u, nDof_u, FunctionClass_p, nDof_p, FunctionClass_x, β_start, F, control, gt_viscosity_type, 
                        sim_time_exp, t_steps_exp)
    conditions = Conditions(camera_matrix=camera_matrix, obj_pose=obj_pose, SIDES=SIDES, filepath=filepath, ANIMATE=false)
    output, gradList, borderPts2DList, splinep, splineq, mdl, pos2D = simulate(model, scene, conditions) # run the simulation

    animate_fields(filepath = string(conditions.filepath,"/Results/images/"), BorderNodes2D=borderPts2DList, pObs=splinexObs[1:(round(Int,steps_exp)+1)], qObs=splineyObs[1:(round(Int,steps_exp)+1)])

    img = load(string(gt_path, "/Results/contour_detection/00000.png"))

    # Plotting the meshes
    p = plot(img; axis=nothing, legend=false)

    display(borderPts2DList[1])
    # Overlay curves
    plot!(p, splinexObs[1], splineyObs[1], label="", lw=2, color=:blue)
    plot!(p, borderPts2DList[1][1,:], borderPts2DList[1][2,:], label="", lw=2, color=:black)

    xlims!(p, 0, 2048)
    ylims!(p, 0, 1536)
    yflip!(p, true)  # if needed to match image coords
    title!(p, "Mesh Initialization")
    xlabel!(p, "X (px)")
    ylabel!(p, "Y (px)")
    display(p)
    savefig(p, "/home/soshala/SMEAR-PhD/temp/mesh_initialization.pdf")
end


initialize_meshes()