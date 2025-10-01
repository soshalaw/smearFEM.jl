using Plots
using smearFEM
using FileIO, Images

function initialize_meshes()
    # Initialization code here

    gt_path = "/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/physical_experiments/ground_truth/5"

    t_obs = read_perception_data(string(gt_path, "/Results/sequence.hdf5"))

    ObsDataList, splinexObs, splineyObs = read_csv(string(gt_path, "/data/img_data/contour_data"))

    # test case 
    r::Float64 = 25  # radius of the cylinder in mm
    h::Float64 = 40  # height of the cylinder in mm
    ne = 20
    ndim = 3
    FunctionClass = "Q2"
    camera_matrix = [[2.39642674e+03, 0.0, 1.00429248e+03] [0.0, 2.40565353e+03, 7.57028161e+02] [0.0, 0.0, 1.0]]'
    obj_pose = Float64.(t_obs)
    # obj_pose = Float64.([-1.0 0.0 0.0 0.0; 0.0 0.0 -1.0 20.0; 0.0 -1.0 0.0 150.0; 0.0 0.0 0.0 1.0])

    filepath = "/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/sim_experiments/visualization_test/"
    set_file(filepath)
    SIDES = false

    obj_pose = get_pose(t_obs[:,:,1:200])  # Example usage of get_pose function
    obj_pose_ = zeros(Float64, 4,4)
    obj_pose_[1,1] = -1.0
    obj_pose_[2,3] = -1.0
    obj_pose_[3,2] = -1.0
    obj_pose_[1:3,4] = obj_pose[1:3,4]

    filter_frames(ObsDataList)

    borderPts2DList, pos2D, splinep, splineq = initialize_mesh(r, h, ne, FunctionClass, camera_matrix, obj_pose_, filepath, SIDES)

    img = load(string(gt_path, "/Results/contour_detection/00000.png"))

    # Plotting the meshes
    p = plot(img; axis=nothing, legend=false)

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