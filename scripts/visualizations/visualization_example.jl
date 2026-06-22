using smearFEM

function main()

    # test case 
    r::Float64 = 25.0  # radius of the cylinder in mm
    h::Float64 = 40.0  # height of the cylinder in mm
    ne = 3 # number of elements in the mesh for the ground truth
    ndim = 3
    element_shape = :Hex
    basis_order = 2
    camera_matrix = [[8*2048/7.07, 0.0, 2048/2] [0.0, 8*1536/5.3, 1536/2] [0.0, 0.0, 1.0]]'
    # obj_pose = Float64.([-1.0 0.0 0.0 0.0; 0.0 0.0 -1.0 20.0; 0.0 -1.0 0.0 150; 0.0 0.0 0.0 1.0])
    obj_pose = [150.0, 0.0, 20.0]

    filepath = "/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/sim_experiments/visualization_test/"
    set_file(filepath)
    geometry = :cube # :cylinder or :cube
    edge_radius = 3.0
    rot_angle_list = [0.0, 30.0, 60.0]

    initialize_mesh(r, h, ne, element_shape, basis_order, camera_matrix, obj_pose, rot_angle_list, geometry=geometry, filepath=filepath, edge_radius=edge_radius)
end    

main()