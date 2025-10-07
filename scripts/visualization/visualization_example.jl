using smearFEM

function main()

    # test case 
    scale = 50
    r::Float64 = 0.5*scale  # radius of the cylinder in mm
    h::Float64 = 1*scale  # height of the cylinder in mm
    ne = 8
    ndim = 3
    FunctionClass = "Q2"
    camera_matrix = [[8*2048/7.07, 0.0, 2048/2] [0.0, 8*1536/5.3, 1536/2] [0.0, 0.0, 1.0]]'
    obj_pose = scale*[0 -0.5 2]'   # camera position in mm

    filepath = "/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/sim_experiments/visualization_test/"
    set_file(filepath)
    SIDES = false

    initialize_mesh(r, h, ne, FunctionClass, camera_matrix, obj_pose, filepath, SIDES)
end    

main()