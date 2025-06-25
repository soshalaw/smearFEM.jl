using smearFEM

function main()

    # test case 
    scale = 100
    r = 0.2*scale  # radius of the cylinder in mm
    h = 0.5*scale  # height of the cylinder in mm
    ne = 4
    ndim = 3
    FunctionClass = "Q2"
    camera_matrix = [[8*2048/7.07, 0.0, 2048/2] [0.0, 8*1536/5.3, 1536/2] [0.0, 0.0, 1.0]]'
    camera_pose = scale*[0 -0.25 2]'   # camera position in mm

    filepath = "/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/sim_experiments/visualization_test/"
    set_file(filepath)
    SIDES = false

    initialize_mesh(r, h, ne, FunctionClass, camera_matrix, camera_pose, filepath, SIDES)
end    

main()