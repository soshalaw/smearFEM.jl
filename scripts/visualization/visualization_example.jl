using smearFEM

function main()

    # test case 
    r = 3.0
    h = 10.0
    ne = 4
    ndim = 3
    FunctionClass = "Q2"
    camera_matrix = [[8*2048/7.07, 0.0, 2048/2] [0.0, 8*1536/5.3, 1536/2] [0.0, 0.0, 1.0]]'
    camera_pose = 10*[0 -0.5 3]'   # camera position in mm

    filepath = "/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/sim_experiments/visualization_test/"
    set_file(filepath)
    SIDES = false

    initialize_mesh(r, h, ne, FunctionClass, camera_matrix, camera_pose, filepath, SIDES)
end    

main()