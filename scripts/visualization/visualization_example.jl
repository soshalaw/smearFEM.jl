using smearFEM

function main()

    # test case 
    r = 0.5
    h = 1
    ne = 4
    ndim = 3
    FunctionClass = "Q2"
    CameraMatrix = [[8*2048/7.07, 0.0, 2048/2] [0.0, 8*1536/5.3, 1536/2] [0.0, 0.0, 1.0]]'

    filepath = "/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/sim_experiments/Lame/test_experiment/"
    set_file(filepath)
    SIDES = false

    initialize_mesh_test(r, h, ne, FunctionClass, CameraMatrix, filepath, SIDES)
end    

main()