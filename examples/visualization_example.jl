using smearFEM


function main()

    # test case 
    x0 = 0
    x1 = 1
    y0 = 0
    y1 = 1
    z0 = 0
    z1 = 1
    ne = 8
    ndim = 3
    FunctionClass = "Q2"
    nDof = ndim  # number of degree of freedom per node
    β = 100
    CameraMatrix = [[8*2048/7.07, 0.0, 2048/2] [0.0, 8*1536/5.3, 1536/2] [0.0, 0.0, 1.0]]
    endTime = 30
    tSteps = 45
    
    Control = "displacement" # "force" or "displacement"

    filepath = "/home/soshala/SMEAR-PhD/SMEAR/Data/sim_experiments/Lame/test_experiment/"

    mode = "standard" # "standard" or "lame"

    initialize_mesh_test(x0, x1, y0, y1, z0, z1, ne, ndim, FunctionClass, CameraMatrix, filepath)
end    

main()