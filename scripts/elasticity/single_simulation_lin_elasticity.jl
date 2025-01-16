using smearFEM
using Dates

function main()
    # test case 
    x0 = 0
    x1 = 1
    y0 = 0
    y1 = 1
    z0 = 0
    z1 = 1
    ne = 4
    ndim = 3
    FunctionClass = "Q2"
    nDof = ndim  # number of degree of freedom per node
    β = 100
    CameraMatrix = [[8*2048/7.07, 0.0, 2048/2] [0.0, 8*1536/5.3, 1536/2] [0.0, 0.0, 1.0]]
    endTime = 15
    tSteps = 45
    Young = 30
    ν = 0.3
    Control = "displacement" # "force" or "displacement"
    mode = "standard" # "standard" or "lame"
    dateTime = Dates.now()

    filepath = string("/home/soshala/SMEAR-PhD/SMEAR/Data/sim_experiments/single_simulation/fem_runs/Linear_Elasticity/",mode,"/",Control,"/",Date(dateTime),"/",Time(dateTime),"/")
    
    write_sim_data(x0, x1, y0, y1, z0, z1, ne, Young, ν, ndim, FunctionClass, nDof, β, CameraMatrix, endTime, tSteps, Control,filepath)
end

main()