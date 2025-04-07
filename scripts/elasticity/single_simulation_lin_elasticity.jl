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
    ne = 3
    ndim = 3
    FunctionClass = "Q2"
    nDof = ndim  # number of degree of freedom per node
    β = 100
    CameraMatrix = [[8*2048/7.07, 0.0, 2048/2] [0.0, 8*1536/5.3, 1536/2] [0.0, 0.0, 1.0]]'


    endTime = 15
    steps = 15
    tSteps = endTime/steps

    Youngtst = 30
    νtst = 0.4
    λ = Youngtst*νtst/((1+νtst)*(1-2*νtst))
    μ = Youngtst/(2*(1+νtst))
    Control = "force" # "force" or "displacement"
    mode = "lame" # "standard" or "lame"
    dateTime = Dates.now()

    # filepath = string("/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/sim_experiments/single_simulation/fem_runs/Linear_Elasticity/",mode,"/",Control,"/",Date(dateTime),"/",Time(dateTime),"/")
    # filepath = string("/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/sim_experiments/single_simulation/fem_runs/Linear_Elasticity/",mode,"/",Control,"/",Date(dateTime),"/09:41:12.027/")
    filepath = string("/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/sim_experiments/single_simulation/fem_runs/Linear_Elasticity/standard/displacement/2025-03-07/09:41:12.027")
    write_sim_data(x0, x1, y0, y1, z0, z1, ne, λ, μ, ndim, FunctionClass, nDof, β, CameraMatrix, endTime, tSteps, Control,filepath)
end

main()