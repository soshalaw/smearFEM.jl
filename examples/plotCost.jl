using LinearAlgebra
using ProgressMeter
using SparseArrays

using smearFEM
using Plots

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
    endTime = 30
    tSteps = 45
    
    Control = "displacement" # "force" or "displacement"

    filepath = "/home/soshala/SMEAR-PhD/SMEAR/Data/sim_experiments/Lame/"

    mode = "standard" # "standard" or "lame"

    Youngtst = 30
    νtst = 0.4

    filepathi = string(filepath,"experiment_single_run")
    write_sim_data(x0, x1, y0, y1, z0, z1, ne, Youngtst, νtst, ndim, FunctionClass, nDof, β, CameraMatrix, endTime, tSteps, Control,filepathi, mode=mode)
    Young = Youngtst-5
    ν = νtst - 0.1
    hcost, cpCost = test(x0, x1, y0, y1, z0, z1, ne, Young, ν, ndim, FunctionClass, nDof, β, CameraMatrix, endTime, tSteps, Control, filepath=filepathi, mode=mode)
    
    # λ = 34.0
    # μ = 14.0
    # write_sim_data(x0, x1, y0, y1, z0, z1, ne, λ, μ, ndim, FunctionClass, nDof, β, CameraMatrix, endTime, tSteps, Control,filepathi, mode=mode)
    # λ = λ - 1
    # μ = μ - 1
    # hcost, cpCost = test(x0, x1, y0, y1, z0, z1, ne, λ, μ, ndim, FunctionClass, nDof, β, CameraMatrix, endTime, tSteps, Control, filepath=filepathi, mode=mode)
  
    plot( cpCost, label="Height Cost")
    xlabel!("Time steps")
    ylabel!("Cost")
    savefig(string(filepathi,"/Results/cost/cost_cp.png"))

    plot(hcost, label="Closest Point Cost")
    xlabel!("Time steps")
    ylabel!("Cost")
    savefig(string(filepathi,"/Results/cost/cost_height.png"))
end

main()