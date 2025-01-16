using LinearAlgebra
using ProgressMeter
using SparseArrays

using smearFEM
using Plots

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
    endTime = 30
    tSteps = 45
    noiseProfile = 1
    noiseLevel = 4
    plot_matches = true
    Control = "displacement" # "force" or "displacement"
    dev = 0.5
    mode = "lame" # "standard" or "lame"
    dateTime = Dates.now()
    filepathi = string("/home/soshala/SMEAR-PhD/SMEAR/Data/sim_experiments/cost_function_test/robusteness/",dateTime,"/")

    Youngtst = 30
    νtst = 0.4

    # Derived Lame constants from Young's modulus and Poisson ratio
    lambdatst = round(Youngtst*νtst/((νtst+1).*(-2*νtst+1)))
    mutst = round(Youngtst/(2*(νtst+1)))
    
    write_sim_data(x0, x1, y0, y1, z0, z1, ne, lambdatst, mutst, ndim, FunctionClass, nDof, β, CameraMatrix, endTime, tSteps, Control,
                    filepathi, mode=mode)
    
    ObsDataList, splinexObs, splineyObs = readData(filepathi, NOISE=false, nProfile=noiseProfile, nFactor=noiseLevel)
    ObsData = [ObsDataList, splinexObs, splineyObs]
    
    Young = lambdatst*(1-dev)
    ν = mutst*(1-dev)

    hcost, cpCost = compare(x0, x1, y0, y1, z0, z1, ne, Young, ν, ndim, FunctionClass, nDof, β, CameraMatrix, endTime, tSteps, Control, "lame", 
                    ObsData, plot_matches, filepathi)

    plot(cpCost, label="Height Cost")
    xlabel!("Time steps")
    ylabel!("Cost")
    savefig(string(filepathi,"/Results/cost/cost_cp.png"))

    plot(hcost, label="Closest Point Cost")
    xlabel!("Time steps")
    ylabel!("Cost")
    savefig(string(filepathi,"/Results/cost/cost_height.png"))

    params = Dict("Paramter type" => "Lame", "λ" => lambdatst, "μ" => mutst, "β" => β, "Control" => Control, "Noise Level" => noiseLevel, 
                  "Noise Profile" => noiseProfile)

    write_json(string(filepathi,"/Results/params"), params) 
end

main()