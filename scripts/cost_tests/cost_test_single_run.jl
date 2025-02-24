using LinearAlgebra
using ProgressMeter
using SparseArrays

using smearFEM
using StatsPlots
using Distributions
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
    endTime = 30
    tSteps = 45
    noiseLevel = 1
    plot_matches = true
    sides_only = false
    Control = "displacement" # "force" or "displacement"
    dev = 0.4
    mode = "lame" # "standard" or "lame"
    dateTime = Dates.now()
    filepathi = string("/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/sim_experiments/cost_function_test/robusteness/",Date(dateTime),"/",Time(dateTime),"/")

    Youngtst = 30
    νtst = 0.4
    nSamples = 10

    # Derived Lame constants from Young's modulus and Poisson ratio
    lambdatst = round(Youngtst*νtst/((νtst+1).*(-2*νtst+1)))
    mutst = round(Youngtst/(2*(νtst+1)))
    
    write_sim_data(x0, x1, y0, y1, z0, z1, ne, lambdatst, mutst, ndim, FunctionClass, nDof, β, CameraMatrix, endTime, tSteps, Control,
                    filepathi, mode=mode)
    
    ObsDataList, splinexObs, splineyObs = read_csv(string(filepathi,"/Results/contour_data"))  
    
    λ = lambdatst*(1-dev)
    μ = mutst*(1-dev)
    cSample = zeros(tSteps+1,nSamples)
    for n = 1:nSamples

        nScene, nSplinex, nSpliney, pd = add_noise(ObsDataList, nFactor=noiseLevel)
        ObsData = [nScene, nSplinex, nSpliney]

        if n == 1
            plot(x->pdf(pd, x))
            savefig(string(filepathi,"/Results"))
        end

        hcost, cpCost = compare(x0, x1, y0, y1, z0, z1, ne, λ, μ, ndim, FunctionClass, nDof, β, CameraMatrix, endTime, tSteps, Control, "lame", 
                    ObsData, sides_only, plot_matches, filepathi)

        cSample[:,n] = cpCost 
    end

    if nSamples == 1
        cost = cSample/nSamples
        plot(cost, label="Cost") 
    else
        errorline(cSample, errorstyle=:ribbon, label="Height Cost")
    end
    xlabel!("Time steps")
    ylabel!("Cost")
    savefig(string(filepathi,"/Results/cost/cost_cp.png"))

    params = Dict("Paramter type" => "Lame", "λ" => lambdatst, "μ" => mutst, "β" => β, "Control" => Control, "Noise Level" => noiseLevel, 
                    "Sample_no" => nSamples)

    write_json(string(filepathi,"/Results/params"), params) 
end

main()