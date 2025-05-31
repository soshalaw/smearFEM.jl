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
    ne = 4
    ndim = 3
    FunctionClass_u = "Q2"
    nDof_u = ndim  # number of degree of freedom per node
    FunctionClass_p = "Q1"
    nDof_p = 1  # number of degree of freedom per node

    CameraMatrix = [[8*2048/7.07, 0.0, 2048/2] [0.0, 8*1536/5.3, 1536/2] [0.0, 0.0, 1.0]]'

    endTime = 10
    steps = 10
    tSteps = endTime/steps
    time = 0:tSteps:endTime

    noiseLevel = 0
    PLOT_MATCHES = true
    SIDES = false
    Control = "velocity" # "force" or "displacement"
    dev = 0.4
    dateTime = Dates.now()
    homeDir = homedir()
    filepathi = string(homeDir,"/SMEAR-PhD/SMEAR-DataFiles/Data/sim_experiments/cost_function_test/robusteness/",Date(dateTime),"/",Time(dateTime),"/")

    βtst = 100
    ηtst = 20
    nSamples = 1
    
    write_sim_data_stokes(x0, x1, y0, y1, z0, z1, ne, ηtst, ndim, FunctionClass_u, nDof_u, FunctionClass_p, nDof_p, βtst, CameraMatrix, endTime, tSteps, Control,
    filepathi)
    
    ObsDataList, splinexObs, splineyObs = read_csv(string(filepathi,"/Results/contour_data"))  
    nScene, nSplinex, nSpliney, pd = add_noise(ObsDataList, nFactor=noiseLevel)
    
    ObsData = [nScene, nSplinex, nSpliney]
    
    obsBorderPts = ObsData[1]
    
    η = ηtst*(1-dev)
    β = βtst*(1-dev)
    cSample = zeros(steps+1,nSamples)

    for n = 1:nSamples
        nScene, nSplinex, nSpliney, pd = add_noise(ObsDataList, nFactor=noiseLevel)
        ObsData = [nScene, nSplinex, nSpliney]

        if n == 1
            plot(x->pdf(pd, x))
            savefig(string(filepathi,"/Results"))
        end

        # hcost, cpCost = compare_stokes(x0, x1, y0, y1, z0, z1, ne, η, ndim, FunctionClass_u, nDof_u, FunctionClass_p, nDof_p, β, CameraMatrix, endTime, tSteps, Control, 
        #             ObsData, SIDES, PLOT_MATCHES, filepathi)

        μ_list, gradList, simBorderPts, splinex, spliney, mdl = test_stokes(x0, x1, y0, y1, z0, z1, ne, η, ndim, FunctionClass_u, nDof_u, FunctionClass_p, 
                                                                            nDof_p, β, CameraMatrix, endTime, tSteps, Control, SIDES=SIDES)

        # test the closest point function
        d, ∂d, ∂2d, pairs = closest_point(simBorderPts, obsBorderPts, gradList)

        println("sum ", sum(d))
        cSample[:,n] = d 
    end

    if nSamples == 1
        cost = cSample/nSamples
        plot(time, cost, label="Cost") 
    else
        errorline(cSample, errorstyle=:ribbon, label="Cost")
    end
    xlabel!("Time steps")
    ylabel!("Cost")
    savefig(string(filepathi,"/Results/cost/cost_cp.png"))

    params = Dict("Paramter type" => "Lame", "η" => ηtst, "β" => β, "Control" => Control, "Noise Level" => noiseLevel, 
                    "Sample_no" => nSamples)

    write_json(string(filepathi,"/Results/params"), params) 
end

main()