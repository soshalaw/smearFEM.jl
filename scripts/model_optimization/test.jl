using LinearAlgebra
using ProgressMeter
using SparseArrays

using smearFEM
using StatsPlots
using Distributions
using Dates

function fit()
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
    CameraMatrix = [[8*2048/7.07, 0.0, 2048/2] [0.0, 8*1536/5.3, 1536/2] [0.0, 0.0, 1.0]]'
    endTime = 15
    tSteps = 45
    dateTime = Dates.now()
    sampleNo = 21
    dev = 0.3

    noiseLevel = 0
    SIDES = false
    Control = "displacement" # "force" or "displacement"
    mode = "lame" # "standard" or "lame"
    filepathi = string("/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/sim_experiments/cost_function_test/optimization/test")

    β = 100
    Youngtst = 30
    νtst = 0.4

    # Derived Lame constants from Young's modulus and Poisson ratio
    lambdatst = round(Youngtst*νtst/((νtst+1)*(-2*νtst+1)))
    mutst = round(Youngtst/(2*(νtst+1)))

    println(lambdatst)

    # write_sim_data(x0, x1, y0, y1, z0, z1, ne, lambdatst, mutst, ndim, FunctionClass, nDof, β, CameraMatrix, endTime, tSteps, Control,
    #                 filepathi, mode=mode)
    
    ObsDataList, splinexObs, splineyObs = read_csv(string(filepathi,"/Results/contour_data"))  
    nScene, nSplinex, nSpliney, pd = add_noise(ObsDataList, nFactor=noiseLevel)
    
    ObsData = [nScene, nSplinex, nSpliney]
    
    obsBorderPts = ObsData[1]

    dev_λ = lambdatst*dev
    lambdatstStart = lambdatst-dev_λ
    λ = lambdatstStart

    μ_list, gradList, simBorderPts, splinex, spliney, mdl = test(x0, x1, y0, y1, z0, z1, ne, lambdatstStart, mutst, ndim, FunctionClass, nDof, β, CameraMatrix, 
                                                                endTime, tSteps, Control, mode=mode, SIDES=SIDES)
    
    # test the closest point function
    d, ∂d, ∂2d, pairs = closest_point(simBorderPts, obsBorderPts, gradList)

    sdinit = sum(d)
    println("cost: ",sdinit)
    ratio = 1
    lambdanList = Float64[λ]
    costList = Float64[sdinit]
    while ratio > 0.01

        t∂2d = 0
        t∂d = 0
        
        sd = 1:length(∂d)

        for i in sd
            t∂2d = t∂2d + ∂2d[i][1,1]
            t∂d = t∂d + ∂d[i][1,1]
        end

        p = t∂d/t∂2d
        α = 0.5

        println("old lambda: ", λ)

        λ = λ + α*p

        println("new lambda: ", λ)

        μ_list, gradList, simBorderPts, splinex, spliney, mdl = test(x0, x1, y0, y1, z0, z1, ne, λ, mutst, ndim, FunctionClass, nDof, β, CameraMatrix, 
                                                                endTime, tSteps, Control, mode=mode, SIDES=SIDES)

        d, ∂d, ∂2d, pairs = closest_point(simBorderPts, obsBorderPts, gradList)

        sd = sum(d)
        ratio = sd/sdinit

        push!(lambdanList,λ)
        push!(costList,sd)
        println("cost ratios: ", ratio)
        println("cost: ",sd)
    end

    sampleNo = 21
    lambdaList = collect(range(lambdatstStart, stop=lambdatst+dev_λ, length=sampleNo))
    println(lambdaList)
    cpCostListλ = zeros(size(lambdaList))
        
    λ_iter = 1:size(lambdaList,1)
    for i in λ_iter
        μ_list, gradList, simBorderPts, splinex, spliney, mdl = test(x0, x1, y0, y1, z0, z1, ne, lambdaList[i], mutst, ndim, FunctionClass, nDof, β, CameraMatrix, 
                                                                    endTime, tSteps, Control, mode=mode, SIDES=SIDES)

        # test the closest point function
        d_cp, pairs = closest_point(simBorderPts, obsBorderPts) 

        cpCostListλ[i] = sum(d_cp)
        println(cpCostListλ[i])
    end

    Plots.plot(lambdaList, cpCostListλ, label="Closest Point Cost", marker=1, dpi=400)
    Plots.plot!(lambdanList,costList,label="Update steps", marker=1, dpi=400)
    Plots.vline!([lambdatst],linestyle=:dash,linecolor=:grey)
    Plots.xlabel!("λ")
    Plots.ylabel!("Error")
    Plots.savefig(string(filepathi,"/Results/cost/cost_cp_lbd.png"))
end

fit()

