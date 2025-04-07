using LinearAlgebra
using ProgressMeter
using SparseArrays

using smearFEM
using StatsPlots
using Distributions
using Dates
using Plots

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
    endTime = 10
    steps = 10
    tSteps = endTime/steps
    dateTime = Dates.now()
    sampleNo = 21
    dev = 0.3

    noiseLevel = 0
    SIDES = false
    Control = "force" # "force" or "displacement"
    mode = "lame" # "standard" or "lame"
    filepathi = string("/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/sim_experiments/cost_function_test/optimization/",Control,"/test1")

    β = 100
    Youngtst = 30
    νtst = 0.4

    # Derived Lame constants from Young's modulus and Poisson ratio
    lambdatst = round(Youngtst*νtst/((νtst+1)*(-2*νtst+1)))
    mutst = round(Youngtst/(2*(νtst+1)))

    println(lambdatst)

    write_sim_data(x0, x1, y0, y1, z0, z1, ne, lambdatst, mutst, ndim, FunctionClass, nDof, β, CameraMatrix, endTime, tSteps, Control,
                    filepathi, mode=mode)
    
    ObsDataList, splinexObs, splineyObs = read_csv(string(filepathi,"/Results/contour_data"))  
    nScene, nSplinex, nSpliney, pd = add_noise(ObsDataList, nFactor=noiseLevel)
    
    ObsData = [nScene, nSplinex, nSpliney]
    
    obsBorderPts = ObsData[1]

    dev_λ = lambdatst*dev
    dev_β = β*dev
    
    lambdatstStart = lambdatst-dev_λ
    βStart = β-dev_β

    θ = [lambdatstStart; βStart]

    μ_list, gradList, simBorderPts, splinex, spliney, mdl = test(x0, x1, y0, y1, z0, z1, ne, θ[1], mutst, ndim, FunctionClass, nDof, θ[2], CameraMatrix, 
                                                                endTime, tSteps, Control, mode=mode, SIDES=SIDES)
    
    # test the closest point function
    d, ∂d, ∂2d, pairs = closest_point(simBorderPts, obsBorderPts, gradList)
    
    λ = θ[1]
    β = θ[2]

    totdinit = sum(d)
    println("Initial cost: ", totdinit)
    println("initial lambda: ", λ, " initial beta: ", β)
    ratio = 1

    lambdapList = Float64[θ[1]]
    βpList = Float64[θ[2]]
    costList = Float64[totdinit]
    iterList = Int64[0]

    iter = 0

    while ratio > 0.005
        t∂2d = zeros(size(∂2d[1]))
        t∂d = zeros(size(∂d[1]))
        
        szd = 1:length(d)

        for i in szd
            t∂2d = t∂2d + ∂2d[i]
            t∂d = t∂d + ∂d[i]
        end
        # println("grad d: ", t∂d)
        # println("grad2 d: ", t∂2d)

        p = t∂2d\t∂d
        α = 1
        
        println("step: ", p)
        θ = θ - α*p

        λ = θ[1]
        β = θ[2]

        println("new lambda: ", λ, " new beta: ", β)

        μ_list, gradList, simBorderPts, splinex, spliney, mdl = test(x0, x1, y0, y1, z0, z1, ne, λ, mutst, ndim, FunctionClass, nDof, β, CameraMatrix, 
                                                                endTime, tSteps, Control, mode=mode, SIDES=SIDES)

        d, ∂d, ∂2d, pairs = closest_point(simBorderPts, obsBorderPts, gradList)

        totd = sum(d)
        ratio = totd/totdinit
        totdinit = totd

        iter = iter + 1
        
        push!(lambdapList,λ)
        push!(βpList,β)
        push!(costList,totd)
        push!(iterList,iter)
    end

    sampleNo = 20
    lambdaList = collect(range(lambdatstStart, stop=lambdatst+dev_λ, length=sampleNo))
    βList = collect(range(βStart, stop=β+dev_β, length=sampleNo))
    CostMat = zeros(size(lambdaList,1),size(βList,1))
        
    λ_iter = 1:size(lambdaList,1)
    β_iter = 1:size(βList,1)

    for i in λ_iter
        for j in β_iter
            θ = [lambdaList[i]; βList[j]]
            λ = θ[1]
            β = θ[2]

            μ_list, gradList, simBorderPts, splinex, spliney, mdl = test(x0, x1, y0, y1, z0, z1, ne, λ, mutst, ndim, FunctionClass, nDof, β, CameraMatrix, 
                                                                    endTime, tSteps, Control, mode=mode, SIDES=SIDES)

            # test the closest point function
            d_cp, pairs = closest_point(simBorderPts, obsBorderPts) 

            CostMat[i,j] = sum(d_cp)
        end
    end

    # Plot the cost function with iterations
    Plots.plot(iterList, costList, marker=1, dpi=400)
    Plots.xlabel!("Iterations")
    Plots.ylabel!("Cost")
    Plots.savefig(string(filepathi,"/Results/cost/cost_steps.png"))
    
    # Plot the cost function surface
    Plots.contour(lambdaList, βList, CostMat, color=:turbo, fill=false, levels=50, xlabel="λ", ylabel="β", dpi=400)
    Plots.plot!(lambdapList, βpList, label="Estimations", marker=1)
    Plots.xlabel!("λ")
    Plots.ylabel!("β") 
    Plots.savefig(string(filepathi,"/Results/cost/cost_surface.png"))

end

fit()

