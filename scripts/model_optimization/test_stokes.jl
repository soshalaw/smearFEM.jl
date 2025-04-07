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

    dateTime = Dates.now()
    sampleNo = 21
    dev = 0.3

    Control = "force" # "force" or "velocity"
    noiseLevel = 0
    SIDES = false
    filepathi = string("/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/sim_experiments/cost_function_test/optimization/Stokes/",Control,"/test1")
    dateTime = Dates.now()

    β = 100
    η = 20

    println("Ground truth: η :", η, " β :", β)
    write_sim_data_stokes(x0, x1, y0, y1, z0, z1, ne, η, ndim, FunctionClass_u, nDof_u, FunctionClass_p, nDof_p, β, CameraMatrix, endTime, tSteps, Control,
                    filepathi)
    
    ObsDataList, splinexObs, splineyObs = read_csv(string(filepathi,"/Results/contour_data"))  
    nScene, nSplinex, nSpliney, pd = add_noise(ObsDataList, nFactor=noiseLevel)
    
    ObsData = [nScene, nSplinex, nSpliney]
    
    obsBorderPts = ObsData[1]

    dev_η = η*dev
    dev_β = β*dev
    
    ηStart = η-dev_η
    βStart = β-dev_β

    # θ = [ηStart, βStart]
    θ = [ηStart; βStart]

    # μ_list, gradList, simBorderPts, splinex, spliney, mdl = test_stokes(x0, x1, y0, y1, z0, z1, ne, θ[1], ndim, FunctionClass_u, nDof_u, FunctionClass_p, 
    #                                                             nDof_p, θ[2], CameraMatrix, endTime, tSteps, Control, SIDES=SIDES)

    μ_list, gradList, simBorderPts, splinex, spliney, mdl = test_stokes(x0, x1, y0, y1, z0, z1, ne, θ[1], ndim, FunctionClass_u, nDof_u, FunctionClass_p, 
                                                                                    nDof_p, θ[2], CameraMatrix, endTime, tSteps, Control, SIDES=SIDES)
    
    # test the closest point function
    d, ∂d, ∂2d, pairs = closest_point(simBorderPts, obsBorderPts, gradList)

    totdinit = sum(d)
    println("initial η: ", θ[1], " initial β: ", θ[2])
    println("Initial cost: ", totdinit)
    ratio = 1

    ηpList = Float64[θ[1]]
    βpList = Float64[θ[1]]
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
        println("grad d: ", t∂d)
        println("grad2 d: ", t∂2d)

        p = t∂2d\t∂d
        α = 1
        
        println("step: ", p)
        θ = θ - α*p

        println("new η: ", θ[1], " new beta: ", θ[2])

        μ_list, gradList, simBorderPts, splinex, spliney, mdl = test_stokes(x0, x1, y0, y1, z0, z1, ne, θ[1], ndim, FunctionClass_u, nDof_u, FunctionClass_p, 
                                                                    nDof_p, θ[2], CameraMatrix, endTime, tSteps, Control, SIDES=SIDES)

        d, ∂d, ∂2d, pairs = closest_point(simBorderPts, obsBorderPts, gradList)

        totd = sum(d)
        ratio = totd/totdinit
        totdinit = totd

        iter = iter + 1
        
        push!(ηpList,θ[1]) 
        push!(βpList,θ[2])
        push!(costList,totd)
        push!(iterList,iter)
    end

    sampleNo = 11
    ηList = collect(range(ηStart, stop=η+dev_η, length=sampleNo))
    βList = collect(range(βStart, stop=β+dev_β, length=sampleNo))
    CostMat = zeros(size(ηList,1),size(βList,1))
    costη = zeros(size(ηList,1))
    costβ = zeros(size(βList,1))

    η_iter = 1:size(ηList,1)
    β_iter = 1:size(βList,1)

    for i in η_iter
        η = ηList[i]
        for j in β_iter
            β = βList[j]

            μ_list, gradList, simBorderPts, splinex, spliney, mdl = test_stokes(x0, x1, y0, y1, z0, z1, ne, η, ndim, FunctionClass_u, nDof_u, FunctionClass_p, 
                                                                        nDof_p, β, CameraMatrix, endTime, tSteps, Control, SIDES=SIDES)

            # test the closest point function
            d_cp, pairs = closest_point(simBorderPts, obsBorderPts) 

            # # if η == 820
                # costβ[j] = sum(d_cp)
            # # elseif β == 130
            #     costη[i] = sum(d_cp)
            # end
            # plot_matches(simBorderPts, splinex, spliney, splinexObs, splineyObs, pairs, string(filepath,"/Results/images/matches","matches_eta_", η, "_beta_", β))
            CostMat[i,j] = sum(d_cp)

            # plot(time, cost, label="Cost") 
            # xlabel!("Time steps")
            # ylabel!("Cost")
            # savefig(string(filepathi,"/Results/cost/cost_cp.png"))
        end
    end

    # Plot the cost function with iterations
    Plots.plot(iterList, costList, label="Cost", marker=1, dpi=400, yscale=:log10)
    # Plots.semilogy(iterList, costList, label="Cost", marker=1, dpi=400)
    # PLots.semilogy
    Plots.xlabel!("Iterations")
    Plots.ylabel!("Error")
    Plots.savefig(string(filepathi,"/Results/cost/cost_steps.png"))
    
    # Plot the cost function surface
    Plots.contour(ηList, βList, CostMat, color=:turbo, fill=false, levels=50, xlabel="η", ylabel="β", dpi=400)
    # Plots.contourf(ηList, βList, CostMat, color=:turbo, fill=false, levels=100, xlabel="η", ylabel="β", dpi=400)
    Plots.plot!(ηpList, βpList, label="Estimations", marker=1)
    Plots.xlabel!("η")
    Plots.ylabel!("β")
    Plots.savefig(string(filepathi,"/Results/cost/cost_surface.png"))

    # Plots.plot(ηList, costη, label="Cost η", marker=1, dpi=400)
    # Plots.plot!(etapList, costList, label="Estimations", marker=1)
    # Plots.xlabel!("η") 
    # Plots.ylabel!("Cost")
    # Plots.savefig(string(filepathi,"/Results/cost/cost_eta.png"))

    # Plots.plot(βList, costList, label="Cost β", marker=1, dpi=400)
    # Plots.plot!(βpList, costList, label="Estimations", marker=1)
    # Plots.xlabel!("β")
    # Plots.ylabel!("Cost")
    # Plots.savefig(string(filepathi,"/Results/cost/cost_beta.png"))

end

fit()

