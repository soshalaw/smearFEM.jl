using LinearAlgebra
using ProgressMeter
using SparseArrays

using smearFEM
using StatsPlots
using Distributions
using Dates
using Plots

function test_opt()
    # test case 
    r::Float64 = 1.0
    h::Float64 = 1.0
    ne::Int = 4
    ndim::Int = 3
    FunctionClass_u::String = "Q2"
    nDof_u::Int = ndim  # number of degree of freedom per node
    FunctionClass_p::String = "Q1"
    nDof_p::Int = 1  # number of degree of freedom per node

    CameraMatrix = [[8*2048/7.07, 0.0, 2048/2] [0.0, 8*1536/5.3, 1536/2] [0.0, 0.0, 1.0]]'

    dateTime = Dates.now()
    dev::Float64 = 0.1

    control::String = "force" # "force" or "velocity"
    viscosity_type::String = "bulk_viscosity" # "constant" or "bulk_viscosity"
    noiseLevel::Float64 = 0
    SIDES::Bool = false
    filepathi::String = string("/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/sim_experiments/cost_function_test/optimization/Stokes/",control,"/test1")

    # simulation parameters for the ground truth
    gt_sim_time::Float64 = 90.0
    gt_steps::Float64 = 90.0
    gt_t_steps::Float64 = gt_sim_time/gt_steps
  
    β::Float64 = 100.0
    η::Float64 = 40.0
    F::Float64 = 1.0

    # Write the ground truth
    gt_model, gt_scene = def_problem(r, h, ne, η, ndim, FunctionClass_u, nDof_u, FunctionClass_p, nDof_p, β, F, control, viscosity_type, gt_sim_time, gt_t_steps)
    # write_sim_data(gt_model, gt_scene, CameraMatrix, filepathi)
    
    # Read the gt data
    ObsDataList, splinexObs, splineyObs = read_csv(string(filepathi,"/Results/contour_data"))  
    nScene, nSplinex, nSpliney, pd = add_noise(ObsDataList, nFactor=noiseLevel)
    ObsData = [nScene, nSplinex, nSpliney]
    obsBorderPts = ObsData[1]

    sim_time::Float64 = 15.0
    steps::Float64 = 15.0
    t_steps::Float64 = sim_time/steps

    viscosity_type = "constant"
    conditions = Conditions(CameraMatrix=CameraMatrix)
    model, scene = def_problem(r, h, ne, η, ndim, FunctionClass_u, nDof_u, FunctionClass_p, nDof_p, β, F, control, viscosity_type, sim_time, t_steps)

    dev_η::Float64 = η*dev
    dev_β::Float64 = β*dev

    ηStart::Float64 = η - dev_η
    βStart::Float64 = β - dev_β

    θ::Vector{Float64} = [ηStart, βStart]
    
    gt_time_frame::Int = round(Int,gt_sim_time)
    sim_time_frame::Int = round(Int,sim_time)
    windows::Int = gt_time_frame/sim_time_frame

    ηpList = Vector{Float64}(undef,gt_time_frame)
    βpList = Vector{Float64}(undef,gt_time_frame)


    println("Number of time windows: ", windows)
    titer::Int = 1
    for ti::Int in 1:windows

        printstyled("Time window: $(ti)\n"; color = :green)
        gt_η = gt_model.η[round(Int,(sim_time*(titer-1))+1):round(Int,sim_time*(titer))]
        av_η = mean(gt_η)
        printstyled("Average ground truth η in the window: $(av_η), ground truth β: $(β)\n"; color = :green)

        stats = fit(model, scene, conditions, obsBorderPts[(sim_time_frame*(titer-1)+1):(sim_time_frame*titer)+1], θ)

        ηpList[(sim_time_frame*(titer-1)+1):(sim_time_frame*titer)] .= stats["η"]
        βpList[(sim_time_frame*(titer-1)+1):(sim_time_frame*titer)] .= stats["β"]
        titer = titer + 1
        update_model(model)
    end

    t_windows = collect(range(start=gt_t_steps, stop=gt_sim_time, step=gt_t_steps))
    Plots.plot(gt_model.η, label="Ground η", marker=1, dpi=400)
    Plots.plot!(t_windows, ηpList, label="Estimated η", marker=1)
    Plots.xlabel!("time (s)")
    Plots.ylabel!("η")
    Plots.savefig(string(filepathi,"/Results/cost/η.png"))

    # if maximum(ηpList) > η+dev_η
    #     ηStop = maximum(ηpList)*1.1
    # else
    #     ηStop = η+dev_η
    # end

    # if minimum(ηpList) < η-dev_η
    #     ηStart = minimum(ηpList)*0.9
    # else
    #     ηStart = η-dev_η
    # end

    # if maximum(βpList) > β+dev_β
    #     βStop = maximum(βpList)*1.1
    # else
    #     βStop = β+dev_β
    # end

    # if minimum(βpList) < β-dev_β
    #     βStart = minimum(βpList)*0.9
    # else
    #     βStart = β-dev_β
    # end

    # sampleNo = 11
    # ηList = collect(range(ηStart, stop=ηStop, length=sampleNo))
    # βList = collect(range(βStart, stop=βStop, length=sampleNo))
    # CostMat = zeros(size(ηList,1),size(βList,1))
    # costη = zeros(size(ηList,1))
    # costβ = zeros(size(βList,1))

    # η_iter = 1:size(ηList,1)
    # β_iter = 1:size(βList,1)

    # for i in η_iter
    #     η = ηList[i]
    #     for j in β_iter
    #         β = βList[j]

    #         μ_list, gradList, simBorderPts, splinex, spliney, mdl = test_stokes(x0, x1, y0, y1, z0, z1, ne, η, ndim, FunctionClass_u, nDof_u, FunctionClass_p, 
    #                                                                     nDof_p, β, CameraMatrix, sim_time, tSteps, control, SIDES=SIDES)

    #         # test the closest point function
    #         d_cp, pairs = closest_point(simBorderPts, obsBorderPts) 
    #         CostMat[i,j] = sum(d_cp)
    #     end
    # end

    # # Plot the cost function with iterations
    # Plots.plot(iterList, costList, label="Cost", marker=1, dpi=400, yscale=:log10)
    # Plots.xlabel!("Iterations")
    # Plots.ylabel!("Error")
    # Plots.savefig(string(filepathi,"/Results/cost/cost_steps.png"))
    
    # # Plot the cost function surface
    # Plots.contour(ηList, βList, CostMat, color=:turbo, fill=false, levels=50, xlabel="η", ylabel="β", dpi=400)
    # Plots.plot!(ηpList, βpList, label="Estimations", marker=1)
    # Plots.xlabel!("η")
    # Plots.ylabel!("β")
    # Plots.savefig(string(filepathi,"/Results/cost/cost_surface_iter.png"))

    # Plots.contourf(ηList, βList, CostMat, color=:turbo, fill=false, levels=100, xlabel="η", ylabel="β", dpi=400)
    # Plots.xlabel!("η")
    # Plots.ylabel!("β")
    # Plots.savefig(string(filepathi,"/Results/cost/cost_surface.png"))

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

function fit(model::Stokes, scene::SqueezeFlow, conditions::Conditions, obsBorderPts::Vector{AbstractArray}, θ::Vector{Float64})

    model.η = [θ[1]]
    scene.β = [θ[2]]

    μ_list, gradList, simBorderPts, splinex, spliney, pos2D = simulate(model, scene, conditions)
    d, ∂d, ∂2d, pairs = closest_point(simBorderPts, obsBorderPts, gradList)
    totdinit::Float64 = sum(d)
    
    ηpList = Vector{Float64}(undef,0)
    βpList = Vector{Float64}(undef,0)
    costList = Vector{Float64}(undef,0)
    iterList = Vector{Float64}(undef,0)

    iter::Int = 0
    c_grad::Float64 = 1.0
    while c_grad > 0.005
        reset_model(model)

        t∂2d = zeros(size(∂2d[1]))
        t∂d = zeros(size(∂d[1]))
        
        szd = 1:length(d)

        for i in szd
            t∂2d = t∂2d + ∂2d[i]
            t∂d = t∂d + ∂d[i]
        end

        p = t∂2d\t∂d
        α = 1
        
        println("step: ", p)
        θ = θ - α*p

        println("new η: ", θ[1], " new beta: ", θ[2])

        model.η = [θ[1]]
        scene.β = [θ[2]]
        μ_list, gradList, simBorderPts, splinex, spliney, pos2D = simulate(model, scene, conditions)

        # test the closest point function
        d, ∂d, ∂2d, pairs = closest_point(simBorderPts, obsBorderPts, gradList)

        totd = sum(d)
        c_grad = abs(totdinit - totd)
        totdinit = totd

        iter = iter + 1
        
        push!(ηpList,θ[1]) 
        push!(βpList,θ[2])
        push!(costList,totd)
        push!(iterList,iter)
        println("iteration: ", iter, " ratio: ", c_grad)
    end

    stats = Dict("η" => θ[1],
                 "β" => θ[2],
                 "ηList" => ηpList, 
                 "βList" => βpList,
                 "costList" => costList,
                 "iterList" => iterList,
                 "η" => θ[1],
                 "β" => θ[2])
    return stats
end
test_opt()

