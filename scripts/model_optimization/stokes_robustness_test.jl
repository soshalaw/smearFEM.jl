using LinearAlgebra
using ProgressMeter
using SparseArrays

using smearFEM
using StatsPlots
using Distributions

using Dates
function main()

    # test case 
    r = 0.5
    h = 1.0
    ne = 4
    ndim = 3
    FunctionClass_u = "Q2"
    nDof_u = ndim  # number of degree of freedom per node
    FunctionClass_p = "Q1"
    nDof_p = 1  # number of degree of freedom per node

    CameraMatrix = [[8*2048/7.07, 0.0, 2048/2] [0.0, 8*1536/5.3, 1536/2] [0.0, 0.0, 1.0]]'

    dateTime = Dates.now()
    sampleNo = 21
    dev = 0.3

    control = "force" # "force" or "velocity"
    viscosity_type = "constant" # "constant" or "bulk_viscosity"

    filepath = string("/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/sim_experiments/cost_function_test/optimization/Stokes/",control,"/test3")

    # simulation parameters for the ground truth
    sim_time = 10.0
    steps = 10.0
    t_steps = sim_time/steps

    ηLst = [40 100 500 1000 2000 5000 10000]
    βLst = [100 1000 10000]

    noiseLevelLst = [0.0]
    sideList = [false]

    # βLst = [100.0]
    # ηLst = [40.0]
    # noiseLevelLst = [0.0]
    # sideList = [false]
    F = 3.0

    exp_iter::Int = 1
    for β_gt::Float64 in βLst
        dev_β::Float64 = dev*β_gt
        βStart::Float64 = β_gt - dev_β
        for η_gt::Float64 in ηLst
            dev_η::Float64 = dev*η_gt
            ηStart::Float64 = η_gt - dev_η

            println("Ground truth: η :", η_gt, " β :", β_gt)
            model, scene = def_problem(r, h, ne, η_gt, ndim, FunctionClass_u, nDof_u, FunctionClass_p, nDof_p, β_gt, F, control, viscosity_type, sim_time, t_steps)
            filepath_gt = string(filepath,"/experiment_$(exp_iter)/ground_truth")
            write_sim_data(model, scene, CameraMatrix, filepath_gt)
            
            for noiseLevel::Float64 in noiseLevelLst
                ObsDataList, splinexObs, splineyObs = read_csv(string(filepath_gt,"/Results/contour_data"))  
                nScene, nSplinex, nSpliney, pd = add_noise(ObsDataList, nFactor=noiseLevel)
                ObsData = [nScene, nSplinex, nSpliney]
                obsBorderPts = ObsData[1]
                
                for sides::Bool in sideList
                    filepathi = string(filepath,"/experiment_$(exp_iter)/trials/noise_$(noiseLevel)/sides_$(sides)")
                    conditions = Conditions(CameraMatrix=CameraMatrix,SIDES=sides,filepath=filepathi)

                    θ = [ηStart, βStart]
                    stats = fit_model(model, scene, conditions, obsBorderPts, θ)
                                        
                    iterList = stats["iterList"]
                    costList = stats["costList"]
                    ηpList = stats["ηList"]
                    βpList = stats["βList"]
                    η = stats["η"]
                    β = stats["β"]

                    η_accuracy = abs(η - η_gt)/η_gt*100
                    β_accuracy = abs(β - β_gt)/β_gt*100
                    printstyled("η accuracy: $(η_accuracy) %\n"; color = :green)
                    printstyled("β accuracy: $(β_accuracy) %\n"; color = :green)

                    params = Dict("η" => η,
                              "β" => β,
                              "η_accuracy" => η_accuracy,
                              "β_accuracy" => β_accuracy)

                    write_json(string(filepathi,"/Results/params"), params)

                    if maximum(ηpList) > η_gt+dev_η
                        ηStop = maximum(ηpList)*1.1
                    else
                        ηStop = η_gt+dev_η
                    end

                    if minimum(ηpList) < η_gt-dev_η
                        ηStart = minimum(ηpList)*0.9
                    else
                        ηStart = η_gt-dev_η
                    end

                    if maximum(βpList) > β_gt+dev_β
                        βStop = maximum(βpList)*1.1
                    else
                        βStop = β_gt+dev_β
                    end

                    if minimum(βpList) < β_gt-dev_β
                        βStart = minimum(βpList)*0.9
                    else
                        βStart = β_gt-dev_β
                    end

                    sample_space = 3
                    ηList = collect(range(ηStart, stop=ηStop, length=sample_space))
                    βList = collect(range(βStart, stop=βStop, length=sample_space))
                    CostMat = zeros(size(ηList,1),size(βList,1))

                    η_iter = 1:size(ηList,1)
                    β_iter = 1:size(βList,1)

                    for i::Int in η_iter
                        η = ηList[i]
                        for j::Int in β_iter
                            β = βList[j]
                            reset_model(model)
                            model.η = [η]
                            scene.β = [β]
                            μ_list, gradList, simBorderPts, splinex, spliney, pos2D = simulate(model, scene, conditions)

                            # test the closest point function
                            d_cp, pairs = closest_point(simBorderPts, obsBorderPts) 

                            CostMat[i,j] = sum(d_cp)
                        end
                    end

                    # Plot the cost function with iterations
                    Plots.plot(iterList, costList, label="Cost", marker=1, dpi=400, yscale=:log10)
                    Plots.xlabel!("Iterations")
                    Plots.ylabel!("Error")
                    Plots.savefig(string(filepathi,"/Results/cost/cost_steps.png"))
                    
                    # Plot the cost function surface
                    Plots.contour(ηList, βList, CostMat, color=:turbo, fill=false, levels=100, xlabel="η", ylabel="β", dpi=400)
                    Plots.plot!(ηpList, βpList, label="Estimations", marker=1)
                    Plots.plot!([η_gt], [β_gt], label="Ground truth", marker=:star5, markersize=8, markercolor=:green)
                    Plots.xlabel!("η")
                    Plots.ylabel!("β")
                    Plots.savefig(string(filepathi,"/Results/cost/cost_surface_iter.png"))

                end
            end
        end
    end
end

main()