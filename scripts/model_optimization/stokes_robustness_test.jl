
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
    dev = 0.1

    control = "force" # "force" or "velocity"
    viscosity_type = "constant" # "constant" or "bulk_viscosity"

    filepath = string("/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/sim_experiments/cost_function_test/optimization/Stokes/",control,"/test3")

    # simulation parameters for the ground truth
    sim_time = 10.0
    steps = 10.0
    t_steps = sim_time/steps

    F = 3.0

    ηLst = [30 40 50 60 70 80 90]
    βLst = [0.00001 100 1000 10000]

    noiseLevelLst = [0 0.25 0.5 1 1.5 2 3 4 5]
    sideList = [false]

    exp_iter::Int = 1
    for β::Float64 in βLst
        dev_β::Float64 = dev*β
        βStart::Float64 = β - dev_β
        # βList = collect(range(βtst-dev_β, stop=βtst+dev_β, length=sampleNo))
        for η::Float64 in ηLst
            dev_η::Float64 = dev*η
            ηStart::Float64 = η - dev_η
            # ηList = collect(range(ηtst-dev_η, stop=ηtst+dev_η, length=sampleNo))

            println("Ground truth: η :", η, " β :", β)
            model, scene = def_problem(r, h, ne, η, ndim, FunctionClass_u, nDof_u, FunctionClass_p, nDof_p, β, F, control, viscosity_type, sim_time, t_steps)
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

                    if maximum(ηpList) > η+dev_η
                        ηStop = maximum(ηpList)*1.1
                    else
                        ηStop = η+dev_η
                    end

                    if minimum(ηpList) < η-dev_η
                        ηStart = minimum(ηpList)*0.9
                    else
                        ηStart = η-dev_η
                    end

                    if maximum(βpList) > β+dev_β
                        βStop = maximum(βpList)*1.1
                    else
                        βStop = β+dev_β
                    end

                    if minimum(βpList) < β-dev_β
                        βStart = minimum(βpList)*0.9
                    else
                        βStart = β-dev_β
                    end

                    sample_space = 21
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
                            model.η = [θ[1]]
                            scene.β = [θ[2]]
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
                    Plots.contour(ηList, βList, CostMat, color=:turbo, fill=false, levels=50, xlabel="η", ylabel="β", dpi=400)
                    Plots.plot!(ηpList, βpList, label="Estimations", marker=1)
                    Plots.xlabel!("η")
                    Plots.ylabel!("β")
                    Plots.savefig(string(filepathi,"/Results/cost/cost_surface_iter.png"))

                    # Plots.contourf(ηList, βList, CostMat, color=:turbo, fill=false, levels=100, xlabel="η", ylabel="β", dpi=400)
                    # Plots.xlabel!("η")
                    # Plots.ylabel!("β")
                    # Plots.savefig(string(filepathi,"/Results/cost/cost_surface.png"))
                end
            end
        end
    end
end

main()