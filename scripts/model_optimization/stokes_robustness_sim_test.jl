using LinearAlgebra
using ProgressMeter
using SparseArrays

using smearFEM
using StatsPlots
using Distributions
using Measures
using Dates
using DelimitedFiles

function run(noiseLevelLst, file_path)

    filepath = string("/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/sim_experiments/Synthetic_data/")
    dev::Float64 = 0.3

    # simulation parameters for the experiments
    sim_time::Float64 = 20.0# simulation time in seconds
    steps::Float64 = 20.0 # number of time steps
    t_steps::Float64 = sim_time/steps

    sampleNo = 21
    sideList = [false]

    dirs = readdir(filepath)

    for dir in dirs
        filepath_gt = string(filepath, dir)
        params = read_json(string(filepath_gt,"/data/params.json"))

        r = params["r"]
        h = params["h"]
        η_gt = params["η"]
        β_gt = params["β"]
        camera_matrix = params["camera_matrix"]
        camera_pose = params["camera_pose"]
        control = params["control_type"]
        viscosity_type = params["viscosity_type"]
        F = params["cParam"]

        sim_time_gt = params["simulation_time"]
        t_steps_gt = params["time_steps"]

        model_gt, scene_gt = def_problem(r, h, ne_exp, η_gt[1], ndim, FunctionClass_u, nDof_u, FunctionClass_p, nDof_p, β_gt[1], F, control, viscosity_type, sim_time_gt, t_steps_gt)
        printstyled("Ground truth η: $(η_gt), ground truth β: $(β_gt)\n"; color = :green)

        dev_η::Float64 = dev*η_gt
        ηStart::Float64 = η_gt - dev_η

        dev_β::Float64 = dev*β_gt
        βStart::Float64 = β_gt - dev_β

        gt_h = readdlm(string(filepath_gt,"/Results/data/h.csv"), ',', Float64, '\n', header=false)
        # Write the ground truth
        params = Dict("gt_η" => η_gt, "gt_β" => β_gt)
        write_json(string(filepath_gt,"/params"), params)

        ne_exp = 4
        model, scene = def_problem(r, h, ne_exp, η_gt, ndim, FunctionClass_u, nDof_u, FunctionClass_p, nDof_p, β_gt, F, control, viscosity_type, sim_time, t_steps)
        for noiseLevel::Float64 in noiseLevelLst
            ObsDataList, splinexObs, splineyObs = read_csv(string(filepath_gt,"/Results/contour_data"))  
            ObsDataList = ObsDataList[1:(round(Int,sim_time/t_steps)+1)]
            if noiseLevel == 0.0
                # Read the gt data 
                obsBorderPts, nSplinex, nSpliney, pd = add_noise(ObsDataList, nFactor=noiseLevel)
                
                for sides::Bool in sideList
                    filepathi = string(filepath_gt,"/trials/noise_$(noiseLevel)")
                    conditions = Conditions(camera_matrix=camera_matrix, camera_pose=camera_pose, SIDES=sides, filepath=filepathi, ANIMATE=false)

                    θ::Vector{Float64} = [ηStart, βStart]
                    stats = fit_model(model, scene, conditions, obsBorderPts, θ)
                                        
                    iterList = stats["iterList"]
                    costList = stats["costList"]
                    ηpList = stats["ηList"]
                    βpList = stats["βList"]
                    η = stats["η"]
                    β = stats["β"]

                    η_accuracy = abs(1 - abs(η-η_gt)/η_gt)*100
                    β_accuracy = abs(1 - abs(β-β_gt)/β_gt)*100
                    printstyled("η accuracy: $(η_accuracy) %\n"; color = :green)
                    printstyled("β accuracy: $(β_accuracy) %\n"; color = :green)
                    
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

                    sample_space = 41
                    ηList = collect(range(ηStart, stop=ηStop, length=sample_space))
                    βList = collect(range(βStart, stop=βStop, length=sample_space))
                    CostMat = zeros(size(ηList,1),size(βList,1))

                    η_iter = 1:size(ηList,1)
                    β_iter = 1:size(βList,1)

                    for i::Int in η_iter
                        η = ηList[i]
                        for j::Int in β_iter
                            β = βList[j]
                            reset_model!(model)
                            model.η = [η]
                            scene.β = [β]
                            μ_list, gradList, simBorderPts, splinex, spliney, pos2D = simulate(model, scene, conditions)

                            # test the closest point function
                            d_cp, pairs = closest_point(simBorderPts, obsBorderPts) 

                            CostMat[i,j] = sum(d_cp)
                        end
                    end

                    model_tst, scene_tst = def_problem(r, h, ne_exp, η, ndim, FunctionClass_u, nDof_u, FunctionClass_p, nDof_p, β, F, control, viscosity_type, sim_time_gt, t_steps_gt)
                    # simulate the model with the estimated parameters
                    est_μ_list, others = simulate(model_tst, scene_tst, conditions)

                    est_h = get_height(est_μ_list, h)

                    write_csv(string(filepathi,"/Results/data/cost_surface"), CostMat)
                    write_csv(string(filepathi,"/Results/data/cost_surface_η"), ηList)
                    write_csv(string(filepathi,"/Results/data/cost_surface_β"), βList)

                    write_csv(string(filepathi,"/Results/data/cost_surface_ηp"), ηpList)
                    write_csv(string(filepathi,"/Results/data/cost_surface_βp"), βpList)

                    write_csv(string(filepathi,"/Results/data/cost_steps"), costList)
                    write_csv(string(filepathi,"/Results/data/cost_steps_iter"), iterList)

                    write_csv(string(filepathi,"/Results/data/est_h"), est_h)

                    write_csv(string(filepathi,"/Results/data/border"), obsBorderPts[round(Int,sim_time/2)])

                    params = Dict("gt_η" => η_gt,
                    "gt_β" => β_gt,
                    "η" => η,
                    "β" => β,
                    "η_accuracy" => η_accuracy,
                    "β_accuracy" => β_accuracy)

                    write_json(string(filepathi,"/Results/data/stats"), params)

                    set_file(string(filepathi,"/Results/plots"))
                    Plots.plot(est_h, label="Estimated height", dpi=400)
                    Plots.plot!(gt_h, label="Ground truth height", dpi=400)
                    Plots.xlabel!("Time (s)")
                    Plots.ylabel!("Height")
                    Plots.savefig(string(filepathi,"/Results/plots/h_est.pdf"))

                    Plots.plot(abs.(est_h-gt_h), label="Height estimation error", dpi=400)
                    Plots.xlabel!("Time (s)")
                    Plots.ylabel!("Error")
                    Plots.savefig(string(filepathi,"/Results/plots/h_est_error.pdf"))

                    # Plot the cost function with iterations
                    Plots.plot(iterList, costList, label="Cost", marker=1, dpi=400, yscale=:log10)
                    Plots.xlabel!("Iterations")
                    Plots.ylabel!("Error")
                    Plots.savefig(string(filepathi,"/Results/plots/cost_steps.pdf"))
                    
                    # Plot the cost function surface
                    Plots.contour(ηList, βList, CostMat, color=:turbo, fill=false, levels=200, xlabel="η", ylabel="β", dpi=400)
                    Plots.plot!(ηpList, βpList, label="Estimations", marker=1)
                    Plots.plot!([η_gt], [β_gt], label="Ground truth", marker=:star5, markersize=5, markercolor=:green)
                    Plots.xlabel!("η")
                    Plots.ylabel!("β")
                    Plots.savefig(string(filepathi,"/Results/plots/cost_surface_iter.pdf"))
                end
            else
                n_samples = 10
                η_pred = zeros(Float64, n_samples)
                β_pred = zeros(Float64, n_samples)
                costnList = Vector{AbstractVector}(undef, n_samples)
                iternList = Vector{AbstractVector}(undef, n_samples)
                filepathi = string(filepath_gt,"/trials/noise_$(noiseLevel)")
                est_h_list = zeros(Float64, round(Int,(steps_gt+1)), n_samples)

                for n::Int in 1:n_samples
                    obsBorderPts, nSplinex, nSpliney, pd = add_noise(ObsDataList, nFactor=noiseLevel)

                    conditions = Conditions(camera_matrix=camera_matrix,filepath=filepathi)

                    θ = [ηStart, βStart]
                    stats = fit_model(model, scene, conditions, obsBorderPts, θ)
                                        
                    iterList = stats["iterList"]
                    costList = stats["costList"]
                    ηpList = stats["ηList"]
                    βpList = stats["βList"]
                    η = stats["η"]
                    β = stats["β"]

                    model_tst, scene_tst = def_problem(r, h, ne_exp, η, ndim, FunctionClass_u, nDof_u, FunctionClass_p, nDof_p, β, F, control, viscosity_type, sim_time_gt, t_steps_gt)
                    
                    # simulate the model with the estimated parameters
                    est_μ_list, others = simulate(model_tst, scene_tst, conditions)

                    est_h_list[:,n] = get_height(est_μ_list, h)

                    η_accuracy = abs(1 - abs(η-η_gt)/η_gt)*100
                    β_accuracy = abs(1 - abs(β-β_gt)/β_gt)*100
                    printstyled("η accuracy: $(η_accuracy) %\n"; color = :green)
                    printstyled("β accuracy: $(β_accuracy) %\n"; color = :green)

                    η_pred[n] = η
                    β_pred[n] = β
                    push!(costnList, costList)
                    push!(iternList, iterList)

                    params = Dict("gt_η" => η_gt,
                    "gt_β" => β_gt,
                    "η" => η,
                    "β" => β,
                    "η_accuracy" => η_accuracy,
                    "β_accuracy" => β_accuracy)

                    write_json(string(filepathi,"/Results/data/stats"), params)

                    write_csv(string(filepathi,"/Results/data/cost_steps/run_$n"), costList)
                    write_csv(string(filepathi,"/Results/data/cost_steps_iter/run_$n"), iterList)
                    write_csv(string(filepathi,"/Results/data/border"), obsBorderPts[round(Int,sim_time/2)])
                end

                write_csv(string(filepathi,"/Results/data/η_est"), η_pred)
                write_csv(string(filepathi,"/Results/data/β_est"), β_pred)
                write_csv(string(filepathi,"/Results/data/h_est"), est_h_list)

                plot_covariance(η_pred, β_pred, string(filepathi,"/Results/plots/"))
                
                plt_error = set_plot(22, xlabel="Time (s)", ylabel="Height")
                StatsPlots.errorline!(est_h_list, label="Estimated height", dpi=400)
                Plots.plot!(gt_h, label="Ground truth height", dpi=400)
                Plots.xlabel!("Time (s)")
                Plots.ylabel!("Height")
                Plots.savefig(string(filepathi,"/Results/plots/h_est.pdf"))

                set_file(string(filepathi,"/Results/plots"))
                set_plot(22)
                plot!(x->pdf(pd, x),  size=(200,150), label="")
                savefig(string(filepathi,"/Results/plots/obs_pdf_$(noiseLevel).pdf"))

                set_plot(22)
                Plots.plot!(nSplinex[10], nSpliney[10], label="", size=(250,300), xaxis=false, yaxis=false) 
                Plots.xlims!(1100,1350)
                Plots.ylims!(420,1100)
                Plots.savefig(string(filepathi,"/Results/plots/noise_contour_$(noiseLevel).pdf"))
            end
        end
    end
end

# save the data to a file and post process

noiseLevelLst = [0.01 0.05 0.1 0.5 1.0]

run(noiseLevelLst, file_path)
plot_height_vs_slip(ηLst, βLst, file_path)
plot_field_at_height(ηLst, βLst, file_path)
plot_noise_covariance(ηLst, βLst, noiseLevelLst, file_path)
plot_data(ηLst, βLst, noiseLevelLst, file_path, n=10)

