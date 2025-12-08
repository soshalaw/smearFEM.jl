using LinearAlgebra
using ProgressMeter
using SparseArrays

using smearFEM
using StatsPlots
using Distributions
using Dates
using Plots
using LaTeXStrings
using DelimitedFiles

using ArgCheck
using Base.Threads
using Random

using Dates
# using JLD2
using Plots.PlotMeasures

global fs = 32  # font size for plots

function optimize(exp_params::Dict)
    η_gt::Float64 = 0.0
    β_gt::Float64 = 0.0
    F_ext::Float64 = 0.0
    sim_time_gt::Float64 = 0.0
    t_steps_gt::Float64 = 0.0
    steps_exp::Float64 = 0.0
    outlier_frames::Vector{Int} = Int[]

    ndim::Int = 3
    nDof_p::Int = 1  # number of degree of freedom per node
    nDof_u::Int = ndim  # number of degree of freedom per node

    FunctionClass_u::String = exp_params["FunctionClass_u"]
    FunctionClass_p::String = exp_params["FunctionClass_p"]    
    FunctionClass_x::String = exp_params["FunctionClass_x"]

    iterList::Vector{Float64} = Vector{Float64}()
    costList::Vector{Float64}  = Vector{Float64}()
    ηpList::Vector{Float64} = Vector{Float64}()
    βpList::Vector{Float64} = Vector{Float64}()
    
    # simulation parameters for the ground truth
    
    filepath_gt::String = exp_params["filepath_gt"]
    filepath_res::String = exp_params["filepath_res"]
    
    obj_pose::AbstractArray = zeros(Float64, 4,4) # initial pose of the object
    camera_matrix::AbstractArray = exp_params["camera_matrix"]
    
    sim_time_exp::Float64 = exp_params["sim_time_exp"]
    if haskey(exp_params, "steps_exp")
        steps_exp = exp_params["steps_exp"]
    else
        steps_exp = sim_time_exp*10.0 # assuming 30 fps for the experiments
    end
    t_steps_exp::Float64 = sim_time_exp/steps_exp
    
    ne_exp::Int = exp_params["ne_exp"] # number of elements in the mesh for the experiment
    
    data_type ::String = exp_params["data_type"] # "simulated" or "physical" or "real"
    noiseLevel::Float64 = exp_params["noise_level"]
    SIDES::Bool = false
    
    if data_type == "simulated" || data_type == "synthetic"
        
        WRITE_GT = exp_params["WRITE_GT"] 
        outlier_frames = Int[]
        
        if WRITE_GT == true # write the ground truth data
            @info "Writing ground truth gt data to with $ne_gt elements to $filepath_res"
            write_gt_data(exp_params)
        end
        
        params = read_json(string(filepath_gt,"/data/sim_params.json"))
        
        r = params["r"]
        h = params["h"]
        
        η_gt = params["η"][1]
        β_gt = params["β"][1]
        gt_viscosity_type = params["viscosity_type"]
        F = Array(float.(params["cParam"]))
        println("Ground truth force: $F")
        
        sim_time_gt = params["simulation_time"]
        t_steps_gt = params["time_steps"]
        
        camera_matrix = reshape(Array(params["camera_matrix"]), 3, 3)
        obj_pose = reshape(Array(params["obj_pose"])*1.0,4,4)
        
        control = params["control_type"]
        
        model_gt, scene_gt = def_problem(r, h, ne_exp, η_gt, ndim, FunctionClass_u, nDof_u, FunctionClass_p, nDof_p, FunctionClass_x, β_gt, F, control, gt_viscosity_type, 
        sim_time_gt, t_steps_gt)
        
        if data_type == "synthetic"
            ObsDataList, splinexObs, splineyObs = read_csv(string(filepath_gt,"/data/img_data/contour_data"))  
        elseif data_type == "simulated"
            ObsDataList, splinexObs, splineyObs = read_csv(string(filepath_gt,"/data/sim_data/contour_data"))  
        end
        
        printstyled("Ground truth η: $(η_gt), ground truth β: $(β_gt)\n"; color = :green)
        
    elseif data_type == "physical"
        r::Float64 = exp_params["r"]  # radius of the cylinder in mm
        h::Float64 = exp_params["h"]  # height of the cylinder in mm
        
        control::String = exp_params["control"]
        gt_viscosity_type::String = exp_params["viscosity_type"] # "constant" or "bulk_viscosity"

        t_obs = read_perception_data(string(filepath_gt, "/Results/sequence.hdf5"))

        _ObsDataList, _splinexObs, _splineyObs = read_csv(string(filepath_gt, "/data/sim_data/contour_data"))
        meta_data = read_json(string(filepath_gt, "/data/video_metadata.json"))

        frame_rate = meta_data["frame_rate"]
        frame_width = meta_data["frame_width"]
        frame_height = meta_data["frame_height"]
        compression_frames = meta_data["compressed_frames"]

        sim_time_gt = length(compression_frames)/frame_rate  # seconds
        steps_exp = sim_time_exp*frame_rate
        t_steps_exp = 1/frame_rate
        t_steps_gt = t_steps_exp

        F_ext = exp_params["F_ext"]
        F = -F_ext*ones(Float64, round(Int, steps_exp)) # force applied to the cylinder in N
        println("Applied force: $F_ext N")

        # t_obs = Float64.(t_obs)

        obj_pose = get_pose(t_obs)  # Example usage of get_pose function
        obj_pose_ = zeros(Float64, 4,4)
        obj_pose_[1,1] = -1.0
        obj_pose_[2,3] = -1.0
        obj_pose_[3,2] = -1.0
        obj_pose_[1:3,4] = obj_pose[1:3,4]

        ObsDataList = _ObsDataList[compression_frames]
        splinexObs = _splinexObs[compression_frames]
        splineyObs = _splineyObs[compression_frames]

        valid_frames, outlier_frames = detect_outlier_observations(ObsDataList)
    else
        error("data_type should be either simulated or physical")
    end

    if sim_time_gt < sim_time_exp
        @warn "Ground truth simulation time $sim_time_gt is less than experimental simulation time $sim_time_exp , switching to ground truth simulation time"
        sim_time_exp = sim_time_gt
    end
    exp_path = string(filepath_res)
    set_file(exp_path)
    write_json(string(exp_path,"/Results/data/experiment_parameters"), exp_params)

    if t_steps_exp > t_steps_gt
        @warn "time resolution of the ground truth $t_steps_gt is larger than the experimental $t_steps_exp, switching to ground truth resolution"
        t_steps_exp = t_steps_gt
    end
    
    if data_type == "physical"  
        η_start = exp_params["η_start"]
        β_start = exp_params["β_start"]
    else
        dev::Float64 = 0.3

        dev_η::Float64 = dev*η_gt
        η_start::Float64 = abs(η_gt - dev_η)

        dev_β::Float64 = dev*β_gt
        β_start::Float64 = abs(β_gt - dev_β)
    end

    θ::Vector{Float64} = [η_start, β_start]

    if gt_viscosity_type == "constant"
        
        _range = 1:(round(Int,sim_time_exp/t_steps_exp)+1)
        @info "Considering from frame $(first(_range)) to frame $(last(_range)) in the observations"
        ObsDataList = ObsDataList[_range] # align the observation points with the simulation time

        # Read the gt data
        
        model, scene = def_problem(r, h, ne_exp, η_gt, ndim, FunctionClass_u, nDof_u, FunctionClass_p, nDof_p, FunctionClass_x, β_gt, F, control, gt_viscosity_type, 
        sim_time_exp, t_steps_exp)
        conditions = Conditions(camera_matrix=camera_matrix, obj_pose=obj_pose, SIDES=SIDES, filepath=exp_path, ANIMATE=false)
        
        gt_h_ = readdlm(string(filepath_gt,"/data/h.csv"), ',', Float64)
        gt_h = gt_h_[1:(round(Int,sim_time_exp/t_steps_exp)+1)]

        time = collect(Float64, range(start=0, stop=sim_time_exp, step=t_steps_exp))
        if noiseLevel == 0.0
                
            obsBorderPts, nSplinex, nSpliney, pd = add_noise(ObsDataList, nFactor=0.0)
            stats = fit_model(model, scene, conditions, obsBorderPts, θ, outliers=outlier_frames)

            iterList = stats["iterList"]
            costList = stats["cost_list"]
            ηpList = stats["ηList"]
            βpList = stats["βList"]

            η = stats["η"]
            β = stats["β"]

            printstyled("Estimated η : $(η), estimated β: $(β)\n"; color = :green)

            η_accuracy = (1-abs((η_gt-η)/η_gt))*100
            β_accuracy = (1-abs((β_gt-β)/β_gt))*100
            printstyled("η accuracy: $(η_accuracy) %\n"; color = :green)
            printstyled("β accuracy: $(β_accuracy) %\n"; color = :green)

            reset_model!(model)
            model.η = [η]
            scene.β = [β]

            # simulate the gt model with the estimated parameters
            est_μ_list, gradList, borderPts2DList, fields, pos3D, pos2D, splinep, splineq = simulate(model, scene, conditions)

            animate_fields(filepath=string(exp_path,"/Results/plots"), p=splinep, q=splineq, pObs=nSplinex, qObs=nSpliney)

            est_h = get_height(est_μ_list, h)

            set_file(string(exp_path,"/Results/plots"))
            set_plot(fs, sz=(1650, 1250))
            Plots.plot!(time, est_h, label="Estimated height", dpi=400, lw=3)
            Plots.plot!(time, gt_h, label="Ground truth height", dpi=400, lw=3)
            # Plots.xticks!(0:1:round(Int,sim_time_exp))
            Plots.xlabel!(L"\mathrm{Time\;(s)}")
            Plots.ylabel!(L"\mathrm{Height\;(mm)}")
            Plots.savefig(string(exp_path,"/Results/plots/h_est.pdf"))

            set_plot(fs, sz=(1650, 1250))
            Plots.plot!(time, abs.(est_h-gt_h), label="Height estimation error", dpi=400, lw=3)
            Plots.xlabel!(L"\mathrm{Time\;(s)}")
            Plots.ylabel!("Height Error (mm)")
            # Plots.xticks!(0:1:round(Int,sim_time_exp))
            Plots.savefig(string(exp_path,"/Results/plots/h_est_error.pdf"))

            # Plot the cost function with iterations
            set_plot(fs, sz=(1650, 1250))
            Plots.plot!(iterList, costList, label="Cost", marker=2, dpi=400, yscale=:log10, xminorgrid = :false, lw=3)
            Plots.xlabel!(L"\mathrm{Iterations}")
            Plots.ylabel!(L"\mathrm{Cost\;(px)}")
            # Plots.xticks!(minimum(iterList):2:maximum(iterList))
            Plots.savefig(string(exp_path,"/Results/plots/cost_steps.pdf"))

            # Plot the cost function with iterations
            set_plot(fs, sz=(1650, 1250))
            Plots.plot!(iterList, costList, label="Cost", marker=2, dpi=400, yscale=:log10, xscale=:log10, lw=3)
            Plots.xlabel!(L"\mathrm{Iterations}")
            Plots.ylabel!(L"\mathrm{Cost\;(px)}")
            Plots.savefig(string(exp_path,"/Results/plots/cost_steps_log.pdf"))

            write_json(string(exp_path,"/Results/data/stats"), stats)
            write_csv(string(exp_path,"/Results/data/η"), ηpList)
            write_csv(string(exp_path,"/Results/data/β"), βpList)
            write_csv(string(exp_path,"/Results/data/est_h"), est_h)
            write_csv(string(exp_path,"/Results/data/gt_h"), gt_h)
            write_csv(string(exp_path,"/Results/data/cost_iter"), costList)

            write_data(string(exp_path,"/Results/data/sim_data/2D_surface_points"), pos2D)
            write_data(string(exp_path,"/Results/data/sim_data/3D_points"), pos3D)
            write_data(string(exp_path,"/Results/data/sim_data/motion_fields "), fields)
            write_data(string(exp_path,"/Results/data/sim_data/2D_border_points"), borderPts2DList)

            if maximum(ηpList) > η+dev_η
                ηStop = maximum(ηpList)*1.1
            else
                ηStop = η+dev_η
            end

            if minimum(ηpList) < η-dev_η
                η_start = minimum(ηpList)*0.9
            else
                η_start = η-dev_η
            end

            if maximum(βpList) > β+dev_β
                βStop = maximum(βpList)*1.1
            else
                βStop = β+dev_β
            end

            if minimum(βpList) < β-dev_β
                β_start = minimum(βpList)*0.9
            else
                β_start = β-dev_β
            end

            # sampleNo = 10
            # ηList = collect(range(η_start, stop=ηStop, length=sampleNo))
            # βList = collect(range(β_start, stop=βStop, length=sampleNo))
            # CostMat = zeros(size(ηList,1),size(βList,1))
            # ∂CostMat = zeros(2,size(ηList,1),size(βList,1))
            # ∂2CostMat = zeros(2,2,size(ηList,1),size(βList,1))

            # costη = zeros(size(ηList,1))
            # costβ = zeros(size(βList,1))

            # η_iter = 1:size(ηList,1)
            # β_iter = 1:size(βList,1)

            # for i::Int in η_iter
            #     η = ηList[i]
            #     for j::Int in β_iter
            #         β = βList[j]
            #         reset_model!(model)
            #         model.η = [η]
            #         scene.β = [β]
            #         μ_list, gradList, simBorderPts, fields_, pos3D_, pos2D_, splinex_, spliney_ = simulate(model, scene, conditions)
                    
            #         # test the closest point function
            #         d, pairs = closest_point(simBorderPts, obsBorderPts)
                    
            #         CostMat[i,j] = sum(d)/length(d)
            #     end
            # end
            
            # # Plot the cost function surface (interactive GLMakie)
            # set_plot(fs, sz=(1650, 1250))
            # # CostMat is constructed as CostMat[i_eta, j_beta]; Makie expects Z with
            # # size (length(beta), length(eta)) => transpose CostMat for plotting.
            # # First attempt: save an interactive PlotlyJS HTML (avoids creating Makie figures)
            # try
            #     htmlpath = string(exp_path, "/Results/plots/cost_surface_iter_interactive.html")
            #     # compute overlay z-values for estimator path
            #     Z_for_interp = CostMat'
            #     zs_est = [interp_z_at(a, b, ηList, βList, Z_for_interp) for (a,b) in zip(ηpList, βpList)]
            #     # compute ground-truth z for legend/marker overlay
            #     z_gt = interp_z_at(η_gt, β_gt, ηList, βList, Z_for_interp)
            #     saved = save_plotly_surface_html(htmlpath, ηList, βList, CostMat'; xs=ηpList, ys=βpList, zs=zs_est, title="Cost surface (iter)", gt_x=η_gt, gt_y=β_gt, gt_z=z_gt, x_label = "\$\\eta\\;\\mathrm{(KPa\\cdot s)}\$", y_label = "\$\\beta\\;\\mathrm{(L/mm^{-1})}\$", font_size = fs, latex_labels=true)
            #     if saved
            #         @info "Saved interactive PlotlyJS HTML: $htmlpath"
            #     else
            #         @warn "PlotlyJS HTML not created; skipping interactive output. To enable interactive PNGs with GLMakie re-enable GLMakie manually."
            #     end
            # catch err
            #     @warn "PlotlyJS path failed; skipping interactive output. Error: $err"
            # end

            # # Also produce static PDF outputs with Plots (preserve previous behavior)
            # try
            #     plt = set_plot(fs, sz=(1650, 1250))
            #     Plots.contour!(plt, ηList, βList, CostMat, color=:turbo, fill=false, levels=100, dpi=400)
            #     Plots.plot!(plt, ηpList, βpList, label="Estimations", ms=:4, m=:x, color=:red, lw=3)
            #     Plots.plot!(plt, [η_gt], [β_gt], label="Ground truth", ms=:8, m=:star5, color=:indianred2, lw=3)
            #     Plots.xlabel!(plt, L"\eta\;\mathrm{(KPa\cdot s)}")
            #     Plots.ylabel!(plt, L"\beta\;\mathrm{(L/mm^{-1})}")
            #     Plots.savefig(plt, string(exp_path,"/Results/plots/cost_surface_iter.pdf"))

            #     plt2 = set_plot(fs, sz=(1650, 1250))
            #     Plots.contourf!(plt2, ηList, βList, CostMat, color=:turbo, fill=false, levels=100, dpi=400)
            #     Plots.plot!(plt2, [η_gt], [β_gt], label="Ground truth", ms=:8, m=:star5, color=:indianred2, lw=3)
            #     Plots.xlabel!(plt2, L"\eta\;\mathrm{(KPa\cdot s)}")
            #     Plots.ylabel!(plt2, L"\beta\;\mathrm{(L/mm^{-1})}")
            #     Plots.savefig(plt2, string(exp_path,"/Results/plots/cost_surface.pdf"))
            # catch err
            #     @warn "Failed to produce static PDF contour outputs: $err"
            # end

            # # Write the results to files
            # contour_plot_params = Dict("η_list" => ηList, "β_list" => βList, "cost_mat" => CostMat)

            # write_json(string(exp_path,"/Results/data/contour_plot_params"), contour_plot_params)

        # @save string(exp_path,"/Results/data/sim_data/Cost_Matrices.jld2") ηList, βList, CostMat, ∂CostMat, ∂2CostMat
        else
            n_samples = 10
            η_pred = zeros(Float64, n_samples)
            β_pred = zeros(Float64, n_samples)
            costnList = Vector{Vector{Float64}}(undef, n_samples)
            iternList = Vector{AbstractVector}(undef, n_samples)
            est_h_list = Matrix{Float64}(undef, n_samples, round(Int,sim_time_exp/t_steps_exp)+1)
            ANIMATED = false

            set_file(string(exp_path,"/Results/plots"))
            for n::Int in 1:n_samples
                obsBorderPts, nSplinex, nSpliney, pd = add_noise(ObsDataList, nFactor=noiseLevel)

                stats = fit_model(model, scene, conditions, obsBorderPts, θ, outliers=outlier_frames)
                                    
                iterList = stats["iterList"]
                costList = stats["cost_list"]
                ηpList = stats["ηList"]
                βpList = stats["βList"]

                η = stats["η"]
                β = stats["β"]
                
                reset_model!(model)
                model.η = [η]
                scene.β = [β]
                # simulate the model with the estimated parameters
                est_μ_list, gradList, borderPts2DList, fields, pos3D, pos2D, splinep, splineq = simulate(model, scene, conditions)

                if ANIMATED == false
                    ANIMATED = true
                    animate_fields(filepath=string(exp_path,"/Results/plots"), p=splinep, q=splineq, pObs=nSplinex, qObs=nSpliney)
                    set_plot(fs, sz=(1650, 1250))
                    Plots.xlims!(-5,5)
                    plot!(x->pdf(pd, x), label="", dpi=400, lw=3)
                    savefig(string(exp_path,"/Results/plots/obs_pdf.pdf"))
                end
                
                est_h = get_height(est_μ_list, h)
                
                η_accuracy = (1-abs((η_gt-η)/η_gt))*100
                β_accuracy = (1-abs((β_gt-β)/β_gt))*100
                printstyled("η accuracy: $(η_accuracy) %\n"; color = :green)
                printstyled("β accuracy: $(β_accuracy) %\n"; color = :green)
                
                η_pred[n] = η
                β_pred[n] = β
                est_h_list[n,:] = est_h
                costnList[n] = costList
                iternList[n] = iterList

                params = Dict("gt_η" => η_gt,
                "gt_β" => β_gt,
                "η" => η,
                "β" => β,
                "η_accuracy" => η_accuracy,
                "β_accuracy" => β_accuracy)

                write_json(string(exp_path,"/Results/data/stats/run_$n"), params)

                write_csv(string(exp_path,"/Results/data/cost_steps/run_$n"), costList)
                write_csv(string(exp_path,"/Results/data/cost_steps_iter/run_$n"), iterList)
                write_csv(string(exp_path,"/Results/data/border/run_$n"), obsBorderPts)
            end

            write_csv(string(exp_path,"/Results/data/eta_est"), η_pred)
            write_csv(string(exp_path,"/Results/data/beta_est"), β_pred)
            write_csv(string(exp_path,"/Results/data/h_est"), est_h_list)

            # plot_covariance(η_pred, β_pred, string(exp_path,"/Results/plots/"))

            plt_error = set_plot(fs, sz=(1650, 1250))
            StatsPlots.errorline!(est_h_list, label="Estimated height", dpi=400)
            Plots.plot!(gt_h, label="Ground truth height", dpi=400)
            Plots.xlabel!("Time (s)")
            Plots.ylabel!("Height")
            Plots.savefig(string(exp_path,"/Results/plots/h_est.pdf"))

            set_plot(fs, sz=(1650, 1250))
            Plots.plot!(nSplinex[1], nSpliney[1], label="", xaxis=false, yaxis=false) 
            Plots.xlims!(1100,1350)
            Plots.ylims!(420,1100)
            Plots.savefig(string(exp_path,"/Results/plots/noise_contour.pdf"))
        end

    elseif gt_viscosity_type == "bulk_viscosity"
        
        av_η::Float64 = 0.0
        obsBorderPts, nSplinex, nSpliney, pd = add_noise(ObsDataList, nFactor=0.0)
        mode::String = exp_params["mode"] # "single_window" or "multiple_window" or "full_time"
        
        viscosity_type = "constant"
        
        conditions = Conditions(camera_matrix=camera_matrix, obj_pose=obj_pose)
        model, scene = def_problem(r, h, ne_exp, η_gt, ndim, FunctionClass_u, nDof_u, FunctionClass_p, nDof_p, FunctionClass_x, β_gt, F, control, viscosity_type, 
        sim_time_exp, t_steps_exp)
        
        gt_time_frame::Int = round(Int,sim_time_gt/t_steps_gt,)
        sim_time_frame::Int = round(Int,sim_time_exp/t_steps_exp)
        window::Int = round(Int,gt_time_frame/sim_time_frame)
        
        set_file(string(exp_path,"/Results/plots"))
        println("Number of time windows: $window")
        
        time_windows, windows, data_ranges_, t_windows = set_time_window(1/t_steps_exp, obsBorderPts, method="fixed", window_size=sim_time_exp)
        _, splinexObs_win, _, _ = set_time_window(1/t_steps_exp, splinexObs, method="fixed", window_size=sim_time_exp)
        _, splineyObs_win, _, _ = set_time_window(1/t_steps_exp, splineyObs, method="fixed", window_size=sim_time_exp)

        println("Time windows: $(time_windows)")
        obs_time = sum(time_windows)

        if obs_time < sim_time_gt
            @warn "Observation time frame $obs_time is less than preset ground truth time frame $sim_time_gt, switching to observation time frame"
            gt_time_frame = obs_time
        end
        
        est_ηpList = Vector{Float64}(undef,gt_time_frame)
        avg_ηList = Vector{Float64}(undef,gt_time_frame)
        est_βpList = Vector{Float64}(undef,gt_time_frame)

        if mode == "single_window"
            @info "Optimizing over a single time window"
            ti = 1
            data_range_ = data_ranges_[ti]
            scene.sim_time = time_windows[ti]
            scene.cParam = F[data_range_]
            @info "Time window $(time_windows[ti])"

            println("Data frame : $(data_range_)")
            println("Time frame : $(scene.sim_time)")

            printstyled("Time window: $(ti), time frames: $(scene.sim_time)\n"; color = :blue)

            # est_μ_list, gradList, borderPts2DList, fields, pos3D, pos2D, splinex, spliney = simulate(model, scene, conditions)
            # println("Simulation completed for time window $(size(borderPts2DList))")
            # animate_fields(filepath=string(exp_path,"/Results/plots/init"), BorderNodes2D=borderPts2DList, pObs=splinexObs_win[ti], qObs=splineyObs_win[ti])

            obsBorderPts_t = windows[ti] # align the observation points with the simulation time

            println("observation window size: $(size(obsBorderPts_t,1))")

            println("Time frame : $data_range_")
            printstyled("Time window: $(ti)\n"; color = :green)

            if data_type != "physical"
                η_gt_ = model_gt.η[data_range_]
                av_η = mean(η_gt_)
                avg_ηList[data_range_] .= av_η
                printstyled("Average ground truth η in the window: $(av_η), ground truth β: $(β_gt)\n"; color = :green)
            end

            stats = fit_model(model, scene, conditions, obsBorderPts_t, θ, outliers=outlier_frames)

            est_ηpList[data_range_] .= stats["η"]
            est_βpList[data_range_] .= stats["β"]

            θ[1] = stats["η"]
            θ[2] = stats["β"]

            update_model!(model)
      
            iterList = stats["iterList"]
            costList = stats["cost_list"]
            ηpList = stats["ηList"]
            βpList = stats["βList"]
    
            # Plot the cost function with iterations
            set_plot(fs, sz=(1650, 1250))
            Plots.plot!(iterList, costList, label="Cost", marker=2, dpi=400, yscale=:log10, xminorgrid = :false, lw=3)
            Plots.xlabel!(L"\mathrm{Iterations}")
            Plots.ylabel!(L"\mathrm{Cost\;(px)}")
            Plots.xticks!(minimum(iterList):2:maximum(iterList))
            Plots.savefig(string(exp_path,"/Results/plots/cost_steps.pdf"))
    
            # Plot the cost function with iterations
            set_plot(fs, sz=(1650, 1250))
            Plots.plot!(iterList, costList, label="Cost", marker=2, dpi=400, yscale=:log10, xminorgrid = :false, lw=3)
            Plots.xlabel!(L"\mathrm{Iterations}")
            Plots.ylabel!(L"\mathrm{Cost\;(px)}")
            Plots.savefig(string(exp_path,"/Results/plots/cost_steps_log.pdf"))

            # Plot the cost function with iterations
            set_plot(fs, sz=(1650, 1250))
            Plots.plot!(iterList, ηpList, label="Cost", marker=2, dpi=400, yscale=:log10, xminorgrid = :false, lw=3)
            Plots.xlabel!(L"\mathrm{Iterations}")
            Plots.ylabel!(L"\mathrm{Cost\;(px)}")
            Plots.savefig(string(exp_path,"/Results/plots/eta_steps.pdf"))

            set_plot(fs, sz=(1650, 1250))
            Plots.plot!(iterList, βpList, label="Cost", marker=2, dpi=400, yscale=:log10, xminorgrid = :false, lw=3)
            Plots.xlabel!(L"\mathrm{Iterations}")
            Plots.ylabel!(L"\mathrm{Cost\;(px)}")
            Plots.xticks!(minimum(iterList):2:maximum(iterList))
            Plots.savefig(string(exp_path,"/Results/plots/beta_steps.pdf"))

            viscosity_type = "bulk_viscosity"
            est_model, est_scene = def_problem(r, h, ne_exp, η_gt, ndim, FunctionClass_u, nDof_u, FunctionClass_p, nDof_p, FunctionClass_x, β_gt, F[data_range_], control, viscosity_type, 
                                sim_time_exp, t_steps_exp)
            est_model.η = est_ηpList[data_range_] 
            est_scene.β = est_βpList[data_range_]
            est_μ_list, gradList, borderPts2DList, fields, pos3D, pos2D, splinex, spliney = simulate(est_model, est_scene, conditions)

            animate_fields(filepath=string(exp_path,"/Results/plots"), BorderNodes2D=borderPts2DList, pObs=splinexObs_win[ti], qObs=splineyObs_win[ti])
            est_h_list = get_height(est_μ_list, h)

            plt_η = set_plot(fs, sz=(1650, 1250))
            t = collect(range(start=t_steps_gt, stop=sim_time_gt, step=t_steps_gt))
            if data_type != "physical"
                Plots.plot!(plt_η, t, model_gt.η, label="Ground truth η(t)", dpi=400, lw=3)
                Plots.plot!(plt_η, t[data_range_], avg_ηList[data_range_], label="Average GT η in window", lw=3)
            end
            Plots.plot!(plt_η, t[data_range_], est_ηpList[data_range_], label="Estimated η(t)", lw=3)
            for t in t_windows
                Plots.vline!(plt_η, [t], color=:gray, lw=3, linestyle=:dash, label=false)
            end
            Plots.xlabel!(L"\mathrm{Time\;(s)}")
            Plots.ylabel!(L"\eta\mathrm{(t)\;(KPa\cdot s)}")
            Plots.savefig(plt_η, string(exp_path,"/Results/plots/η.pdf"))
        else
            println("Number of time windows: $(length(windows))")

            for ti::Int in 1:length(windows)
                data_range_ = data_ranges_[ti]
                scene.sim_time = time_windows[ti]
                scene.cParam = F[data_range_]
                obsBorderPts_t = windows[ti] # align the observation points with the simulation time

                println("Data frame : $(data_range_)")
                println("Time frame : $(scene.sim_time)")
                @info "Time window $(t_windows[ti])"
                printstyled("Time window: $(ti), time frames: $(scene.sim_time)\n"; color = :blue)
                println("observation window size: $(size(obsBorderPts_t,1))")

                if data_type != "physical"
                    η_gt_ = model_gt.η[data_range_]
                    av_η = mean(η_gt_)
                    avg_ηList[data_range_] .= av_η
                    printstyled("Average ground truth η in the window: $(av_η), ground truth β: $(β_gt)\n"; color = :green)
                end

                @info "Fitting model in time window $(ti)..."
                stats = fit_model(model, scene, conditions, obsBorderPts_t, θ, outliers=outlier_frames)

                est_ηpList[data_range_] .= stats["η"]
                est_βpList[data_range_] .= stats["β"]

                if ti == 1 && data_type == "physical"
                    θ[1] = stats["η"]*100
                else
                    θ[1] = stats["η"]
                end
                θ[2] = stats["β"]

                update_model!(model)
            end

            @info "Completed all time windows."

            plt_η = set_plot(fs, sz=(1650, 1250))
            t_full = collect(range(start=t_steps_gt, stop=sim_time_gt, step=t_steps_gt))
            if data_type != "physical"
                Plots.plot!(plt_η, t_full, model_gt.η, label="Ground truth η(t)", dpi=400, lw=3)
            end
            for ti::Int in 1:size(data_ranges_, 1)
                t = t_windows[ti]
                Plots.vline!(plt_η, [t], color=:gray, lw=3, linestyle=:dash, label=false)
            end
            Plots.xlabel!(L"\mathrm{Time\;(s)}")
            Plots.ylabel!(L"\eta\mathrm{(t)\;(KPa\cdot s)}")
            Plots.savefig(plt_η, string(exp_path,"/Results/plots/η_gt.pdf"))
            t_prev = 0.1
            for ti::Int in 1:length(data_ranges_)
                t = t_windows[ti]
                data_range_ = data_ranges_[ti]
                t_win = collect(range(start=t_prev, stop=t, step=t_steps_gt))
                if ti == 1
                    Plots.plot!(plt_η, t_win, est_ηpList[data_range_], label="Estimated η(t)", lw=3, color=:orange)
                    if data_type != "physical"
                        Plots.plot!(plt_η, t_win, avg_ηList[data_range_], label="Average GT η in window", lw=3, color=:gray)
                    end 
                else
                    Plots.plot!(plt_η, t_win, est_ηpList[data_range_], lw=3, color=:orange, label=false)
                    if data_type != "physical"
                        Plots.plot!(plt_η, t_win, avg_ηList[data_range_], lw=3, color=:gray, label=false)
                    end 
                end
                t_prev = t+t_steps_gt
            end
            Plots.savefig(plt_η, string(exp_path,"/Results/plots/η.pdf"))
            
            plt_β = set_plot(fs, sz=(1650, 1250))
            if data_type != "physical"
                Plots.hline!(plt_β, [β_gt], label="Ground truth β(t)", dpi=400, lw=3)
            end
            t_prev = 0.1
            for ti::Int in 1:length(data_ranges_)
                t = t_windows[ti]
                Plots.vline!(plt_β, [t], color=:gray, lw=3, linestyle=:dash, label=false)
                data_range_ = data_ranges_[ti]
                t_win = collect(range(start=t_prev, stop=t, step=t_steps_gt))
                if ti == 1
                    Plots.plot!(plt_β, t_win, est_βpList[data_range_], label="Estimated β(t)", lw=3, color=:orange)
                else
                    Plots.plot!(plt_β, t_win, est_βpList[data_range_], lw=3, color=:orange, label=false)
                end
                t_prev = t+t_steps_gt
            end
            Plots.xlabel!(L"\mathrm{Time\;(s)}")
            Plots.ylabel!(L"\beta\mathrm{(t)\;(mm^{-1})}")
            Plots.savefig(plt_β, string(exp_path,"/Results/plots/β.pdf"))
            
            viscosity_type = "bulk_viscosity"
            est_model, est_scene = def_problem(r, h, ne_exp, η_gt, ndim, FunctionClass_u, nDof_u, FunctionClass_p, nDof_p, FunctionClass_x, β_gt, F, control, viscosity_type, 
            sim_time_gt, t_steps_gt)
            est_model.η = est_ηpList
            
            animate_fields(filepath=string(exp_path,"/Results/plots"), p=splinex, q=spliney, pObs=splinexObs, qObs=splineyObs)
            
            gt_μ_list, gradList, borderPts2DList, fields, pos3D, pos2D, splinex, spliney = simulate(model_gt, scene_gt, conditions)
            est_μ_list, gradList, simBorderPts, splinex, spliney, pos2D = simulate(est_model, est_scene, conditions)

            if data_type != "physical"

                est_h_list = get_height(est_μ_list, h)
                gt_h_list = get_height(gt_μ_list, h)

                plt_h = set_plot(fs, sz=(1650, 1250))
                Plots.plot!(gt_h_list, label="Ground truth height", dpi=400, lw=3)
                Plots.plot!(est_h_list, label="Estimated height")
                Plots.xlabel!(L"\mathrm{Time\;(s)}")
                Plots.ylabel!(L"\mathrm{Height\;(mm)}")
                Plots.savefig(string(exp_path,"/Results/plots/h.pdf"))

                plt_error = set_plot(fs, sz=(1650, 1250))
                Plots.plot!(abs.(est_h_list-gt_h_list), label="Height estimation error", dpi=400, lw=3)
                Plots.savefig(string(exp_path,"/Results/plots/h_est_error.pdf"))
                Plots.xlabel!(L"\mathrm{Time\;(s)}")
                Plots.ylabel!(L"\mathrm{Height\;Error\;(px)}")
                Plots.savefig(string(exp_path,"/Results/plots/error.pdf"))

                write_csv(string(exp_path,"/Results/data/η_gt"), model_gt.η)
                write_csv(string(exp_path,"/Results/data/β_gt"), β_gt)
                write_csv(string(exp_path,"/Results/data/gt_h"), gt_h_list)
            end
        end

        write_csv(string(exp_path,"/Results/data/est_η"), est_ηpList)
        write_csv(string(exp_path,"/Results/data/est_β"), est_βpList)
        write_csv(string(exp_path,"/Results/data/est_h"), est_h_list)
        write_csv(string(exp_path,"/Results/data/avg_η"), avg_ηList)
        write_csv(string(exp_path,"/Results/data/window_data/time_windows"), time_windows)
        write_csv(string(exp_path,"/Results/data/window_data/t_windows"), t_windows)
        write_csv(string(exp_path,"/Results/data/window_data/data_ranges"), data_ranges_)
        write_csv(string(exp_path,"/Results/data/window_data/windows_sizes"), windows)
        

        write_data(string(exp_path,"/Results/data/sim_data/2D_surface_points"), pos2D)
        write_data(string(exp_path,"/Results/data/sim_data/3D_points"), pos3D)
        write_data(string(exp_path,"/Results/data/sim_data/motion_fields "), fields)
        write_data(string(exp_path,"/Results/data/sim_data/2D_border_points"), borderPts2DList)
    end
end

function compare_stats(filepath,filepath_gt)

    path_1 = []
    path_2 = []
    path_3 = []
    path_4 = []
    path_5 = []
    run_filepath = readdir(string(filepath,"/runs/"))

    for exp_path_ in run_filepath

        ne = exp_path_[4]

        if ne == '2'
            push!(path_1,exp_path_)
        elseif ne == '4'
            push!(path_2,exp_path_)
        elseif ne == '6'
            push!(path_3,exp_path_)
        elseif ne == '8'
            push!(path_4,exp_path_)
        elseif ne == "16"
            push!(path_5,exp_path_)
        end
    end

    if length(path_1) != 0
        cost_list_list = []
        est_η_list = []
        est_β_list = []
        h_error_list = []
        iter_list = []
        η_gt = 0
        β_gt = 0

        plot_η = set_plot(fs, sz=(1650, 1250))
        plot_β = set_plot(fs, sz=(1650, 1250))
        plot_cost = set_plot(fs, sz=(1650, 1250))
        plot_error = set_plot(fs, sz=(1650, 1250))

        for file in path_1
            exp_path = string(filepath,"/runs/",file)
            sim_params = read_json(string(exp_path,"/Results/data/sim_params.json")) 

            sim_time = 10.0 # simulation time in seconds 
            steps = 25.0 # number of time steps
            t_steps = sim_time/steps

            time = collect(Float64, range(start=0, stop=sim_time, step=t_steps))

            viscosity_type = sim_params["viscosity_type"]

            if viscosity_type == "constant"
                η_gt = sim_params["η_gt"]
                β_gt = sim_params["β_gt"]

                est_η = readdlm(string(exp_path,"/Results/data/η.csv"), ',', Float64)
                est_β = readdlm(string(exp_path,"/Results/data/β.csv"), ',', Float64)
                est_h = readdlm(string(exp_path,"/Results/data/est_h.csv"), ',', Float64)
                gt_h = readdlm(string(exp_path,"/Results/data/gt_h.csv"), ',', Float64)

                h_error = abs.(est_h-gt_h)

                stats = read_json(string(exp_path,"/Results/data/stats.json")) 

                costList = stats["cost_list"]
                iterList = stats["iterList"]

                label = " "
                if file[1:2] == "Q2"
                    label = "Lagrange basis"
                elseif file[1:2] == "S2"
                    label = "NURBS basis"
                end

                plot!(plot_η, iterList, est_η, label=label, dpi=400, lw=3)
                xlabel!(plot_η, L"\mathrm{Iterations}")
                ylabel!(plot_η, L"\eta\;\mathrm{(KPa\cdot s)}")

                plot!(plot_β, iterList, est_β, label=label, dpi=400, lw=3)
                xlabel!(plot_β, L"\mathrm{Iterations}")
                ylabel!(plot_β, L"\beta\;\mathrm{mm^{-1}}")

                plot!(plot_cost, iterList, costList, label=label, dpi=400, lw=3)
                xlabel!(plot_cost, L"\mathrm{Iterations}")
                ylabel!(plot_cost, L"\mathrm{Cost\;(px)}")

                plot!(plot_cost_log, iterList, costList, label=label, dpi=400, lw=3, yscale=:log10)
                xlabel!(plot_cost_log, L"\mathrm{Iterations}")
                ylabel!(plot_cost_log, L"\mathrm{Cost\;(px)}")

                plot!(plot_error, time, h_error, label=label, dpi=400, lw=3)
                xlabel!(plot_error, L"\mathrm{Time\;(s)}")
                ylabel!(plot_error, L"\mathrm{Height\;Error\;(mm)}")

                plot!(plot_h, time, est_h, label=label, dpi=400, lw=3)
                xlabel!(plot_h, L"\mathrm{Time\;(s)}")
                ylabel!(plot_h, L"\mathrm{Height\;(mm)}")

                push!(cost_list_list,costList)
                push!(est_η_list,est_η)
                push!(est_β_list,est_β)
                push!(h_error_list,h_error)
                push!(iter_list,iterList)
            end
        end

        res_path = string(filepath,"/runs/post_analysis/ne_2")
        set_file(res_path)

        Plots.hline!(plot_η, [η_gt], label="Ground truth η", lw=3)
        savefig(plot_η, string(res_path,"/eta.pdf"))

        Plots.hline!(plot_β, [β_gt], label="Ground truth β", lw=3)
        savefig(plot_β, string(res_path,"/beta.pdf"))

        savefig(plot_cost, string(res_path,"/convergence.pdf"))
        savefig(plot_cost_log, string(res_path,"/convergence_log.pdf"))
        savefig(plot_error, string(res_path,"/error.pdf"))

        plot!(plot_h, time, gt_h, label="Ground truth height", dpi=400, lw=3)
        savefig(plot_h, string(res_path,"/h.pdf"))
    end

    if length(path_2) != 0
        println(path_2)
        cost_list_list = []
        est_η_list = []
        est_β_list = []
        h_error_list = []
        iter_list = []
        η_gt = 0
        β_gt = 0

        x = []
        y = []

        plot_η = set_plot(fs, sz=(1650, 1250))
        plot_β = set_plot(fs, sz=(1650, 1250))
        plot_cost = set_plot(fs, sz=(1650, 1250))
        plot_cost_log = set_plot(fs, sz=(1650, 1250))
        plot_error = set_plot(fs, sz=(1650, 1250))
        plot_h = set_plot(fs, sz=(1650, 1250))

        for file in path_2
            exp_path = string(filepath,"/runs/",file)
            println(exp_path)
            sim_params = read_json(string(exp_path,"/Results/data/sim_params.json")) 

            ObsDataList, splinex, spliney = read_csv(string(exp_path,"/Results/data/sim_data/2D_border_points/"))  

            push!(x, splinex)
            push!(y, spliney)

            sim_time = 10.0 # simulation time in seconds 
            steps = 25.0 # number of time steps
            t_steps = sim_time/steps

            time = collect(Float64, range(start=0, stop=sim_time, step=t_steps))

            viscosity_type = sim_params["viscosity_type"]

            if viscosity_type == "constant"
                η_gt = sim_params["η_gt"]
                β_gt = sim_params["β_gt"]

                est_η = readdlm(string(exp_path,"/Results/data/η.csv"), ',', Float64)
                est_β = readdlm(string(exp_path,"/Results/data/β.csv"), ',', Float64)
                est_h = readdlm(string(exp_path,"/Results/data/est_h.csv"), ',', Float64)
                gt_h = readdlm(string(exp_path,"/Results/data/gt_h.csv"), ',', Float64)

                h_error = abs.(est_h-gt_h)

                stats = read_json(string(exp_path,"/Results/data/stats.json")) 

                costList = stats["cost_list"]
                iterList = stats["iterList"]

                label = " "
                if file[1:2] == "Q2"
                    label = "Lagrange basis"
                elseif file[1:2] == "S2"
                    label = "NURBS basis"
                end

                plot!(plot_η, iterList, est_η, label=label, dpi=400, lw=3)
                xlabel!(plot_η, L"\mathrm{Iterations}")
                ylabel!(plot_η, L"\eta\;\mathrm{(KPa\cdot s)}")

                plot!(plot_β, iterList, est_β, label=label, dpi=400, lw=3)
                xlabel!(plot_β, L"\mathrm{Iterations}")
                ylabel!(plot_β,L"\beta\;\mathrm{mm^{-1}}")

                plot!(plot_cost, iterList, costList, label=label, dpi=400, lw=3)
                xlabel!(plot_cost, L"\mathrm{Iterations}")
                ylabel!(plot_cost, L"\mathrm{Cost\;(px)}")

                plot!(plot_cost_log, iterList, costList, label=label, dpi=400, lw=3, yscale=:log10)
                xlabel!(plot_cost_log, L"\mathrm{Iterations}")
                ylabel!(plot_cost_log, L"\mathrm{Cost\;(px)}")

                plot!(plot_error, time, h_error, label=label, dpi=400, lw=3)
                xlabel!(plot_error, L"\mathrm{Time\;(s)}")
                ylabel!(plot_error, L"\mathrm{Height\;Error\;(mm)}")

                plot!(plot_h, time, est_h, label=label, dpi=400, lw=3)
                xlabel!(plot_h, L"\mathrm{Time\;(s)}")
                ylabel!(plot_h, L"\mathrm{Height\;(mm)}")

                push!(cost_list_list,costList)
                push!(est_η_list,est_η)
                push!(est_β_list,est_β)
                push!(h_error_list,h_error)
                push!(iter_list,iterList)
            end
        end

        ObsDataList, splinex_gt, spliney_gt = read_csv(string(filepath_gt,"/data/sim_data/2D_border_points"))  

        res_path = string(filepath,"/runs/post_analysis/ne_4")
        set_file(res_path)

        animate_fields(filepath=string(res_path,"/plots"), p=x[1], q=y[1], pObs=x[2], qObs=y[2], pgt= splinex_gt, qgt=spliney_gt)

        Plots.hline!(plot_η, [η_gt], label="Ground truth η", lw=3)
        savefig(plot_η, string(res_path,"/eta.pdf"))

        Plots.hline!(plot_β, [β_gt], label="Ground truth β", lw=3)
        savefig(plot_β, string(res_path,"/beta.pdf"))

        savefig(plot_cost, string(res_path,"/convergence.pdf"))
        savefig(plot_cost_log, string(res_path,"/convergence_log.pdf"))
        savefig(plot_error, string(res_path,"/error.pdf"))

        plot!(plot_h, time, gt_h, label="Ground truth height", dpi=400, lw=3)
        savefig(plot_h, string(res_path,"/h.pdf"))

    end

    if length(path_3) != 0
       cost_list_list = []
        est_η_list = []
        est_β_list = []
        h_error_list = []
        iter_list = []
        η_gt = 0
        β_gt = 0
        gt_h = 0

        plot_η = set_plot(fs, sz=(1650, 1250))
        plot_β = set_plot(fs, sz=(1650, 1250))
        plot_cost = set_plot(fs, sz=(1650, 1250))
        plot_cost_log = set_plot(fs, sz=(1650, 1250))
        plot_error = set_plot(fs, sz=(1650, 1250))
        plot_h = set_plot(fs, sz=(1650, 1250))

        sim_time = 10.0 # simulation time in seconds 
        steps = 25.0 # number of time steps
        t_steps = sim_time/steps

        time = collect(Float64, range(start=0, stop=sim_time, step=t_steps))

        for file in path_3
            exp_path = string(filepath,"/runs/",file)
            sim_params = read_json(string(exp_path,"/Results/data/sim_params.json")) 

            viscosity_type = sim_params["viscosity_type"]

            if viscosity_type == "constant"
                η_gt = sim_params["η_gt"]
                β_gt = sim_params["β_gt"]

                est_η = readdlm(string(exp_path,"/Results/data/η.csv"), ',', Float64)
                est_β = readdlm(string(exp_path,"/Results/data/β.csv"), ',', Float64)
                est_h = readdlm(string(exp_path,"/Results/data/est_h.csv"), ',', Float64)
                gt_h = readdlm(string(exp_path,"/Results/data/gt_h.csv"), ',', Float64)

                h_error = abs.(est_h-gt_h)

                stats = read_json(string(exp_path,"/Results/data/stats.json")) 

                costList = stats["cost_list"]
                iterList = stats["iterList"]

                label = " "
                if file[1:2] == "Q2"
                    label = "Lagrange basis"
                elseif file[1:2] == "S2"
                    label = "NURBS basis"
                end
                plot!(plot_η, iterList, est_η, label=label, dpi=400, lw=3)
                xlabel!(plot_η, L"\mathrm{Iterations}")
                ylabel!(plot_η, L"\eta\;\mathrm{(KPa\cdot s)}")

                plot!(plot_β, iterList, est_β, label=label, dpi=400, lw=3)
                xlabel!(plot_β, L"\mathrm{Iterations}")
                ylabel!(plot_β,L"\beta\;\mathrm{mm^{-1}}")

                plot!(plot_cost, iterList, costList, label=label, dpi=400, lw=3)
                xlabel!(plot_cost, L"\mathrm{Iterations}")
                ylabel!(plot_cost, L"\mathrm{Cost\;(px)}")

                plot!(plot_cost_log, iterList, costList, label=label, dpi=400, lw=3, yscale=:log10)
                xlabel!(plot_cost_log, L"\mathrm{Iterations}")
                ylabel!(plot_cost_log, L"\mathrm{Cost\;(px)}")

                plot!(plot_error, time, h_error, label=label, dpi=400, lw=3)
                xlabel!(plot_error, L"\mathrm{Time\;(s)}")
                ylabel!(plot_error, "Height Error (mm)")

                plot!(plot_h, time, est_h, label=label, dpi=400, lw=3)
                xlabel!(plot_h, L"\mathrm{Time\;(s)}")
                ylabel!(plot_h, L"\mathrm{Height\;(mm)}")

                push!(cost_list_list,costList)
                push!(est_η_list,est_η)
                push!(est_β_list,est_β)
                push!(h_error_list,h_error) 
                push!(iter_list,iterList)
            end
        end

        res_path = string(filepath,"/runs/post_analysis/ne_8")
        set_file(res_path)

        Plots.hline!(plot_η, [η_gt], label="Ground truth η", lw=3)
        savefig(plot_η, string(res_path,"/eta.pdf"))

        Plots.hline!(plot_β, [β_gt], label="Ground truth β", lw=3)
        savefig(plot_β, string(res_path,"/beta.pdf"))

        savefig(plot_cost, string(res_path,"/convergence.pdf"))
        savefig(plot_cost_log, string(res_path,"/convergence_log.pdf"))
        savefig(plot_error, string(res_path,"/error.pdf"))

        plot!(plot_h, time, gt_h, label="Ground truth height", dpi=400, lw=3)
        savefig(plot_h, string(res_path,"/h.pdf"))

    end

    if length(path_4) != 0
        cost_list_list = []
        est_η_list = []
        est_β_list = []
        h_error_list = []
        iter_list = []
        η_gt = 0
        β_gt = 0
        gt_h = 0

        plot_η = set_plot(fs, sz=(1650, 1250))
        plot_β = set_plot(fs, sz=(1650, 1250))
        plot_cost = set_plot(fs, sz=(1650, 1250))
        plot_cost_log = set_plot(fs, sz=(1650, 1250))
        plot_error = set_plot(fs, sz=(1650, 1250))
        plot_h = set_plot(fs, sz=(1650, 1250))

        sim_time = 10.0 # simulation time in seconds 
        steps = 25.0 # number of time steps
        t_steps = sim_time/steps

        time = collect(Float64, range(start=0, stop=sim_time, step=t_steps))

        for file in path_4
            exp_path = string(filepath,"/runs/",file)
            sim_params = read_json(string(exp_path,"/Results/data/sim_params.json")) 

            viscosity_type = sim_params["viscosity_type"]

            if viscosity_type == "constant"
                η_gt = sim_params["η_gt"]
                β_gt = sim_params["β_gt"]

                est_η = readdlm(string(exp_path,"/Results/data/η.csv"), ',', Float64)
                est_β = readdlm(string(exp_path,"/Results/data/β.csv"), ',', Float64)
                est_h = readdlm(string(exp_path,"/Results/data/est_h.csv"), ',', Float64)
                gt_h = readdlm(string(exp_path,"/Results/data/gt_h.csv"), ',', Float64)

                h_error = abs.(est_h-gt_h)

                stats = read_json(string(exp_path,"/Results/data/stats.json")) 

                costList = stats["cost_list"]
                iterList = stats["iterList"]

                label = " "
                if file[1:2] == "Q2"
                    label = "Lagrange basis"
                elseif file[1:2] == "S2"
                    label = "NURBS basis"
                end
                plot!(plot_η, iterList, est_η, label=label, dpi=400, lw=3)
                xlabel!(plot_η, L"\mathrm{Iterations}")
                ylabel!(plot_η, L"\eta\;\mathrm{(KPa\cdot s)}")

                plot!(plot_β, iterList, est_β, label=label, dpi=400, lw=3)
                xlabel!(plot_β, L"\mathrm{Iterations}")
                ylabel!(plot_β,L"\beta\;\mathrm{mm^{-1}}")

                plot!(plot_cost, iterList, costList, label=label, dpi=400, lw=3)
                xlabel!(plot_cost, L"\mathrm{Iterations}")
                ylabel!(plot_cost, L"\mathrm{Cost\;(px)}")

                plot!(plot_cost_log, iterList, costList, label=label, dpi=400, lw=3, yscale=:log10)
                xlabel!(plot_cost_log, L"\mathrm{Iterations}")
                ylabel!(plot_cost_log, L"\mathrm{Cost\;(px)}")

                plot!(plot_error, time, h_error, label=label, dpi=400, lw=3)
                xlabel!(plot_error, L"\mathrm{Time\;(s)}")
                ylabel!(plot_error, "Height Error (mm)")

                plot!(plot_h, time, est_h, label=label, dpi=400, lw=3)
                xlabel!(plot_h, L"\mathrm{Time\;(s)}")
                ylabel!(plot_h, L"\mathrm{Height\;(mm)}")

                push!(cost_list_list,costList)
                push!(est_η_list,est_η)
                push!(est_β_list,est_β)
                push!(h_error_list,h_error)
                push!(iter_list,iterList)
            end
        end

        res_path = string(filepath,"/runs/post_analysis/ne_8")
        set_file(res_path)

        Plots.hline!(plot_η, [η_gt], label="Ground truth η", lw=3)
        savefig(plot_η, string(res_path,"/eta.pdf"))

        Plots.hline!(plot_β, [β_gt], label="Ground truth β", lw=3)
        savefig(plot_β, string(res_path,"/beta.pdf"))

        savefig(plot_cost, string(res_path,"/convergence.pdf"))
        savefig(plot_cost_log, string(res_path,"/convergence_log.pdf"))
        savefig(plot_error, string(res_path,"/error.pdf"))

        plot!(plot_h, time, gt_h, label="Ground truth height", dpi=400, lw=3)
        savefig(plot_h, string(res_path,"/h.pdf"))

    end

    if length(path_5) != 0
        cost_list_list = []
        est_η_list = []
        est_β_list = []
        h_error_list = []
        iter_list = []
        η_gt = 0
        β_gt = 0

        plot_η = set_plot(fs, sz=(1650, 1250))
        plot_β = set_plot(fs, sz=(1650, 1250))
        plot_cost = set_plot(fs, sz=(1650, 1250))
        plot_error = set_plot(fs, sz=(1650, 1250))

        for file in path_5
            exp_path = string(filepath,"/runs/",file)
            sim_params = read_json(string(exp_path,"/Results/data/sim_params.json")) 

            sim_time = 10.0 # simulation time in seconds 
            steps = 25.0 # number of time steps
            t_steps = sim_time/steps

            time = collect(Float64, range(start=0, stop=sim_time, step=t_steps))

            viscosity_type = sim_params["viscosity_type"]

            if viscosity_type == "constant"
                η_gt = sim_params["η_gt"]
                β_gt = sim_params["β_gt"]

                est_η = readdlm(string(exp_path,"/Results/data/η.csv"), ',', Float64)
                est_β = readdlm(string(exp_path,"/Results/data/β.csv"), ',', Float64)
                est_h = readdlm(string(exp_path,"/Results/data/est_h.csv"), ',', Float64)
                gt_h = readdlm(string(exp_path,"/Results/data/gt_h.csv"), ',', Float64)

                h_error = abs.(est_h-gt_h)

                stats = read_json(string(exp_path,"/Results/data/stats.json")) 

                costList = stats["cost_list"]
                iterList = stats["iterList"]

                label = " "
                if file[1:2] == "Q2"
                    label = "Lagrange basis"
                elseif file[1:2] == "S2"
                    label = "NURBS basis"
                end

                plot!(plot_η, iterList, est_η, label=label, dpi=400, lw=3)
                xlabel!(plot_η, L"\mathrm{Iterations}")
                ylabel!(plot_η, L"\eta\;\mathrm{(KPa\cdot s)}")

                plot!(plot_β, iterList, est_β, label=label, dpi=400, lw=3)
                xlabel!(plot_β, L"\mathrm{Iterations}")
                ylabel!(plot_β, L"\beta\;\mathrm{mm^{-1}}")

                plot!(plot_cost, iterList, costList, label=label, dpi=400, lw=3)
                xlabel!(plot_cost, L"\mathrm{Iterations}")
                ylabel!(plot_cost, L"\mathrm{Cost\;(px)}")

                plot!(plot_cost_log, iterList, costList, label=label, dpi=400, lw=3, yscale=:log10)
                xlabel!(plot_cost_log, L"\mathrm{Iterations}")
                ylabel!(plot_cost_log, L"\mathrm{Cost\;(px)}")

                plot!(plot_error, time, h_error, label=label, dpi=400, lw=3)
                xlabel!(plot_error, L"\mathrm{Time\;(s)}")
                ylabel!(plot_error, L"\mathrm{Height\;Error\;(mm)}")

                plot!(plot_h, time, est_h, label=label, dpi=400, lw=3)
                xlabel!(plot_h, L"\mathrm{Time\;(s)}")
                ylabel!(plot_h, L"\mathrm{Height\;(mm)}")

                push!(cost_list_list,costList)
                push!(est_η_list,est_η)
                push!(est_β_list,est_β)
                push!(h_error_list,h_error)
                push!(iter_list,iterList)
            end
        end

        res_path = string(filepath,"/runs/post_analysis/ne_2")
        set_file(res_path)

        Plots.hline!(plot_η, [η_gt], label="Ground truth η", lw=3)
        savefig(plot_η, string(res_path,"/eta.pdf"))

        Plots.hline!(plot_β, [β_gt], label="Ground truth β", lw=3)
        savefig(plot_β, string(res_path,"/beta.pdf"))

        savefig(plot_cost, string(res_path,"/convergence.pdf"))
        savefig(plot_cost_log, string(res_path,"/convergence_log.pdf"))
        savefig(plot_error, string(res_path,"/error.pdf"))

        plot!(plot_h, time, gt_h, label="Ground truth height", dpi=400, lw=3)
        savefig(plot_h, string(res_path,"/h.pdf"))
    end
end

function replot(filepath, filepath_gt)
    
    elem_size_folders = readdir(string(filepath))
    
    sim_params = read_json(string(filepath_gt,"/data/sim_params.json")) 

    ObsDataList, splinexObs, splineyObs = read_csv(string(filepath_gt,"/data/img_data/contour_data"))  
    
    sim_time = sim_params["simulation_time"]    
    t_steps = sim_params["time_steps"]
    viscosity_type = sim_params["viscosity_type"]
    time = collect(Float64, range(start=0, stop=sim_time, step=t_steps))

    for elem_size_folder in elem_size_folders
        if elem_size_folder == "post_analysis"
            continue
        end
        
        sim_time_folders = readdir(string(filepath,"/",elem_size_folder,"/"))

        # Iterate over simulation time folders
        for sim_time_folder in sim_time_folders
            if sim_time_folder == "post_analysis_time" || sim_time_folder == "Results"
                continue
            end

            noise_folders = readdir(string(filepath,"/",elem_size_folder,"/",sim_time_folder,"/"))
            for noise_folder in noise_folders
                if noise_folder == "post_analysis_noise" || noise_folder == "Results"
                    continue
                end
                exp_path = string(filepath, elem_size_folder, "/", sim_time_folder, "/", noise_folder)
                println("Comparing experiments in: $exp_path")

                exp_params = read_json(string(exp_path,"/Results/data/experiment_parameters.json"))
                noise_level = exp_params["noise_level"]

                obsBorderPts, nSplinex, nSpliney, pd = add_noise(ObsDataList, nFactor=0.0)

                if viscosity_type == "constant"
                    
                    if noise_level == 0.0
                        η_gt = sim_params["η"]
                        β_gt = sim_params["β"]

                        est_η = readdlm(string(exp_path,"/Results/data/η.csv"), ',', Float64)
                        est_β = readdlm(string(exp_path,"/Results/data/β.csv"), ',', Float64)
                        est_h = readdlm(string(exp_path,"/Results/data/est_h.csv"), ',', Float64)
                        gt_h = readdlm(string(exp_path,"/Results/data/gt_h.csv"), ',', Float64)
                        cost_iter = readdlm(string(exp_path,"/Results/data/cost_iter.csv"), ',', Float64)

                        stats = read_json(string(exp_path,"/Results/data/stats.json")) 
                        contour_plot_params = read_json(string(exp_path,"/Results/data/contour_plot_params.json")) 
                        exp_params = read_json(string(exp_path,"/Results/data/experiment_parameters.json"))

                        sim_time_exp = exp_params["sim_time_exp"]

                        symBorderPts, splinex, spliney = read_csv(string(exp_path,"/Results/data/sim_data/2D_border_points"))

                        if sim_time_exp != sim_time
                            @warn "Simulation time in experiment parameters ($sim_time_exp) does not match that in ground truth ($sim_time)."
                            if sim_time_exp < sim_time
                                @warn "Truncating obsBorderPts to match sim_time_exp."
                                obsBorderPts = obsBorderPts[1:Int(sim_time_exp/ t_steps)+1, :]
                            else
                                @warn " Truncating simBorderPts to match gt sim time."
                                symBorderPts = symBorderPts[1:Int(sim_time/ t_steps)+1, :]
                            end
                            println("obsBorderPts length: ", length(obsBorderPts))
                            println("symBorderPts length: ", length(symBorderPts))
                        end
                        
        
                        d, pairs = closest_point(symBorderPts, obsBorderPts)
        
                        plt_cnt_error = set_plot(fs, sz=(1650, 1250))
                        Plots.plot!(plt_cnt_error, d, label="Closest point distance error", dpi=400, lw=3)
                        Plots.xlabel!(L"\mathrm{Time\;(s)}")
                        Plots.ylabel!(L"\mathrm{Closest\;Point\;Distance\;(px)}")
                        Plots.savefig(plt_cnt_error, string(exp_path,"/Results/plots/closest_point_distance_error.pdf"))

                        # Coerce contour parameters to concrete arrays (be defensive against scalar values)
                        function _to_vector(x)
                            if isa(x, Number)
                                return Float64[x]
                            elseif isa(x, AbstractVector)
                                return Float64.(collect(x))
                            else
                                try
                                    return Float64.(collect(x))
                                catch err
                                    error("Cannot coerce input to vector: $err")
                                end
                            end
                        end

                        ηList = _to_vector(contour_plot_params["η_list"])
                        βList = _to_vector(contour_plot_params["β_list"])
                        # ensure CostMat is a dense Float64 matrix
                        CostMat = try
                            Array(contour_plot_params["cost_mat"]) |> x -> Float64.(x)
                        catch err
                            error("Failed to coerce cost_mat to dense Float64 matrix: $err")
                        end

                        # --- Post-process: compute Hessian-based principal directions from CostMat ---
                        try
                            # Ensure CostMat shape matches (length(ηList), length(βList))
                            nx = length(ηList); ny = length(βList)
                            sz = size(CostMat)
                            if ndims(CostMat) == 1 && length(CostMat) == nx*ny
                                CostMat = reshape(CostMat, nx, ny)
                            elseif ndims(CostMat) == 2
                                if sz[1] == ny && sz[2] == nx
                                    # CostMat is (ny, nx) => transpose to (nx, ny)
                                    CostMat = CostMat'
                                elseif sz[1] == nx && sz[2] == ny
                                    # already in expected orientation
                                else
                                    @warn "CostMat shape $(sz) does not match (nx,ny)=($(nx),$(ny)). Attempting to reshape if possible."
                                    if prod(sz) == nx*ny
                                        CostMat = reshape(vec(CostMat), nx, ny)
                                    end
                                end
                            end
                            # Prepare interpolation grid (rows->β, cols->η)
                            Z_for_interp = CostMat'

                            # locate minimum in the cost matrix (η-index, β-index)
                            minval, minidx = findmin(CostMat)
                            i0 = Int(minidx[1]); j0 = Int(minidx[2])
                            # guard against out-of-range indices
                            ni = length(ηList); nj = length(βList)
                            if i0 < 1 || i0 > ni || j0 < 1 || j0 > nj
                                error("Minimum index from CostMat (i0=$i0, j0=$j0) is out of bounds for ηList/βList lengths (", ni, ",", nj, ")")
                            end
                            η0 = ηList[i0]; β0 = βList[j0]

                            # choose interior index for finite-difference Hessian computation
                            ii = clamp(i0, 2, max(2, ni-1))
                            jj = clamp(j0, 2, max(2, nj-1))
                            if ii != i0 || jj != j0
                                @warn "Minimum on boundary; using nearby interior index (ii,jj)=($(ii),$(jj)) for Hessian"
                            end

                            # grid spacings (assume near-uniform)
                            dη = (ni > 1) ? (ηList[min(2,ni)] - ηList[1]) : 1.0
                            dβ = (nj > 1) ? (βList[min(2,nj)] - βList[1]) : 1.0

                            C = CostMat
                            # centered finite differences for second derivatives
                            d2C_dη2 = (C[ii+1,jj] - 2C[ii,jj] + C[ii-1,jj]) / (dη^2)
                            d2C_dβ2 = (C[ii,jj+1] - 2C[ii,jj] + C[ii,jj-1]) / (dβ^2)
                            d2C_dηdβ = (C[ii+1,jj+1] - C[ii+1,jj-1] - C[ii-1,jj+1] + C[ii-1,jj-1]) / (4.0*dη*dβ)

                            H = [d2C_dη2 d2C_dηdβ; d2C_dηdβ d2C_dβ2]
                            ev = eigen(Symmetric(H))
                            evals = ev.values
                            evecs = ev.vectors

                            # steepest = largest eigenvalue, flattest = smallest
                            idx_max = argmax(evals); idx_min = argmin(evals)
                            v_steep = evecs[:, idx_max] / norm(evecs[:, idx_max])
                            v_flat  = evecs[:, idx_min] / norm(evecs[:, idx_min])

                            # param-space sampling along directions (normalized) centered at minimum
                            span_η = maximum(ηList) - minimum(ηList)
                            span_β = maximum(βList) - minimum(βList)
                            diag_span = sqrt(span_η^2 + span_β^2)
                            tvec = collect(range(-0.6*diag_span, stop=0.6*diag_span, length=401))

                            etas_steep = η0 .+ tvec .* v_steep[1]
                            betas_steep = β0 .+ tvec .* v_steep[2]
                            etas_flat = η0 .+ tvec .* v_flat[1]
                            betas_flat = β0 .+ tvec .* v_flat[2]

                            # keep samples inside the original grid extents
                            minη, maxη = minimum(ηList), maximum(ηList)
                            minβ, maxβ = minimum(βList), maximum(βList)
                            mask_steep = (etas_steep .>= minη) .& (etas_steep .<= maxη) .& (betas_steep .>= minβ) .& (betas_steep .<= maxβ)
                            mask_flat  = (etas_flat  .>= minη) .& (etas_flat  .<= maxη) .& (betas_flat  .>= minβ) .& (betas_flat  .<= maxβ)

                            t_steep = tvec[mask_steep]; t_flat = tvec[mask_flat]
                            etas_steep = etas_steep[mask_steep]; betas_steep = betas_steep[mask_steep]
                            etas_flat  = etas_flat[mask_flat];  betas_flat  = betas_flat[mask_flat]

                            # interpolate cost along these samples
                            zs_steep = interp_z_at(etas_steep, betas_steep, ηList, βList, Z_for_interp)
                            zs_flat  = interp_z_at(etas_flat,  betas_flat,  ηList, βList, Z_for_interp)

                            # Contour with overlaid directions and minimum (use est_η/est_β as iteration path)
                            # ensure estimation paths are 1D Float64 vectors
                            est_η = vec(Float64.(collect(est_η)))
                            est_β = vec(Float64.(collect(est_β)))

                            @info "Post-process grid sizes: len(ηList)=$(length(ηList)), len(βList)=$(length(βList)), CostMat_size=$(size(CostMat))"
                            plt_dirs = set_plot(fs, sz=(1650, 1250))
                            Plots.contour!(plt_dirs, ηList, βList, CostMat', color=:turbo, fill=false, levels=100, dpi=400, legend=:outerbottom, legend_column=3, bottom_margin = -15mm)
                            Plots.plot!(plt_dirs, est_η, est_β, label="Estimations", ms=:4, m=:x, color=:red, lw=3, legend=:outerbottom, legend_column=3, bottom_margin = -15mm)
                            Plots.plot!(plt_dirs, etas_steep, betas_steep, label = "Steepest dir", lw=3, color=:black, legend=:outerbottom, legend_column=3, bottom_margin = -15mm)
                            Plots.plot!(plt_dirs, etas_flat,  betas_flat,  label = "Flattest dir",  lw=3, color=:magenta, legend=:outerbottom, legend_column=3, bottom_margin = -15mm)
                            Plots.scatter!(plt_dirs, [η0], [β0], label="Minimum Cost", ms=15, m=:star5, color=:black, legend=:outerbottom, legend_column=3, bottom_margin = -15mm)
                            Plots.scatter!(plt_dirs, [η_gt], [β_gt], label="Ground truth", ms=:10, m=:circle, color=:indianred2, lw=3, legend=:outerbottom, legend_column=3, bottom_margin = -15mm)
                            Plots.xlabel!(plt_dirs, L"\eta\;\mathrm{(KPa\cdot s)}")
                            Plots.ylabel!(plt_dirs, L"\beta\;\mathrm{(L/mm^{-1})}")
                            Plots.savefig(plt_dirs, string(exp_path, "/Results/plots/cost_surface_with_directions.pdf"))
                            
                            plt_cont = set_plot(fs, sz=(1650, 1250))
                            Plots.contour!(plt_cont, ηList, βList, CostMat', color=:turbo, fill=false, levels=100, dpi=400, legend=:outerbottom, legend_column=3, bottom_margin = -15mm)
                            Plots.plot!(plt_cont, est_η, est_β, label="Estimations", ms=:4, m=:x, color=:red, lw=3, legend=:outerbottom, legend_column=3, bottom_margin = -15mm)
                            Plots.scatter!(plt_cont, [η0], [β0], label="Minimum Cost", ms=15, m=:star5, color=:black, legend=:outerbottom, legend_column=3, bottom_margin = -15mm)
                            Plots.scatter!(plt_cont, [η_gt], [β_gt], label="Ground truth", ms=:10, m=:circle, color=:indianred2, lw=3, legend=:outerbottom, legend_column=3, bottom_margin = -15mm)
                            Plots.xlabel!(plt_cont, L"\eta\;\mathrm{(KPa\cdot s)}")
                            Plots.ylabel!(plt_cont, L"\beta\;\mathrm{(L/mm^{-1})}")
                            Plots.savefig(plt_cont, string(exp_path, "/Results/plots/cost_surface.pdf"))

                            plt_cont = set_plot(fs, sz=(1650, 1250))
                            Plots.contourf!(plt_cont, ηList, βList, CostMat', color=:turbo, fill=false, levels=100, dpi=400, legend=:outerbottom, legend_column=3, bottom_margin = -15mm)
                            Plots.plot!(plt_cont, est_η, est_β, label="Estimations", ms=:4, m=:x, color=:red, lw=3, legend=:outerbottom, legend_column=3, bottom_margin = -15mm)
                            Plots.scatter!(plt_cont, [η0], [β0], label="Minimum Cost", ms=15, m=:star5, color=:black, legend=:outerbottom, legend_column=3, bottom_margin = -15mm)
                            Plots.scatter!(plt_cont, [η_gt], [β_gt], label="Ground truth", ms=:10, m=:circle, color=:indianred2, lw=3, legend=:outerbottom, legend_column=3, bottom_margin = -15mm)
                            Plots.xlabel!(plt_cont, L"\eta\;\mathrm{(KPa\cdot s)}")
                            Plots.ylabel!(plt_cont, L"\beta\;\mathrm{(L/mm^{-1})}")
                            Plots.savefig(plt_cont, string(exp_path, "/Results/plots/cost_surface_iter.pdf"))

                            # 2D slices: cost vs distance along the two directions
                            plt_slices = set_plot(fs, sz=(1650, 1250))
                            if length(t_steep) > 0 && length(zs_steep) == length(t_steep)
                                Plots.plot!(plt_slices, t_steep, zs_steep, label = "Steepest direction", lw=3, color=:black, legend=:outerbottom, legend_column=2, bottom_margin = -15mm)
                            else
                                @warn "Skipping steep slice plot: empty or mismatched lengths: $(length(t_steep)) vs $(length(zs_steep))"
                            end
                            if length(t_flat) > 0 && length(zs_flat) == length(t_flat)
                                Plots.plot!(plt_slices, t_flat,  zs_flat,  label = "Flattest direction",  lw=3, color=:gray, legend=:outerbottom, legend_column=2, bottom_margin = -15mm)
                            else
                                @warn "Skipping flat slice plot: empty or mismatched lengths: $(length(t_flat)) vs $(length(zs_flat))"
                            end
                            Plots.vline!(plt_slices, [0.0], color=:blue, linestyle=:dash, label="Minimum")
                            Plots.xlabel!(plt_slices, L"\mathrm{Distance\;from\;minimum\;(px)}")
                            Plots.ylabel!(plt_slices, L"\mathrm{Cost}")
                            Plots.savefig(plt_slices, string(exp_path, "/Results/plots/cost_slices_along_directions.pdf"))

                            Plots.ylims!(plt_slices, 0, 10)
                            Plots.savefig(plt_slices, string(exp_path, "/Results/plots/cost_slices_along_directions_zoomed.pdf"))
                    

                            # save analysis metadata
                            dir_info = Dict("η_min"=>η0, "β_min"=>β0, "Hessian"=>H, "evals"=>evals, "v_steep"=>v_steep, "v_flat"=>v_flat)
                            slice_data = Dict(
                                "steep"=>Dict("t"=>t_steep, "etas"=>etas_steep, "betas"=>betas_steep, "zs"=>zs_steep),
                                "flat"=>Dict("t"=>t_flat, "etas"=>etas_flat, "betas"=>betas_flat, "zs"=>zs_flat)
                            )
                            write_json(string(exp_path, "/Results/data/direction_analysis"), dir_info)
                            write_json(string(exp_path, "/Results/data/slice_data"), slice_data)

                        catch err
                            @warn "Post-process Hessian/direction analysis failed: $err"
                        end

                        costList = stats["cost_list"]
                        iterList = stats["iterList"]

                        # Plot the estimated and ground truth height
                        set_plot(fs, sz=(1650, 1250))
                        est_h = vec(Float64.(collect(est_h)))
                        gt_h = vec(Float64.(collect(gt_h)))
                        n_time = min(length(time), length(est_h), length(gt_h))
                        if n_time < length(time) || n_time < length(est_h) || n_time < length(gt_h)
                            @warn "Time and height vectors have mismatched lengths: time=$(length(time)), est_h=$(length(est_h)), gt_h=$(length(gt_h)). Truncating to $n_time samples for plotting."
                        end
                        Plots.plot!(time[1:n_time], est_h[1:n_time], label="Estimated height", dpi=400, lw=3, legend=:outerbottom, legend_column=2, bottom_margin = -15mm)
                        Plots.plot!(time[1:n_time], gt_h[1:n_time], label="Ground truth height", dpi=400, lw=3, legend=:outerbottom, legend_column=2, bottom_margin = -15mm)
                        Plots.xlabel!(L"\mathrm{Time\;(s)}")
                        Plots.ylabel!(L"\mathrm{Height\;(mm)}")
                        Plots.savefig(string(exp_path,"/Results/plots/h_est.pdf"))

                        # Plot the height estimation error
                        set_plot(fs, sz=(1650, 1250))
                        # ensure vectors are aligned for error plot
                        n_err = min(length(time), length(est_h), length(gt_h))
                        if n_err < max(length(time), length(est_h), length(gt_h))
                            @warn "Truncating arrays for height error plot to $n_err samples"
                        end
                        Plots.plot!(time[1:n_err], abs.(est_h[1:n_err] .- gt_h[1:n_err]), label="Height estimation error", dpi=400, lw=3, legend=:outerbottom, legend_column=1, bottom_margin = -15mm)
                        Plots.xlabel!(L"\mathrm{Time\;(s)}")
                        Plots.ylabel!("Height Error (mm)")
                        Plots.savefig(string(exp_path,"/Results/plots/h_est_error.pdf"))

                        # Plot the estimated and ground truth parameters
                        set_plot(fs, sz=(1650, 1250))
                        Plots.plot!(est_η, label="Estimated η", marker=2, dpi=400, lw=3, legend=:outerbottom, legend_column=2, bottom_margin = -15mm)
                        Plots.hline!([η_gt], label="Ground truth η", lw=3, legend=:outerbottom, legend_column=2, bottom_margin = -15mm)
                        Plots.xlabel!(L"\mathrm{Iterations}")
                        Plots.ylabel!(L"\eta\;\mathrm{(KPa\cdot s)}")
                        Plots.savefig(string(exp_path,"/Results/plots/η.pdf"))
                        
                        set_plot(fs, sz=(1650, 1250))
                        Plots.plot!(est_β, label="Estimated β", marker=2, dpi=400, lw=3, legend=:outerbottom, legend_column=2, bottom_margin = -15mm)
                        Plots.hline!([β_gt], label="Ground truth β", lw=3,legend=:outerbottom, legend_column=2, bottom_margin = -15mm)
                        Plots.xlabel!(L"\mathrm{Iterations}")
                        Plots.ylabel!(L"\beta\;\mathrm{mm^{-1}}")
                        Plots.savefig(string(exp_path,"/Results/plots/β.pdf"))

                        # Plot the cost function with iterations
                        set_plot(fs, sz=(1650, 1250))
                        Plots.plot!(iterList, costList, label="Cost", marker=2, dpi=400, yscale=:log10, xminorgrid = :false, lw=3,legend=:outerbottom, legend_column=1, bottom_margin = -15mm)
                        Plots.xlabel!(L"\mathrm{Iterations}")
                        Plots.ylabel!(L"\mathrm{Cost\;(px)}")
                        # Plots.xticks!(minimum(iterList):2:maximum(iterList))
                        Plots.savefig(string(exp_path,"/Results/plots/cost_steps.pdf"))

                        # Plot the cost function with iterations
                        set_plot(fs, sz=(1650, 1250))
                        Plots.plot!(iterList, costList, label="Cost", marker=2, dpi=400, yscale=:log10, xscale=:log10, lw=3, legend=:outerbottom, legend_column=1, bottom_margin = -15mm)
                        Plots.xlabel!(L"\mathrm{Iterations}")
                        Plots.ylabel!(L"\mathrm{Cost\;(px)}")
                        Plots.savefig(string(exp_path,"/Results/plots/cost_steps_log.pdf"))


                        # Plot the cost function surface: prefer PlotlyJS HTML export to avoid
                        # creating a GLMakie Figure (which may require a compatible WebIO setup).
                        # set_plot(fs, sz=(1650, 1250))
                        # try

                        #     htmlpath = string(exp_path, "/Results/plots/cost_surface_interactive.html")
                        #     # compute z for ground-truth on the provided grid
                        #     Z_for_interp = CostMat'
                        #     z_gt = interp_z_at(η_gt, β_gt, ηList, βList, Z_for_interp)
                        #     saved = save_plotly_surface_html(htmlpath, ηList, βList, CostMat'; xs=Float64.(collect(est_η)), ys=Float64.(collect(est_β)), zs=Float64.(collect(cost_iter)), title="Cost surface", gt_x=η_gt, gt_y=β_gt, gt_z=z_gt, x_label = "\$\\eta\\;\\mathrm{(KPa\\cdot s)}\$", y_label = "\$\\beta\\;\\mathrm{(L/mm^{-1})}\$", font_size = fs, latex_labels=true)
                        #     if saved
                        #         @info "Saved interactive PlotlyJS HTML: $htmlpath"
                        #     else
                        #         @warn "PlotlyJS HTML not created; skipping interactive output. To enable interactive PNGs with GLMakie re-enable GLMakie manually."
                        #     end
                        # catch err
                        #     @warn "PlotlyJS path failed; skipping interactive output. Error: $err"
                        # end
                        # Also preserve the static PDFs via Plots
                        # try
                        #     plt = set_plot(fs, sz=(1650, 1250))
                        #     Plots.contour!(plt, ηList, βList, CostMat', color=:turbo, fill=false, levels=100, dpi=400, legend=:outerbottom, legend_column=3, bottom_margin = -15mm)
                        #     Plots.plot!(plt, est_η, est_β, label="Estimations", ms=:4, m=:x, color=:red, lw=3, legend=:outerbottom, legend_column=3, bottom_margin = -15mm)
                        #     Plots.scatter!(plt, [η_gt], [β_gt], label="Ground truth", ms=:8, m=:circle, color=:indianred2, lw=3, legend=:outerbottom, legend_column=3, bottom_margin = -15mm)
                        #     Plots.scatter!(plt_dirs, [η0], [β0], label="Minimum Cost", ms=15, m=:star5, color=:black, legend=:outerbottom, legend_column=3, bottom_margin = -15mm)
                        #     Plots.xlabel!(plt, L"\eta\;\mathrm{(KPa\cdot s)}")
                        #     Plots.ylabel!(plt, L"\beta\;\mathrm{mm^{-1}}")
                        #     Plots.savefig(plt, string(exp_path,"/Results/plots/cost_surface_iter.pdf"))

                        #     plt2 = set_plot(fs, sz=(1650, 1250))
                        #     Plots.contourf!(plt2, ηList, βList, CostMat', color=:turbo, fill=false, levels=100, dpi=400, legend=:outerbottom, legend_column=3, bottom_margin = -15mm)
                        #     Plots.plot!(plt2, est_η, est_β, label="Estimations", ms=:4, m=:x, color=:red, lw=3, legend=:outerbottom, legend_column=3, bottom_margin = -15mm)
                        #     Plots.scatter!(plt2, [η_gt], [β_gt], label="Ground truth", ms=:8, m=:circle, color=:indianred2, lw=3, legend=:outerbottom, legend_column=3, bottom_margin = -15mm)
                        #     Plots.xlabel!(plt2, L"\eta\;\mathrm{(KPa\cdot s)}")
                        #     Plots.ylabel!(plt2, L"\beta\;\mathrm{mm^{-1}}")
                        #     Plots.savefig(plt2, string(exp_path,"/Results/plots/cost_surface.pdf"))
                        # catch err
                        #     @warn "Failed to produce static PDF contour outputs: $err"
                        # end
                    end
                
                elseif viscosity_type == "bulk_viscosity"
                    window_dirs = readdir(string(exp_path))
                    for window_dir in window_dirs
                        if window_dir == "Results" || window_dir == "post_analysis_window" || window_dir == "single_window"
                            continue
                        end
                        win_exp_path = string(filepath, elem_size_folder, "/", sim_time_folder, "/", noise_folder, "/", window_dir)
                        
                        println("Processing window: $win_exp_path")
                        exp_params = read_json(string(win_exp_path ,"/Results/data/experiment_parameters.json"))
                        sim_time_exp = exp_params["sim_time_exp"]
                        data_type = exp_params["data_type"]

                        sim_time_gt = 40.0 # simulation time in seconds
                        steps_gt = 400.0 # number of time steps
                        t_steps_gt = sim_time_gt/steps_gt

                        est_ηpList = readdlm(string(win_exp_path,"/Results/data/est_η.csv"), ',', Float64)
                        est_βpList = readdlm(string(win_exp_path,"/Results/data/est_β.csv"), ',', Float64)
                        avg_ηList = readdlm(string(win_exp_path,"/Results/data/avg_η.csv"), ',', Float64)
                        η_gt = readdlm(string(win_exp_path,"/Results/data/η_gt.csv"), ',', Float64)
                        β_gt = readdlm(string(win_exp_path,"/Results/data/β_gt.csv"), ',', Float64)
                        est_h_list = readdlm(string(win_exp_path,"/Results/data/est_h.csv"), ',', Float64)
                        gt_h_list = readdlm(string(win_exp_path,"/Results/data/gt_h.csv"), ',', Float64)
                        
                        data_ranges_ = readdlm(string(win_exp_path,"/Results/data/window_data/data_ranges.csv"), ',', Int)
                        t_windows = readdlm(string(win_exp_path,"/Results/data/window_data/t_windows.csv"), ',', Float64)
                        
                        symBorderPts, splinex, spliney = read_csv(string(win_exp_path,"/Results/data/sim_data/2D_border_points"))

                        d, pairs = closest_point(symBorderPts, obsBorderPts)

                        plt_cnt_error = set_plot(fs, sz=(1650, 1250))
                        Plots.plot!(plt_cnt_error, d, label="Closest point distance error", dpi=400, lw=3)
                        Plots.xlabel!(L"\mathrm{Time\;(s)}")
                        Plots.ylabel!(L"\mathrm{Closest\;Point\;Distance\;(px)}")
                        Plots.savefig(plt_cnt_error, string(win_exp_path,"/Results/plots/closest_point_distance_error.pdf"))

                        plt_η = set_plot(fs, sz=(1650, 1250))
                        t_full = collect(range(start=t_steps_gt, stop=sim_time_gt, step=t_steps_gt))
                        if data_type != "physical"
                            Plots.plot!(plt_η, t_full, η_gt, label="Ground truth η(t)", dpi=400, lw=3)
                        end
                        for ti::Int in 1:size(data_ranges_, 1)
                            t = t_windows[ti]
                            Plots.vline!(plt_η, [t], color=:gray, lw=3, linestyle=:dash, label=false)
                        end
                        Plots.xlabel!(L"\mathrm{Time\;(s)}")
                        Plots.ylabel!(L"\eta\mathrm{(t)\;(KPa\cdot s)}")
                        Plots.savefig(plt_η, string(win_exp_path,"/Results/plots/η_gt.pdf"))
                        t_prev = 0.1
                        for ti::Int in 1:size(data_ranges_, 1)
                            t = t_windows[ti]
                            data_range_ = data_ranges_[ti,:]
                            t_win = collect(range(start=t_prev, stop=t, step=t_steps_gt))
                            if ti == 1
                                Plots.plot!(plt_η, t_win, est_ηpList[data_range_], label="Estimated η(t)", lw=3, color=:orange)
                                if data_type != "physical"
                                    Plots.plot!(plt_η, t_win, avg_ηList[data_range_], label="Average GT η in window", lw=3, color=:gray)
                                end 
                            else
                                Plots.plot!(plt_η, t_win, est_ηpList[data_range_], lw=3, color=:orange, label=false)
                                if data_type != "physical"
                                    Plots.plot!(plt_η, t_win, avg_ηList[data_range_], lw=3, color=:gray, label=false)
                                end 
                            end
                            t_prev = t+t_steps_gt
                        end
                        Plots.savefig(plt_η, string(win_exp_path,"/Results/plots/η.pdf"))
                        
                        plt_β = set_plot(fs, sz=(1650, 1250))
                        println(size(β_gt))
                        if data_type != "physical"
                            Plots.hline!(plt_β, β_gt, label="Ground truth β(t)", dpi=400, lw=3)
                        end
                        t_prev = 0.1
                        for ti::Int in 1:size(data_ranges_, 1)
                            t = t_windows[ti]
                            Plots.vline!(plt_β, [t], color=:gray, lw=3, linestyle=:dash, label=false)
                            data_range_ = data_ranges_[ti,:]
                            t_win = collect(range(start=t_prev, stop=t, step=t_steps_gt))
                            if ti == 1
                                Plots.plot!(plt_β, t_win, est_βpList[data_range_], label="Estimated β(t)", lw=3, color=:orange)
                            else
                                Plots.plot!(plt_β, t_win, est_βpList[data_range_], lw=3, color=:orange, label=false)
                            end
                            t_prev = t+t_steps_gt
                        end
                        Plots.xlabel!(L"\mathrm{Time\;(s)}")
                        Plots.ylabel!(L"\beta\mathrm{(t)\;(mm^{-1})}")
                        Plots.savefig(plt_β, string(win_exp_path,"/Results/plots/β.pdf"))
                        
                        # animate_fields(filepath=string(win_exp_path,"/Results/plots"), p=splinex, q=spliney, pObs=nSplinex, qObs=nSpliney)

                        if data_type != "physical"

                            plt_h = set_plot(fs, sz=(1650, 1250))
                            Plots.plot!(gt_h_list, label="Ground truth height", dpi=400, lw=3)
                            Plots.plot!(est_h_list, label="Estimated height")
                            Plots.xlabel!(L"\mathrm{Time\;(s)}")
                            Plots.ylabel!(L"\mathrm{Height\;(mm)}")
                            Plots.savefig(string(win_exp_path,"/Results/plots/h.pdf"))

                            plt_error = set_plot(fs, sz=(1650, 1250))
                            Plots.plot!(abs.(est_h_list-gt_h_list), label="Height estimation error", dpi=400, lw=3)
                            Plots.savefig(string(win_exp_path,"/Results/plots/h_est_error.pdf"))
                            Plots.xlabel!(L"\mathrm{Time\;(s)}")
                            Plots.ylabel!(L"\mathrm{Height\;Error\;(px)}")
                            Plots.savefig(string(win_exp_path,"/Results/plots/error.pdf"))
                        end
                    end

                end
            end
        end
    end

    return 
end

function post_analysis(filepath_gt_::String, filepath::String, avoid_list)
    # ηList = readdlm(string(filepath_gt,"/data/η.csv"), ',', Float64)
    # βList = readdlm(string(filepath_gt,"/data/β.csv"), ',', Float64)

    dir_list = readdir(filepath)

    # plot for convergence per slip case
    plot_conv_10 = set_plot(fs, sz=(1650, 1250))

    plot_conv_20 = set_plot(fs, sz=(1650, 1250))

    plot_conv_30 = set_plot(fs, sz=(1650, 1250))

    plot_conv_log_10 = set_plot(fs, sz=(1650, 1250))

    plot_conv_log_20 = set_plot(fs, sz=(1650, 1250))

    plot_conv_log_30 = set_plot(fs, sz=(1650, 1250))

    for dir in dir_list

        if dir in avoid_list || dir == "post_analysis_global"
            println("Skipping directory: ", dir)
            continue
        end
        filepath_dir = string(filepath,dir)
        filepath_gt = string(filepath_gt_,"/",dir)

        printstyled("Processing directory: $(filepath_dir)\n", color=:green)
        params = read_json(string(filepath_gt,"/data/sim_params.json"))
        η_gt = params["η"]
        β_gt = params["β"]  

        println("Ground truth η: ", η_gt[1])
        η_gt_list = η_gt[1]*ones(40,1)
        β_gt_list = β_gt[1]*ones(40,1)

        elem_size_folders = readdir(filepath_dir)

        # plots for element vise comparison
        fig1 = set_plot(fs, sz=(1650, 1250))
        # Plots.plot!(fig1, η_gt_list, label="Ground truth η", dpi=400, lw=3)
        Plots.hline!(fig1, [η_gt[1]], label="Ground truth η", lw=3)
        Plots.xlabel!(fig1,L"\mathrm{Iterations}")
        Plots.ylabel!(fig1,L"\eta\;\mathrm{(KPa\cdot s)}")

        fig2 = set_plot(fs, sz=(1650, 1250))
        # Plots.plot!(fig2, β_gt_list, label="Ground truth β", dpi=400, lw=3)
        Plots.hline!(fig2, [β_gt[1]], label="Ground truth β", lw=3)
        Plots.xlabel!(fig2, L"\mathrm{Iterations}")
        Plots.ylabel!(fig2, L"\beta\;\mathrm{mm^{-1}}")
        
        for elem_size_folder_ in elem_size_folders
            if elem_size_folder_ == "post_analysis"
                continue
            end

            elem_size_folder = string(filepath_dir, "/", elem_size_folder_)
            printstyled("Processing element size folder: $(elem_size_folder)\n", color=:blue)  
            sim_time_folders = readdir(elem_size_folder)

            # figures for Simulation window vise camparison
            fig3 = set_plot(fs, sz=(1650, 1250))
            Plots.plot!(fig3, η_gt_list, label="Ground truth η", dpi=400, lw=3)
            Plots.xlabel!(fig3,L"\mathrm{Iterations}")
            Plots.ylabel!(fig3,L"\eta\;\mathrm{(KPa\cdot s)}")

            fig4 = set_plot(fs, sz=(1650, 1250))
            Plots.plot!(fig4, β_gt_list, label="Ground truth β", dpi=400, lw=3)
            Plots.xlabel!(fig4, L"\mathrm{Iterations}")
            Plots.ylabel!(fig4, L"\beta\;\mathrm{mm^{-1}}")

            plt_slices = set_plot(fs, sz=(1800, 750))
            Plots.vline!(plt_slices, [0.0], color=:blue, linestyle=:dash, label="Minimum", bottom_margin = 15mm, left_margin=12mm, lw=3)
            Plots.xlabel!(plt_slices, L"\mathrm{Distance\;from\;minimum\;(px)}")
            Plots.ylabel!(plt_slices, L"\mathrm{Cost}")
            Plots.ylims!(plt_slices, 0, 50)

            for sim_time_folder_ in sim_time_folders

                if sim_time_folder_ == "post_analysis_time" || sim_time_folder_ == "Results" 
                    continue
                end

                noise_folders = readdir(string(elem_size_folder,"/",sim_time_folder_))
                printstyled("Processing simulation time folder: $(string(elem_size_folder,"/",sim_time_folder_))\n", color=:cyan)
                for noise_folder_ in noise_folders

                    if noise_folder_ == "post_analysis_noise" || noise_folder_ == "Results" 
                        continue
                    end

                    exp_folder = string(filepath_dir,"/",elem_size_folder_,"/",sim_time_folder_,"/",noise_folder_)

                    printstyled("Processing experiment folder: $(exp_folder)\n", color=:magenta)
                    params = read_json(string(exp_folder,"/Results/data/experiment_parameters.json"))

                    noise_level = params["noise_level"]

                    FunctionClass_x = params["FunctionClass_x"]
                    ne = params["ne_exp"]
                    sim_time = params["sim_time_exp"]
            
                    
                    if noise_level == 0.0
                        printstyled("Processing for zeros noise level: $(exp_folder)\n", color=:yellow)
                        η = readdlm(string(exp_folder,"/Results/data/η.csv"), ',', Float64)
                        β = readdlm(string(exp_folder,"/Results/data/β.csv"), ',', Float64)
                        
                        cost_list = readdlm(string(exp_folder,"/Results/data/cost_iter.csv"), ',', Float64)
                
                        Plots.plot!(fig3, η, label=string("Window size - $(sim_time)s"," - ne: ",ne), marker=2, dpi=400, lw=3, legend=:outerbottom, legend_column=2, bottom_margin = -15mm)
                        Plots.xlabel!(fig3, L"\mathrm{Iterations}")
                        Plots.ylabel!(fig3, L"\eta\;\mathrm{(KPa\cdot s)}")
                
                        Plots.plot!(fig4, β, label=string("Window size - $(sim_time)s"," - ne: ",ne), marker=2, dpi=400, lw=3, legend=:outerbottom, legendcolumn=2, bottom_margin = -15mm)
                        Plots.xlabel!(fig4, L"\mathrm{Iterations}")
                        Plots.ylabel!(fig4, L"\beta\;\mathrm{mm^{-1}}")

                        slice_data = read_json(string(exp_folder,"/Results/data/slice_data.json"))

                        t_steep = Float64.(collect(slice_data["steep"]["t"]))
                        zs_steep = Float64.(collect(slice_data["steep"]["zs"]))
                        t_flat = Float64.(collect(slice_data["flat"]["t"]))
                        zs_flat = Float64.(collect(slice_data["flat"]["zs"]))

                        # 2D slices: cost vs distance along the two directions
                        
                        if length(t_steep) > 0 && length(zs_steep) == length(t_steep)
                            Plots.plot!(plt_slices, t_steep, zs_steep, lw=3, legend = :outerright, label = "Steepest direction; window size = $(sim_time)s")
                        else
                            @warn "Skipping steep slice plot: empty or mismatched lengths: $(length(t_steep)) vs $(length(zs_steep))"
                        end

                        if length(t_flat) > 0 && length(zs_flat) == length(t_flat)
                            Plots.plot!(plt_slices, t_flat,  zs_flat, lw=3, legend = :outerright, label = "Flattest direction; window size = $(sim_time)s", legendfontsize=20)
                        else
                            @warn "Skipping flat slice plot: empty or mismatched lengths: $(length(t_flat)) vs $(length(zs_flat))"
                        end

                        if sim_time == 10.0
                            if ne == 6
                                Plots.plot!(plot_conv_10, cost_list, label=string(L"\beta:\;",β_gt), marker=2, dpi=400, lw=3, legend=:outerbottom, legendcolumn=5, bottom_margin = -15mm)
                                Plots.xlabel!(plot_conv_10, L"\mathrm{Iterations}")
                                Plots.ylabel!(plot_conv_10, L"mathrm{Cost\;(px)}")

                                Plots.plot!(plot_conv_log_10, cost_list, label=string(L"\beta:\;",β_gt), marker=2, dpi=400, lw=3, legend=:outerbottom, legendcolumn=5, bottom_margin = -15mm, xscale=:log10, yscale=:log10)
                                Plots.xlabel!(plot_conv_log_10, L"\mathrm{Iterations}")
                                Plots.ylabel!(plot_conv_log_10, L"mathrm{Cost\;(px)}")
                            end

                        elseif sim_time == 20.0
                            if ne == 6
                                Plots.plot!(plot_conv_20, cost_list, label=string(L"\beta:\;",β_gt), marker=2, dpi=400, lw=3, legend=:outerbottom, legendcolumn=5, bottom_margin = -15mm)
                                Plots.xlabel!(plot_conv_20, L"\mathrm{Iterations}")
                                Plots.ylabel!(plot_conv_20, L"mathrm{Cost\;(px)}")

                                Plots.plot!(plot_conv_log_20, cost_list, label=string(L"\beta:\;",β_gt), marker=2, dpi=400, lw=3, legend=:outerbottom, legendcolumn=5, bottom_margin = -15mm, xscale=:log10, yscale=:log10)
                                Plots.xlabel!(plot_conv_log_20, L"\mathrm{Iterations}")
                                Plots.ylabel!(plot_conv_log_20, L"mathrm{Cost\;(px)}")
                            end

                        elseif sim_time == 30.0
                            if ne == 6
                                Plots.plot!(plot_conv_30, cost_list, label=string(L"\beta:\;",β_gt), marker=2, dpi=400, lw=3, legend=:outerbottom, legendcolumn=5, bottom_margin = -15mm)
                                Plots.xlabel!(plot_conv_30, L"\mathrm{Iterations}")
                                Plots.ylabel!(plot_conv_30, L"mathrm{Cost\;(px)}")

                                Plots.plot!(plot_conv_log_30, cost_list, label=string(L"\beta:\;",β_gt), marker=2, dpi=400, lw=3, legend=:outerbottom, legendcolumn=5, bottom_margin = -15mm, xscale=:log10, yscale=:log10)
                                Plots.xlabel!(plot_conv_log_30, L"\mathrm{Iterations}")
                                Plots.ylabel!(plot_conv_log_30, L"mathrm{Cost\;(px)}")
                            end

                            # Plots.plot!(fig1, η, label=string("Basis - ",FunctionClass_x," - ne: ",ne), marker=2, dpi=400, lw=3, legend=:outerbottom, legendcolumn=2, bottom_margin = -15mm)
                            Plots.plot!(fig1, η, label=string("Number of elements: ",ne), marker=2, dpi=400, lw=3, legend=:outerbottom, legendcolumn=2, bottom_margin = -15mm)
                            Plots.xlabel!(fig1, L"\mathrm{Iterations}")
                            Plots.ylabel!(fig1, L"\eta\;\mathrm{(KPa\cdot s)}")

                            # Plots.plot!(fig2, β, label=string("Basis - ",FunctionClass_x," - ne: ",ne), marker=2, dpi=400, lw=3, legend=:outerbottom, legendcolumn=2, bottom_margin = -15mm)
                            Plots.plot!(fig2, β, label=string("Number of elements: ",ne), marker=2, dpi=400, lw=3, legend=:outerbottom, legendcolumn=2, bottom_margin = -15mm)
                            Plots.xlabel!(fig2, L"\mathrm{Iterations}")
                            Plots.ylabel!(fig2, L"\beta\;\mathrm{mm^{-1}}")

                        end
                    end
                end
            end
            
            plot_path_sim_time = string(elem_size_folder,"/post_analysis_time/plots")
            set_file(plot_path_sim_time)
            
            Plots.savefig(fig3, string(plot_path_sim_time,"/η.pdf"))
            Plots.savefig(fig4, string(plot_path_sim_time,"/β.pdf"))
            Plots.savefig(plt_slices, string(plot_path_sim_time,"/cost_slices_along_directions.pdf"))
        end

        plot_path_elems = string(filepath_dir,"/post_analysis/plots")
        set_file(plot_path_elems)

        @info "Saving plots to $plot_path_elems"
        Plots.savefig(fig1, string(plot_path_elems,"/η_20.pdf"))
        Plots.savefig(fig2, string(plot_path_elems,"/β_20.pdf"))
    end
    plot_path_global = string(filepath,"/post_analysis_global/plots")
    set_file(plot_path_global)

    @info "Saving plots to $plot_path_global"
    Plots.savefig(plot_conv_10, string(plot_path_global,"/conv_10.pdf"))
    Plots.savefig(plot_conv_log_10, string(plot_path_global,"/conv_log_10.pdf"))

    Plots.savefig(plot_conv_20, string(plot_path_global,"/conv_20.pdf"))
    Plots.savefig(plot_conv_log_20, string(plot_path_global,"/conv_log_20.pdf"))

    Plots.savefig(plot_conv_30, string(plot_path_global,"/conv_30.pdf"))
    Plots.savefig(plot_conv_log_30, string(plot_path_global,"/conv_log_30.pdf"))

end

# Note: GLMakie / Makie support was removed from automated paths; this script
# now prefers PlotlyJS for saved interactive HTML output and keeps Plots for
# static PDF output. If you need Makie-based interactive windows, re-enable
# GLMakie imports and the helper function manually.

# Helper: ensure xg, yg are Float64 vectors and Z is a dense Float64 matrix
# shaped (ny, nx) where ny = length(yg), nx = length(xg). Returns (xg_arr, yg_arr, Zmat)
function _ensure_grid_matrix(xg, yg, Z)
    xg_arr = Float64.(collect(xg))
    yg_arr = Float64.(collect(yg))
    nx = length(xg_arr); ny = length(yg_arr)
    Zmat = try
        Array(Z) |> x -> Float64.(x)
    catch err
        error("_ensure_grid_matrix: Failed to coerce Z to dense Float64 matrix: $err")
    end
    if ndims(Zmat) == 1
        if length(Zmat) == nx*ny
            Zmat = reshape(Zmat, ny, nx)
        else
            error("_ensure_grid_matrix: Z is 1D with length=$(length(Zmat)) which is not nx*ny=$(nx*ny)")
        end
    elseif ndims(Zmat) == 2
        if size(Zmat,1) == ny && size(Zmat,2) == nx
            # ok
        elseif size(Zmat,1) == nx && size(Zmat,2) == ny
            Zmat = Zmat'
        elseif prod(size(Zmat)) == nx*ny
            Zmat = reshape(vec(Zmat), ny, nx)
        else
            @warn "_ensure_grid_matrix: Z matrix shape $(size(Zmat)) does not match (ny,nx)=($(ny),$(nx)). Attempting to continue but result may be wrong."
        end
    else
        error("_ensure_grid_matrix: Z has unexpected dims: $(ndims(Zmat))")
    end
    return xg_arr, yg_arr, Zmat
end

# PlotlyJS exporter: create a Plotly surface and optional scatter overlays,
# then save as a self-contained HTML file.
function save_plotly_surface_html(filename::AbstractString, xg, yg, Z; xs=AbstractVector[], ys=AbstractVector[], zs=AbstractVector[],
                                  title="", colormap="Viridis",
                                  surface_label::AbstractString = "Cost surface",
                                  path_label::AbstractString = "Estimations",
                                  gt_label::AbstractString = "Ground truth",
                                  path_color::AbstractString = "red",
                                  gt_color::AbstractString = "indianred2",
                                  path_marker_size::Int = 14,
                                  gt_marker_size::Int = 18,
                                  gt_x = nothing, gt_y = nothing, gt_z = nothing,
                                  x_label::AbstractString = "",
                                  y_label::AbstractString = "",
                                  font_size::Int = 14,
                                  latex_labels::Bool = false,
                                  z_offset::Real = 0.0,
                                  surface_opacity::Real = 1.0)
    try
        @eval using PlotlyJS
    catch err
        @warn "PlotlyJS not available; cannot save interactive HTML: $err"
        return false
    end
    try
        xg_arr, yg_arr, Zmat = _ensure_grid_matrix(xg, yg, Z)
    # Plotly expects z as 2D array with rows corresponding to y
    # Surface trace: do NOT include the surface in the legend so the HTML
    # matches the PDF (PDF legend only shows Estimations and Ground truth).
    # Set transparency so overlay lines/markers remain visible. Use a lower
    # opacity to reduce occlusion of the estimation path.
    surf_trace = PlotlyJS.surface(x = xg_arr, y = yg_arr, z = Zmat, colorscale = colormap,
                      name = surface_label, showlegend = false, opacity = surface_opacity, showscale = true)
    # compute a tiny z offset (if z_offset==0.0 we pick a conservative
    # default proportional to the Z range) so overlays sit slightly above
    # the surface and avoid z-fighting/occlusion in the browser.
    zrange = maximum(Zmat) - minimum(Zmat)
    # use a slightly larger default offset so overlay traces are visible
    # in WebGL renderers; user can override via z_offset kwarg
    eff_offset = Float64(z_offset == 0.0 ? 1e-3 * (zrange + eps(Float64)) : z_offset)

    # Build traces list and create Plot without relying on add_traces! helper
        traces = [surf_trace]
        # add overlays if provided
        if length(xs) > 0 && length(ys) > 0
            xsv = Float64.(collect(xs))
            ysv = Float64.(collect(ys))
            # If zs not provided or wrong length, try to compute via interpolation
            if length(zs) == length(xs)
                zsv = Float64.(collect(zs))
            else
                try
                    zsv = Float64.(interp_z_at(xsv, ysv, xg_arr, yg_arr, Zmat))
                catch err
                    @warn "Estimator zs length mismatch and interpolation failed: $err; skipping estimation trace"
                    zsv = Float64[]
                end
            end
            if length(zsv) == length(xsv) && length(zsv) > 0
                @info "Adding estimation trace: $(length(xsv)) points"
                # Plot estimation path as line+markers so the trajectory is visible.
                # Use a star marker and strong black outline so it stands out above the surface.
                # lift the points slightly above the surface to avoid z-fighting
                zsv = zsv .+ eff_offset
                scatter_trace = PlotlyJS.scatter3d(x = xsv, y = ysv, z = zsv, mode = "lines+markers",
                                                   name = path_label,
                                                   marker = PlotlyJS.attr(size = path_marker_size, color = path_color, opacity = 1.0, symbol = "star", line = Dict(:width => 1.5, :color => "black")),
                                                   line = PlotlyJS.attr(width = 6, color = path_color))
                push!(traces, scatter_trace)
                # Add a non-legend projected overlay slightly above the surface
                # as a visual aid so the path is always visible in the HTML viewer.
                try
                    proj_offset = 0.02 * (zrange + eps(Float64))
                    zsv_proj = zsv .+ proj_offset
                    proj_marker = PlotlyJS.attr(size = max(path_marker_size+6, 22), color = path_color, opacity = 1.0, symbol = "circle", line = Dict(:width => 2, :color => "black"))
                    proj_trace = PlotlyJS.scatter3d(x = xsv, y = ysv, z = zsv_proj, mode = "markers",
                                                    name = "", showlegend = false,
                                                    marker = proj_marker)
                    push!(traces, proj_trace)
                catch err
                    @warn "Failed to add projected overlay trace: $err"
                end
            else
                @warn "Estimation trace not added: zs length $(length(zsv)) does not match xs length $(length(xsv))"
            end
        end
        # Optional ground-truth point (single marker)
        if gt_x !== nothing && gt_y !== nothing && gt_z !== nothing
            gx = _ensure_scalar_float(gt_x)
            gy = _ensure_scalar_float(gt_y)
            gz = _ensure_scalar_float(gt_z)
            # lift ground-truth marker slightly as well so it is visible above the surface
            gz = gz + eff_offset
            gt_trace = PlotlyJS.scatter3d(x = [gx], y = [gy], z = [gz], mode = "markers",
                                          name = gt_label,
                                          marker = PlotlyJS.attr(size = gt_marker_size, color = gt_color, symbol = "star"))
            push!(traces, gt_trace)
        end
    # Construct the Plot from the traces vector (some PlotlyJS versions
    # expect a single array argument rather than varargs).
        # Utility: convert a small subset of LaTeX expressions to a plain
        # unicode string for axis titles so they render correctly in Plotly's
        # SVG text when MathJax cannot typeset SVG contents.
        function _tex_to_plain(s::AbstractString)
            mp = Dict("\\eta"=>"η", "\\beta"=>"β", "\\cdot"=>"·")
            out = s
            # unwrap \mathrm{...} -> inner text
            out = replace(out, r"\\mathrm\{([^}]*)\}" => s"\1")
            # replace common macros
            for (k,v) in mp
                out = replace(out, k => v)
            end
            # remove dollar signs and braces
            out = replace(out, "\$" => "")
            out = replace(out, '{' => "")
            out = replace(out, '}' => "")
            return out
        end

        # Layout: place legend inside the plotting area (top-left) so it does
        # not overlap the surface colorbar on the right. Provide a translucent
        # white background for readability and a thin border.
        legend_cfg = Dict(:x => 0.02,
                          :y => 0.98,
                          :xanchor => "left",
                          :yanchor => "top",
                          :bgcolor => "rgba(255,255,255,0.85)",
                          :bordercolor => "black",
                          :borderwidth => 1,
                          :font => Dict(:size => Int(round(font_size*0.6))))
    # reasonable margins to accommodate colorbar and axes labels
    margin_cfg = Dict(:l => 60, :r => 120, :t => 80, :b => 60)

    # Axis label configs to align with PDF labels (2D axes). Use a smaller
    # tick font so numbers don't appear oversized in the HTML.
    tick_font_size = max(8, Int(round(font_size*0.36)))
    xaxis_cfg = Dict(:title => x_label, :titlefont => Dict(:size => font_size), :tickfont => Dict(:size => tick_font_size))
    yaxis_cfg = Dict(:title => y_label, :titlefont => Dict(:size => font_size), :tickfont => Dict(:size => tick_font_size))
    # Scene (3D) axis labels — ensure 3D plot axes show the same titles.
    # Use plain (unicode) versions for SVG rendering; MathJax may still
    # typeset the LaTeX string if available.
    plain_x = _tex_to_plain(x_label)
    plain_y = _tex_to_plain(y_label)
    scene_cfg = Dict(:xaxis => Dict(:title => plain_x, :titlefont => Dict(:size => font_size), :tickfont => Dict(:size => tick_font_size)),
             :yaxis => Dict(:title => plain_y, :titlefont => Dict(:size => font_size), :tickfont => Dict(:size => tick_font_size)),
             :zaxis => Dict(:title => "Cost", :titlefont => Dict(:size => font_size), :tickfont => Dict(:size => tick_font_size)))

        # Try to construct the Plot with a Layout first; if the PlotlyJS
        # version here doesn't accept that signature, fall back to creating
        # the Plot and calling relayout!.
        # Provide a default camera view that looks slightly down at the surface
        # so overlaid paths are visible by default.
        camera_cfg = Dict(:eye => Dict(:x => 1.2, :y => 1.2, :z => 0.9))
        scene_with_camera = merge(scene_cfg, Dict(:camera => camera_cfg))

        layout_obj = try
            PlotlyJS.Layout(title = title, legend = legend_cfg, margin = margin_cfg, autosize = true,
                            xaxis = xaxis_cfg, yaxis = yaxis_cfg, scene = scene_with_camera)
        catch _
          # Older versions may expect a Dict for layout
            Dict(:title => title, :legend => legend_cfg, :margin => margin_cfg, :autosize => true,
                 :xaxis => xaxis_cfg, :yaxis => yaxis_cfg, :scene => scene_with_camera)
        end

        plt = try
            PlotlyJS.Plot(traces, layout_obj)
        catch _
            plt = PlotlyJS.Plot(traces)
            try
                PlotlyJS.relayout!(plt, layout_obj)
            catch _
                # If relayout also fails, at least set the title if available
                if !isempty(title)
                    try
                        PlotlyJS.relayout!(plt, Dict(:title => title))
                    catch _
                    end
                end
            end
        end
        # Render the plot to a HTML string using the MIME show method which is
        # available across PlotlyJS/PlotlyBase versions, then write the string.
        html_str = sprint(io -> show(io, MIME("text/html"), plt))
        # If LaTeX labels requested, inject MathJax into the HTML so TeX is rendered
        if latex_labels
            mj_cfg = raw"""<script>window.MathJax = {tex: {inlineMath: [['$','$'], ['\(','\)']]}};</script>""" * "\n"
            mj_src = "<script src=\"https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js\"></script>\n"
            mj_run = "<script>if(window.MathJax && MathJax.typesetPromise) MathJax.typesetPromise();</script>\n"
            inject = mj_cfg * mj_src * mj_run
            # Try to insert into <head> if present, else before <body>, else prepend
            if occursin("</head>", html_str)
                html_str = replace(html_str, "</head>" => inject * "</head>")
            elseif occursin("<body>", html_str)
                html_str = replace(html_str, "<body>" => "<head>" * inject * "</head><body>")
            else
                html_str = inject * html_str
            end
        end
        open(filename, "w") do io
            write(io, html_str)
        end
        return true
    catch err
        @warn "Failed to create/save PlotlyJS HTML: $err"
        return false
    end
end

# GLMakie helper removed: this project now uses PlotlyJS for saved
# interactive HTML output. If you previously relied on the Makie helper,
# use `save_plotly_surface_html` above to write interactive HTML files.

## Bilinear interpolation helper: given vectors xg (length nx), yg (length ny)
## and a Z matrix of size (ny, nx) (rows->y, cols->x), return interpolated
## z value at (x, y). Values outside the grid are clamped to the edges.
function interp_z_at(x::Real, y::Real, xg::AbstractVector, yg::AbstractVector, Z::AbstractMatrix)
    # Coerce grid vectors and Z to concrete Float64 arrays/matrix and
    # normalize shapes. This avoids errors when Z is an Adjoint, 1D Vec,
    # or contains Union element types (JSON-backed arrays).
    xg_arr = Float64.(collect(xg))
    yg_arr = Float64.(collect(yg))
    nx = length(xg_arr); ny = length(yg_arr)

    Zmat = try
        Array(Z) |> x -> Float64.(x)
    catch err
        error("Failed to coerce Z for interpolation: $err")
    end

    # If Zmat ended up 1D (or different shape), try to reshape/transpose to
    # match (ny, nx) where rows -> y and cols -> x.
    if ndims(Zmat) == 1
        if length(Zmat) == nx*ny
            Zmat = reshape(Zmat, ny, nx)
        else
            error("interp_z_at: Z is 1D with length=$(length(Zmat)) which is not nx*ny=$(nx*ny)")
        end
    elseif ndims(Zmat) == 2
        if size(Zmat,1) == ny && size(Zmat,2) == nx
            # ok
        elseif size(Zmat,1) == nx && size(Zmat,2) == ny
            Zmat = Zmat'
        elseif prod(size(Zmat)) == nx*ny
            Zmat = reshape(vec(Zmat), ny, nx)
        else
            @warn "interp_z_at: Z matrix shape $(size(Zmat)) does not match (ny,nx)=($(ny),$(nx)). Attempting to continue but interpolation may be wrong."
        end
    else
        error("interp_z_at: Z has unexpected number of dimensions: $(ndims(Zmat))")
    end

    # clamp x,y to grid extents
    xcl = clamp(x, first(xg_arr), last(xg_arr))
    ycl = clamp(y, first(yg_arr), last(yg_arr))

    # find surrounding indices
    ix = searchsortedfirst(xg_arr, xcl)
    if ix == 1
        i0 = 1; i1 = 1; tx = 0.0
    elseif ix > nx
        i0 = nx; i1 = nx; tx = 0.0
    else
        i1 = ix; i0 = max(1, ix-1)
        x0 = xg_arr[i0]; x1 = xg_arr[i1]
        tx = x1==x0 ? 0.0 : (xcl - x0)/(x1 - x0)
    end

    iy = searchsortedfirst(yg, ycl)
    if iy == 1
        j0 = 1; j1 = 1; ty = 0.0
    elseif iy > ny
        j0 = ny; j1 = ny; ty = 0.0
    else
        j1 = iy; j0 = max(1, iy-1)
        y0 = yg_arr[j0]; y1 = yg_arr[j1]
        ty = y1==y0 ? 0.0 : (ycl - y0)/(y1 - y0)
    end

    # values at corners: note Z rows->y (j index), cols->x (i index)
    z00 = float(Zmat[j0, i0]); z10 = float(Zmat[j0, i1])
    z01 = float(Zmat[j1, i0]); z11 = float(Zmat[j1, i1])

    # bilinear interpolation
    z0 = (1-tx)*z00 + tx*z10
    z1 = (1-tx)*z01 + tx*z11
    z = (1-ty)*z0 + ty*z1
    return z
end

# Overload to accept vector inputs (single-element vectors or paired vectors)
function interp_z_at(xs::AbstractVector{<:Real}, ys::AbstractVector{<:Real}, xg::AbstractVector, yg::AbstractVector, Z::AbstractMatrix)
    # If single-element vectors, delegate to scalar version
    if length(xs) == 1 && length(ys) == 1
        return interp_z_at(xs[1], ys[1], xg, yg, Z)
    end
    # If paired vectors, compute element-wise
    if length(xs) == length(ys)
        return [interp_z_at(xi, yi, xg, yg, Z) for (xi, yi) in zip(xs, ys)]
    end
    error("interp_z_at: expecting scalar x,y or vectors of equal length")
end

# Ensure a numeric scalar Float64 is returned for numbers or single-element containers
function _ensure_scalar_float(x)
    if isa(x, Number)
        return Float64(x)
    end
    # otherwise try to collect and take the first element
    try
        c = collect(x)
        if length(c) >= 1
            return Float64(c[1])
        else
            error("_ensure_scalar_float: container is empty")
        end
    catch err
        error("_ensure_scalar_float: cannot convert input to scalar Float64: $err")
    end
end

function plot_results()
    ne_gt::Int = 8 # number of elements in the mesh for the ground truth
    ne_exp::Int = 2 # number of elements in the mesh for the experiment 

    β_gt_list = [10.0, 50.0, 100.0, 1e3]
    η_gt_list = [60.0]
    FunctionClass_x_List = ["S2", "Q2"]
    # refine_list = [1, 2, 3] # refinement levels, ne = ne_exp^refine
    refine_list = [2] # refinement levels, ne = ne_exp^refine
    control = "force" # "force" or "velocity"

    viscosity_type_list = ["constant"] #,"bulk_viscosity"]
    FunctionClass_x_gt_list = ["S2"] # Function space for the ground truth

    exp_size = size(FunctionClass_x_List,1)*size(refine_list,1)
    η_mat = zeros(30, exp_size)
    β_mat = zeros(30, exp_size)

    WRITE_GT = true
    for viscosity_type in viscosity_type_list
        for FunctionClass_x_gt in FunctionClass_x_gt_list
            run_id = 1
            for β_gt in β_gt_list
                for η_gt in η_gt_list
                    filepath = string("/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/sim_experiments/cost_function_test/optimization/Stokes/$control/model_validation/$viscosity_type/$(FunctionClass_x_gt)_$(ne_gt)/$run_id")
                    filepath_gt = string("/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/sim_experiments/ground_truth/fem/Stokes/$control/$viscosity_type/$(FunctionClass_x_gt)_$(ne_gt)/$run_id")
                    for ref in refine_list
                        ne = ne_exp^ref            
                        @info "reading optimization with ne = $ne"
                        for FunctionClass_x in FunctionClass_x_List
                            
                        end
                        replot(filepath, filepath_gt)
                        compare_stats(filepath,filepath_gt)
                    end
                    post_analysis(filepath_gt, filepath)
                    run_id = run_id + 1
                end
            end
        end
    end
end

function set_time_window(step_len::Float64, data::AbstractArray; method::String="fixed", window_size::Float64=10.0)
    windows::Vector{AbstractArray} = Vector{AbstractArray}()
    time_windows::Vector{Float64} = Vector{Float64}()
    data_ranges::Vector{AbstractArray} = Vector{AbstractArray}()
    t_windows::Vector{Float64} = Vector{Float64}()

    function get_t_window(window_size::Float64, step_len::Float64, iter::Int, method::String)::Float64
        t_window = 0.0
        if method == "fixed"
            t_window = round(window_size*iter, digits=1)
        elseif method == "exponential"
            t_window = round(step_len*exp(3*(iter-1)), digits=1)
        end
        return t_window
    end

    iter::Int = 1
    start_point::Int = 1
    t_window_prev::Float64 = 0.0
    # t_window::Float64 = round(0.5*exp(3*(iter-1)), digits=1) 
    t_window_end::Float64 = get_t_window(window_size, step_len, iter, method)
    end_point::Int = round(Int,t_window_end*step_len)+1

    END_FLAG::Bool = false
    end_clause::String = ""
    while true
            # If the computed end_point reaches or exceeds the data length, adjust.
            # Keep the original requested value for clearer messages.
            if end_point >= size(data, 1)
                if end_point == size(data, 1)
                    @info "Reached end of data at point $end_point."
                else
                    requested_end = end_point
                    # adjust the time window end to the last available sample
                    t_window_end = round((size(data, 1)-1)/step_len, digits=1)
                    end_point = size(data, 1)
                    @warn "Requested end point $requested_end exceeds data size $(size(data,1)); adjusting to end of data (end_point=$end_point). Adjusted end time to $t_window_end seconds."
                end
                END_FLAG = true
            end

        data_range = start_point:end_point
        data_range_ = start_point:(end_point-1)
        
        println("Data frame : $data_range_")
        println("time windows from : $t_window_prev to $t_window_end")
        println("time window size : $(t_window_end - t_window_prev) seconds")
        println("----------")
        t_window_size = round(t_window_end - t_window_prev, digits=1)
        push!(time_windows, t_window_size)
        println("data range size : $(length(data_range_))")
        push!(windows, data[data_range])
        push!(data_ranges, data_range_)
        push!(t_windows, t_window_end)

        if END_FLAG == true
            break
        end
        iter = iter + 1
        start_point = end_point
        t_window_prev = t_window_end
        t_window_end = get_t_window(window_size, step_len, iter, method)
        end_point = round(Int,t_window_end*step_len)+1
    end
    return time_windows, windows, data_ranges, t_windows
end

# helper: run a vector of parameter Dicts with limited concurrent workers
function run_param_list(params_list::Vector{Dict}; max_workers::Int=8, base_seed::Int=12345)

    nparams = length(params_list)
    if nparams == 0
        return
    end

    workers = max(1, min(Threads.nthreads(), max_workers))
    @info "Running $nparams experiments with $workers workers (Threads.nthreads()=$(Threads.nthreads()))"

    ch = Channel{Int}(nparams)
    @sync begin
        # enqueue indices
        for i in 1:nparams
            put!(ch, i)
        end

        # close the channel so workers can exit when it's empty
        close(ch)

        tasks = Vector{Task}(undef, workers)
        for w in 1:workers
            tasks[w] = Threads.@spawn begin
                while true
                    idx = try
                        take!(ch)
                    catch
                        # channel closed and empty -> exit loop
                        break
                    end
                    params = params_list[idx]
                    try
                        # Avoid BLAS thread oversubscription per worker
                        LinearAlgebra.BLAS.set_num_threads(1)
                        optimize(params)
                    catch err
                        # capture backtrace and format similar to native Julia error output
                        bt = catch_backtrace()
                        formatted = sprint(io -> begin
                            try
                                println("Error during optimize for params $params:")
                                showerror(io, err, bt)
                            catch
                                # fallback to simple show if showerror with bt fails
                                println("Error during optimize for params $params:")
                                showerror(io, err)
                            end
                            try
                                println("Error during optimize for params $params:")
                                Base.show_backtrace(io, bt)
                            catch
                            end
                        end)
                        println("Error during optimize for params $params:")
                        @error "optimize failed for params index $idx:\n$formatted"
                        # Attempt to write a detailed error log next to the experiment results
                        try
                            dest_dir = "."
                            if isa(params, AbstractDict) && haskey(params, "filepath_res")
                                dest_dir = string(params["filepath_res"], "/Results/logs")
                            else
                                dest_dir = "./logs"
                            end
                            write_error_log(err, bt; params=params, dest_dir=dest_dir)
                        catch ewrite
                            @error "Failed to write error log: $ewrite"
                        end
                    end
                end
            end
        end

        # wait for all worker tasks to complete
        for t in tasks
            wait(t)
        end
    end

    # signal completion for this batch
    println("All experiments in this run_param_list completed.")
end

function optimize_sim()

    FunctionClass_x_List = ["Q2"]
    # refine_list = [1, 2, 3] # refinement levels, ne = ne_exp^refine
    refine_list = [4, 6, 8] # refinement levels, ne = ne_exp^refine
    noise_level_list = [0.0, 0.5, 1.0]
    control = "force" # "force" or "velocity"
    viscosity_type_list = ["constant"]

    camera_matrix::AbstractArray = [[2.39642674e+03, 0.0, 1.00429248e+03] [0.0, 2.40565353e+03, 7.57028161e+02] [0.0, 0.0, 1.0]]'
    filepath_res::String = ""
    param_list = Vector{Dict}(undef, 0)

    avoid_dirs = ["3_less_noise"]
    for viscosity_type in viscosity_type_list
        _filepath_gt = string("/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/ground_truth/sim_data/Stokes/$control/$viscosity_type/Q2_16")
        dir_list = readdir(_filepath_gt)
        for dir in dir_list
            if dir in avoid_dirs
                continue
                println("Skipping dir $dir")
            end
            @info "Processing ground truth directory: $dir for $viscosity_type"
            filepath_gt = string(_filepath_gt,"/",dir)
            for ne in refine_list
                if ne == 6 && viscosity_type == "constant"
                    noise_level_list = [0.0]
                else
                    noise_level_list = [0.0]
                end
                for noise_level in noise_level_list 
                    if noise_level == 0.0 && viscosity_type == "constant"
                        sim_time_exp_list = [30.0] # simulation time in seconds
                    elseif noise_level != 0.0 && viscosity_type == "constant" && ne == 6
                        sim_time_exp_list = [10.0, 20.0, 30.0] # simulation time in seconds
                    elseif noise_level == 0.0 && viscosity_type == "bulk_viscosity"
                        sim_time_exp_list = [5.0] # simulation time in seconds
                    else
                        sim_time_exp_list = [30.0] # simulation time in seconds
                    end
                    for sim_time_exp::Float16 in sim_time_exp_list
                        # ne = ne_exp^ref
                        @info "Running optimization with ne = $ne and simulation time = $sim_time_exp with noise level = $noise_level."
                        for FunctionClass_x in FunctionClass_x_List

                            start_time = Dates.now()
                            filepath_res = string("/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/experiments/sim_data/optimization/Stokes/$control/$viscosity_type/Q2_16/$dir/$(FunctionClass_x)_$(ne)/simtime_$(sim_time_exp)/noise_$(noise_level)/multi_window")
                            @info "Running optimization with FunctionClass_x = $FunctionClass_x with $ne elements"

                            exp_params = Dict("FunctionClass_x" => FunctionClass_x, "FunctionClass_u" => "Q2", "FunctionClass_p" => "Q1", "ne_exp" => ne, "sim_time_exp" => sim_time_exp, 
                            "filepath_res" => filepath_res, "filepath_gt"=>filepath_gt, "control" => control, "data_type"=>"simulated", "camera_matrix" => camera_matrix, "WRITE_GT"=> false,
                            "noise_level"=>noise_level, "mode"=>"multi_window")

                            optimize(exp_params)

                            end_time = Dates.now()
                            write_time_log(start_time, end_time, exp_params; dest_dir=string(filepath_res,"/Results/logs"))
                        end
                    end
                end
            end
        end
        # run_param_list(param_list; max_workers=11)
    end
end

function optimize_syn()

    FunctionClass_x_List = ["Q2"]
    # refine_list = [1, 2, 3] # refinement levels, ne = ne_exp^refine
    refine_list = [8] # refinement levels, ne = ne_exp^refine
    noise_level_list = [0.0, 0.5, 1.0]
    control = "force" # "force" or "velocity"
    viscosity_type_list = ["bulk_viscosity"]

    camera_matrix::AbstractArray = [[2.39642674e+03, 0.0, 1.00429248e+03] [0.0, 2.40565353e+03, 7.57028161e+02] [0.0, 0.0, 1.0]]'
    filepath_res::String = ""
    param_list = Vector{Dict}(undef, 0)

    avoid_dirs = ["3_less_noise", "2", "3", "4", "5"]
    for viscosity_type in viscosity_type_list
        _filepath_gt = string("/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/ground_truth/sim_data/Stokes/$control/$viscosity_type/Q2_16")
        dir_list = readdir(_filepath_gt)
        for dir in dir_list
            if dir in avoid_dirs
                continue
                println("Skipping dir $dir")
            end
            filepath_gt = string(_filepath_gt,"/",dir)
            @info "Processing ground truth directory: $filepath_gt for $viscosity_type"
            for ne in refine_list
                if ne == 4 && viscosity_type == "constant"
                    noise_level_list = [0.0]
                else
                    noise_level_list = [0.0]
                end
                for noise_level in noise_level_list 
                    if noise_level == 0.0 && viscosity_type == "constant"
                        sim_time_exp_list = [10.0, 20.0, 30.0] # simulation time in seconds
                    elseif noise_level != 0.0 && viscosity_type == "constant" && ne == 6
                        sim_time_exp_list = [10.0, 20.0, 30.0] # simulation time in seconds
                    elseif noise_level == 0.0 && viscosity_type == "bulk_viscosity"
                        sim_time_exp_list = [5.0, 10.0] # simulation time in seconds
                    else
                        sim_time_exp_list = [30.0] # simulation time in seconds
                    end
                    for sim_time_exp::Float16 in sim_time_exp_list
                        # ne = ne_exp^ref
                        @info "Running optimization with ne = $ne and simulation time = $sim_time_exp with noise level = $noise_level"
                        for FunctionClass_x in FunctionClass_x_List

                            start_time = Dates.now()
                            filepath_res = string("/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/experiments/syn_data/optimization/Stokes/$control/$viscosity_type/Q2_16/$dir/$(FunctionClass_x)_$(ne)/simtime_$(sim_time_exp)/noise_$(noise_level)/multi_window")

                            exp_params = Dict("FunctionClass_x" => FunctionClass_x, "FunctionClass_u" => "Q2", "FunctionClass_p" => "Q1", "ne_exp" => ne, "sim_time_exp" => sim_time_exp, 
                            "filepath_res" => filepath_res, "filepath_gt"=>filepath_gt, "control" => control, "data_type"=>"synthetic", "camera_matrix" => camera_matrix, "WRITE_GT"=> false,
                            "noise_level"=>noise_level, "mode"=>"multi_window")

                            optimize(exp_params)

                            end_time = Dates.now()
                            write_time_log(start_time, end_time, exp_params; dest_dir=string(filepath_res,"/Results/logs"))
                        end
                    end
                end
            end
        end
        # run_param_list(param_list; max_workers=11)
    end
end

function plot_syn()

    control = "force" # "force" or "velocity"
    viscosity_type_list = ["constant"] #,"constant"]
    avoid_dirs = ["3_less_noise"]
    for viscosity_type in viscosity_type_list
        filepath_gt = string("/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/ground_truth/sim_data/Stokes/$control/$viscosity_type/Q2_16/")
        filepath_res = string("/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/experiments/syn_data/optimization/Stokes/$control/$viscosity_type/Q2_16/")
        dirs = readdir(filepath_gt)
        for dir in dirs
            if dir in avoid_dirs
                continue
                println("Skipping dir $dir")
            end
            filepath_gt_dir = string(filepath_gt,"/$dir/")
            filepath_res_dir = string(filepath_res,"/$dir/")
            replot(filepath_res_dir, filepath_gt_dir)
        end
        if viscosity_type == "constant"
            post_analysis(filepath_gt, filepath_res, avoid_dirs)
        end
    end
end

function optimize_real()

    ne_exp::Int = 2 # number of elements in the mesh for the experiment 
    FunctionClass_x_List = ["Q2"]
    # refine_list = [1, 2, 3] # refinement levels, ne = ne_exp^refine
    refine_list = [2] # refinement levels, ne = ne_exp^refine
    control = "force" # "force" or "velocity"
    η_start = 1.0
    β_start = 1.0
    viscosity_type = "bulk_viscosity"

    r::Float64 = 25.0  # radius of the cylinder in mm
    h::Float64 = 40.0  # height of the cylinder in mm
    camera_matrix::AbstractArray = [[2.39642674e+03, 0.0, 1.00429248e+03] [0.0, 2.40565353e+03, 7.57028161e+02] [0.0, 0.0, 1.0]]'
    # sim_time_exp::Float64 = 5.0 # simulation time in seconds
    F_ext::Float64 = 1*9.812*1e3 # force applied to the cylinder in N

    sim_time_list = reverse([0.5]) # simulation time in seconds
    for sim_time_exp in sim_time_list
        println("Simulation time: $sim_time_exp seconds")
        file_id = 5
        filepath_gt = string("/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/ground_truth/physical_data/$file_id")
        filepath_res = string("/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/experiments/physical_data/integration_tests/single_window/$sim_time_exp")
        
        for ref in refine_list
            ne = 4
            @info "Running optimization with ne = $ne"
            for FunctionClass_x in FunctionClass_x_List
                @info "Running optimization with FunctionClass_x = $FunctionClass_x with $ne elements"

                exp_params = Dict("FunctionClass_x" => FunctionClass_x, "FunctionClass_u" => "Q2", "FunctionClass_p" => "Q1", "ne_exp" => ne, "sim_time_exp" => sim_time_exp, 
                "η_start" => η_start, "β_start" => β_start, "filepath_res" => filepath_res, "filepath_gt"=>filepath_gt, "control" => control, "viscosity_type"=>viscosity_type, 
                "data_type"=>"physical", "r" => r, "h" => h, "camera_matrix" => camera_matrix, "F_ext" => F_ext, "mode"=>"multi_window")

                # plot_(exp_params)
            end
        end
        # post_analysis(filepath_gt, filepath_res)
        file_id = file_id + 1
    end
end

# main()
# plot_syn()
optimize_sim()
# optimize_syn()
