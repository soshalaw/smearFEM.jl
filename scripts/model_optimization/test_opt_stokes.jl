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

global def_orange = theme_palette(:default)[2]
global end_obs_win = 25.1

# for 1/2 linewidth
# fs::Int = 12
# plt_height::Int = 350
# plt_width ::Int = 477
# plt_lft_margin = 1pt
# plt_right_margin = 5pt
# plt_top_margin = 1pt

# for 1/3 linewidth
global fs::Int = 10
global plt_height::Int = 360
global plt_width::Int = 330
global plt_lft_margin = -6pt
global plt_right_margin = 10pt
global plt_top_margin = -5pt

global y_lims_h_norm = (0.98, 1.002)
global y_lims_rel_error = (-0.05, 2)

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
    start_time = Dates.now()
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
    noiseLevel::Float64 = 0.0
    SIDES::Bool = false
    
    viscosity_model::String = "" # "constant" or "bulk_viscosity"
    if data_type == "simulated" || data_type == "synthetic"
        
        WRITE_GT = exp_params["WRITE_GT"] 
        noiseLevel = exp_params["noise_level"]
        outlier_frames = Int[]
        
        if WRITE_GT == true # write the ground truth data
            @info "Writing ground truth gt data to with $ne_gt elements to $filepath_res"
            write_gt_data(exp_params)
        end
        
        params = read_json(joinpath(filepath_gt,"data","sim_params.json"))
        
        r = params["r"]
        h = params["h"]
        
        gt_viscosity_type = params["viscosity_type"]
        F = Array(float.(params["cParam"]))
        
        if gt_viscosity_type == "bulk_viscosity"
            if haskey(params, "model_type") && params["model_type"] == "carreau"
                viscosity_model = params["model_type"]
            else
                viscosity_model = "power_law"
                η_gt = params["η"][1]
                β_gt = params["β"][1]
            end
        else
            η_gt = params["η"][1]
            β_gt = params["β"][1]
        end

        sim_time_gt = params["simulation_time"]
        t_steps_gt = params["time_steps"]
        
        camera_matrix = reshape(Array(params["camera_matrix"]), 3, 3)
        obj_pose = reshape(Array(params["obj_pose"])*1.0,4,4)
        
        control = params["control_type"]
        
        model_gt, scene_gt = def_problem(r, h, ne_exp, η_gt, ndim, FunctionClass_u, nDof_u, FunctionClass_p, nDof_p, FunctionClass_x, β_gt, F, control, gt_viscosity_type, 
        sim_time_gt, t_steps_gt, viscosity_model=viscosity_model)
        
        if data_type == "synthetic"
            ObsDataList, splinexObs, splineyObs = read_csv(joinpath(filepath_gt,"data","img_data","contour_data"))  
        elseif data_type == "simulated"
            ObsDataList, splinexObs, splineyObs = read_csv(joinpath(filepath_gt,"data","sim_data","contour_data"))  
        end
        
        printstyled("Ground truth η: $(η_gt), ground truth β: $(β_gt)\n"; color = :green)
        
    elseif data_type == "physical"
        r::Float64 = exp_params["r"]  # radius of the cylinder in mm
        h::Float64 = exp_params["h"]  # height of the cylinder in mm
        
        control::String = exp_params["control"]
        gt_viscosity_type::String = exp_params["viscosity_type"] # "constant" or "bulk_viscosity"

        t_obs, t_top_plt, t_btm_plt = read_perception_data(joinpath(filepath_gt, "data", "sequence.hdf5"))

        _ObsDataList, _splinexObs, _splineyObs = read_csv(joinpath(filepath_gt, "data", "img_data", "contour_data"))
        meta_data = read_json(joinpath(filepath_gt, "data", "video_metadata.json"))

        frame_rate::Float64 = round(meta_data["frame_rate"], digits=1)
        frame_width::Int = meta_data["frame_width"]
        frame_height::Int = meta_data["frame_height"]
        compression_frames = meta_data["compressed_frames"]
        
        sim_time_gt = length(compression_frames)/frame_rate  # seconds
        steps_exp = sim_time_exp*frame_rate
        t_steps_exp = 1/frame_rate
        t_steps_gt = t_steps_exp
        
        F_ext = exp_params["F_ext"]
        println("Applied force: $F_ext N")
        
        F = -F_ext*ones(Float64, round(Int, sim_time_gt*frame_rate)) # force applied to the cylinder in N
        
        obj_pose = get_pose(t_obs)  # Example usage of get_pose function
        obj_pose_ = zeros(Float64, 4,4)
        obj_pose_[1,1] = -1.0
        obj_pose_[2,3] = -1.0
        obj_pose_[3,2] = -1.0
        obj_pose_[1:3,4] = obj_pose[1:3,4]
        
        ObsDataList = _ObsDataList[compression_frames]
        splinexObs = _splinexObs[compression_frames]
        splineyObs = _splineyObs[compression_frames]
        
        function get_gt_height(t_top_plt::AbstractArray, t_btm_plt::AbstractArray, h_gt_init::Float64)
            n_frames = size(t_top_plt,2)
            Δh_gt = zeros(Float64, n_frames)
            h_init = abs(t_top_plt[2,1] - t_btm_plt[2,1])
            for i::Int in 1:n_frames
                Δh_gt[i] = h_init - abs(t_top_plt[2,i] - t_btm_plt[2,i])
            end

            h_gt = h_gt_init .- Δh_gt

            return h_gt, Δh_gt
        end

        h_gt , _ = get_gt_height(t_top_plt[1:3,4,compression_frames], t_btm_plt[1:3,4,compression_frames], 38.5)

        sim_data = Dict("camera_matrix" => camera_matrix,
                    "obj_pose" => obj_pose,
                    "time_steps" => t_steps_exp,
                    "simulation_time" => sim_time_gt,
                    "cParam" => F,
                    "r" => r,
                    "h" => h,
                    "control_type" => control,
                    "viscosity_type" => gt_viscosity_type,
                    "data_type" => data_type)

        write_json(joinpath(filepath_gt,"data","sim_params"), sim_data)
        write_csv(joinpath(filepath_gt,"data","h"), h_gt)

        valid_frames, outlier_frames = detect_outlier_observations(ObsDataList)
    else
        error("data_type should be either simulated or physical")
    end

    # if sim_time_gt < sim_time_exp
    #     @warn "Ground truth simulation time $sim_time_gt is less than experimental simulation time $sim_time_exp , switching to ground truth simulation time"
    #     sim_time_exp = sim_time_gt
    # end
    # exp_path = filepath_res
    # set_file(exp_path)
    # write_json(joinpath(exp_path,"Results","data","experiment_parameters"), exp_params)

    # if t_steps_exp > t_steps_gt
    #     @warn "time resolution of the ground truth $t_steps_gt is larger than the experimental $t_steps_exp, switching to ground truth resolution"
    #     t_steps_exp = t_steps_gt
    # end
    
    # if data_type == "physical"  || viscosity_model == "carreau"
    #     η_start = exp_params["η_start"]
    #     β_start = exp_params["β_start"]
    # else
    #     dev::Float64 = 0.3

    #     dev_η::Float64 = dev*η_gt
    #     η_start::Float64 = abs(η_gt - dev_η)

    #     dev_β::Float64 = dev*β_gt
    #     β_start::Float64 = abs(β_gt - dev_β)
    # end

    # θ::Vector{Float64} = [η_start, β_start]

    # if gt_viscosity_type == "constant"
        
    #     _range = 1:(round(Int,sim_time_exp/t_steps_exp)+1)
    #     @info "Considering from frame $(first(_range)) to frame $(last(_range)) in the observations"
    #     ObsDataList = ObsDataList[_range] # align the observation points with the simulation time

    #     # Read the gt data
    #     model, scene = def_problem(r, h, ne_exp, η_gt, ndim, FunctionClass_u, nDof_u, FunctionClass_p, nDof_p, FunctionClass_x, β_gt, F, control, gt_viscosity_type, 
    #     sim_time_exp, t_steps_exp)
    #     conditions = Conditions(camera_matrix=camera_matrix, obj_pose=obj_pose, SIDES=SIDES, filepath=exp_path, ANIMATE=false)
        
    #     gt_h_ = readdlm(joinpath(filepath_gt,"data","h.csv"), ',', Float64)
    #     gt_h = gt_h_[1:(round(Int,sim_time_gt/t_steps_exp)+1)]

    #     # time = collect(Float64, range(start=0, stop=sim_time_exp, step=t_steps_exp))
    #     if noiseLevel == 0.0
                
    #         obsBorderPts, nSplinex, nSpliney, pd = add_noise(ObsDataList, nFactor=0.0)
    #         stats = fit_model(model, scene, conditions, obsBorderPts, θ, outliers=outlier_frames)

    #         iterList = stats["iterList"]
    #         costList = stats["cost_list"]
    #         ηpList = stats["ηList"]
    #         βpList = stats["βList"]

    #         η = stats["η"]
    #         β = stats["β"]

    #         printstyled("Estimated η : $(η), estimated β: $(β)\n"; color = :green)

    #         η_accuracy = (1-abs((η_gt-η)/η_gt))*100
    #         β_accuracy = (1-abs((β_gt-β)/β_gt))*100
    #         printstyled("η accuracy: $(η_accuracy) %\n"; color = :green)
    #         printstyled("β accuracy: $(β_accuracy) %\n"; color = :green)

    #         set_file(joinpath(exp_path,"Results","plots"))
            
    #         write_json(joinpath(exp_path,"Results","data","stats"), stats)
    #         write_csv(joinpath(exp_path,"Results","data","η"), ηpList)
    #         write_csv(joinpath(exp_path,"Results","data","β"), βpList)
    #         write_csv(joinpath(exp_path,"Results","data","gt_h"), gt_h)
    #         write_csv(joinpath(exp_path,"Results","data","cost_iter"), costList)
            
    #         reset_model!(model)
    #         model.η = [η]
    #         scene.β = [β]
    #         scene.sim_time = sim_time_gt
            
    #         # simulate the model with the estimated parameters
    #         est_μ_list, gradList, borderPts2DList, fields, pos3D, pos2D, splinep, splineq = simulate(model, scene, conditions)
    #         est_h = get_height(est_μ_list, h)
            
    #         write_csv(joinpath(exp_path,"Results","data","est_h"), est_h)
    #         write_data(joinpath(exp_path,"Results","data","sim_data","2D_surface_points"), pos2D)
    #         write_data(joinpath(exp_path,"Results","data","sim_data","3D_points"), pos3D)
    #         write_data(joinpath(exp_path,"Results","data","sim_data","motion_fields "), fields)
    #         write_data(joinpath(exp_path,"Results","data","sim_data","2D_border_points"), borderPts2DList)

    #         # if maximum(ηpList) > η+dev_η
    #         #     ηStop = maximum(ηpList)*1.1
    #         # else
    #         #     ηStop = η+dev_η
    #         # end

    #         # if minimum(ηpList) < η-dev_η
    #         #     if minimum(ηpList) <= 0.0
    #         #         η_start = 0.1
    #         #     else
    #         #         η_start = minimum(ηpList)*0.9
    #         #     end
    #         # else
    #         #     η_start = η-dev_η
    #         # end

    #         # if maximum(βpList) > β+dev_β
    #         #     βStop = maximum(βpList)*1.1
    #         # else
    #         #     βStop = β+dev_β
    #         # end

    #         # if minimum(βpList) < β-dev_β
    #         #     β_start = minimum(βpList)*0.9
    #         # else
    #         #     β_start = β-dev_β
    #         # end

    #         # sampleNo = 10
    #         # ηList = collect(range(η_start, stop=ηStop, length=sampleNo))
    #         # βList = collect(range(β_start, stop=βStop, length=sampleNo))
    #         # CostMat = zeros(size(ηList,1),size(βList,1))
    #         # ∂CostMat = zeros(2,size(ηList,1),size(βList,1))
    #         # ∂2CostMat = zeros(2,2,size(ηList,1),size(βList,1))

    #         # costη = zeros(size(ηList,1))
    #         # costβ = zeros(size(βList,1))

    #         # η_iter = 1:size(ηList,1)
    #         # β_iter = 1:size(βList,1)

    #         # for i::Int in η_iter
    #         #     η = ηList[i]
    #         #     for j::Int in β_iter
    #         #         β = βList[j]
    #         #         reset_model!(model)
    #         #         model.η = [η]
    #         #         scene.β = [β]
    #         #         μ_list, gradList, simBorderPts, fields_, pos3D_, pos2D_, splinex_, spliney_ = simulate(model, scene, conditions)
                    
    #         #         # test the closest point function
    #         #         d, pairs = closest_point(simBorderPts, obsBorderPts)
                    
    #         #         CostMat[i,j] = sum(d)/length(d)
    #         #     end
    #         # end
            
    #         # # Plot the cost function surface (interactive GLMakie)
    #         # # set_plot(fs, sz=(plt_width, plt_height))
    #         # # # CostMat is constructed as CostMat[i_eta, j_beta]; Makie expects Z with
    #         # # # size (length(beta), length(eta)) => transpose CostMat for plotting.
    #         # # # First attempt: save an interactive PlotlyJS HTML (avoids creating Makie figures)
    #         # # try
    #         # #     htmlpath = joinpath(exp_path, "Results", "plots", "cost_surface_iter_interactive.html")
    #         # #     # compute overlay z-values for estimator path
    #         # #     Z_for_interp = CostMat'
    #         # #     zs_est = [interp_z_at(a, b, ηList, βList, Z_for_interp) for (a,b) in zip(ηpList, βpList)]
    #         # #     # compute ground-truth z for legend/marker overlay
    #         # #     z_gt = interp_z_at(η_gt, β_gt, ηList, βList, Z_for_interp)
    #         # #     saved = save_plotly_surface_html(htmlpath, ηList, βList, CostMat'; xs=ηpList, ys=βpList, zs=zs_est, title="Cost surface (iter)", gt_x=η_gt, gt_y=β_gt, gt_z=z_gt, x_label = "\$\\eta\\;\\mathrm{(KPa\\cdot s)}\$", y_label = "\$\\beta\\;\\mathrm{(L/mm^{-1})}\$", font_size = fs, latex_labels=true)
    #         # #     if saved
    #         # #         @info "Saved interactive PlotlyJS HTML: $htmlpath"
    #         # #     else
    #         # #         @warn "PlotlyJS HTML not created; skipping interactive output. To enable interactive PNGs with GLMakie re-enable GLMakie manually."
    #         # #     end
    #         # # catch err
    #         # #     @warn "PlotlyJS path failed; skipping interactive output. Error: $err"
    #         # # end

    #         # # # Also produce static PDF outputs with Plots (preserve previous behavior)
    #         # # try
    #         # #     plt = set_plot(fs, sz=(plt_width, plt_height))
    #         # #     Plots.contour!(plt, ηList, βList, CostMat, color=:turbo, fill=false, levels=100, dpi=400)
    #         # #     Plots.plot!(plt, ηpList, βpList, label="Estimations", ms=:4, m=:x, color=:red)
    #         # #     Plots.plot!(plt, [η_gt], [β_gt], label="Ground truth", ms=:8, m=:star5, color=def_orange)
    #         # #     Plots.xlabel!(plt, L"\eta\;\mathrm{(KPa\cdot s)}")
    #         # #     Plots.ylabel!(plt, L"\beta\;\mathrm{(L/mm^{-1})}")
    #         # #     Plots.savefig(plt, joinpath(exp_path,"Results","plots","cost_surface_iter.pdf"))

    #         # #     plt2 = set_plot(fs, sz=(plt_width, plt_height))
    #         # #     Plots.contourf!(plt2, ηList, βList, CostMat, color=:turbo, fill=false, levels=100, dpi=400)
    #         # #     Plots.plot!(plt2, [η_gt], [β_gt], label="Ground truth", ms=:8, m=:star5, color=def_orange)
    #         # #     Plots.xlabel!(plt2, L"\eta\;\mathrm{(KPa\cdot s)}")
    #         # #     Plots.ylabel!(plt2, L"\beta\;\mathrm{(L/mm^{-1})}")
    #         # #     Plots.savefig(plt2, joinpath(exp_path,"Results","plots","cost_surface.pdf"))
    #         # # catch err
    #         # #     @warn "Failed to produce static PDF contour outputs: $err"
    #         # # end

    #         # # Write the results to files
    #         # contour_plot_params = Dict("η_list" => ηList, "β_list" => βList, "cost_mat" => CostMat)

    #         # write_json(joinpath(exp_path,"Results","data","contour_plot_params"), contour_plot_params)

    #     # @save joinpath(exp_path,"Results","data","sim_data","Cost_Matrices.jld2") ηList, βList, CostMat, ∂CostMat, ∂2CostMat
    #     else
    #         n_samples = 10
    #         η_pred = zeros(Float64, n_samples)
    #         β_pred = zeros(Float64, n_samples)
    #         costnList = Vector{Vector{Float64}}(undef, n_samples)
    #         iternList = Vector{AbstractVector}(undef, n_samples)
    #         est_h_list = Matrix{Float64}(undef, n_samples, round(Int,sim_time_exp/t_steps_exp)+1)
    #         ANIMATED = false

    #         set_file(joinpath(exp_path,"Results","plots"))
    #         for n::Int in 1:n_samples
    #             obsBorderPts, nSplinex, nSpliney, pd = add_noise(ObsDataList, nFactor=noiseLevel)

    #             stats = fit_model(model, scene, conditions, obsBorderPts, θ, outliers=outlier_frames)
                                    
    #             iterList = stats["iterList"]
    #             costList = stats["cost_list"]
    #             ηpList = stats["ηList"]
    #             βpList = stats["βList"]

    #             η = stats["η"]
    #             β = stats["β"]
                
    #             reset_model!(model)
    #             model.η = [η]
    #             scene.β = [β]
    #             scene.sim_time = sim_time_gt
    #             # simulate the model with the estimated parameters
    #             est_μ_list, gradList, borderPts2DList, fields, pos3D, pos2D, splinep, splineq = simulate(model, scene, conditions)
                
    #             est_h = get_height(est_μ_list, h)
                
    #             η_accuracy = (1-abs((η_gt-η)/η_gt))*100
    #             β_accuracy = (1-abs((β_gt-β)/β_gt))*100
    #             printstyled("η accuracy: $(η_accuracy) %\n"; color = :green)
    #             printstyled("β accuracy: $(β_accuracy) %\n"; color = :green)
                
    #             η_pred[n] = η
    #             β_pred[n] = β
    #             est_h_list[n,:] = est_h
    #             costnList[n] = costList
    #             iternList[n] = iterList

    #             params = Dict("gt_η" => η_gt,
    #             "gt_β" => β_gt,
    #             "η" => η,
    #             "β" => β,
    #             "η_accuracy" => η_accuracy,
    #             "β_accuracy" => β_accuracy)

    #             write_json(joinpath(exp_path,"Results","data","stats","run_$n"), params)

    #             write_csv(joinpath(exp_path,"Results","data","opt_data","cost_steps","run_$n"), costList)
    #             write_csv(joinpath(exp_path,"Results","data","opt_data","eta_steps","run_$n"), ηpList)
    #             write_csv(joinpath(exp_path,"Results","data","opt_data","beta_steps","run_$n"), βpList)
    #             write_csv(joinpath(exp_path,"Results","data","opt_data","iter","run_$n"), iterList)
    #             write_csv(joinpath(exp_path,"Results","data","opt_data","est_height","run_$n"), est_h)

    #             write_data(joinpath(exp_path,"Results","data","sim_data","2D_points","run_$n"), pos2D)
    #             write_data(joinpath(exp_path,"Results","data","sim_data","3D_points","run_$n"), pos3D)
    #             write_data(joinpath(exp_path,"Results","data","sim_data","motion_fields ","run_$n"), fields)
    #             write_data(joinpath(exp_path,"Results","data","sim_data","2D_border_points","run_$n"), borderPts2DList)
    #             write_data(joinpath(exp_path,"Results","data","sim_data","spline_p","run_$n"), splinep)
    #             write_data(joinpath(exp_path,"Results","data","sim_data","spline_q","run_$n"), splineq)

    #         end

    #         write_csv(joinpath(exp_path,"Results","data","eta_est"), η_pred)
    #         write_csv(joinpath(exp_path,"Results","data","beta_est"), β_pred)
    #         write_csv(joinpath(exp_path,"Results","data","h_est"), est_h_list)

    #         @info "Data writing completed. Results saved to $exp_path"
    #     end

    # elseif gt_viscosity_type == "bulk_viscosity"
        
    #     av_η::Float64 = 0.0
    #     obsBorderPts, nSplinex, nSpliney, pd = add_noise(ObsDataList, nFactor=0.0)
    #     mode::String = exp_params["mode"] # "single_window" or "multiple_window" or "full_time"
    #     sim_time_window::Float64 = 30.0 # time window size for optimization in seconds
        
    #     viscosity_type = "constant"
        
    #     conditions = Conditions(camera_matrix=camera_matrix, obj_pose=obj_pose)
    #     model, scene = def_problem(r, h, ne_exp, η_gt, ndim, FunctionClass_u, nDof_u, FunctionClass_p, nDof_p, FunctionClass_x, β_gt, F, control, viscosity_type, 
    #     sim_time_exp, t_steps_exp)
        
    #     set_file(joinpath(exp_path,"Results","plots"))
        
    #     time_windows, windows, data_ranges_, t_windows = set_time_window(1/t_steps_exp, obsBorderPts, method="quadratic", window_size=sim_time_exp)
    #     _, splinexObs_win, _, _ = set_time_window(1/t_steps_exp, splinexObs, method="quadratic", window_size=sim_time_exp)
    #     _, splineyObs_win, _, _ = set_time_window(1/t_steps_exp, splineyObs, method="quadratic", window_size=sim_time_exp)
    #     println("Time windows: $(time_windows)")
    #     obs_time = sum(time_windows)

    #     if obs_time < sim_time_gt
    #         @warn "Observation time frame $obs_time is less than preset ground truth time frame $sim_time_gt, switching to observation time frame"
    #         sim_time_gt = obs_time
    #     end

    #     if obs_time < sim_time_exp
    #         @warn "Observation time frame $obs_time is less than experimental simulation time frame $sim_time_exp, switching to observation time frame"
    #         sim_time_exp = obs_time
    #     end
        
    #     data_pt_len = round(Int,obs_time/t_steps_exp)
    #     est_ηpList = Vector{Float64}(undef,data_pt_len)
    #     avg_ηList = Vector{Float64}(undef,data_pt_len)
    #     est_βpList = Vector{Float64}(undef,data_pt_len)

    #     if mode == "single_window"
    #         @info "Optimizing over a single time window"
    #         ti = 1
    #         data_range_ = data_ranges_[ti]
    #         scene.sim_time = time_windows[ti]
    #         if data_type == "physical"
    #             _F = -F_ext*ones(Float64, round(Int, scene.sim_time*frame_rate)) # force applied to the cylinder in N
    #             scene.cParam = _F
    #         else
    #             scene.cParam = F[data_range_]
    #         end
    #         @info "Time window $(time_windows[ti])"

    #         println("Data frame : $(data_range_)")
    #         println("Time frame : $(scene.sim_time)")

    #         printstyled("Time window: $(ti), time frames: $(scene.sim_time)\n"; color = :blue)

    #         obsBorderPts_t = windows[ti] # align the observation points with the simulation time
    #         splinexObs_t = splinexObs_win[ti]
    #         splineyObs_t = splineyObs_win[ti]

    #         println("observation Window size: $(size(obsBorderPts_t,1)) seconds")

    #         println("Time frame : $data_range_")
    #         printstyled("Time window: $(ti)\n"; color = :green)

    #         if data_type != "physical"
    #             η_gt_ = model_gt.η[data_range_]
    #             av_η = mean(η_gt_)
    #             avg_ηList[data_range_] .= av_η
    #             printstyled("Average ground truth η in the window: $(av_η), ground truth β: $(β_gt)\n"; color = :green)
    #         end

    #         stats = fit_model(model, scene, conditions, obsBorderPts_t, θ, outliers=outlier_frames)

    #         est_ηpList[data_range_] .= stats["η"]
    #         est_βpList[data_range_] .= stats["β"]

    #         θ[1] = stats["η"]
    #         θ[2] = stats["β"]

    #         update_model!(model)
      
    #         iterList = stats["iterList"]
    #         costList = stats["cost_list"]
    #         ηpList = stats["ηList"]
    #         βpList = stats["βList"]

    #         viscosity_type = "bulk_viscosity"
    #         est_model, est_scene = def_problem(r, h, ne_exp, η_gt, ndim, FunctionClass_u, nDof_u, FunctionClass_p, nDof_p, FunctionClass_x, β_gt, F[data_range_], control, viscosity_type, 
    #                             sim_time_exp, t_steps_exp)
    #         est_model.η = est_ηpList[data_range_] 
    #         est_scene.β = est_βpList[data_range_]
    #         est_μ_list, gradList, borderPts2DList, fields_est, pos3D_est, pos2D_est, splinex_est, spliney_est = simulate(est_model, est_scene, conditions)

    #         est_h_list = get_height(est_μ_list, h)

    #         animate_fields(filepath=joinpath(exp_path,"Results","plots"), p=splinex_est, q=spliney_est, pObs=splinexObs_t, qObs=splineyObs_t)
    #     else
    #         println("Number of time windows: $(length(windows))")

    #         for ti::Int in 1:length(windows)
    #             data_range_ = data_ranges_[ti]
    #             scene.sim_time = time_windows[ti]
    #             if data_type == "physical"
    #                 _F = -F_ext*ones(Float64, round(Int, scene.sim_time*frame_rate)) # force applied to the cylinder in N
    #                 scene.cParam = _F
    #             else
    #                 scene.cParam = F[data_range_]
    #             end
    #             obsBorderPts_t = windows[ti] # align the observation points with the simulation time

    #             printstyled("Time window: $(ti), time frames: $(scene.sim_time)\n"; color = :blue)
    #             println("Data frame : $(data_range_)")
    #             println("Time frame : $(scene.sim_time)")
    #             @info "Time window $(t_windows[ti])"
    #             println("observation Window size: $(size(obsBorderPts_t,1)) seconds")

    #             if data_type != "physical" && viscosity_model != "carreau"
    #                 η_gt_ = model_gt.η[data_range_]
    #                 av_η = mean(η_gt_)
    #                 avg_ηList[data_range_] .= av_η
    #                 printstyled("Average ground truth η in the window: $(av_η), ground truth β: $(β_gt)\n"; color = :green)
    #             end

    #             @info "Fitting model in time window $(ti)..."
    #             stats = fit_model(model, scene, conditions, obsBorderPts_t, θ, outliers=outlier_frames)

    #             est_ηpList[data_range_] .= stats["η"]
    #             est_βpList[data_range_] .= stats["β"]

    #             if data_type == "physical" 
    #                 θ[1] = stats["η"]
    #             else
    #                 θ[1] = stats["η"]
    #             end
    #             θ[2] = stats["β"]

    #             update_model!(model)
    #         end

    #         @info "Completed all time windows."
            
    #         viscosity_type = "bulk_viscosity"
    #         est_model, est_scene = def_problem(r, h, ne_exp, η_start, ndim, FunctionClass_u, nDof_u, FunctionClass_p, nDof_p, FunctionClass_x, β_start, F, control, viscosity_type, 
    #         sim_time_gt, t_steps_gt, viscosity_model=viscosity_model)
    #         est_model.η = est_ηpList
    #         est_scene.β = est_βpList
            
    #         write_csv(joinpath(exp_path,"Results","data","est_η"), est_ηpList)
    #         write_csv(joinpath(exp_path,"Results","data","est_β"), est_βpList)
    #         write_csv(joinpath(exp_path,"Results","data","avg_η"), avg_ηList)
    #         write_csv(joinpath(exp_path,"Results","data","window_data","time_windows"), time_windows)
    #         write_csv(joinpath(exp_path,"Results","data","window_data","t_windows"), t_windows)
    #         write_csv(joinpath(exp_path,"Results","data","window_data","data_ranges"), data_ranges_)
    #         write_csv(joinpath(exp_path,"Results","data","window_data","windows_sizes"), windows)
            
    #         est_μ_list, gradList, simBorderPts, fields_est, pos3D_est, pos2D_est, splinex_est, spliney_est = simulate(est_model, est_scene, conditions)
    #         est_h_list = get_height(est_μ_list, h)
            
    #         if data_type != "physical" && viscosity_model != "carreau"
    #             gt_μ_list, gradList, borderPts2DList, fields_gt, pos3D_gt, pos2D_gt, splinex_gt, spliney_gt = simulate(model_gt, scene_gt, conditions)
    #             gt_h_list = get_height(gt_μ_list, h)
    #             write_csv(joinpath(exp_path,"Results","data","η_gt"), model_gt.η)
    #             write_csv(joinpath(exp_path,"Results","data","β_gt"), β_gt)
    #             write_csv(joinpath(exp_path,"Results","data","gt_h"), gt_h_list)
                
    #             write_data(joinpath(exp_path,"Results","data","sim_data","2D_surface_points_gt"), pos2D_gt)
    #             write_data(joinpath(exp_path,"Results","data","sim_data","3D_points_gt"), pos3D_gt)
    #             write_data(joinpath(exp_path,"Results","data","sim_data","motion_fields_gt "), fields_gt)
    #             write_data(joinpath(exp_path,"Results","data","sim_data","2D_border_points_gt"), borderPts2DList)
    #         end
    #     end

    #     write_csv(joinpath(exp_path,"Results","data","est_h"), est_h_list)

    #     write_data(joinpath(exp_path,"Results","data","sim_data","2D_surface_points"), pos2D_est)
    #     write_data(joinpath(exp_path,"Results","data","sim_data","3D_points"), pos3D_est)
    #     write_data(joinpath(exp_path,"Results","data","sim_data","motion_fields "), fields_est)
    #     write_data(joinpath(exp_path,"Results","data","sim_data","2D_border_points"), simBorderPts)
    # end
    # end_time = Dates.now()
    # write_time_log(start_time, end_time, exp_params; dest_dir=joinpath(exp_path,"Results","logs"))
    # @info "Data writing completed. Results saved to $exp_path"
end

function compare_stats(filepath,filepath_gt)

    path_1 = []
    path_2 = []
    path_3 = []
    path_4 = []
    path_5 = []
    run_filepath = readdir(joinpath(filepath,"runs"))

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

        plot_η = set_plot(fs, sz=(plt_width, plt_height))
        plot_β = set_plot(fs, sz=(plt_width, plt_height))
        plot_cost = set_plot(fs, sz=(plt_width, plt_height))
        plot_error = set_plot(fs, sz=(plt_width, plt_height))

        for file in path_1
            exp_path = joinpath(filepath,"runs",file)
            sim_params = read_json(joinpath(exp_path,"Results","data","sim_params.json")) 

            sim_time = 10.0 # simulation time in seconds 
            steps = 25.0 # number of time steps
            t_steps = sim_time/steps

            time = collect(Float64, range(start=0, stop=sim_time, step=t_steps))

            viscosity_type = sim_params["viscosity_type"]

            if viscosity_type == "constant"
                η_gt = sim_params["η_gt"]
                β_gt = sim_params["β_gt"]

                est_η = readdlm(joinpath(exp_path,"Results","data","η.csv"), ',', Float64)
                    est_β = readdlm(joinpath(exp_path,"Results","data","β.csv"), ',', Float64)
                    est_h = readdlm(joinpath(exp_path,"Results","data","est_h.csv"), ',', Float64)
                    gt_h = readdlm(joinpath(exp_path,"Results","data","gt_h.csv"), ',', Float64)

                h_error = abs.(est_h-gt_h)

                stats = read_json(joinpath(exp_path,"Results","data","stats.json")) 

                costList = stats["cost_list"]
                iterList = stats["iterList"]

                label = " "
                if file[1:2] == "Q2"
                    label = "Lagrange basis"
                elseif file[1:2] == "S2"
                    label = "NURBS basis"
                end

                plot!(plot_η, iterList, est_η, label=label)
                xlabel!(plot_η, L"\mathrm{Iterations}")
                ylabel!(plot_η, L"\eta\;\mathrm{(KPa\cdot s)}")

                plot!(plot_β, iterList, est_β, label=label)
                xlabel!(plot_β, L"\mathrm{Iterations}")
                ylabel!(plot_β, L"\beta\;\mathrm{(mm^{-1})}")

                plot!(plot_cost, iterList, costList, label=label)
                xlabel!(plot_cost, L"\mathrm{Iterations}")
                ylabel!(plot_cost, L"\mathrm{Cost\;(px)}")

                plot!(plot_cost_log, iterList, costList, label=label, yscale=:log10)
                xlabel!(plot_cost_log, L"\mathrm{Iterations}")
                ylabel!(plot_cost_log, L"\mathrm{Cost\;(px)}")

                plot!(plot_error, time, h_error, label=label)
                xlabel!(plot_error, L"\mathrm{Time\;(s)}")
                ylabel!(plot_error, L"\mathrm{Height\;Error\;(mm)}")

                plot!(plot_h, time, est_h, label=label)
                xlabel!(plot_h, L"\mathrm{Time\;(s)}")
                ylabel!(plot_h, L"\mathrm{Height\;(mm)}")

                push!(cost_list_list,costList)
                push!(est_η_list,est_η)
                push!(est_β_list,est_β)
                push!(h_error_list,h_error)
                push!(iter_list,iterList)
            end
        end

        res_path = joinpath(filepath,"runs","post_analysis","ne_2")
        set_file(res_path)

        Plots.hline!(plot_η, [η_gt], label="Ground truth η")
        savefig(plot_η, joinpath(res_path,"eta.pdf"))

        Plots.hline!(plot_β, [β_gt], label="Ground truth β")
        savefig(plot_β, joinpath(res_path,"beta.pdf"))

        savefig(plot_cost, joinpath(res_path,"convergence.pdf"))
        savefig(plot_cost_log, joinpath(res_path,"convergence_log.pdf"))
        savefig(plot_error, joinpath(res_path,"error.pdf"))

        plot!(plot_h, time, gt_h, label=L"h_{gt}(t)")
        savefig(plot_h, joinpath(res_path,"h.pdf"))
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

        plot_η = set_plot(fs, sz=(plt_width, plt_height))
        plot_β = set_plot(fs, sz=(plt_width, plt_height))
        plot_cost = set_plot(fs, sz=(plt_width, plt_height))
        plot_cost_log = set_plot(fs, sz=(plt_width, plt_height))
        plot_error = set_plot(fs, sz=(plt_width, plt_height))
        plot_h = set_plot(fs, sz=(plt_width, plt_height))

        for file in path_2
            exp_path = joinpath(filepath,"runs",file)
            println(exp_path)
            sim_params = read_json(joinpath(exp_path,"Results","data","sim_params.json")) 

            ObsDataList, splinex, spliney = read_csv(joinpath(exp_path,"Results","data","sim_data","2D_border_points"))

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

                est_η = readdlm(joinpath(exp_path,"Results","data","η.csv"), ',', Float64)
                est_β = readdlm(joinpath(exp_path,"Results","data","β.csv"), ',', Float64)
                est_h = readdlm(joinpath(exp_path,"Results","data","est_h.csv"), ',', Float64)
                gt_h = readdlm(joinpath(exp_path,"Results","data","gt_h.csv"), ',', Float64)

                h_error = abs.(est_h-gt_h)

                stats = read_json(joinpath(exp_path,"Results","data","stats.json")) 

                costList = stats["cost_list"]
                iterList = stats["iterList"]

                label = " "
                if file[1:2] == "Q2"
                    label = "Lagrange basis"
                elseif file[1:2] == "S2"
                    label = "NURBS basis"
                end

                plot!(plot_η, iterList, est_η, label=label)
                xlabel!(plot_η, L"\mathrm{Iterations}")
                ylabel!(plot_η, L"\eta\;\mathrm{(KPa\cdot s)}")

                plot!(plot_β, iterList, est_β, label=label)
                xlabel!(plot_β, L"\mathrm{Iterations}")
                ylabel!(plot_β,L"\beta\;\mathrm{(mm^{-1})}")

                plot!(plot_cost, iterList, costList, label=label)
                xlabel!(plot_cost, L"\mathrm{Iterations}")
                ylabel!(plot_cost, L"\mathrm{Cost\;(px)}")

                plot!(plot_cost_log, iterList, costList, label=label, yscale=:log10)
                xlabel!(plot_cost_log, L"\mathrm{Iterations}")
                ylabel!(plot_cost_log, L"\mathrm{Cost\;(px)}")

                plot!(plot_error, time, h_error, label=label)
                xlabel!(plot_error, L"\mathrm{Time\;(s)}")
                ylabel!(plot_error, L"\mathrm{Height\;Error\;(mm)}")

                plot!(plot_h, time, est_h, label=label)
                xlabel!(plot_h, L"\mathrm{Time\;(s)}")
                ylabel!(plot_h, L"\mathrm{Height\;(mm)}")

                push!(cost_list_list,costList)
                push!(est_η_list,est_η)
                push!(est_β_list,est_β)
                push!(h_error_list,h_error)
                push!(iter_list,iterList)
            end
        end

        ObsDataList, splinex_gt, spliney_gt = read_csv(joinpath(filepath_gt,"data","sim_data","2D_border_points"))

        res_path = joinpath(filepath,"runs","post_analysis","ne_4")
        set_file(res_path)

        animate_fields(filepath=joinpath(res_path,"plots"), p=x[1], q=y[1], pObs=x[2], qObs=y[2], pgt= splinex_gt, qgt=spliney_gt)

        Plots.hline!(plot_η, [η_gt], label="Ground truth η")
        savefig(plot_η, joinpath(res_path,"eta.pdf"))

        Plots.hline!(plot_β, [β_gt], label="Ground truth β")
        savefig(plot_β, joinpath(res_path,"beta.pdf"))

        savefig(plot_cost, joinpath(res_path,"convergence.pdf"))
        savefig(plot_cost_log, joinpath(res_path,"convergence_log.pdf"))
        savefig(plot_error, joinpath(res_path,"error.pdf"))

        plot!(plot_h, time, gt_h, label=L"h_{gt}(t)")
        savefig(plot_h, joinpath(res_path,"h.pdf"))

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

        plot_η = set_plot(fs, sz=(plt_width, plt_height))
        plot_β = set_plot(fs, sz=(plt_width, plt_height))
        plot_cost = set_plot(fs, sz=(plt_width, plt_height))
        plot_cost_log = set_plot(fs, sz=(plt_width, plt_height))
        plot_error = set_plot(fs, sz=(plt_width, plt_height))
        plot_h = set_plot(fs, sz=(plt_width, plt_height))

        sim_time = 10.0 # simulation time in seconds 
        steps = 25.0 # number of time steps
        t_steps = sim_time/steps

        time = collect(Float64, range(start=0, stop=sim_time, step=t_steps))

        for file in path_3
            exp_path = joinpath(filepath,"runs",file)
            sim_params = read_json(joinpath(exp_path,"Results","data","sim_params.json")) 

            viscosity_type = sim_params["viscosity_type"]

            if viscosity_type == "constant"
                η_gt = sim_params["η_gt"]
                β_gt = sim_params["β_gt"]

                est_η = readdlm(joinpath(exp_path,"Results","data","η.csv"), ',', Float64)
                est_β = readdlm(joinpath(exp_path,"Results","data","β.csv"), ',', Float64)
                est_h = readdlm(joinpath(exp_path,"Results","data","est_h.csv"), ',', Float64)
                gt_h = readdlm(joinpath(exp_path,"Results","data","gt_h.csv"), ',', Float64)

                h_error = abs.(est_h-gt_h)

                stats = read_json(joinpath(exp_path,"Results","data","stats.json")) 

                costList = stats["cost_list"]
                iterList = stats["iterList"]

                label = " "
                if file[1:2] == "Q2"
                    label = "Lagrange basis"
                elseif file[1:2] == "S2"
                    label = "NURBS basis"
                end
                plot!(plot_η, iterList, est_η, label=label)
                xlabel!(plot_η, L"\mathrm{Iterations}")
                ylabel!(plot_η, L"\eta\;\mathrm{(KPa\cdot s)}")

                plot!(plot_β, iterList, est_β, label=label)
                xlabel!(plot_β, L"\mathrm{Iterations}")
                ylabel!(plot_β,L"\beta\;\mathrm{(mm^{-1})}")

                plot!(plot_cost, iterList, costList, label=label)
                xlabel!(plot_cost, L"\mathrm{Iterations}")
                ylabel!(plot_cost, L"\mathrm{Cost\;(px)}")

                plot!(plot_cost_log, iterList, costList, label=label, yscale=:log10)
                xlabel!(plot_cost_log, L"\mathrm{Iterations}")
                ylabel!(plot_cost_log, L"\mathrm{Cost\;(px)}")

                plot!(plot_error, time, h_error, label=label)
                xlabel!(plot_error, L"\mathrm{Time\;(s)}")
                ylabel!(plot_error, "Height Error (mm)")

                plot!(plot_h, time, est_h, label=label)
                xlabel!(plot_h, L"\mathrm{Time\;(s)}")
                ylabel!(plot_h, L"\mathrm{Height\;(mm)}")

                push!(cost_list_list,costList)
                push!(est_η_list,est_η)
                push!(est_β_list,est_β)
                push!(h_error_list,h_error) 
                push!(iter_list,iterList)
            end
        end

        res_path = joinpath(filepath,"runs","post_analysis","ne_8")
        set_file(res_path)

        Plots.hline!(plot_η, [η_gt], label="Ground truth η")
        savefig(plot_η, joinpath(res_path,"eta.pdf"))

        Plots.hline!(plot_β, [β_gt], label="Ground truth β")
        savefig(plot_β, joinpath(res_path,"beta.pdf"))

        savefig(plot_cost, joinpath(res_path,"convergence.pdf"))
        savefig(plot_cost_log, joinpath(res_path,"convergence_log.pdf"))
        savefig(plot_error, joinpath(res_path,"error.pdf"))

        plot!(plot_h, time, gt_h, label=L"h_{gt}(t)")
        savefig(plot_h, joinpath(res_path,"h.pdf"))

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

        plot_η = set_plot(fs, sz=(plt_width, plt_height))
        plot_β = set_plot(fs, sz=(plt_width, plt_height))
        plot_cost = set_plot(fs, sz=(plt_width, plt_height))
        plot_cost_log = set_plot(fs, sz=(plt_width, plt_height))
        plot_error = set_plot(fs, sz=(plt_width, plt_height))
        plot_h = set_plot(fs, sz=(plt_width, plt_height))

        sim_time = 10.0 # simulation time in seconds 
        steps = 25.0 # number of time steps
        t_steps = sim_time/steps

        time = collect(Float64, range(start=0, stop=sim_time, step=t_steps))

        for file in path_4
            exp_path = joinpath(filepath,"runs",file)
            sim_params = read_json(joinpath(exp_path,"Results","data","sim_params.json")) 

            viscosity_type = sim_params["viscosity_type"]

            if viscosity_type == "constant"
                η_gt = sim_params["η_gt"]
                β_gt = sim_params["β_gt"]

                est_η = readdlm(joinpath(exp_path,"Results","data","η.csv"), ',', Float64)
                est_β = readdlm(joinpath(exp_path,"Results","data","β.csv"), ',', Float64)
                est_h = readdlm(joinpath(exp_path,"Results","data","est_h.csv"), ',', Float64)
                gt_h = readdlm(joinpath(exp_path,"Results","data","gt_h.csv"), ',', Float64)

                h_error = abs.(est_h-gt_h)

                stats = read_json(joinpath(exp_path,"Results","data","stats.json")) 

                costList = stats["cost_list"]
                iterList = stats["iterList"]

                label = " "
                if file[1:2] == "Q2"
                    label = "Lagrange basis"
                elseif file[1:2] == "S2"
                    label = "NURBS basis"
                end
                plot!(plot_η, iterList, est_η, label=label)
                xlabel!(plot_η, L"\mathrm{Iterations}")
                ylabel!(plot_η, L"\eta\;\mathrm{(KPa\cdot s)}")

                plot!(plot_β, iterList, est_β, label=label)
                xlabel!(plot_β, L"\mathrm{Iterations}")
                ylabel!(plot_β,L"\beta\;\mathrm{(mm^{-1})}")

                plot!(plot_cost, iterList, costList, label=label)
                xlabel!(plot_cost, L"\mathrm{Iterations}")
                ylabel!(plot_cost, L"\mathrm{Cost\;(px)}")

                plot!(plot_cost_log, iterList, costList, label=label, yscale=:log10)
                xlabel!(plot_cost_log, L"\mathrm{Iterations}")
                ylabel!(plot_cost_log, L"\mathrm{Cost\;(px)}")

                plot!(plot_error, time, h_error, label=label)
                xlabel!(plot_error, L"\mathrm{Time\;(s)}")
                ylabel!(plot_error, "Height Error (mm)")

                plot!(plot_h, time, est_h, label=label)
                xlabel!(plot_h, L"\mathrm{Time\;(s)}")
                ylabel!(plot_h, L"\mathrm{Height\;(mm)}")

                push!(cost_list_list,costList)
                push!(est_η_list,est_η)
                push!(est_β_list,est_β)
                push!(h_error_list,h_error)
                push!(iter_list,iterList)
            end
        end

        res_path = joinpath(filepath,"runs","post_analysis","ne_8")
        set_file(res_path)

        Plots.hline!(plot_η, [η_gt], label="Ground truth η")
        savefig(plot_η, joinpath(res_path,"eta.pdf"))

        Plots.hline!(plot_β, [β_gt], label="Ground truth β")
        savefig(plot_β, joinpath(res_path,"beta.pdf"))

        savefig(plot_cost, joinpath(res_path,"convergence.pdf"))
        savefig(plot_cost_log, joinpath(res_path,"convergence_log.pdf"))
        savefig(plot_error, joinpath(res_path,"error.pdf"))

        plot!(plot_h, time, gt_h, label=L"h_{gt}(t)")
        savefig(plot_h, joinpath(res_path,"h.pdf"))

    end

    if length(path_5) != 0
        cost_list_list = []
        est_η_list = []
        est_β_list = []
        h_error_list = []
        iter_list = []
        η_gt = 0
        β_gt = 0

        plot_η = set_plot(fs, sz=(plt_width, plt_height))
        plot_β = set_plot(fs, sz=(plt_width, plt_height))
        plot_cost = set_plot(fs, sz=(plt_width, plt_height))
        plot_error = set_plot(fs, sz=(plt_width, plt_height))

        for file in path_5
            exp_path = joinpath(filepath,"runs",file)
            sim_params = read_json(joinpath(exp_path,"Results","data","sim_params.json")) 

            sim_time = 10.0 # simulation time in seconds 
            steps = 25.0 # number of time steps
            t_steps = sim_time/steps

            time = collect(Float64, range(start=0, stop=sim_time, step=t_steps))

            viscosity_type = sim_params["viscosity_type"]

            if viscosity_type == "constant"
                η_gt = sim_params["η_gt"]
                β_gt = sim_params["β_gt"]

                est_η = readdlm(joinpath(exp_path,"Results","data","η.csv"), ',', Float64)
                est_β = readdlm(joinpath(exp_path,"Results","data","β.csv"), ',', Float64)
                est_h = readdlm(joinpath(exp_path,"Results","data","est_h.csv"), ',', Float64)
                gt_h = readdlm(joinpath(exp_path,"Results","data","gt_h.csv"), ',', Float64)

                h_error = abs.(est_h-gt_h)

                stats = read_json(joinpath(exp_path,"Results","data","stats.json")) 

                costList = stats["cost_list"]
                iterList = stats["iterList"]

                label = " "
                if file[1:2] == "Q2"
                    label = "Lagrange basis"
                elseif file[1:2] == "S2"
                    label = "NURBS basis"
                end

                plot!(plot_η, iterList, est_η, label=label)
                xlabel!(plot_η, L"\mathrm{Iterations}")
                ylabel!(plot_η, L"\eta\;\mathrm{(KPa\cdot s)}")

                plot!(plot_β, iterList, est_β, label=label)
                xlabel!(plot_β, L"\mathrm{Iterations}")
                ylabel!(plot_β, L"\beta\;\mathrm{(mm^{-1})}")

                plot!(plot_cost, iterList, costList, label=label)
                xlabel!(plot_cost, L"\mathrm{Iterations}")
                ylabel!(plot_cost, L"\mathrm{Cost\;(px)}")

                plot!(plot_cost_log, iterList, costList, label=label, yscale=:log10)
                xlabel!(plot_cost_log, L"\mathrm{Iterations}")
                ylabel!(plot_cost_log, L"\mathrm{Cost\;(px)}")

                plot!(plot_error, time, h_error, label=label)
                xlabel!(plot_error, L"\mathrm{Time\;(s)}")
                ylabel!(plot_error, L"\mathrm{Height\;Error\;(mm)}")

                plot!(plot_h, time, est_h, label=label)
                xlabel!(plot_h, L"\mathrm{Time\;(s)}")
                ylabel!(plot_h, L"\mathrm{Height\;(mm)}")

                push!(cost_list_list,costList)
                push!(est_η_list,est_η)
                push!(est_β_list,est_β)
                push!(h_error_list,h_error)
                push!(iter_list,iterList)
            end
        end

        res_path = joinpath(filepath,"runs","post_analysis","ne_2")
        set_file(res_path)

        Plots.hline!(plot_η, [η_gt], label="Ground truth η")
        savefig(plot_η, joinpath(res_path,"eta.pdf"))

        Plots.hline!(plot_β, [β_gt], label="Ground truth β")
        savefig(plot_β, joinpath(res_path,"beta.pdf"))

        savefig(plot_cost, joinpath(res_path,"convergence.pdf"))
        savefig(plot_cost_log, joinpath(res_path,"convergence_log.pdf"))
        savefig(plot_error, joinpath(res_path,"error.pdf"))

        plot!(plot_h, time, gt_h, label=L"h_{gt}(t)")
        savefig(plot_h, joinpath(res_path,"h.pdf"))
    end
end

function replot(filepath, filepath_gt)
    
    elem_size_folders = readdir(joinpath(filepath))
    
    sim_params = read_json(joinpath(filepath_gt,"data","sim_params.json")) 

    r::Float64 = sim_params["r"]
    h::Float64 = sim_params["h"]
    
    viscosity_type::String = sim_params["viscosity_type"]
    F = Array(float.(sim_params["cParam"]))

    sim_time::Float64 = sim_params["simulation_time"]
    t_steps::Float64 = sim_params["time_steps"]
    
    camera_matrix::AbstractArray = reshape(Array(sim_params["camera_matrix"]), 3, 3)
    obj_pose::AbstractArray = reshape(Array(sim_params["obj_pose"])*1.0,4,4)
    
    control::String = sim_params["control_type"]

    viscosity_model::String = ""
    if viscosity_type == "bulk_viscosity" && haskey(sim_params, "model_type") && sim_params["model_type"] == "carreau"
            viscosity_model = sim_params["model_type"]
    end

    ndim::Int = 3
    nDof_p::Int = 1  # number of degree of freedom per node
    nDof_u::Int = ndim  # number of degree of freedom per node

    gt_h_::Matrix{Float64} = Matrix{Float64}(undef,0,0)

     # Read ground truth height data

    for elem_size_folder in elem_size_folders
        if elem_size_folder == "post_analysis"
            continue
        end
        
        sim_time_folders = readdir(joinpath(filepath, elem_size_folder))

        # Iterate over simulation time folders
        for sim_time_folder in sim_time_folders
            if sim_time_folder == "post_analysis_time" || sim_time_folder == "Results"  || sim_time_folder == "simtime_2.0"  || sim_time_folder == "simtime_20.0" || sim_time_folder == "simtime_10.0" || sim_time_folder == "simtime_30.0"
                continue
            end

            noise_folders = readdir(joinpath(filepath, elem_size_folder, sim_time_folder))
            for noise_folder in noise_folders
                if noise_folder == "post_analysis_noise" || noise_folder == "Results"
                    continue
                end
                exp_path = joinpath(filepath, elem_size_folder, sim_time_folder, noise_folder)
                
                if viscosity_type == "constant"
                    println("Comparing experiments in: $exp_path")
                    
                    η_gt = float(sim_params["η"][1])
                    β_gt = float(sim_params["β"][1])

                    exp_params = read_json(joinpath(exp_path,"Results","data","experiment_parameters.json"))

                    noise_level = exp_params["noise_level"]
                    FunctionClass_p = exp_params["FunctionClass_p"]
                    FunctionClass_u = exp_params["FunctionClass_u"]
                    FunctionClass_x = exp_params["FunctionClass_x"]

                    sim_time_exp = exp_params["sim_time_exp"]   
                    ne_exp = exp_params["ne_exp"]
                    data_type = exp_params["data_type"]
                    control = exp_params["control"]

                    if data_type != "physical"
                        gt_h_ = readdlm(joinpath(filepath_gt,"data","h.csv"), ',', Float64)
                    end
                    if data_type  == "synthetic"
                        ObsDataList, splinexObs, splineyObs = read_csv(joinpath(filepath_gt,"data","img_data","contour_data"))  
                        @info "Data type $data_type Reading synthetic contour data of $(length(ObsDataList)) time steps"
                    elseif data_type == "simulated"
                        ObsDataList, splinexObs, splineyObs = read_csv(joinpath(filepath_gt,"data","sim_data","contour_data"))  
                        @info "Data type $data_type Reading simulated contour data from $(joinpath(filepath_gt,"data","sim_data","contour_data")) of $(length(ObsDataList)) time steps"
                    else
                        error("Unknown data type: $data_type")
                    end

                    time = collect(Float64, range(start=0, stop=sim_time, step=t_steps))
                    gt_h = gt_h_[1:(round(Int, sim_time/ t_steps))+1]

                    obsBorderPts, symBorderPts, nSplinex, nSpliney, splinex, spliney = _get_borders(data_type, filepath_gt, exp_path)
                    conditions = Conditions(camera_matrix=camera_matrix, obj_pose=obj_pose)
                    if noise_level == 0.0

                        est_η = readdlm(joinpath(exp_path,"Results","data","η.csv"), ',', Float64)
                        est_β = readdlm(joinpath(exp_path,"Results","data","β.csv"), ',', Float64)
                        est_h = readdlm(joinpath(exp_path,"Results","data","est_h.csv"), ',', Float64)
                        cost_iter = readdlm(joinpath(exp_path,"Results","data","cost_iter.csv"), ',', Float64)
                        
                        stats = read_json(joinpath(exp_path,"Results","data","stats.json")) 
                        η = stats["η"]
                        β = stats["β"]
                        
                        d, pairs = closest_point(symBorderPts, obsBorderPts)

                        if sim_time_exp == 5.0
                            cont_y_min = 350
                        elseif sim_time_exp == 10.0
                            cont_y_min = 400
                        elseif sim_time_exp == 20.0
                            cont_y_min = 420
                        else
                            cont_y_min = 480
                        end
                        
                        contour_plt = set_plot(fs, sz=(plt_width, plt_height))
                        Plots.plot!(contour_plt, [], label=false, legend=:outerbottom, legend_column=2, aspect_ratio = :equal)
                        # Plots.scatter!(contour_plt, symBorderPts[end][1,:], symBorderPts[end][2,:], label="Simulated contour", ms=:10, mc=:royalblue, ma=:0.7, markerstrokewidth=0.2)
                        Plots.plot!(contour_plt, nSplinex[end], nSpliney[end], label="Ground truth contour", color=:red)
                        Plots.yflip!(true)
                        Plots.xlims!(contour_plt, 480, 1520)
                        Plots.ylims!(contour_plt, cont_y_min, 1200)
                        Plots.xlabel!(contour_plt, L"x\;\mathrm{(px)}")
                        Plots.ylabel!(contour_plt, L"y\;\mathrm{(px)}")
                        Plots.savefig(contour_plt, joinpath(exp_path,"Results","plots","contour_comparison.pdf"))

                        contour_plt_zoom = set_plot(fs, sz=(1200, plt_height))
                        Plots.plot!(contour_plt_zoom, [], label=false, legend=:outerbottom, legend_column=2, aspect_ratio = :equal)
                        # Plots.scatter!(contour_plt_zoom, symBorderPts[end][1,:], symBorderPts[end][2,:], label="Simulated contour", ms=:10, mc=:royalblue, ma=:0.7, markerstrokewidth=0.2)
                        Plots.plot!(contour_plt_zoom, nSplinex[end], nSpliney[end], label=false, color=:red)
                        Plots.xlims!(contour_plt_zoom, 1000, 1520)
                        Plots.ylims!(contour_plt_zoom, cont_y_min, 1200)
                        Plots.xticks!(contour_plt_zoom, 1100:200:1520)
                        Plots.yflip!(true)
                        Plots.xlabel!(contour_plt_zoom, L"x\;\mathrm{(px)}")
                        Plots.ylabel!(contour_plt_zoom, L"y\;\mathrm{(px)}")
                        Plots.savefig(contour_plt_zoom, joinpath(exp_path,"Results","plots","contour_comparison_zoom.pdf"))

                        # animate_fields(filepath=joinpath(exp_path,"Results","plots"), BorderNodes2D=symBorderPts, pObs=nSplinex, qObs=nSpliney)
                        # animate_fields(filepath=joinpath(exp_path,"Results","plots"), pObs=nSplinex, qObs=nSpliney)
                    
                        costList = stats["cost_list"]
                        iterList = stats["iterList"]

                        plt_cnt_error = set_plot(fs, sz=(plt_width, plt_height))
                        Plots.plot!(plt_cnt_error, time[1:length(d)], d, label="Closest point distance error", legend=:outerbottom, legend_column=2)
                        Plots.xlabel!(plt_cnt_error, L"\mathrm{Time\;(s)}")
                        Plots.ylabel!(plt_cnt_error, L"\mathrm{Closest\;Point\;Distance\;(px)}")
                        Plots.xlims!(plt_cnt_error, 0, end_obs_win)
                        Plots.savefig(plt_cnt_error, joinpath(exp_path,"Results","plots","closest_point_distance_error.pdf"))
                        
                        # Plot the estimated and ground truth height with inset zoom
                        h_plt = set_plot(fs, sz=(plt_width, plt_height))
                        est_h = vec(Float64.(collect(est_h)))
                        gt_h = vec(Float64.(collect(gt_h)))
                        
                        n_time = min(length(time), length(est_h), length(gt_h)) #, (end_obs_win/t_steps + 1))
                        if n_time < length(time) || n_time < length(est_h) || n_time < length(gt_h)
                            @warn "Time and height vectors have mismatched lengths: time=$(length(time)), est_h=$(length(est_h)), gt_h=$(length(gt_h)). Truncating to $n_time samples for plotting."
                        end

                        Plots.plot!(h_plt, time[1:n_time], est_h[1:n_time], label=L"h_{est}(t)", legend=:outerbottom, legend_column=2)
                        Plots.plot!(h_plt, time[1:n_time], gt_h[1:n_time], label=L"h_{gt}(t)", legend=:outerbottom, legend_column=2)
                        Plots.xlabel!(h_plt, L"\mathrm{Time\;(s)}")
                        Plots.ylabel!(h_plt, L"\mathrm{Height\;(mm)}")
                        Plots.xlims!(h_plt, 0, end_obs_win)
                        Plots.savefig(h_plt, joinpath(exp_path,"Results","plots","h_est.pdf"))

                        # Plot the height estimation error
                        error_plt = set_plot(fs, sz=(plt_width, plt_height))
                        # ensure vectors are aligned for error plot
                        Plots.plot!(error_plt, time[1:n_time], abs.(est_h[1:n_time] .- gt_h[1:n_time]), label="Height estimation error", legend=:outerbottom, legend_column=1)
                        Plots.xlabel!(error_plt, L"\mathrm{Time\;(s)}")
                        Plots.ylabel!(error_plt, "Height Error (mm)")
                        Plots.xlims!(error_plt, 0, end_obs_win)
                        Plots.savefig(error_plt, joinpath(exp_path,"Results","plots","h_est_error.pdf"))

                        # Plot the estimated and ground truth parameters
                        η_plt = set_plot(fs, sz=(plt_width, plt_height))
                        Plots.plot!(η_plt, est_η, label="Estimated η", marker=1, legend=:outerbottom, legend_column=2)
                        Plots.hline!([η_gt], label="Ground truth η", legend=:outerbottom, legend_column=2)
                        Plots.xlabel!(η_plt, L"\mathrm{Iterations}")
                        Plots.ylabel!(η_plt, L"\eta\;\mathrm{(KPa\cdot s)}")
                        Plots.savefig(η_plt, joinpath(exp_path,"Results","plots","η.pdf"))
                        
                        β_plt = set_plot(fs, sz=(plt_width, plt_height))
                        Plots.plot!(β_plt, est_β, label="Estimated β", marker=1, legend=:outerbottom, legend_column=2)
                        Plots.hline!(β_plt, [β_gt], label="Ground truth β",legend=:outerbottom, legend_column=2)
                        Plots.xlabel!(β_plt, L"\mathrm{Iterations}")
                        Plots.ylabel!(β_plt, L"\beta\;\mathrm{(mm^{-1})}")
                        Plots.savefig(β_plt, joinpath(exp_path,"Results","plots","β.pdf"))

                        # Plot the cost function with iterations
                        cost_plt = set_plot(fs, sz=(plt_width, plt_height))
                        Plots.plot!(cost_plt, iterList, costList, label="Cost", marker=1, yscale=:log10, xminorgrid = :false,legend=:outerbottom, legend_column=1)
                        Plots.xlabel!(cost_plt, L"\mathrm{Iterations}")
                        Plots.ylabel!(cost_plt, L"\mathrm{Cost\;(px)}")
                        # Plots.xticks!(minimum(iterList):2:maximum(iterList))
                        Plots.savefig(cost_plt, joinpath(exp_path,"Results","plots","cost_steps.pdf"))

                        # Plot the cost function with iterations
                        cost_log_plt = set_plot(fs, sz=(plt_width, plt_height))
                        Plots.plot!(cost_log_plt, iterList, costList, label="Cost", marker=1, yscale=:log10, xscale=:log10, legend=:outerbottom, legend_column=1)
                        Plots.xlabel!(cost_log_plt, L"\mathrm{Iterations}")
                        Plots.ylabel!(cost_log_plt, L"\mathrm{Cost\;(px)}")
                        Plots.savefig(cost_log_plt, joinpath(exp_path,"Results","plots","cost_steps_log.pdf"))

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
                        
                        try 
                            contour_plot_params = read_json(joinpath(exp_path,"Results","data","contour_plot_params.json")) 
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
                                # Compute adaptive clims from CostMat to represent all values
                                cost_min, cost_max = extrema(CostMat)
                                cost_max = min(cost_max, 1e3) # avoid zero max
                                cost_clims = (0, cost_max)
                                plt_dirs = set_plot(fs, sz=(plt_width, plt_height))
                                Plots.contour!(plt_dirs, ηList, βList, CostMat', color=:turbo, fill=false, levels=100, legend=:outerbottom, legend_column=3, clims=cost_clims)
                                Plots.plot!(plt_dirs, est_η, est_β, label="Estimations", ms=:4, m=:x, color=:red, legend=:outerbottom, legend_column=3, markerstrokewidth=0.1)
                                Plots.plot!(plt_dirs, etas_steep, betas_steep, label = "Steepest dir", color=:black, legend=:outerbottom, legend_column=3)
                                Plots.plot!(plt_dirs, etas_flat,  betas_flat,  label = "Flattest dir",  lw=1, color=:magenta, legend=:outerbottom, legend_column=3)
                                Plots.scatter!(plt_dirs, [η0], [β0], label="Minimum Cost", ms=15, m=:star5, color=:black, legend=:outerbottom, legend_column=3, markerstrokewidth=0.1)
                                Plots.scatter!(plt_dirs, [η_gt], [β_gt], label="Ground truth", ms=:15, m=:star5, color=def_orange, legend=:outerbottom, legend_column=3, markerstrokewidth=0.1)
                                Plots.xlabel!(plt_dirs, L"\eta\;\mathrm{(KPa\cdot s)}")
                                Plots.ylabel!(plt_dirs, L"\beta\;\mathrm{(mm^{-1})}")
                                Plots.savefig(plt_dirs, joinpath(exp_path, "Results", "plots", "cost_surface_with_directions.pdf"))
                                
                                plt_cont = set_plot(fs, sz=(plt_width, plt_height))
                                Plots.contour!(plt_cont, ηList, βList, CostMat', color=:turbo, fill=false, levels=100, legend=:outerbottom, legend_column=3, clims=cost_clims)
                                Plots.plot!(plt_cont, est_η, est_β, label="Estimations", ms=:4, m=:x, color=:red, legend=:outerbottom, legend_column=3, markerstrokewidth=0.1)
                                Plots.scatter!(plt_cont, [η0], [β0], label="Minimum Cost", ms=15, m=:star5, color=:black, legend=:outerbottom, legend_column=3, markerstrokewidth=0.1)
                                Plots.scatter!(plt_cont, [η_gt], [β_gt], label="Ground truth", ms=:15, m=:star5, color=def_orange, legend=:outerbottom, legend_column=3, markerstrokewidth=0.1)
                                Plots.xlabel!(plt_cont, L"\eta\;\mathrm{(KPa\cdot s)}")
                                Plots.ylabel!(plt_cont, L"\beta\;\mathrm{(mm^{-1})}")
                                Plots.savefig(plt_cont, joinpath(exp_path, "Results", "plots", "cost_surface.pdf"))

                                plt_cont = set_plot(fs, sz=(plt_width, plt_height))
                                Plots.contourf!(plt_cont, ηList, βList, CostMat', color=:turbo, fill=false, levels=100, legend=:outerbottom, legend_column=3, clims=cost_clims)
                                Plots.plot!(plt_cont, est_η, est_β, label="Estimations", ms=:4, m=:x, color=:red, legend=:outerbottom, legend_column=3, markerstrokewidth=0.1)
                                Plots.scatter!(plt_cont, [η0], [β0], label="Minimum Cost", ms=15, m=:star5, color=:black, legend=:outerbottom, legend_column=3, markerstrokewidth=0.1)
                                Plots.scatter!(plt_cont, [η_gt], [β_gt], label="Ground truth", ms=15, m=:star5, color=def_orange, legend=:outerbottom, legend_column=3, markerstrokewidth=0.1)
                                Plots.xlabel!(plt_cont, L"\eta\;\mathrm{(KPa\cdot s)}")
                                Plots.ylabel!(plt_cont, L"\beta\;\mathrm{(mm^{-1})}")
                                Plots.savefig(plt_cont, joinpath(exp_path, "Results", "plots", "cost_surface_iter.pdf"))

                                # 2D slices: cost vs distance along the two directions
                                plt_slices = set_plot(fs, sz=(plt_width, plt_height))
                                if length(t_steep) > 0 && length(zs_steep) == length(t_steep)
                                    Plots.plot!(plt_slices, t_steep, zs_steep, label = "Steepest direction", color=:black, legend=:outerbottom, legend_column=2)
                                else
                                    @warn "Skipping steep slice plot: empty or mismatched lengths: $(length(t_steep)) vs $(length(zs_steep))"
                                end
                                if length(t_flat) > 0 && length(zs_flat) == length(t_flat)
                                    Plots.plot!(plt_slices, t_flat,  zs_flat,  label = "Flattest direction",  lw=1, color=:gray, legend=:outerbottom, legend_column=2)
                                else
                                    @warn "Skipping flat slice plot: empty or mismatched lengths: $(length(t_flat)) vs $(length(zs_flat))"
                                end
                                Plots.vline!(plt_slices, [0.0], color=:blue, linestyle=:dash, label="Minimum")
                                Plots.xlabel!(plt_slices, L"\mathrm{Distance\;from\;minimum\;(px)}")
                                Plots.ylabel!(plt_slices, L"\mathrm{Cost}")
                                Plots.savefig(plt_slices, joinpath(exp_path, "Results", "plots", "cost_slices_along_directions.pdf"))

                                Plots.ylims!(plt_slices, 0, 10)
                                Plots.savefig(plt_slices, joinpath(exp_path, "Results", "plots", "cost_slices_along_directions_zoomed.pdf"))
                        
                                # save analysis metadata
                                dir_info = Dict("η_min"=>η0, "β_min"=>β0, "Hessian"=>H, "evals"=>evals, "v_steep"=>v_steep, "v_flat"=>v_flat)
                                slice_data = Dict(
                                    "steep"=>Dict("t"=>t_steep, "etas"=>etas_steep, "betas"=>betas_steep, "zs"=>zs_steep),
                                    "flat"=>Dict("t"=>t_flat, "etas"=>etas_flat, "betas"=>betas_flat, "zs"=>zs_flat)
                                )
                                write_json(joinpath(exp_path, "Results", "data", "direction_analysis"), dir_info)
                                write_json(joinpath(exp_path, "Results", "data", "slice_data"), slice_data)

                            catch err
                                @warn "Post-process Hessian/direction analysis failed: $err"
                            end

                            # Plot the cost function surface: prefer PlotlyJS HTML export to avoid
                            # creating a GLMakie Figure (which may require a compatible WebIO setup).
                            # set_plot(fs, sz=(plt_width, plt_height))
                            # try

                            #     htmlpath = joinpath(exp_path, "Results", "plots", "cost_surface_interactive.html")
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
                            #     plt = set_plot(fs, sz=(plt_width, plt_height))
                            #     Plots.contour!(plt, ηList, βList, CostMat', color=:turbo, fill=false, levels=100, legend=:outerbottom, legend_column=3)
                            #     Plots.plot!(plt, est_η, est_β, label="Estimations", ms=:4, m=:x, color=:red, legend=:outerbottom, legend_column=3)
                            #     Plots.scatter!(plt, [η_gt], [β_gt], label="Ground truth", ms=:8, m=:circle, color=def_orange, legend=:outerbottom, legend_column=3)
                            #     Plots.scatter!(plt_dirs, [η0], [β0], label="Minimum Cost", ms=15, m=:star5, color=:black, legend=:outerbottom, legend_column=3)
                            #     Plots.xlabel!(plt, L"\eta\;\mathrm{(KPa\cdot s)}")
                            #     Plots.ylabel!(plt, L"\beta\;\mathrm{(mm^{-1})}")
                            #     Plots.savefig(plt, joinpath(exp_path,"Results","plots","cost_surface_iter.pdf"))

                            #     plt2 = set_plot(fs, sz=(plt_width, plt_height))
                            #     Plots.contourf!(plt2, ηList, βList, CostMat', color=:turbo, fill=false, levels=100, legend=:outerbottom, legend_column=3)
                            #     Plots.plot!(plt2, est_η, est_β, label="Estimations", ms=:4, m=:x, color=:red, legend=:outerbottom, legend_column=3)
                            #     Plots.scatter!(plt2, [η_gt], [β_gt], label="Ground truth", ms=:8, m=:circle, color=def_orange, legend=:outerbottom, legend_column=3)
                            #     Plots.xlabel!(plt2, L"\eta\;\mathrm{(KPa\cdot s)}")
                            #     Plots.ylabel!(plt2, L"\beta\;\mathrm{(mm^{-1})}")
                            #     Plots.savefig(plt2, joinpath(exp_path,"Results","plots","cost_surface.pdf"))
                            # catch err
                            #     @warn "Failed to produce static PDF contour outputs: $err"
                            # end
                    catch err
                        @warn "Failed to post-process contour parameters and cost surface: $err"
                    end
                    else
                        # try
                            η_pred = readdlm(joinpath(exp_path,"Results","data","eta_est.csv"), ',', Float64) # estimated η values per sample [n x n_iter]
                            β_pred = readdlm(joinpath(exp_path,"Results","data","beta_est.csv"), ',', Float64) # estimated β values per sample [n x n_iter]
                            h_pred = readdlm(joinpath(exp_path,"Results","data","h_est.csv"), ',', Float64) # estimated height values per sample [n x sim_time]

                            exp_params = read_json(joinpath(exp_path,"Results","data","experiment_parameters.json"))
                            
                            sim_time_exp = exp_params["sim_time_exp"]
                            obsBorderPts, nSplinex, nSpliney, pd = add_noise(ObsDataList, nFactor=noise_level)
                            # symBorderPts, splinex, spliney = read_csv(joinpath(exp_path,"Results","data","sim_data","2D_border_points"))
                            
                            # if sim_time_exp != sim_time
                            #     @warn "Simulation time in experiment parameters ($sim_time_exp) does not match that in ground truth ($sim_time)."
                            #     @warn "Resimulating to match ground truth sim time."
                            #     n_samples = size(η_pred, 1)
                            #     h_est_lst = zeros(n_samples, Int(sim_time/ t_steps)+1)
                            #     for n in 1:n_samples
                            #         η = η_pred[n]
                            #         β = β_pred[n]
                                    
                            #         @info "Resimulating for sample $n with η=$η, β=$β"
                            #         est_model, est_scene = def_problem(r, h, ne_exp, η, ndim, FunctionClass_u, nDof_u, FunctionClass_p, nDof_p, FunctionClass_x, β, F, control, viscosity_type, 
                            #         sim_time, t_steps)
                            #         est_μ_list, gradList, borderPts2DList, fields, pos3D, pos2D, splinex, spliney = simulate(est_model, est_scene, conditions)
                                    
                            #         est_h = get_height(est_μ_list, h)
                            #         h_est_lst[n, :] = est_h
                
                            #         write_csv(joinpath(exp_path,"Results","data","opt_data","est_height","run_$n"), est_h)

                            #         write_data(joinpath(exp_path,"Results","data","sim_data","2D_points","run_$n"), pos2D)
                            #         write_data(joinpath(exp_path,"Results","data","sim_data","3D_points","run_$n"), pos3D)
                            #         write_data(joinpath(exp_path,"Results","data","sim_data","motion_fields ","run_$n"), fields)
                            #         write_data(joinpath(exp_path,"Results","data","sim_data","2D_border_points","run_$n"), borderPts2DList)
                            #         write_data(joinpath(exp_path,"Results","data","sim_data","spline_p","run_$n"), splinex)
                            #         write_data(joinpath(exp_path,"Results","data","sim_data","spline_q","run_$n"), spliney)
                            #     end
                            #     write_csv(joinpath(exp_path,"Results","data","h_est"), h_est_lst)
                                # if sim_time_exp < sim_time
                                #     @warn "Truncating obsBorderPts to match sim_time_exp."
                                #     obsBorderPts = obsBorderPts[1:Int(sim_time_exp/ t_steps)+1, :]
                                #     # symBorderPts = symBorderPts[1:Int(sim_time_exp/ t_steps)+1, :]

                                #     nSplinex = nSplinex[1:Int(sim_time_exp/ t_steps)+1, :]
                                #     nSpliney = nSpliney[1:Int(sim_time_exp/ t_steps)+1, :]

                                #     # splinex = splinex[1:Int(sim_time_exp/ t_steps)+1, :]
                                #     # spliney = spliney[1:Int(sim_time_exp/ t_steps)+1, :]

                                # else
                                #     @warn " Truncating simBorderPts to match gt sim time."
                                #     # symBorderPts = symBorderPts[1:Int(sim_time/ t_steps)+1, :]
                                #     obsBorderPts = obsBorderPts[1:Int(sim_time/ t_steps)+1, :]

                                #     nSplinex = nSplinex[1:Int(sim_time/ t_steps)+1, :]
                                #     nSpliney = nSpliney[1:Int(sim_time/ t_steps)+1, :]

                                #     # splinex = splinex[1:Int(sim_time/ t_steps)+1, :]
                                #     # spliney = spliney[1:Int(sim_time/ t_steps)+1, :]
                                # end
                            # elseif length(obsBorderPts) != length(symBorderPts)
                            #     n_time = min(length(obsBorderPts), length(symBorderPts))
                            #     @warn "Mismatched number of time steps between observed border points ($(length(obsBorderPts))) and simulated border points ($(length(symBorderPts))). Truncating to $n_time time steps."
                            #     obsBorderPts = obsBorderPts[1:n_time, :]
                            #     # symBorderPts = symBorderPts[1:n_time, :]

                            #     # nSplinex = nSplinex[1:n_time, :]
                            #     # nSpliney = nSpliney[1:n_time, :]

                            #     splinex = splinex[1:n_time, :]
                            #     spliney = spliney[1:n_time, :]
                            # end
                            
                            # d, pairs = closest_point(symBorderPts, obsBorderPts)
                            
                            # contour_plt = set_plot(fs, sz=(plt_width, plt_height))
                            # Plots.plot!(contour_plt, [], label=false, legend=:outerbottom, legend_column=2, aspect_ratio = :equal)
                            # Plots.scatter!(contour_plt, symBorderPts[end][1,:], symBorderPts[end][2,:], label="Simulated contour", ms=:10, mc=:royalblue, ma=:0.7, markerstrokewidth=0.2)
                            # Plots.plot!(contour_plt, nSplinex[end], nSpliney[end], label="Ground truth contour", color=:red)
                            # Plots.yflip!(true)
                            # Plots.xlims!(contour_plt, 480, 1520)
                            # Plots.ylims!(contour_plt, 480, 1200)
                            # Plots.xlabel!(contour_plt, L"\mathrm{x\;(px)}")
                            # Plots.ylabel!(contour_plt, L"\mathrm{y\;(px)}")
                            # Plots.savefig(contour_plt, joinpath(exp_path,"Results","plots","contour_comparison.pdf"))

                            # contour_plt_zoom = set_plot(fs, sz=(350, 350))
                            # Plots.plot!(contour_plt_zoom, [], label=false, legend=:outerbottom, legend_column=2, aspect_ratio = :equal)
                            # Plots.scatter!(contour_plt_zoom, symBorderPts[end][1,:], symBorderPts[end][2,:], label="Simulated contour", ms=:10, mc=:royalblue, ma=:0.7, markerstrokewidth=0.2)
                            # Plots.plot!(contour_plt_zoom, nSplinex[end], nSpliney[end], label="Ground truth contour", color=:red)
                            # Plots.xlims!(contour_plt_zoom, 1000, 1520)
                            # Plots.ylims!(contour_plt_zoom, 480, 1200)
                            # Plots.yflip!(true)
                            # Plots.xlabel!(contour_plt_zoom, L"\mathrm{x\;(px)}")
                            # Plots.ylabel!(contour_plt_zoom, L"\mathrm{y\;(px)}")
                            # Plots.savefig(contour_plt_zoom, joinpath(exp_path,"Results","plots","contour_comparison_zoom.pdf"))

                            # animate_fields(filepath=joinpath(exp_path,"Results","plots"), pObs=nSplinex, qObs=nSpliney)
                        
                            # plt_cnt_error = set_plot(fs, sz=(plt_width, plt_height))
                            # Plots.plot!(plt_cnt_error, d, label="Closest point distance error", legend=:outerbottom, legend_column=2)
                            # Plots.xlabel!(plt_cnt_error, L"\mathrm{Time\;(s)}")
                            # Plots.ylabel!(plt_cnt_error, L"\mathrm{Closest\;Point\;Distance\;(px)}")
                            # Plots.savefig(plt_cnt_error, joinpath(exp_path,"Results","plots","closest_point_distance_error.pdf"))

                            covarience_plt = plot_covariance(η_pred[:,1], β_pred[:,1], legend_column=2, fs=fs)
                            Plots.xlabel!(covarience_plt, L"\eta\;\mathrm{(KPa\cdot s)}")
                            Plots.ylabel!(covarience_plt, L"\beta\;\mathrm{(mm^{-1})}")
                            Plots.savefig(covarience_plt, joinpath(exp_path,"Results","plots","covariance.pdf"))

                            n_time = min(length(time), size(h_pred, 2), length(gt_h)) #, (end_obs_win/t_steps + 1))
                            if n_time < length(time) || n_time < size(h_pred, 2) || n_time < length(gt_h)
                                @warn "Time and height vectors have mismatched lengths: time=$(length(time)), est_h=$(size(h_pred, 2)), gt_h=$(length(gt_h)). Truncating to $n_time samples for plotting."
                            end

                            h_plot = set_plot(fs, sz=(plt_width, plt_height))
                            Plots.plot!(h_plot, time[1:n_time], gt_h[1:n_time], label=L"h_{gt}(t)", legend=false)
                            StatsPlots.errorline!(h_plot, time[1:n_time], h_pred[:,1:n_time], label=L"h_{est}(t)", legend=false)
                            Plots.xlabel!(L"\mathrm{Time\;(s)}")
                            Plots.ylabel!(L"\mathrm{Height\;(mm)}")
                            Plots.xlims!(h_plot, 0, end_obs_win)
                            Plots.savefig(joinpath(exp_path,"Results","plots","h_est_noisy.pdf"))
                            
                            # Add inset subplot for zoomed region (10-20 seconds)
                            t_min, t_max = 5.0, 5.1
                            mask_zoom = (time .>= t_min) .& (time .<= t_max)
                            if sum(mask_zoom) > 0
                                idx_zoom = findall(mask_zoom)
                                # Build a single combined figure with two subplots and one shared legend
                                # Use a 2-row layout: top row 2 subplots, bottom row small height for a centered legend
                                plt_layout = @layout [a b; c{0.12h}]

                                # Build a single combined figure and draw into subplots 1,2 and 3
                                plt_combined = set_subplot(fs, sz=(3500, 350), layout=plt_layout)

                                # Full time series on left subplot (subplot=1)
                                Plots.plot!(plt_combined, time, gt_h, subplot=1, label="", color=:blue)
                                try
                                    StatsPlots.errorline!(plt_combined, time[1:n_time], h_pred[:,1:n_time], subplot=1, label="", color=def_orange)
                                catch
                                    StatsPlots.errorline!(plt_combined, time[1:n_time], h_pred[:,1:n_time]', subplot=1, label="", color=def_orange)
                                end
                                Plots.xlabel!(plt_combined, L"\mathrm{Time\;(s)}", subplot=1)
                                Plots.ylabel!(plt_combined, L"\mathrm{Height\;(mm)}", subplot=1)

                                # Zoomed region on right subplot (subplot=2)
                                Plots.plot!(plt_combined, time[idx_zoom], gt_h[idx_zoom], subplot=2, label="", color=:blue)
                                try
                                    StatsPlots.errorline!(plt_combined, time[idx_zoom], h_pred[:, idx_zoom], subplot=2, label="", color=def_orange)
                                catch
                                    StatsPlots.errorline!(plt_combined, time[idx_zoom], h_pred[idx_zoom], subplot=2, label="", color=def_orange)
                                end
                                Plots.xlabel!(plt_combined, L"\mathrm{Time\;(s)}", subplot=2)
                                Plots.xlims!(plt_combined, t_min, t_max, subplot=2)

                                # Bottom subplot (subplot=3) — draw a single centered annotation as deterministic legend
                                Plots.plot!(plt_combined, [], [], subplot=3, label=L"h_{gt}(t)", framestyle=:none, legend=:outerbottom, color=:blue, legend_column=2, background_color=:transparent)
                                Plots.plot!(plt_combined, [], [], subplot=3, label=L"h_{est}(t)", framestyle=:none, legend=:outerbottom, color=def_orange, background_color=:transparent)
                                Plots.xlims!(plt_combined, 0.0, 1.0, subplot=3)
                                Plots.ylims!(plt_combined, 0.0, 1.0, subplot=3)

                                # Save combined figure
                                Plots.savefig(plt_combined, joinpath(exp_path, "Results", "plots", "h_est_noisy_zoomed.pdf"))
                            end

                            error = abs.(h_pred[:,1:n_time]' .- gt_h[1:n_time]) ./ gt_h[1:n_time]
                            h_error_plot = set_plot(fs, sz=(plt_width, plt_height))
                            StatsPlots.errorline!(h_error_plot, time[1:n_time], error', label="Height estimation error", legend=:outerbottom, legend_column=1)
                            Plots.xlabel!(L"\mathrm{Time\;(s)}")
                            Plots.ylabel!(L"\mathrm{Height\;(mm)}")
                            Plots.xlims!(h_error_plot, 0, end_obs_win)
                            Plots.savefig(joinpath(exp_path,"Results","plots","h_rel_error_noisy.pdf"))

                            h_norm = h_pred[:,1:n_time]' ./ gt_h[1:n_time]
                            h_normalized_plot = set_plot(fs, sz=(plt_width, plt_height))
                            StatsPlots.errorline!(h_normalized_plot, time[1:n_time], h_norm', label="Normalized height estimation error", legend=:outerbottom, legend_column=1)
                            Plots.xlabel!(L"\mathrm{Time\;(s)}")
                            Plots.ylabel!(L"\mathrm{Normalized\;Height\;(mm)}")
                            Plots.xlims!(h_normalized_plot, 0, end_obs_win)
                            Plots.savefig(joinpath(exp_path,"Results","plots","h_normalized_rel_error_noisy.pdf"))
                        # catch err
                        #     @error "Replotting for noisy data with noise $noise_level failed: $err"
                        # end
                    end
                elseif viscosity_type == "bulk_viscosity"
                    local η_gt = Vector{Float64}(undef, 0)
                    local β_gt = Vector{Float64}(undef, 0)

                    # for 1/2 linewidth
                    # fs::Int = 12
                    # bulk_plt_height::Int = 350
                    # bulk_plt_width ::Int = 477
                    # bulk_lft_margin = 1pt
                    # bulk_right_margin = 5pt
                    # bulk_top_margin = 1pt

                    # for 1/3 linewidth
                    # fs::Int = 12
                    # plt_height::Int = 255
                    # plt_width::Int = 318
                    # lft_margin = -5pt
                    # right_margin = -1pt
                    # top_margin = -5pt

                    window_dirs = readdir(exp_path)
                    for window_dir in window_dirs
                        if window_dir == "Results" || window_dir == "post_analysis_window" || window_dir == "single_window" || window_dir == "post_analysis_noise"
                            println("Skipping directory: $window_dir")
                            continue
                        end
                        win_exp_path = joinpath(filepath, elem_size_folder, sim_time_folder, noise_folder, window_dir)
                        
                        println("Processing window: $win_exp_path")
                        exp_params = read_json(joinpath(win_exp_path ,"Results","data","experiment_parameters.json"))
                        sim_time_exp = exp_params["sim_time_exp"]
                        data_type = exp_params["data_type"]
                        ne = exp_params["ne_exp"]
                        
                        if data_type != "physical"
                            noise_level = exp_params["noise_level"]
                            gt_h_ = readdlm(joinpath(filepath_gt,"data","h.csv"), ',', Float64)
                        else
                            noise_level = 0.0
                        end

                        # try 
                            if data_type  == "synthetic" || data_type == "physical"
                                ObsDataList, splinexObs, splineyObs = read_csv(joinpath(filepath_gt,"data","img_data","contour_data"))  
                                @info "Data type $data_type Reading synthetic contour data of $(length(ObsDataList)) time steps"
                            elseif data_type == "simulated"
                                ObsDataList, splinexObs, splineyObs = read_csv(joinpath(filepath_gt,"data","sim_data","contour_data"))  
                                @info "Data type $data_type Reading simulated contour data from $(joinpath(filepath_gt,"data","sim_data","contour_data")) of $(length(ObsDataList)) time steps"
                            else
                                error("Unknown data type: $data_type")
                            end

                            obsBorderPts, nSplinex, nSpliney, pd = add_noise(ObsDataList, nFactor=0.0)

                            est_ηpList = readdlm(joinpath(win_exp_path,"Results","data","est_η.csv"), ',', Float64)
                            est_βpList = readdlm(joinpath(win_exp_path,"Results","data","est_β.csv"), ',', Float64)
                            est_h_list = readdlm(joinpath(win_exp_path,"Results","data","est_h.csv"), ',', Float64)
                            ratio_gt = similar(est_ηpList)
                            ratio_opt = similar(est_ηpList)
                            avg_ηList = similar(est_ηpList)
                            
                            data_ranges_ = get_time_windows(joinpath(win_exp_path,"Results","data","window_data","data_ranges.csv"))
                            t_windows = readdlm(joinpath(win_exp_path,"Results","data","window_data","t_windows.csv"),',',Float64)
                            time_windows = readdlm(joinpath(win_exp_path,"Results","data","window_data","time_windows.csv"),',',Float64)

                            symBorderPts, splinex, spliney = read_csv(joinpath(win_exp_path,"Results","data","sim_data","2D_border_points"))
                            
                            println("Time windows: $(time_windows)")
                            obs_time = sum(time_windows)
                            
                            if obs_time != sim_time
                                @warn "Observation time frame $obs_time is less than preset ground truth time frame $sim_time, updating time frame"
                                min_time = min(obs_time, sim_time)
                                sim_time = min_time
                                obs_time = min_time
                            end
                            
                            if obs_time < sim_time_exp
                                @warn "Observation time frame $obs_time is less than experimental simulation time frame $sim_time_exp, switching to observation time frame"
                                sim_time_exp = obs_time
                            end
                                                
                            data_point_len = round(Int, obs_time/t_steps)
                            if data_type != "physical"
                                if viscosity_model != "carreau"
                                    η_gt = float.(sim_params["η"])
                                    β_gt = float.(sim_params["β"])
                                    η_gt = η_gt[1:data_point_len]
                                    avg_ηList = readdlm(joinpath(win_exp_path,"Results","data","avg_η.csv"), ',', Float64)
                                    ratio_opt = est_ηpList ./ est_βpList
                                    ratio_gt = η_gt ./ β_gt
                                end
                                gt_h = gt_h_[1:(data_point_len+1)]
                            end
                            symBorderPts = symBorderPts[1:data_point_len+1, :]
                            obsBorderPts = obsBorderPts[1:data_point_len+1, :]
                            est_h_list = est_h_list[1:data_point_len+1, :]

                            d, pairs = closest_point(symBorderPts, obsBorderPts)
                            
                            t_full_h = collect(range(start=0, stop=sim_time, step=t_steps))
                            
                            plt_cnt_error = set_plot(fs, sz=(plt_width, plt_height))
                            Plots.plot!(plt_cnt_error, [], label=false, legend=:outerbottom, left_margin = plt_lft_margin, right_margin = plt_right_margin, top_margin = plt_top_margin)
                            Plots.plot!(plt_cnt_error, t_full_h, d, label="Closest point distance error", legend=:outerbottom, legend_column=2)
                            for ti::Int in 1:(size(data_ranges_, 1)-1)
                                t = t_windows[ti]
                                Plots.vline!(plt_cnt_error, [t], color=:gray, linestyle=:dash, label=false)
                            end
                            Plots.xlabel!(L"\mathrm{Time\;(s)}")
                            Plots.ylabel!(L"\mathrm{Closest\;Point\;Distance\;(px)}")
                            Plots.xlims!(plt_cnt_error, 0, end_obs_win)
                            Plots.savefig(plt_cnt_error, joinpath(win_exp_path,"Results","plots","closest_point_distance_error.pdf"))
                            
                            t_full = collect(range(start=t_steps, stop=sim_time, step=t_steps))
                            if data_type != "physical" && viscosity_model != "carreau"
                                plot_ratios = set_plot(fs, sz=(plt_width, plt_height))
                                Plots.plot!(plot_ratios, [], label=false, legend=:outerbottom, left_margin = plt_lft_margin, right_margin = plt_right_margin, top_margin = plt_top_margin)
                                Plots.plot!(plot_ratios, t_full, ratio_gt, label=L"\eta_{gt}/\beta_{gt}", legend=:outerbottom, legend_column=2)
                                # Plots.plot!(plot_ratios, ratio_opt, label="Est η/β", legend=:outerbottom, legend_column=2)
                                t_prev = 0.1
                                for ti::Int in 1:(size(data_ranges_, 1)-1)
                                    t = t_windows[ti]
                                    Plots.vline!(plot_ratios, [t], color=:gray, linestyle=:dash, label=false)
                                end
                                for ti::Int in 1:(size(data_ranges_, 1)-1)
                                    t = t_windows[ti]
                                    data_range_ = data_ranges_[ti]
                                    t_win = collect(range(start=t_prev, stop=t, step=t_steps))
                                    if ti == 1
                                        Plots.plot!(plot_ratios, t_win, ratio_opt[data_range_], label=L"\eta_{est}/\beta_{est}", color=def_orange)
                                    else
                                        Plots.plot!(plot_ratios, t_win, ratio_opt[data_range_], color=def_orange, label=false)
                                    end
                                    t_prev = t+t_steps
                                end
                                Plots.xlabel!(L"\mathrm{Time\;(s)}")
                                Plots.ylabel!(latexstring("\$\\eta/\\beta\\;(KPa\\cdot s \\cdot mm^{-1})\$"))
                                Plots.xlims!(plot_ratios, 0, end_obs_win)
                                Plots.savefig(plot_ratios, joinpath(win_exp_path,"Results","plots","η_over_β.pdf"))  
                            end

                            plt_η = set_plot(fs, sz=(plt_width, plt_height), legend_column=2)
                            Plots.plot!(plt_η, [], label=false, legend=:outerbottom, left_margin = plt_lft_margin, right_margin = plt_right_margin, top_margin = plt_top_margin)
                            if data_type != "physical" && viscosity_model != "carreau"
                                Plots.plot!(plt_η, t_full, η_gt, label=L"\eta_{gt}(t)")
                            end
                            for ti::Int in 1:(size(data_ranges_, 1)-1)
                                t = t_windows[ti]
                                Plots.vline!(plt_η, [t], color=:gray, linestyle=:dash, label=false)
                            end
                            Plots.plot(plt_η, [], label=false, legend=:outerbottom, legend_column=2)
                            # Plots.savefig(plt_η, joinpath(win_exp_path,"Results","plots","η_gt.pdf"))
                            t_prev = 0.1
                            for ti::Int in 1:(size(data_ranges_, 1)-1)
                                t = t_windows[ti]
                                data_range_ = data_ranges_[ti]
                                t_win = collect(range(start=t_prev, stop=t, step=t_steps))
                                if ti == 1
                                    if data_type != "physical" && viscosity_model != "carreau"
                                        Plots.plot!(plt_η, t_win, avg_ηList[data_range_], label=L"\eta_{avg}(t)", color=:gray, legend_column=3)
                                    end 
                                    Plots.plot!(plt_η, t_win, est_ηpList[data_range_], label=L"\eta_{est}(t)", color=def_orange)
                                else
                                    Plots.plot!(plt_η, t_win, est_ηpList[data_range_], color=def_orange, label=false)
                                    if data_type != "physical" && viscosity_model != "carreau"
                                        Plots.plot!(plt_η, t_win, avg_ηList[data_range_], color=:gray, label=false)
                                    end 
                                end
                                t_prev = t+t_steps
                            end
                            Plots.xlabel!(L"\mathrm{Time\;(s)}")
                            Plots.ylabel!(latexstring("\$\\eta_{est}(t)\$ (KPa \$\\cdot\$ s)"))
                            Plots.xlims!(plt_η, 0, end_obs_win)
                            Plots.savefig(plt_η, joinpath(win_exp_path,"Results","plots","η.pdf"))
                            
                            plt_β = set_plot(fs, sz=(plt_width, plt_height), legend_column=2)
                            Plots.plot!(plt_β, [], label=false, legend=:outerbottom, left_margin = plt_lft_margin, right_margin = plt_right_margin, top_margin = plt_top_margin)
                            if data_type != "physical" && viscosity_model != "carreau"
                                Plots.hline!(plt_β, β_gt, label=L"\beta_{gt}")
                            end
                            t_prev = 0.1
                            for ti::Int in 1:(size(data_ranges_, 1)-1)
                                t = t_windows[ti]
                                Plots.vline!(plt_β, [t], color=:gray, linestyle=:dash, label=false)
                                data_range_ = data_ranges_[ti]
                                t_win = collect(range(start=t_prev, stop=t, step=t_steps))
                                if ti == 1
                                    Plots.plot!(plt_β, t_win, est_βpList[data_range_], label=L"\beta_{est}(t)", color=def_orange)
                                else
                                    Plots.plot!(plt_β, t_win, est_βpList[data_range_], color=def_orange, label=false)
                                end
                                t_prev = t+t_steps
                            end
                            Plots.xlabel!(L"\mathrm{Time\;(s)}")
                            Plots.ylabel!(latexstring("\$\\beta_{est}(t)\$ (mm\$^{-1}\$)"))
                            Plots.xlims!(plt_β, 0, end_obs_win)
                            Plots.savefig(plt_β, joinpath(win_exp_path,"Results","plots","β.pdf"))
                            
                            # animate_fields(filepath=joinpath(win_exp_path,"Results","plots"), p=splinex, q=spliney, pObs=nSplinex, qObs=nSpliney)

                        
                            if data_type != "physical"
                                h_plt = set_plot(fs, sz=(plt_width, plt_height), legend_column=2)
                                Plots.plot!(h_plt, t_full_h, gt_h, label=L"h_{gt}", left_margin = plt_lft_margin, right_margin = plt_right_margin, top_margin = plt_top_margin)
                                Plots.plot!(h_plt, t_full_h, est_h_list, label=L"h_{est}")
                                for ti::Int in 1:(size(data_ranges_, 1)-1)
                                    t = t_windows[ti]
                                    Plots.vline!(h_plt, [t], color=:gray, linestyle=:dash, label=false)
                                end
                                Plots.xlabel!(h_plt, L"\mathrm{Time\;(s)}")
                                Plots.ylabel!(h_plt, L"\mathrm{Height\;(mm)}")
                                Plots.xlims!(h_plt, 0, end_obs_win)
                                Plots.savefig(joinpath(win_exp_path,"Results","plots","h.pdf"))

                                h_normalized_plt = set_plot(fs, sz=(plt_width, plt_height))
                                Plots.plot!(h_normalized_plt, t_full_h, est_h_list./gt_h, label=L"h_{est}/h_{gt}", left_margin = plt_lft_margin, right_margin = plt_right_margin, top_margin = plt_top_margin)
                                for ti::Int in 1:(size(data_ranges_, 1)-1)
                                    t = t_windows[ti]
                                    Plots.vline!(h_normalized_plt, [t], color=:gray, linestyle=:dash, label=false)
                                end
                                Plots.xlabel!(h_normalized_plt, L"\mathrm{Time\;(s)}")
                                Plots.ylabel!(h_normalized_plt, L"h_{est}/h_{gt}")
                                Plots.xlims!(h_normalized_plt, 0, end_obs_win)
                                Plots.savefig(joinpath(win_exp_path,"Results","plots","h_normalized.pdf"))

                                error_plt = set_plot(fs, sz=(plt_width, plt_height))
                                Plots.plot!(error_plt, t_full_h, abs.(est_h_list-gt_h), label="Height estimation error", left_margin = plt_lft_margin, right_margin = plt_right_margin, top_margin = plt_top_margin)
                                for ti::Int in 1:(size(data_ranges_, 1)-1)
                                    t = t_windows[ti]
                                    Plots.vline!(error_plt, [t], color=:gray, linestyle=:dash, label=false)
                                end
                                Plots.savefig(joinpath(win_exp_path,"Results","plots","h_est_error.pdf"))
                                Plots.xlabel!(error_plt, L"\mathrm{Time\;(s)}")
                                Plots.ylabel!(error_plt, L"\mathrm{Height\;Error\;(mm)}")
                                Plots.xlims!(error_plt, 0, end_obs_win)
                                Plots.savefig(joinpath(win_exp_path,"Results","plots","h_est_error.pdf"))

                                rel_error_plt = set_plot(fs, sz=(plt_width, plt_height))
                                Plots.plot!(rel_error_plt, t_full_h, abs.(est_h_list-gt_h)./gt_h*100, label="Relative Height estimation error", left_margin = plt_lft_margin, right_margin = plt_right_margin, top_margin = plt_top_margin)
                                for ti::Int in 1:(size(data_ranges_, 1)-1)
                                    t = t_windows[ti]
                                    Plots.vline!(rel_error_plt, [t], color=:gray, linestyle=:dash, label=false)
                                end
                                Plots.xlabel!(rel_error_plt, L"\mathrm{Time\;(s)}")
                                Plots.ylabel!(rel_error_plt, latexstring("Relative Height Error (\$\\%\$)"))
                                Plots.xlims!(rel_error_plt, 0, end_obs_win)
                                Plots.savefig(joinpath(win_exp_path,"Results","plots","rel_error.pdf"))
                            end
                        # catch err
                        #         @warn "Post-analysis for window $window_dir failed: $err"
                        # end
                    end

                end
            end
        end
    end
    return 
end

function post_analysis_const(filepath_gt_::String, filepath::String, avoid_list)
    # ηList = readdlm(joinpath(filepath_gt,"data","η.csv"), ',', Float64)
    # βList = readdlm(joinpath(filepath_gt,"data","β.csv"), ',', Float64)

    # fs::Int = 12
    # noise_cols::Int = 2
    # plt_height::Int = 255
    # plt_width::Int = 318
    # plt_lft_margin = -5pt
    # plt_right_margin = -1pt
    # plt_top_margin = -5pt

    dir_list = readdir(filepath)

    # plot for convergence per slip case
    plot_conv_2 = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, legend_column=3)
    Plots.plot!(plot_conv_2, [],  label=false)
    Plots.xlabel!(plot_conv_2, L"\mathrm{Iterations}")
    Plots.ylabel!(plot_conv_2, L"\mathrm{Cost\;(px)}")
    
    plot_conv_5 = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, legend_column=3)
    Plots.plot!(plot_conv_5, [],  label=false)
    Plots.xlabel!(plot_conv_5, L"\mathrm{Iterations}")
    Plots.ylabel!(plot_conv_5, L"\mathrm{Cost\;(px)}")
    
    plot_conv_10 = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, legend_column=3)
    Plots.plot!(plot_conv_10, [],  label=false)
    Plots.xlabel!(plot_conv_10, L"\mathrm{Iterations}")
    Plots.ylabel!(plot_conv_10, L"\mathrm{Cost\;(px)}")

    plot_conv_20 = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, legend_column=3)
    Plots.plot!(plot_conv_20, [],  label=false)
    Plots.xlabel!(plot_conv_20, L"\mathrm{Iterations}")
    Plots.ylabel!(plot_conv_20, L"\mathrm{Cost\;(px)}")

    plot_conv_30 = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, legend_column=3)
    Plots.plot!(plot_conv_30, [],  label=false)
    Plots.xlabel!(plot_conv_30, L"\mathrm{Iterations}")
    Plots.ylabel!(plot_conv_30, L"\mathrm{Cost\;(px)}")

    plot_conv_log_2 = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, legend_column=3)
    Plots.plot!(plot_conv_log_2, [],  label=false)
    Plots.xlabel!(plot_conv_log_2, L"\mathrm{Iterations}")
    Plots.ylabel!(plot_conv_log_2, L"\mathrm{Cost\;(px)}")

    plot_conv_log_5 = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, legend_column=3)
    Plots.plot!(plot_conv_log_5, [],  label=false)
    Plots.xlabel!(plot_conv_log_5, L"\mathrm{Iterations}")
    Plots.ylabel!(plot_conv_log_5, L"\mathrm{Cost\;(px)}")

    plot_conv_log_10 = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, legend_column=3)
    Plots.plot!(plot_conv_log_10, [],  label=false)
    Plots.xlabel!(plot_conv_log_10, L"\mathrm{Iterations}")
    Plots.ylabel!(plot_conv_log_10, L"\mathrm{Cost\;(px)}")

    plot_conv_log_20 = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, legend_column=3)
    Plots.plot!(plot_conv_log_20, [],  label=false)
    Plots.xlabel!(plot_conv_log_20, L"\mathrm{Iterations}")
    Plots.ylabel!(plot_conv_log_20, L"\mathrm{Cost\;(px)}")

    plot_conv_log_30 = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, legend_column=3)
    Plots.plot!(plot_conv_log_30, [],  label=false)
    Plots.xlabel!(plot_conv_log_30, L"\mathrm{Iterations}")
    Plots.ylabel!(plot_conv_log_30, L"\mathrm{Cost\;(px)}")

    η_norm_plot_2 = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, legend_column=3)
    Plots.hline!(η_norm_plot_2, [1.0],  linestyle=:dash, label=false, color=:black)
    Plots.xlabel!(η_norm_plot_2,L"\mathrm{Iterations}")
    Plots.ylabel!(η_norm_plot_2,L"\eta_{est}/\eta_{gt}")

    β_norm_plot_2 = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, legend_column=3)
    Plots.hline!(β_norm_plot_2, [1.0],  linestyle=:dash, label=false, color=:black)
    Plots.xlabel!(β_norm_plot_2, L"\mathrm{Iterations}")
    Plots.ylabel!(β_norm_plot_2, L"\beta_{est}/\beta_{gt}")

    η_norm_plot_5 = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, legend_column=3)
    Plots.hline!(η_norm_plot_5, [1.0],  linestyle=:dash, label=false, color=:black)
    Plots.xlabel!(η_norm_plot_5,L"\mathrm{Iterations}")
    Plots.ylabel!(η_norm_plot_5,L"\eta_{est}/\eta_{gt}")

    β_norm_plot_5 = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, legend_column=3)
    Plots.hline!(β_norm_plot_5, [1.0],  linestyle=:dash, label=false, color=:black)
    Plots.xlabel!(β_norm_plot_5, L"\mathrm{Iterations}")
    Plots.ylabel!(β_norm_plot_5, L"\beta_{est}/\beta_{gt}")
    
    η_norm_plot_10 = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, legend_column=3)
    Plots.hline!(η_norm_plot_10, [1.0],  linestyle=:dash, label=false, color=:black)
    Plots.xlabel!(η_norm_plot_10,L"\mathrm{Iterations}")
    Plots.ylabel!(η_norm_plot_10,L"\eta_{est}/\eta_{gt}")

    β_norm_plot_10 = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, legend_column=3)
    Plots.hline!(β_norm_plot_10, [1.0],  linestyle=:dash, label=false, color=:black)
    Plots.xlabel!(β_norm_plot_10, L"\mathrm{Iterations}")
    Plots.ylabel!(β_norm_plot_10, L"\beta_{est}/\beta_{gt}")

    η_norm_plot_20 = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, legend_column=3)
    Plots.hline!(η_norm_plot_20, [1.0],  linestyle=:dash, label=false, color=:black)
    Plots.xlabel!(η_norm_plot_20,L"\mathrm{Iterations}")
    Plots.ylabel!(η_norm_plot_20,L"\eta_{est}/\eta_{gt}")

    β_norm_plot_20 = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, legend_column=3)
    Plots.hline!(β_norm_plot_20, [1.0],  linestyle=:dash, label=false, color=:black)
    Plots.xlabel!(β_norm_plot_20, L"\mathrm{Iterations}")
    Plots.ylabel!(β_norm_plot_20, L"\beta_{est}/\beta_{gt}")

    η_norm_plot_30 = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, legend_column=3)
    Plots.hline!(η_norm_plot_30, [1.0],  linestyle=:dash, label=false, color=:black)
    Plots.xlabel!(η_norm_plot_30,L"\mathrm{Iterations}")
    Plots.ylabel!(η_norm_plot_30,L"\eta_{est}/\eta_{gt}")

    β_norm_plot_30 = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, legend_column=3)
    Plots.hline!(β_norm_plot_30, [1.0],  linestyle=:dash, label=false, color=:black)
    Plots.xlabel!(β_norm_plot_30, L"\mathrm{Iterations}")
    Plots.ylabel!(β_norm_plot_30, L"\beta_{est}/\beta_{gt}")

    # η\β ratio plot
    ratio_plot_2 = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, legend_column=3)
    Plots.hline!(ratio_plot_2, [1.0],  linestyle=:dash, label=false)
    Plots.xlabel!(ratio_plot_2,L"\mathrm{Iterations}")
    Plots.ylabel!(ratio_plot_2,L"\eta_{est}/\eta_{gt}")

    ratio_plot_5 = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, legend_column=3)
    Plots.hline!(ratio_plot_5, [1.0],  linestyle=:dash, label=false)
    Plots.xlabel!(ratio_plot_5,L"\mathrm{Iterations}")
    Plots.ylabel!(ratio_plot_5,L"\eta_{est}/\eta_{gt}")
    
    ratio_plot_10 = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, legend_column=3)
    Plots.hline!(ratio_plot_10, [1.0],  linestyle=:dash, label=false)
    Plots.xlabel!(ratio_plot_10,L"\mathrm{Iterations}")
    Plots.ylabel!(ratio_plot_10,L"\eta_{est}/\eta_{gt}")

    ratio_plot_20 = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, legend_column=3)
    Plots.hline!(ratio_plot_20, [1.0],  linestyle=:dash, label=false)
    Plots.xlabel!(ratio_plot_20, L"\mathrm{Iterations}")
    Plots.ylabel!(ratio_plot_20, L"\beta_{est}/\beta_{gt}")

    ratio_plot_30 = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, legend_column=3)
    Plots.hline!(ratio_plot_30, [1.0],  linestyle=:dash, label=false)
    Plots.xlabel!(ratio_plot_30,L"\mathrm{Iterations}")
    Plots.ylabel!(ratio_plot_30,L"\eta_{est}/\eta_{gt}")

    # height plots
    h_glob_plot_2 = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, legend_column=3)
    Plots.xlabel!(h_glob_plot_2,L"\mathrm{Time\;(s)}")
    Plots.ylabel!(h_glob_plot_2,L"h_{est}\;\mathrm{(mm)}")
    Plots.xlims!(h_glob_plot_2, 0, end_obs_win)

    h_glob_plot_5 = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, legend_column=3)
    Plots.xlabel!(h_glob_plot_5,L"\mathrm{Time\;(s)}")
    Plots.ylabel!(h_glob_plot_5,L"h_{est}\;\mathrm{(mm)}")  
    Plots.xlims!(h_glob_plot_5, 0, end_obs_win) 

    h_glob_plot_10 = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, legend_column=3)
    Plots.xlabel!(h_glob_plot_10,L"\mathrm{Time\;(s)}")
    Plots.ylabel!(h_glob_plot_10,L"h_{est}\;\mathrm{(mm)}") 
    Plots.xlims!(h_glob_plot_10, 0, end_obs_win)

    #normalized height plots

    h_norm_plot_2 = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, legend_column=3)
    Plots.hline!(h_norm_plot_2, [1.0],  left_margin=plt_lft_margin, linestyle=:dash, label=false, color=:black)
    Plots.xlabel!(h_norm_plot_2,L"\mathrm{Time\;(s)}")
    Plots.ylabel!(h_norm_plot_2,L"h_{est}/h_{gt}")
    Plots.xlims!(h_norm_plot_2, 0, end_obs_win)
    Plots.ylims!(h_norm_plot_2, y_lims_h_norm)

    h_norm_plot_5 = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, legend_column=3)
    Plots.hline!(h_norm_plot_5, [1.0],  left_margin=plt_lft_margin, linestyle=:dash, label=false, color=:black)
    Plots.xlabel!(h_norm_plot_5,L"\mathrm{Time\;(s)}")
    Plots.ylabel!(h_norm_plot_5,L"h_{est}/h_{gt}")  
    Plots.xlims!(h_norm_plot_5, 0, end_obs_win)
    Plots.ylims!(h_norm_plot_5, y_lims_h_norm)

    h_norm_plot_10 = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, legend_column=3)
    Plots.hline!(h_norm_plot_10, [1.0],  left_margin=plt_lft_margin, linestyle=:dash, label=false, color=:black)
    Plots.xlabel!(h_norm_plot_10,L"\mathrm{Time\;(s)}")
    Plots.ylabel!(h_norm_plot_10,L"h_{est}/h_{gt}") 
    Plots.xlims!(h_norm_plot_10, 0, end_obs_win)
    Plots.ylims!(h_norm_plot_10, y_lims_h_norm)

    h_norm_plot_20 = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, legend_column=3)
    Plots.hline!(h_norm_plot_20, [1.0],  left_margin=plt_lft_margin, linestyle=:dash, label=false, color=:black)
    Plots.xlabel!(h_norm_plot_20,L"\mathrm{Time\;(s)}")
    Plots.ylabel!(h_norm_plot_20,L"h_{est}/h_{gt}") 
    Plots.xlims!(h_norm_plot_20, 0, end_obs_win)
    Plots.ylims!(h_norm_plot_20, y_lims_h_norm)

    h_norm_plot_30 = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, legend_column=3)
    Plots.hline!(h_norm_plot_30, [1.0],  left_margin=plt_lft_margin, linestyle=:dash, label=false, color=:black)
    Plots.xlabel!(h_norm_plot_30,L"\mathrm{Time\;(s)}")
    Plots.ylabel!(h_norm_plot_30,L"h_{est}/h_{gt}")
    Plots.xlims!(h_norm_plot_30, 0, end_obs_win)
    Plots.ylims!(h_norm_plot_30, y_lims_h_norm)

    # relative height error plots
    rel_height_error_glob_plot_2 = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, legend_column=3)
    Plots.xlabel!(rel_height_error_glob_plot_2, L"\mathrm{Time\;(s)}")
    Plots.ylabel!(rel_height_error_glob_plot_2, latexstring("Relative Height Error (\$\\%\$)"))
    Plots.xlims!(rel_height_error_glob_plot_2, 0, end_obs_win)
    Plots.ylims!(rel_height_error_glob_plot_2, y_lims_rel_error)

    rel_height_error_glob_plot_5 = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, legend_column=3)
    Plots.xlabel!(rel_height_error_glob_plot_5, L"\mathrm{Time\;(s)}")
    Plots.ylabel!(rel_height_error_glob_plot_5, latexstring("Relative Height Error (\$\\%\$)"))
    Plots.xlims!(rel_height_error_glob_plot_5, 0, end_obs_win)
    Plots.ylims!(rel_height_error_glob_plot_5, y_lims_rel_error)

    rel_height_error_glob_plot_10 = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, legend_column=3)
    Plots.xlabel!(rel_height_error_glob_plot_10, L"\mathrm{Time\;(s)}")
    Plots.ylabel!(rel_height_error_glob_plot_10, latexstring("Relative Height Error (\$\\%\$)"))
    Plots.xlims!(rel_height_error_glob_plot_10, 0, end_obs_win)
    Plots.ylims!(rel_height_error_glob_plot_10, y_lims_rel_error)

    rel_height_error_glob_plot_20 = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, legend_column=3)
    Plots.xlabel!(rel_height_error_glob_plot_20, L"\mathrm{Time\;(s)}")
    Plots.ylabel!(rel_height_error_glob_plot_20, latexstring("Relative Height Error (\$\\%\$)"))
    Plots.xlims!(rel_height_error_glob_plot_20, 0, end_obs_win)
    Plots.ylims!(rel_height_error_glob_plot_20, y_lims_rel_error)

    rel_height_error_glob_plot_30 = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, legend_column=3)
    Plots.xlabel!(rel_height_error_glob_plot_30, L"\mathrm{Time\;(s)}")
    Plots.ylabel!(rel_height_error_glob_plot_30, latexstring("Relative Height Error (\$\\%\$)"))
    Plots.xlims!(rel_height_error_glob_plot_30, 0, end_obs_win)
    Plots.ylims!(rel_height_error_glob_plot_30, y_lims_rel_error)

    # η \ β * β_gt \ η_gt plot
    ratio_norm_plot_2 = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, legend_column=3)
    Plots.hline!(ratio_norm_plot_2, [1.0],  linestyle=:dash, label=false, color=:black)
    Plots.xlabel!(ratio_norm_plot_2,L"\mathrm{Iterations}")
    Plots.ylabel!(ratio_norm_plot_2,L"\frac{\eta_{est}/\beta_{est}}{\eta_{gt}/\beta_{gt}}")


    ratio_norm_plot_5 = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, legend_column=3)
    Plots.hline!(ratio_norm_plot_5, [1.0],  linestyle=:dash, label=false, color=:black)
    Plots.xlabel!(ratio_norm_plot_5,L"\mathrm{Iterations}")
    Plots.ylabel!(ratio_norm_plot_5,L"\frac{\eta_{est}/\beta_{est}}{\eta_{gt}/\beta_{gt}}")

    ratio_norm_plot_10 = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, legend_column=3)
    Plots.hline!(ratio_norm_plot_10, [1.0],  linestyle=:dash, label=false, color=:black)
    Plots.xlabel!(ratio_norm_plot_10,L"\mathrm{Iterations}")
    Plots.ylabel!(ratio_norm_plot_10,L"\frac{\eta_{est}/\beta_{est}}{\eta_{gt}/\beta_{gt}}")

    ratio_norm_plot_20 = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, legend_column=3)
    Plots.hline!(ratio_norm_plot_20, [1.0],  linestyle=:dash, label=false, color=:black)
    Plots.xlabel!(ratio_norm_plot_20, L"\mathrm{Iterations}")
    Plots.ylabel!(ratio_norm_plot_20, L"\frac{\eta_{est}/\beta_{est}}{\eta_{gt}/\beta_{gt}}")

    ratio_norm_plot_30 = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, legend_column=3)
    Plots.hline!(ratio_norm_plot_30, [1.0],  linestyle=:dash, label=false, color=:black)
    Plots.xlabel!(ratio_norm_plot_30,L"\mathrm{Iterations}")
    Plots.ylabel!(ratio_norm_plot_30,L"\frac{\eta_{est}/\beta_{est}}{\eta_{gt}/\beta_{gt}}")
    
    # conotour_plots
    cont_y_min = 350
    contour_plt = set_plot(fs, sz=(plt_width, plt_height))
    Plots.plot!(contour_plt, [], label=false, legend=:outerbottom, legend_column=2, aspect_ratio = :equal)
    Plots.yflip!(true)
    Plots.xlims!(contour_plt, 480, 1520)
    Plots.ylims!(contour_plt, cont_y_min, 1200)
    Plots.xlabel!(contour_plt, L"x\;\mathrm{(px)}")
    Plots.ylabel!(contour_plt, L"y\;\mathrm{(px)}")

    contour_plt_zoom = set_plot(fs, sz=(1200, plt_height))
    Plots.plot!(contour_plt_zoom, [], label=false, legend=:outerbottom, legend_column=2, aspect_ratio = :equal)
    Plots.xticks!(contour_plt_zoom, 1100:200:1520)
    Plots.yflip!(true)
    Plots.xlims!(contour_plt_zoom, 1000, 1520)
    Plots.ylims!(contour_plt_zoom, cont_y_min, 1200)
    Plots.xlabel!(contour_plt_zoom, L"x\;\mathrm{(px)}")
    Plots.ylabel!(contour_plt_zoom, L"y\;\mathrm{(px)}")

    max_iter = 1
    for dir in dir_list
    
        if dir in avoid_list || dir == "post_analysis_global"
            println("Skipping directory: ", dir)
            continue
        end
        filepath_dir = joinpath(filepath, dir)
        filepath_gt = joinpath(filepath_gt_, dir)

        printstyled("Processing directory: $(filepath_dir)\n", color=:green)
        sim_params = read_json(joinpath(filepath_gt,"data","sim_params.json"))
        η_gt = sim_params["η"]
        β_gt = sim_params["β"]  
        sim_time = sim_params["simulation_time"]    
        t_steps = sim_params["time_steps"]
        gt_h = readdlm(joinpath(filepath_gt,"data","h.csv"), ',', Float64)

        time = collect(Float64, range(start=0, stop=sim_time, step=t_steps))

        println("Ground truth η: ", η_gt[1])

        elem_size_folders = readdir(filepath_dir)

        # plots for element vise comparison
        # plots for convergence 
        elem_η_plt = set_plot(fs, sz=(plt_width, plt_height))
        Plots.hline!(elem_η_plt, [η_gt[1]], label="Ground truth η",  left_margin=plt_lft_margin)
        Plots.xlabel!(elem_η_plt,L"\mathrm{Iterations}")
        Plots.ylabel!(elem_η_plt,L"\eta\;\mathrm{(KPa\cdot s)}")

        elem_β_plt = set_plot(fs, sz=(plt_width, plt_height))
        Plots.hline!(elem_β_plt, [β_gt[1]], label="Ground truth β",  left_margin=plt_lft_margin)
        Plots.xlabel!(elem_β_plt, L"\mathrm{Iterations}")
        Plots.ylabel!(elem_β_plt, L"\beta\;\mathrm{(mm^{-1})}")

        elem_ratio_plt = set_plot(fs, sz=(plt_width, plt_height))
        Plots.hline!(elem_ratio_plt, [η_gt[1]/β_gt[1]], label="Ground truth η/β",  left_margin=plt_lft_margin)
        Plots.xlabel!(elem_ratio_plt, L"\mathrm{Iterations}")   
        Plots.ylabel!(elem_ratio_plt, L"\eta/\beta\;(KPa\cdot s \cdot mm^{-1})")
        
        # plots for normalized values
        elem_η_norm_plt = set_plot(fs, sz=(plt_width, plt_height))
        Plots.hline!(elem_η_norm_plt, [1.0], label=false,  left_margin=plt_lft_margin, linestyle=:dash, color=:black)
        Plots.xlabel!(elem_η_norm_plt,L"\mathrm{Iterations}")
        Plots.ylabel!(elem_η_norm_plt,L"\eta_{est}/\eta_{gt}")

        elem_β_norm_plt = set_plot(fs, sz=(plt_width, plt_height))
        Plots.hline!(elem_β_norm_plt, [1.0], label=false,  left_margin=plt_lft_margin, linestyle=:dash, color=:black)
        Plots.xlabel!(elem_β_norm_plt, L"\mathrm{Iterations}")
        Plots.ylabel!(elem_β_norm_plt, L"\beta_{est}/\beta_{gt}")

        elem_ratio_norm_plt = set_plot(fs, sz=(plt_width, plt_height))
        Plots.hline!(elem_ratio_norm_plt, [1.0], label=false,  left_margin=plt_lft_margin, linestyle=:dash, color=:black)
        Plots.xlabel!(elem_ratio_norm_plt, L"\mathrm{Iterations}")   
        Plots.ylabel!(elem_ratio_norm_plt, L"\frac{\eta_{est}/\beta_{est}}{\eta_{gt}/\beta_{gt}}")

        # plots for height error
        elem_rel_height_error_plt = set_plot(fs, sz=(plt_width, plt_height))
        Plots.plot!(elem_rel_height_error_plt, [], left_margin=plt_lft_margin, label=false)
        Plots.xlabel!(elem_rel_height_error_plt, L"\mathrm{Time\;(s)}")
        Plots.ylabel!(elem_rel_height_error_plt, latexstring("Relative Height Error (\$\\%\$)"))
        Plots.xlims!(elem_rel_height_error_plt, 0, end_obs_win)

        elem_height_plt = set_plot(fs, sz=(plt_width, plt_height))
        Plots.plot!(elem_height_plt, [], left_margin=plt_lft_margin, label=false)
        Plots.xlabel!(elem_height_plt, L"\mathrm{Time\;(s)}")
        Plots.ylabel!(elem_height_plt, L"\mathrm{Height\;(mm)}")
        Plots.xlims!(elem_height_plt, 0, end_obs_win)
        
        for elem_size_folder_ in elem_size_folders
            if elem_size_folder_ == "post_analysis" || elem_size_folder_ == "Q2_2"
                continue
            end
            
            elem_size_folder = joinpath(filepath_dir, elem_size_folder_)
            printstyled("Processing element size folder: $(elem_size_folder)\n", color=:blue)  
            sim_time_folders = readdir(elem_size_folder)
            
            # figures for Simulation window vise camparison
            sim_window_η_plt = set_plot(fs, sz=(plt_width, plt_height))
            Plots.hline!(sim_window_η_plt, [η_gt[1]], label="Ground truth η")
            Plots.xlabel!(sim_window_η_plt,L"\mathrm{Iterations}")
            Plots.ylabel!(sim_window_η_plt,L"\eta\;\mathrm{(KPa\cdot s)}")
            
            sim_window_β_plt = set_plot(fs, sz=(plt_width, plt_height))
            Plots.hline!(sim_window_β_plt, [β_gt[1]], label="Ground truth β")
            Plots.xlabel!(sim_window_β_plt, L"\mathrm{Iterations}")
            Plots.ylabel!(sim_window_β_plt, L"\beta\;\mathrm{(mm^{-1})}")

            sim_window_ratio_plt = set_plot(fs, sz=(plt_width, plt_height))
            Plots.hline!(sim_window_ratio_plt, [η_gt[1]/β_gt[1]], label="Ground truth η/β",  left_margin=plt_lft_margin)
            Plots.xlabel!(sim_window_ratio_plt, L"\mathrm{Iterations}")   
            Plots.ylabel!(sim_window_ratio_plt, L"\eta/\beta\;(KPa\cdot s \cdot mm^{-1})")

            # plots for normalized values
            sim_window_η_norm_plt = set_plot(fs, sz=(plt_width, plt_height))
            Plots.hline!(sim_window_η_norm_plt, [1.0], label=false,  left_margin=plt_lft_margin, linestyle=:dash, color=:black)
            Plots.xlabel!(sim_window_η_norm_plt,L"\mathrm{Iterations}")
            Plots.ylabel!(sim_window_η_norm_plt,L"\eta_{est}/\eta_{gt}")

            sim_window_β_norm_plt = set_plot(fs, sz=(plt_width, plt_height))
            Plots.hline!(sim_window_β_norm_plt, [1.0], label=false,  left_margin=plt_lft_margin, linestyle=:dash, color=:black)
            Plots.xlabel!(sim_window_β_norm_plt, L"\mathrm{Iterations}")
            Plots.ylabel!(sim_window_β_norm_plt, L"\beta_{est}/\beta_{gt}")

            sim_window_ratio_norm_plt = set_plot(fs, sz=(plt_width, plt_height))
            Plots.hline!(sim_window_ratio_norm_plt, [1.0], label=false,  left_margin=plt_lft_margin, linestyle=:dash, color=:black)
            Plots.xlabel!(sim_window_ratio_norm_plt, L"\mathrm{Iterations}")   
            Plots.ylabel!(sim_window_ratio_norm_plt, L"\frac{\eta_{est}/\beta_{est}}{\eta_{gt}/\beta_{gt}}")

            # plots for height error
            sim_window_rel_height_error_plt = set_plot(fs, sz=(plt_width, plt_height))
            Plots.plot!(sim_window_rel_height_error_plt, [], left_margin=plt_lft_margin, label=false)
            Plots.xlabel!(sim_window_rel_height_error_plt, L"\mathrm{Time\;(s)}")
            Plots.ylabel!(sim_window_rel_height_error_plt, latexstring("Relative Height Error (\$\\%\$)"))
            Plots.xlims!(sim_window_rel_height_error_plt, 0, end_obs_win)

            sim_window_height_plt = set_plot(fs, sz=(plt_width, plt_height))
            Plots.plot!(sim_window_height_plt, [], left_margin=plt_lft_margin, label=false)
            Plots.xlabel!(sim_window_height_plt, L"\mathrm{Time\;(s)}")
            Plots.ylabel!(sim_window_height_plt, L"\mathrm{Height\;(mm)}")
            Plots.xlims!(sim_window_height_plt, 0, end_obs_win)

            plt_slices = set_plot(fs, sz=(350, 750))
            Plots.vline!(plt_slices, [0.0], color=:blue, linestyle=:dash, label="Minimum",  lw=1)
            Plots.xlabel!(plt_slices, L"\mathrm{Distance\;from\;minimum\;(px)}")
            Plots.ylabel!(plt_slices, L"\mathrm{Cost}")
            Plots.ylims!(plt_slices, 0, 50)
                
                
            for sim_time_folder_ in reverse(sort(sim_time_folders))

                if sim_time_folder_ == "post_analysis_time" || sim_time_folder_ == "Results" || sim_time_folder_ == "simtime_20.0" || sim_time_folder_ == "simtime_10.0" || sim_time_folder_ == "simtime_30.0"
                    continue
                end
                
                height_error_plt = set_plot(fs, sz=(plt_width, plt_height))
                Plots.plot!(height_error_plt, [], label=false)
                Plots.xlabel!(height_error_plt, L"\mathrm{Time\;(s)}")
                Plots.ylabel!(height_error_plt, L"\mathrm{Height\;Error\;(mm)}")
                Plots.xlims!(height_error_plt, 0, end_obs_win)
                
                rel_height_error_plt = set_plot(fs, sz=(plt_width, plt_height))
                Plots.plot!(rel_height_error_plt, [], label=false)
                Plots.xlabel!(rel_height_error_plt, L"\mathrm{Time\;(s)}")
                Plots.ylabel!(rel_height_error_plt, latexstring("Relative Height Error (\$\\%\$)"))
                Plots.xlims!(rel_height_error_plt, 0, end_obs_win)
                
                η_β_ratio_plt = set_plot(fs, sz=(plt_width, plt_height))
                Plots.plot!(η_β_ratio_plt, [], label=false)
                Plots.hline!(η_β_ratio_plt, [1], label=false, linestyle=:dash)
                Plots.xlabel!(η_β_ratio_plt, L"\mathrm{Iterations}")
                Plots.ylabel!(η_β_ratio_plt, L"\eta/\beta\;(KPa\cdot s \cdot mm^{-1})")
                
                # noise analysis plots
                noise_cols::Int = 2
                covarience_plt = set_plot(fs, sz=(plt_width, plt_height))
                Plots.plot!(covarience_plt, [], label=false, left_margin=plt_lft_margin, right_margin=plt_right_margin, top_margin=plt_top_margin)
                Plots.xlabel!(covarience_plt, L"\eta/\eta_{gt}")
                Plots.ylabel!(covarience_plt, L"\beta/\beta_{gt}")

                height_noise_plt = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, top_margin=plt_top_margin, legend_column=3)
                Plots.xlabel!(height_noise_plt, L"\mathrm{Time\;(s)}")
                Plots.ylabel!(height_noise_plt, L"\mathrm{Height\;(mm)}")
                Plots.xlims!(height_noise_plt, 0, end_obs_win)

                height_error_noise_plt = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, top_margin=plt_top_margin, legend_column=noise_cols)
                Plots.xlabel!(height_error_noise_plt, L"\mathrm{Time\;(s)}")
                Plots.ylabel!(height_error_noise_plt, L"\mathrm{Height\;Error\;(mm)}")
                Plots.xlims!(height_error_noise_plt, 0, end_obs_win)

                rel_height_error_noise_plt = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, top_margin=plt_top_margin, legend_column=noise_cols)
                Plots.plot!(rel_height_error_noise_plt, [], label=false)
                Plots.xlabel!(rel_height_error_noise_plt, L"\mathrm{Time\;(s)}")
                Plots.ylabel!(rel_height_error_noise_plt, latexstring("Relative Height Error (\$\\%\$)"))
                Plots.xlims!(rel_height_error_noise_plt, 0, end_obs_win)

                normalized_height_noise_plt = set_plot(fs, sz=(plt_width, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, top_margin=plt_top_margin, legend_column=noise_cols)
                Plots.plot!(normalized_height_noise_plt, [1.0], linestyle=:dash, color=:black, label=false)
                Plots.ylabel!(normalized_height_noise_plt, L"h_{est}/h_{gt}")
                Plots.xlabel!(normalized_height_noise_plt, L"\mathrm{Time\;(s)}")
                Plots.xlims!(normalized_height_noise_plt, 0, end_obs_win)
                
                sim_time_folder = joinpath(elem_size_folder, sim_time_folder_)
                noise_folders = readdir(sim_time_folder)

                printstyled("Processing simulation time folder: $(sim_time_folder)\n", color=:cyan)
                for noise_folder_ in reverse(sort(noise_folders))
                    
                    if noise_folder_ == "post_analysis_noise" || noise_folder_ == "Results" 
                        continue
                    end
                    
                    exp_path = joinpath(filepath_dir, elem_size_folder_, sim_time_folder_, noise_folder_)

                    printstyled("Processing experiment folder: $(exp_path)\n", color=:magenta)
                    exp_params = read_json(joinpath(exp_path,"Results","data","experiment_parameters.json"))

                    noise_level = exp_params["noise_level"]
                    
                    FunctionClass_x = exp_params["FunctionClass_x"]
                    ne = exp_params["ne_exp"]
                    sim_time = exp_params["sim_time_exp"]
                    data_type = exp_params["data_type"]
                    exp_path_n0 = replace(exp_path, "$noise_folder_" => "noise_0.0")
                    exp_points::Int = round(Int,sim_time/t_steps)
                    exp_points = sim_time/t_steps

                    obsBorderPts, symBorderPts, nSplinex, nSpliney, splinex, spliney = _get_borders(data_type, filepath_gt, exp_path)

                    printstyled("Processing for noise level: $(noise_level): $(exp_path)\n", color=:yellow)
                    if noise_level == 0.0
                        est_η = readdlm(joinpath(exp_path,"Results","data","η.csv"), ',', Float64)
                        est_β = readdlm(joinpath(exp_path,"Results","data","β.csv"), ',', Float64)
                        est_h = readdlm(joinpath(exp_path,"Results","data","est_h.csv"), ',', Float64)
                        iter = collect(range(start=1, stop=size(est_η,1), step=1))

                        max_iter = max(max_iter, size(iter,1))

                        if length(est_h) != length(gt_h)
                            @warn "Length mismatch between estimated height ($(length(est_h))) and ground truth height ($(length(gt_h))). Adjusting..."
                            min_length = min(length(est_h), length(gt_h))
                            est_h = est_h[1:min_length]
                            gt_h = gt_h[1:min_length]
                        end

                        height_error = abs.(est_h - gt_h)
                        rel_height_error = height_error ./ gt_h * 100.0

                        # normalize by ground-truth                        
                        est_η_norm = est_η./η_gt
                        est_β_norm = est_β./β_gt
                        h_norm = est_h ./ gt_h

                        ratio_est = est_η ./ est_β
                        ratio_gt = η_gt / β_gt
                        normalized_ratio = (est_η ./ est_β) * (β_gt / η_gt)

                        cost_list = readdlm(joinpath(exp_path,"Results","data","cost_iter.csv"), ',', Float64)
                
                        Plots.plot!(sim_window_η_plt, iter, est_η, label=string("Window - $(sim_time)s"), marker=1, legend=:outerbottom, legend_column=2)
                        Plots.plot!(sim_window_β_plt, iter, est_β, label=string("Window - $(sim_time)s"), marker=1, legend=:outerbottom, legend_column=2)
                        Plots.plot!(sim_window_ratio_plt, iter, ratio_est, label=string("Window - $(sim_time)s"), marker=1, legend=:outerbottom, legend_column=2)
                        
                        time_h = copy(time)
                        if length(rel_height_error) != length(time) 
                            @warn "Length mismatch between height plot x-axis ($(length(rel_height_error))) and ground truth height ($(length(gt_h))). Adjusting..."
                            min_length = min(length(rel_height_error), length(time))
                            time_h = time[1:min_length]
                            rel_height_error = rel_height_error[1:min_length]
                            est_h = est_h[1:min_length]
                            gt_h = gt_h[1:min_length]
                        end
                        Plots.plot!(sim_window_rel_height_error_plt, time_h, rel_height_error, label=string("Window - $(sim_time)s"," - ne: ",ne), legend=:outerbottom, legend_column=2)
                        Plots.plot!(sim_window_height_plt, time_h, est_h, label=string("Window - $(sim_time)s"), legend=:outerbottom, legend_column=2)
                        Plots.plot!(sim_window_height_plt, time_h, gt_h, label=L"h_{gt}(t)", legend=:outerbottom, legend_column=2)
                        
                        Plots.plot!(sim_window_η_norm_plt, iter, est_η_norm, label=string("Window - $(sim_time)s"), marker=1, legend=:outerbottom, legend_column=2)
                        Plots.plot!(sim_window_β_norm_plt, iter, est_β_norm, label=string("Window - $(sim_time)s"), marker=1, legend=:outerbottom, legend_column=2)
                        Plots.plot!(sim_window_ratio_norm_plt, iter, normalized_ratio, label=string("Window - $(sim_time)s"), marker=1, legend=:outerbottom, legend_column=2)
                        
                        try
                            slice_data = read_json(joinpath(exp_path,"Results","data","slice_data.json"))
                            t_steep = Float64.(collect(slice_data["steep"]["t"]))
                            zs_steep = Float64.(collect(slice_data["steep"]["zs"]))
                            t_flat = Float64.(collect(slice_data["flat"]["t"]))
                            zs_flat = Float64.(collect(slice_data["flat"]["zs"]))
                            
                            # 2D slices: cost vs distance along the two directions
                            if length(t_steep) > 0 && length(zs_steep) == length(t_steep)
                                Plots.plot!(plt_slices, t_steep, zs_steep, legend = :outerright, label = "Steepest direction; Window = $(sim_time)s")
                            else
                                @warn "Skipping steep slice plot: empty or mismatched lengths: $(length(t_steep)) vs $(length(zs_steep))"
                            end

                            if length(t_flat) > 0 && length(zs_flat) == length(t_flat)
                                Plots.plot!(plt_slices, t_flat,  zs_flat, legend = :outerright, label = "Flattest direction; Window = $(sim_time)s", legendfontsize=20)
                            else
                                @warn "Skipping flat slice plot: empty or mismatched lengths: $(length(t_flat)) vs $(length(zs_flat))"
                            end
                        catch e
                            @warn "Failed to read or plot slice data: $e"
                        end

                        if sim_time == 2.0
                            if ne == 6
                                Plots.plot!(η_norm_plot_2, iter, est_η_norm, label=latexstring("\$\\beta_{gt}:$(β_gt[1])\\mathrm{\\frac{Pa\\cdot s}{m}}\$"), marker=1, color=palette(:Set1)[(parse(Int, dir)-1) % length(palette(:Set1))+1])
                                Plots.xticks!(η_norm_plot_2, 0:2:(max_iter+1))

                                Plots.plot!(β_norm_plot_2, iter, est_β_norm, label=latexstring("\$\\beta_{gt}:$(β_gt[1])\\mathrm{\\frac{Pa\\cdot s}{m}}\$"), marker=1, color=palette(:Set1)[(parse(Int, dir)-1) % length(palette(:Set1))+1])
                                Plots.xticks!(β_norm_plot_2, 0:2:(max_iter+1))

                                Plots.plot!(plot_conv_2, iter, cost_list, label=latexstring("\$\\beta_{gt}:$(β_gt[1])\\mathrm{\\frac{Pa\\cdot s}{m}}\$"), marker=1, color=palette(:Set1)[(parse(Int, dir)-1) % length(palette(:Set1))+1])
                                Plots.xticks!(plot_conv_2, 0:2:(max_iter+1))
                                
                                Plots.plot!(plot_conv_log_2, iter, cost_list, label=latexstring("\$\\beta_{gt}:$(β_gt[1])\\mathrm{\\frac{Pa\\cdot s}{m}}\$"), marker=1, color=palette(:Set1)[(parse(Int, dir)-1) % length(palette(:Set1))+1], xscale=:log10, yscale=:log10)
                                Plots.xticks!(plot_conv_log_2, 1:10:max_iter)
                                
                                Plots.plot!(ratio_plot_2, iter, ratio_est, label=latexstring("\$\\beta_{gt}:$(β_gt[1])\\mathrm{\\frac{Pa\\cdot s}{m}}\$"), marker=1, color=palette(:Set1)[(parse(Int, dir)-1) % length(palette(:Set1))+1], yscale=:log10)
                                Plots.xticks!(ratio_plot_2, 0:2:(max_iter+1))
                                Plots.hline!(ratio_plot_2, [ratio_gt], linestyle=:dash, label=false, color=:black, yscale=:log10)

                                Plots.plot!(ratio_norm_plot_2, iter, normalized_ratio, label=latexstring("\$\\beta_{gt}:$(β_gt[1])\\mathrm{\\frac{Pa\\cdot s}{m}}\$"), marker=1, color=palette(:Set1)[(parse(Int, dir)-1) % length(palette(:Set1))+1])
                                Plots.xticks!(ratio_norm_plot_2, 0:2:(max_iter+1))
                                
                                Plots.plot!(rel_height_error_glob_plot_2, time_h, rel_height_error, label=latexstring("\$\\beta_{gt}:$(β_gt[1])\\mathrm{\\frac{Pa\\cdot s}{m}}\$"), color=palette(:Set1)[(parse(Int, dir)-1) % length(palette(:Set1))+1])
                                Plots.plot!(h_norm_plot_2, h_norm, label=latexstring("\$\\beta_{gt}:$(β_gt[1])\\mathrm{\\frac{Pa\\cdot s}{m}}\$"), color=palette(:Set1)[(parse(Int, dir)-1) % length(palette(:Set1))+1])
                            end
                        elseif sim_time == 5.0
                            if ne == 6
                                Plots.plot!(η_norm_plot_5, iter, est_η_norm, label=latexstring("\$\\beta_{gt}:$(β_gt[1])\\mathrm{\\frac{Pa\\cdot s}{m}}\$"), marker=1, color=palette(:Set1)[(parse(Int, dir)-1) % length(palette(:Set1))+1])
                                Plots.xticks!(η_norm_plot_5, 0:2:(max_iter+1))
                                
                                Plots.plot!(β_norm_plot_5, iter, est_β_norm, label=latexstring("\$\\beta_{gt}:$(β_gt[1])\\mathrm{\\frac{Pa\\cdot s}{m}}\$"), marker=1, color=palette(:Set1)[(parse(Int, dir)-1) % length(palette(:Set1))+1])
                                Plots.xticks!(β_norm_plot_5, 0:2:(max_iter+1))

                                Plots.plot!(plot_conv_5, iter, cost_list, label=latexstring("\$\\beta_{gt}:$(β_gt[1])\\mathrm{\\frac{Pa\\cdot s}{m}}\$"), marker=1, color=palette(:Set1)[(parse(Int, dir)-1) % length(palette(:Set1))+1])
                                Plots.xticks!(plot_conv_5, 0:2:(max_iter+1))

                                Plots.plot!(plot_conv_log_5, iter, cost_list, label=latexstring("\$\\beta_{gt}:$(β_gt[1])\\mathrm{\\frac{Pa\\cdot s}{m}}\$"), marker=1, color=palette(:Set1)[(parse(Int, dir)-1) % length(palette(:Set1))+1], xscale=:log10, yscale=:log10)
                                Plots.xticks!(plot_conv_log_5, 1:10:max_iter)

                                Plots.plot!(ratio_plot_5, iter, ratio_est, label=latexstring("\$\\beta_{gt}:$(β_gt[1])\\mathrm{\\frac{Pa\\cdot s}{m}}\$"), marker=1, color=palette(:Set1)[(parse(Int, dir)-1) % length(palette(:Set1))+1], yscale=:log10)
                                Plots.hline!(ratio_plot_5, [ratio_gt], linestyle=:dash, label=false, color=:black, yscale=:log10)
                                Plots.xticks!(ratio_plot_5, 0:2:(max_iter+1))

                                Plots.plot!(ratio_norm_plot_5, iter, normalized_ratio, label=latexstring("\$\\beta_{gt}:$(β_gt[1])\\mathrm{\\frac{Pa\\cdot s}{m}}\$"), marker=1, color=palette(:Set1)[(parse(Int, dir)-1) % length(palette(:Set1))+1])
                                Plots.xticks!(ratio_norm_plot_5, 0:2:(max_iter+1))

                                Plots.plot!(rel_height_error_glob_plot_5, time_h, rel_height_error, label=latexstring("\$\\beta_{gt}:$(β_gt[1])\\mathrm{\\frac{Pa\\cdot s}{m}}\$"), color=palette(:Set1)[(parse(Int, dir)-1) % length(palette(:Set1))+1])
                                Plots.plot!(h_norm_plot_5, time_h, h_norm, label=latexstring("\$\\beta_{gt}:$(β_gt[1])\\mathrm{\\frac{Pa\\cdot s}{m}}\$"), color=palette(:Set1)[(parse(Int, dir)-1) % length(palette(:Set1))+1])
                                Plots.plot!(h_glob_plot_5, time_h, h_norm, label=latexstring("\$\\beta_{gt}:$(β_gt[1])\\mathrm{\\frac{Pa\\cdot s}{m}}\$"), color=palette(:Set1)[(parse(Int, dir)-1) % length(palette(:Set1))+1])
                                Plots.plot!(h_glob_plot_5, time_h, gt_h, label=false, style=:dash, color=palette(:Set1)[(parse(Int, dir)-1) % length(palette(:Set1))+1])
                            
                                Plots.plot!(contour_plt, nSplinex[end], nSpliney[end], label="Ground truth contour", color=:red)
                                Plots.plot!(contour_plt_zoom, nSplinex[end], nSpliney[end], label=false, color=:red)
                            end
                        elseif sim_time == 10.0
                            if ne == 6
                                Plots.plot!(η_norm_plot_10, iter, est_η_norm, label=latexstring("\$\\beta_{gt}:$(β_gt[1])\\mathrm{\\frac{Pa\\cdot s}{m}}\$"), marker=1, color=palette(:Set1)[(parse(Int, dir)-1) % length(palette(:Set1))+1])
                                Plots.xticks!(η_norm_plot_10, 0:2:(max_iter+1))
                                
                                Plots.plot!(β_norm_plot_10, iter, est_β_norm, label=latexstring("\$\\beta_{gt}:$(β_gt[1])\\mathrm{\\frac{Pa\\cdot s}{m}}\$"), marker=1, color=palette(:Set1)[(parse(Int, dir)-1) % length(palette(:Set1))+1])
                                Plots.xticks!(β_norm_plot_10, 0:2:(max_iter+1))

                                Plots.plot!(plot_conv_10, iter, cost_list, label=latexstring("\$\\beta_{gt}:$(β_gt[1])\\mathrm{\\frac{Pa\\cdot s}{m}}\$"), marker=1, color=palette(:Set1)[(parse(Int, dir)-1) % length(palette(:Set1))+1])
                                Plots.xticks!(plot_conv_10, 0:2:(max_iter+1))
                                
                                Plots.plot!(plot_conv_log_10, iter, cost_list, label=latexstring("\$\\beta_{gt}:$(β_gt[1])\\mathrm{\\frac{Pa\\cdot s}{m}}\$"), marker=1, color=palette(:Set1)[(parse(Int, dir)-1) % length(palette(:Set1))+1], xscale=:log10, yscale=:log10)
                                Plots.xticks!(plot_conv_log_10, 1:10:max_iter)

                                Plots.plot!(ratio_plot_10, iter, ratio_est, label=latexstring("\$\\beta_{gt}:$(β_gt[1])\\mathrm{\\frac{Pa\\cdot s}{m}}\$"), marker=1, color=palette(:Set1)[(parse(Int, dir)-1) % length(palette(:Set1))+1], yscale=:log10)
                                Plots.hline!(ratio_plot_10, [ratio_gt], linestyle=:dash, label=false, color=:black, yscale=:log10)
                                Plots.xticks!(ratio_plot_10, 0:2:(max_iter+1))

                                Plots.plot!(ratio_norm_plot_10, iter, normalized_ratio, label=latexstring("\$\\beta_{gt}:$(β_gt[1])\\mathrm{\\frac{Pa\\cdot s}{m}}\$"), marker=1, color=palette(:Set1)[(parse(Int, dir)-1) % length(palette(:Set1))+1])
                                Plots.xticks!(ratio_norm_plot_10, 0:2:(max_iter+1))

                                Plots.plot!(rel_height_error_glob_plot_10, time_h, rel_height_error, label=latexstring("\$\\beta_{gt}:$(β_gt[1])\\mathrm{\\frac{Pa\\cdot s}{m}}\$"), color=palette(:Set1)[(parse(Int, dir)-1) % length(palette(:Set1))+1])
                                Plots.plot!(h_norm_plot_10, time_h, h_norm, label=latexstring("\$\\beta_{gt}:$(β_gt[1])\\mathrm{\\frac{Pa\\cdot s}{m}}\$"), color=palette(:Set1)[(parse(Int, dir)-1) % length(palette(:Set1))+1])
                                Plots.plot!(h_plot_10, time_h, h_norm, label=latexstring("\$\\beta_{gt}:$(β_gt[1])\\mathrm{\\frac{Pa\\cdot s}{m}}\$"), color=palette(:Set1)[(parse(Int, dir)-1) % length(palette(:Set1))+1])
                                Plots.plot!(h_plot_10, time_h, gt_h, label=false, style=:dash, color=palette(:Set1)[(parse(Int, dir)-1) % length(palette(:Set1))+1])
                            end
                        elseif sim_time == 20.0
                            if ne == 6
                                Plots.plot!(η_norm_plot_20, iter, est_η_norm, label=latexstring("\$\\beta_{gt}:$(β_gt[1])\\mathrm{\\frac{Pa\\cdot s}{m}}\$"), marker=1, color=palette(:Set1)[(parse(Int, dir)-1) % length(palette(:Set1))+1])
                                Plots.xticks!(η_norm_plot_20, 0:2:(max_iter+1))
                                
                                Plots.plot!(β_norm_plot_20, iter, est_β_norm, label=latexstring("\$\\beta_{gt}:$(β_gt[1])\\mathrm{\\frac{Pa\\cdot s}{m}}\$"), marker=1, color=palette(:Set1)[(parse(Int, dir)-1) % length(palette(:Set1))+1])
                                Plots.xticks!(β_norm_plot_20, 0:2:(max_iter+1))

                                Plots.plot!(plot_conv_20, iter, cost_list, label=latexstring("\$\\beta_{gt}:$(β_gt[1])\\mathrm{\\frac{Pa\\cdot s}{m}}\$"), marker=1, color=palette(:Set1)[(parse(Int, dir)-1) % length(palette(:Set1))+1])
                                Plots.xticks!(plot_conv_20, 0:2:(max_iter+1))
                                
                                Plots.plot!(plot_conv_log_20, iter, cost_list, label=latexstring("\$\\beta_{gt}:$(β_gt[1])\\mathrm{\\frac{Pa\\cdot s}{m}}\$"), marker=1, color=palette(:Set1)[(parse(Int, dir)-1) % length(palette(:Set1))+1], xscale=:log10, yscale=:log10)
                                Plots.xticks!(plot_conv_log_20, 1:10:max_iter)

                                Plots.plot!(ratio_plot_20, iter, ratio_est, label=latexstring("\$\\beta_{gt}:$(β_gt[1])\\mathrm{\\frac{Pa\\cdot s}{m}}\$"), marker=1, color=palette(:Set1)[(parse(Int, dir)-1) % length(palette(:Set1))+1], yscale=:log10)
                                Plots.hline!(ratio_plot_20, [ratio_gt], linestyle=:dash, label=false, color=:black, yscale=:log10)
                                Plots.xticks!(ratio_plot_20, 0:2:(max_iter+1))

                                Plots.plot!(ratio_norm_plot_20, iter, normalized_ratio, label=latexstring("\$\\beta_{gt}:$(β_gt[1])\\mathrm{\\frac{Pa\\cdot s}{m}}\$"), marker=1, color=palette(:Set1)[(parse(Int, dir)-1) % length(palette(:Set1))+1])
                                Plots.xticks!(ratio_norm_plot_20, 0:2:(max_iter+1))

                                Plots.plot!(rel_height_error_glob_plot_20, time_h, rel_height_error, label=string("Window - $(sim_time)s"," - ne: ",ne), legend=:outerbottom, legend_column=2)
                                Plots.plot!(h_norm_plot_20, time_h, h_norm, label=latexstring("\$\\beta_{gt}:$(β_gt[1])\\mathrm{\\frac{Pa\\cdot s}{m}}\$"), color=palette(:Set1)[(parse(Int, dir)-1) % length(palette(:Set1))+1])
                            end
                        elseif sim_time == 30.0
                            if ne == 6
                                Plots.plot!(η_norm_plot_30, iter, est_η_norm, label=latexstring("\$\\beta_{gt}:$(β_gt[1])\\mathrm{\\frac{Pa\\cdot s}{m}}\$"), marker=1, color=palette(:Set1)[(parse(Int, dir)-1) % length(palette(:Set1))+1])
                                Plots.xticks!(η_norm_plot_30, 0:2:(max_iter+1))
                                
                                Plots.plot!(β_norm_plot_30, iter, est_β_norm, label=latexstring("\$\\beta_{gt}:$(β_gt[1])\\mathrm{\\frac{Pa\\cdot s}{m}}\$"), marker=1, color=palette(:Set1)[(parse(Int, dir)-1) % length(palette(:Set1))+1])
                                Plots.xticks!(β_norm_plot_30, 0:2:(max_iter+1))

                                Plots.plot!(plot_conv_30, iter, cost_list, label=latexstring("\$\\beta_{gt}:$(β_gt[1])\\mathrm{\\frac{Pa\\cdot s}{m}}\$"), marker=1, color=palette(:Set1)[(parse(Int, dir)-1) % length(palette(:Set1))+1])
                                Plots.xticks!(plot_conv_30, 0:2:(max_iter+1))

                                Plots.plot!(plot_conv_log_30, iter, cost_list, label=latexstring("\$\\beta_{gt}:$(β_gt[1])\\mathrm{\\frac{Pa\\cdot s}{m}}\$"), marker=1, color=palette(:Set1)[(parse(Int, dir)-1) % length(palette(:Set1))+1], xscale=:log10, yscale=:log10)
                                Plots.xticks!(plot_conv_log_30, 1:10:max_iter)

                                Plots.plot!(ratio_plot_30, iter, ratio_est, label=latexstring("\$\\beta_{gt}:$(β_gt[1])\\mathrm{\\frac{Pa\\cdot s}{m}}\$"), marker=1, color=palette(:Set1)[(parse(Int, dir)-1) % length(palette(:Set1))+1], yscale=:log10)
                                Plots.hline!(ratio_plot_30, [ratio_gt], linestyle=:dash, label=false, color=:black, yscale=:log10)
                                Plots.xticks!(ratio_plot_30, 0:2:(max_iter+1))

                                Plots.plot!(ratio_norm_plot_30, iter, normalized_ratio, label=latexstring("\$\\beta_{gt}:$(β_gt[1])\\mathrm{\\frac{Pa\\cdot s}{m}}\$"), marker=1, color=palette(:Set1)[(parse(Int, dir)-1) % length(palette(:Set1))+1])
                                Plots.xticks!(ratio_norm_plot_30, 0:2:(max_iter+1))

                                Plots.plot!(rel_height_error_glob_plot_30, time_h, rel_height_error, label=latexstring("\$\\beta_{gt}:$(β_gt[1])\\mathrm{\\frac{Pa\\cdot s}{m}}\$"), legend=:outerbottom, legend_column=2)
                                Plots.plot!(h_norm_plot_30, time_h, h_norm, label=latexstring("\$\\beta_{gt}:$(β_gt[1])\\mathrm{\\frac{Pa\\cdot s}{m}}\$"), color=palette(:Set1)[(parse(Int, dir)-1) % length(palette(:Set1))+1])
                            end
                            println(size(height_error, 1), " vs ", size(gt_h), " vs ", size(time_h), " vs ", length(rel_height_error))
                            Plots.plot!(height_error_plt, time_h, height_error, label=string("Number of elements: ",ne), marker=1, legend=:outerbottom, legend_column=2)
                            Plots.plot!(rel_height_error_plt, time_h, rel_height_error, label=string("Number of elements: ",ne), marker=1, legend=:outerbottom, legend_column=2)
                            
                            Plots.plot!(elem_η_plt, iter, est_η, label=string("Number of elements: ",ne), marker=1, legend=:outerbottom, legend_column=2)
                            Plots.xticks!(elem_η_plt, 0:2:(max_iter+1))
                            
                            Plots.plot!(elem_β_plt, iter, est_β, label=string("Number of elements: ",ne), marker=1, legend=:outerbottom, legend_column=2)
                            Plots.xticks!(elem_β_plt, 0:2:(max_iter+1))
                            Plots.plot!(elem_ratio_plt, iter, ratio_est, label=string("Number of elements: ",ne), marker=1, legend=:outerbottom, legend_column=2)
                            Plots.xticks!(elem_ratio_plt, 0:2:(max_iter+1))
                            
                            Plots.plot!(elem_rel_height_error_plt, time_h, rel_height_error, label=string("Number of elements: ",ne), legend=:outerbottom, legend_column=2)
                            Plots.plot!(elem_height_plt, time_h, est_h, label=string("Number of elements: ",ne), legend=:outerbottom, legend_column=2)
                            Plots.plot!(elem_height_plt, time_h, gt_h, label=L"h_{gt}(t)", legend=:outerbottom, legend_column=2)
                        
                            Plots.plot!(elem_η_norm_plt, iter, est_η_norm, label=string("Number of elements: ",ne), marker=1, legend=:outerbottom, legend_column=2)
                            Plots.xticks!(elem_η_norm_plt, 0:2:(max_iter+1))
                            Plots.plot!(elem_β_norm_plt, iter, est_β_norm, label=string("Number of elements: ",ne), marker=1, legend=:outerbottom, legend_column=2)
                            Plots.xticks!(elem_β_norm_plt, 0:2:(max_iter+1))
                            Plots.plot!(elem_ratio_norm_plt, iter, normalized_ratio, label=string("Number of elements: ",ne), marker=1, legend=:outerbottom, legend_column=2)
                            Plots.xticks!(elem_ratio_norm_plt, 0:2:(max_iter+1))
                        end
                    else
                        η_pred = readdlm(joinpath(exp_path,"Results","data","eta_est.csv"), ',', Float64) # estimated η values per sample
                        β_pred = readdlm(joinpath(exp_path,"Results","data","beta_est.csv"), ',', Float64) # estimated β values per sample
                        h_pred = readdlm(joinpath(exp_path,"Results","data","h_est.csv"), ',', Float64) # estimated height values per sample

                        n_samples = size(η_pred, 2)
                        η_β_list = []
                        iter_list = []
                        for n in 1:n_samples
                            β = readdlm(joinpath(exp_path,"Results","data","opt_data","beta_steps","run_$n.csv"), ',', Float64)
                            η = readdlm(joinpath(exp_path,"Results","data","opt_data","eta_steps","run_$n.csv"), ',', Float64)
                            iter = readdlm(joinpath(exp_path,"Results","data","opt_data","iter","run_$n.csv"), ',', Float64)
                            ratio = η./β * (β_gt / η_gt)
                            push!(η_β_list, ratio[end])
                            push!(iter_list, iter)
                        end
                        plot_covariance!(covarience_plt, η_pred[:,1]./η_gt, β_pred[:,1]./β_gt, label=string(L"\sigma:\;",(round(noise_level,digits=2))," px  "), legend_column=noise_cols)

                        n_time = min(length(time), size(h_pred, 2), length(gt_h))
                        if n_time < length(time) || n_time < size(h_pred, 2) || n_time < length(gt_h)
                            @warn "Time and height vectors have mismatched lengths: time=$(length(time)), est_h=$(size(h_pred, 2)), gt_h=$(length(gt_h)). Truncating to $n_time samples for plotting."
                        end
                        height_error = abs.(h_pred[:,1:n_time]' .- gt_h[1:n_time])
                        normalized_height = h_pred[:,1:n_time]' ./ gt_h[1:n_time]
                        rel_height_error = abs.(h_pred[:,1:n_time]' .- gt_h[1:n_time]) ./ gt_h[1:n_time] * 100.0 # in percentage

                        time_h = copy(time)
                        if size(rel_height_error,1) != length(time) 
                            @warn "Length mismatch between height plot x-axis ($(size(rel_height_error,1))) and ground truth height ($(length(gt_h))). Adjusting..."
                            min_length = min(size(rel_height_error,1), length(time))
                            time_h = time[1:min_length] 
                            h_pred = h_pred[:,1:min_length]
                            gt_h = gt_h[1:min_length]
                        end
                        Plots.plot!(height_noise_plt, time_h, gt_h, label=L"h_{gt}", style=:dash)
                        StatsPlots.errorline!(height_noise_plt, time_h, h_pred', label=string(L"\sigma:\;",(round(noise_level,digits=2))," px  "))
                        StatsPlots.errorline!(height_error_noise_plt, time_h, height_error', label=string(L"\sigma:\;",(round(noise_level,digits=2))," px  "))
                        StatsPlots.errorline!(rel_height_error_noise_plt, time_h, rel_height_error', label=string(L"\sigma:\;",(round(noise_level,digits=2))," px  "))
                        StatsPlots.errorline!(normalized_height_noise_plt, time_h, normalized_height', label=string(L"\sigma:\;",(round(noise_level,digits=2))," px  "))
                    end
                end

                # Plots.scatter!(covarience_plt, η_gt, β_gt, label="Ground truth", ms=:15, m=:star5, color=def_orange, markerstrokewidth=0.1)
                Plots.hline!(covarience_plt, [1], linestyle=:dash, label=false, color=:black)
                Plots.vline!(covarience_plt, [1], linestyle=:dash, label=false, color=:black)

                plot_path_noise = joinpath(sim_time_folder,"post_analysis_noise","plots")
                set_file(plot_path_noise)
                Plots.savefig(covarience_plt, joinpath(plot_path_noise,"covariance_$dir.pdf"))
                Plots.savefig(height_noise_plt, joinpath(plot_path_noise,"height_$dir.pdf"))
                Plots.savefig(height_error_noise_plt, joinpath(plot_path_noise,"height_error_$dir.pdf"))
                Plots.savefig(rel_height_error_noise_plt, joinpath(plot_path_noise,"relative_height_error_$dir.pdf"))
                Plots.savefig(normalized_height_noise_plt, joinpath(plot_path_noise,"normalized_height_$dir.pdf"))
                @info "Saved plots to $plot_path_noise"
            end
            plot_path_sim_time = joinpath(elem_size_folder,"post_analysis_time","plots")
            set_file(plot_path_sim_time)
            # Plots.plot!(height_noise_plt, time, gt_h, label="Ground truth")

            Plots.savefig(sim_window_η_plt, joinpath(plot_path_sim_time,"η_$dir.pdf"))
            Plots.savefig(sim_window_β_plt, joinpath(plot_path_sim_time,"β_$dir.pdf"))
            Plots.savefig(sim_window_ratio_plt, joinpath(plot_path_sim_time,"η_β_ratio_$dir.pdf"))
            Plots.savefig(sim_window_rel_height_error_plt, joinpath(plot_path_sim_time,"relative_height_error_$dir.pdf"))
            Plots.savefig(sim_window_height_plt, joinpath(plot_path_sim_time,"height_$dir.pdf"))
            Plots.savefig(sim_window_η_norm_plt, joinpath(plot_path_sim_time,"η_normalized_$dir.pdf"))
            Plots.savefig(sim_window_β_norm_plt, joinpath(plot_path_sim_time,"β_normalized_$dir.pdf"))
            Plots.savefig(sim_window_ratio_norm_plt, joinpath(plot_path_sim_time,"η_β_ratio_normalized_$dir.pdf"))
            Plots.savefig(plt_slices, joinpath(plot_path_sim_time,"cost_slices_along_directions_$dir.pdf"))
            @info "Saved plots to $plot_path_sim_time"
        end
        
        plot_path_elems = joinpath(filepath_dir,"post_analysis","plots")
        set_file(plot_path_elems)

        Plots.savefig(elem_ratio_plt, joinpath(plot_path_elems,"η_β_ratio_$dir.pdf"))
        Plots.savefig(elem_η_plt, joinpath(plot_path_elems,"η_$dir.pdf"))
        Plots.savefig(elem_β_plt, joinpath(plot_path_elems,"β_$dir.pdf"))
        Plots.savefig(elem_ratio_norm_plt, joinpath(plot_path_elems,"η_β_ratio_normalized_$dir.pdf"))
        Plots.savefig(elem_η_norm_plt, joinpath(plot_path_elems,"η_normalized_$dir.pdf"))
        Plots.savefig(elem_β_norm_plt, joinpath(plot_path_elems,"β_normalized_$dir.pdf"))
        Plots.savefig(elem_rel_height_error_plt, joinpath(plot_path_elems,"relative_height_error_$dir.pdf"))
        Plots.savefig(elem_height_plt, joinpath(plot_path_elems,"height_$dir.pdf"))
        @info "Saved plots to $plot_path_elems"
    end
    plot_path_global = joinpath(filepath,"post_analysis_global","plots")
    set_file(plot_path_global)

    @info "Saving plots to $plot_path_global"
    Plots.savefig(plot_conv_2, joinpath(plot_path_global,"conv_2.pdf"))
    Plots.savefig(plot_conv_log_2, joinpath(plot_path_global,"conv_log_2.pdf"))
    Plots.savefig(η_norm_plot_2, joinpath(plot_path_global,"η_normalized_2.pdf"))
    Plots.savefig(β_norm_plot_2, joinpath(plot_path_global,"β_normalized_2.pdf"))
    Plots.savefig(ratio_plot_2, joinpath(plot_path_global,"η_β_ratio_2.pdf"))
    Plots.savefig(ratio_norm_plot_2, joinpath(plot_path_global,"η_β_ratio_normalized_2.pdf"))
    Plots.savefig(rel_height_error_glob_plot_2, joinpath(plot_path_global,"relative_height_error_2.pdf"))
    Plots.savefig(h_norm_plot_2, joinpath(plot_path_global,"h_normalized_2.pdf"))
    
    Plots.savefig(plot_conv_5, joinpath(plot_path_global,"conv_5.pdf"))
    Plots.savefig(plot_conv_log_5, joinpath(plot_path_global,"conv_log_5.pdf"))
    Plots.savefig(η_norm_plot_5, joinpath(plot_path_global,"η_normalized_5.pdf"))
    Plots.savefig(β_norm_plot_5, joinpath(plot_path_global,"β_normalized_5.pdf"))
    Plots.savefig(ratio_plot_5, joinpath(plot_path_global,"η_β_ratio_5.pdf"))
    Plots.savefig(ratio_norm_plot_5, joinpath(plot_path_global,"η_β_ratio_normalized_5.pdf"))
    Plots.savefig(rel_height_error_glob_plot_5, joinpath(plot_path_global,"relative_height_error_5.pdf"))
    Plots.savefig(h_norm_plot_5, joinpath(plot_path_global,"h_normalized_5.pdf"))
    Plots.savefig(h_glob_plot_5, joinpath(plot_path_global,"height_comparison_5.pdf"))
    Plots.savefig(contour_plt, joinpath(plot_path_global,"contour_comparison_5.pdf"))
    Plots.savefig(contour_plt_zoom, joinpath(plot_path_global,"contour_comparison_5_zoomed.pdf"))
    
    Plots.savefig(plot_conv_10, joinpath(plot_path_global,"conv_10.pdf"))
    Plots.savefig(plot_conv_log_10, joinpath(plot_path_global,"conv_log_10.pdf"))
    Plots.savefig(η_norm_plot_10, joinpath(plot_path_global,"η_normalized_10.pdf"))
    Plots.savefig(β_norm_plot_10, joinpath(plot_path_global,"β_normalized_10.pdf"))
    Plots.savefig(ratio_plot_10, joinpath(plot_path_global,"η_β_ratio_10.pdf"))
    Plots.savefig(ratio_norm_plot_10, joinpath(plot_path_global,"η_β_ratio_normalized_10.pdf"))
    Plots.savefig(rel_height_error_glob_plot_10, joinpath(plot_path_global,"relative_height_error_10.pdf"))
    Plots.savefig(h_glob_plot_10, joinpath(plot_path_global,"height_comparison_10.pdf"))
    
    Plots.savefig(plot_conv_20, joinpath(plot_path_global,"conv_20.pdf"))
    Plots.savefig(plot_conv_log_20, joinpath(plot_path_global,"conv_log_20.pdf"))
    Plots.savefig(η_norm_plot_20, joinpath(plot_path_global,"η_normalized_20.pdf"))
    Plots.savefig(β_norm_plot_20, joinpath(plot_path_global,"β_normalized_20.pdf"))
    Plots.savefig(ratio_plot_20, joinpath(plot_path_global,"η_β_ratio_20.pdf"))
    Plots.savefig(ratio_norm_plot_20, joinpath(plot_path_global,"η_β_ratio_normalized_20.pdf"))
    Plots.savefig(rel_height_error_glob_plot_20, joinpath(plot_path_global,"relative_height_error_20.pdf"))
    Plots.savefig(h_norm_plot_20, joinpath(plot_path_global,"h_normalized_20.pdf"))
    
    Plots.savefig(plot_conv_30, joinpath(plot_path_global,"conv_30.pdf"))
    Plots.savefig(plot_conv_log_30, joinpath(plot_path_global,"conv_log_30.pdf"))
    Plots.savefig(η_norm_plot_30, joinpath(plot_path_global,"η_normalized_30.pdf"))
    Plots.savefig(β_norm_plot_30, joinpath(plot_path_global,"β_normalized_30.pdf"))
    Plots.savefig(ratio_plot_30, joinpath(plot_path_global,"η_β_ratio_30.pdf"))
    Plots.savefig(ratio_norm_plot_30, joinpath(plot_path_global,"η_β_ratio_normalized_30.pdf"))
    Plots.savefig(rel_height_error_glob_plot_30, joinpath(plot_path_global,"relative_height_error_30.pdf"))
    Plots.savefig(h_norm_plot_30, joinpath(plot_path_global,"h_normalized_30.pdf"))
    
    @info "Saved plots to $plot_path_global"
end

function post_analysis_bulk(filepath_gt_::String, filepath::String, avoid_list)
    dir_list = readdir(filepath)

    # bulk_right_margin = 16pt
    # bulk_lft_margin = -3pt

    # for 1/2 linewidth
    # fs::Int = 12
    # bulk_plt_height::Int = 350
    # bulk_plt_width ::Int = 477
    # bulk_lft_margin = 1pt
    # bulk_right_margin = 5pt
    # bulk_top_margin = 1pt

    # for 1/3 linewidth
    # fs::Int = 12
    # plt_height::Int = 255
    # plt_width::Int = 318
    # plt_lft_margin = -5pt
    # plt_right_margin = -1pt
    # plt_top_margin = -5pt
    
    # figures to compare convergence with slip levels
    η_norm_plot_5 = set_plot(fs, sz=(plt_width, plt_height), legend_column=3, right_margin=plt_right_margin, left_margin=plt_lft_margin)
    Plots.hline!(η_norm_plot_5, [1.0],  linestyle=:dash, label=false, color=:black)
    Plots.xlabel!(η_norm_plot_5,L"\mathrm{Iterations}")
    Plots.ylabel!(η_norm_plot_5,L"\eta_{est}/\eta_{avg}")

    β_norm_plot_5 = set_plot(fs, sz=(plt_width, plt_height), legend_column=3, right_margin=plt_right_margin, left_margin=plt_lft_margin)
    Plots.hline!(β_norm_plot_5, [1.0],  linestyle=:dash, label=false, color=:black)
    Plots.xlabel!(β_norm_plot_5, L"\mathrm{Iterations}")
    Plots.ylabel!(β_norm_plot_5, L"\beta_{est}/\beta_{gt}")

    η_norm_plot_10 = set_plot(fs, sz=(plt_width, plt_height), legend_column=3, right_margin=plt_right_margin, left_margin=plt_lft_margin)
    Plots.hline!(η_norm_plot_10, [1.0],  linestyle=:dash, label=false, color=:black)
    Plots.xlabel!(η_norm_plot_10,L"\mathrm{Iterations}")
    Plots.ylabel!(η_norm_plot_10,L"\eta_{est}/\eta_{avg}")

    β_norm_plot_10 = set_plot(fs, sz=(plt_width, plt_height), legend_column=3, right_margin=plt_right_margin, left_margin=plt_lft_margin)
    Plots.hline!(β_norm_plot_10, [1.0],  linestyle=:dash, label=false, color=:black)
    Plots.xlabel!(β_norm_plot_10, L"\mathrm{Iterations}")
    Plots.ylabel!(β_norm_plot_10, L"\beta_{est}/\beta_{gt}")

    # height plots
    h_plot_5 = set_plot(fs, sz=(plt_width, plt_height), legend_column=3, right_margin=plt_right_margin, left_margin=plt_lft_margin)
    Plots.xlabel!(h_plot_5,L"\mathrm{Time\;(s)}")
    Plots.ylabel!(h_plot_5,L"h_{est}\;\mathrm{(mm)}")  
    Plots.xlims!(h_plot_5, 0, end_obs_win)

    h_plot_10 = set_plot(fs, sz=(plt_width, plt_height), legend_column=3, right_margin=plt_right_margin, left_margin=plt_lft_margin)
    Plots.xlabel!(h_plot_10,L"\mathrm{Time\;(s)}")
    Plots.ylabel!(h_plot_10,L"h_{est}\;\mathrm{(mm)}") 
    Plots.xlims!(h_plot_10, 0, end_obs_win)

    # normalised height plots
    h_norm_plot_5 = set_plot(fs, sz=(plt_width, plt_height), legend_column=3, right_margin=plt_right_margin, left_margin=plt_lft_margin)
    Plots.hline!(h_norm_plot_5, [1.0],  linestyle=:dash, label=false, color=:black)
    Plots.xlabel!(h_norm_plot_5,L"\mathrm{Time\;(s)}")
    Plots.ylabel!(h_norm_plot_5,L"h_{est}/h_{gt}")  
    Plots.xlims!(h_norm_plot_5, 0, end_obs_win)
    Plots.ylims!(h_norm_plot_5, y_lims_h_norm)

    h_norm_plot_10 = set_plot(fs, sz=(plt_width, plt_height), legend_column=3, right_margin=plt_right_margin, left_margin=plt_lft_margin)
    Plots.hline!(h_norm_plot_10, [1.0],  linestyle=:dash, label=false, color=:black)
    Plots.xlabel!(h_norm_plot_10,L"\mathrm{Time\;(s)}")
    Plots.ylabel!(h_norm_plot_10,L"h_{est}/h_{gt}") 
    Plots.xlims!(h_norm_plot_10, 0, end_obs_win)
    Plots.ylims!(h_norm_plot_10, y_lims_h_norm)

    rel_height_error_glob_plot_5 = set_plot(fs, sz=(plt_width, plt_height), legend_column=3, right_margin=plt_right_margin, left_margin=plt_lft_margin)
    Plots.xlabel!(rel_height_error_glob_plot_5, L"\mathrm{Time\;(s)}")
    Plots.ylabel!(rel_height_error_glob_plot_5, latexstring("Relative Height Error (\$\\%\$)"))
    Plots.xlims!(rel_height_error_glob_plot_5, 0, end_obs_win)
    Plots.ylims!(rel_height_error_glob_plot_5, y_lims_rel_error)

    rel_height_error_glob_plot_10 = set_plot(fs, sz=(plt_width, plt_height), legend_column=3, right_margin=plt_right_margin, left_margin=plt_lft_margin)
    Plots.xlabel!(rel_height_error_glob_plot_10, L"\mathrm{Time\;(s)}")
    Plots.ylabel!(rel_height_error_glob_plot_10, latexstring("Relative Height Error (\$\\%\$)"))
    Plots.xlims!(rel_height_error_glob_plot_10, 0, end_obs_win)
    Plots.ylims!(rel_height_error_glob_plot_10, y_lims_rel_error)
    
    # η\β ratio plot
    ratio_plot_5 = set_plot(fs, sz=(plt_width, plt_height), legend_column=3, right_margin=plt_right_margin, left_margin=plt_lft_margin)
    Plots.hline!(ratio_plot_5, [1.0],  linestyle=:dash, label=false)
    Plots.xlabel!(ratio_plot_5,L"\mathrm{Iterations}")
    Plots.ylabel!(ratio_plot_5,L"\eta_{est}/\eta_{gt}")

    ratio_plot_10 = set_plot(fs, sz=(plt_width, plt_height), legend_column=3, right_margin=plt_right_margin, left_margin=plt_lft_margin)
    Plots.hline!(ratio_plot_10, [1.0],  linestyle=:dash, label=false)
    Plots.xlabel!(ratio_plot_10,L"\mathrm{Iterations}")
    Plots.ylabel!(ratio_plot_10,L"\eta_{est}/\eta_{gt}")

    # η \ β * β_gt \ η_gt plot
    ratio_norm_plot_5 = set_plot(fs, sz=(plt_width, plt_height), legend_column=3, right_margin=plt_right_margin, left_margin=plt_lft_margin)
    Plots.hline!(ratio_norm_plot_5, [1.0],  linestyle=:dash, label=false)
    Plots.xlabel!(ratio_norm_plot_5,L"\mathrm{Iterations}")
    Plots.ylabel!(ratio_norm_plot_5,L"\frac{\eta_{est}/\beta_{est}}{\eta_{gt}/\beta_{gt}}")

    ratio_norm_plot_10 = set_plot(fs, sz=(plt_width, plt_height), legend_column=3, right_margin=plt_right_margin, left_margin=plt_lft_margin)
    Plots.hline!(ratio_norm_plot_10, [1.0],  linestyle=:dash, label=false)
    Plots.xlabel!(ratio_norm_plot_10,L"\mathrm{Iterations}")
    Plots.ylabel!(ratio_norm_plot_10,L"\frac{\eta_{est}/\beta_{est}}{\eta_{gt}/\beta_{gt}}")
    
    for dir in dir_list
        if dir in avoid_list || dir == "post_analysis_global"
            println("Skipping directory: ", dir)
            continue
        end
        filepath_dir = joinpath(filepath, dir)
        filepath_gt = joinpath(filepath_gt_, dir)

        printstyled("Processing directory: $(filepath_dir)\n", color=:green)
        sim_params = read_json(joinpath(filepath_gt,"data","sim_params.json"))

        η_gt = sim_params["η"]
        β_gt = sim_params["β"]  
        gt_h_ = readdlm(joinpath(filepath_gt,"data","h.csv"), ',', Float64)

        sim_time = sim_params["simulation_time"]    
        t_steps = sim_params["time_steps"]

        time = collect(Float64, range(start=0, stop=sim_time, step=t_steps))

        elem_size_folders = readdir(filepath_dir)

        # plots for element vise comparison
        # plots for convergence 
        elem_η_plt = set_plot(fs, sz=(plt_width, plt_height))
        Plots.hline!(elem_η_plt, [η_gt[1]], label="Ground truth η",  left_margin=plt_lft_margin)
        Plots.xlabel!(elem_η_plt,L"\mathrm{Iterations}")
        Plots.ylabel!(elem_η_plt,L"\eta\;\mathrm{(KPa\cdot s)}")

        elem_β_plt = set_plot(fs, sz=(plt_width, plt_height))
        Plots.hline!(elem_β_plt, [β_gt[1]], label="Ground truth β",  left_margin=plt_lft_margin)
        Plots.xlabel!(elem_β_plt, L"\mathrm{Iterations}")
        Plots.ylabel!(elem_β_plt, L"\beta\;\mathrm{(mm^{-1})}")

        elem_ratio_plt = set_plot(fs, sz=(plt_width, plt_height))
        Plots.hline!(elem_ratio_plt, [η_gt[1]/β_gt[1]], label="Ground truth η/β",  left_margin=plt_lft_margin)
        Plots.xlabel!(elem_ratio_plt, L"\mathrm{Iterations}")   
        Plots.ylabel!(elem_ratio_plt, L"\eta/\beta\;(KPa\cdot s \cdot mm^{-1})")
        
        # plots for normalized values
        elem_η_norm_plt = set_plot(fs, sz=(plt_width, plt_height))
        Plots.hline!(elem_η_norm_plt, [1.0], label=false,  left_margin=plt_lft_margin, linestyle=:dash)
        Plots.xlabel!(elem_η_norm_plt,L"\mathrm{Iterations}")
        Plots.ylabel!(elem_η_norm_plt,L"\eta_{est}/\eta_{gt}")

        elem_β_norm_plt = set_plot(fs, sz=(plt_width, plt_height))
        Plots.hline!(elem_β_norm_plt, [1.0], label=false,  left_margin=plt_lft_margin, linestyle=:dash)
        Plots.xlabel!(elem_β_norm_plt, L"\mathrm{Iterations}")
        Plots.ylabel!(elem_β_norm_plt, L"\beta_{est}/\beta_{gt}")

        elem_ratio_norm_plt = set_plot(fs, sz=(plt_width, plt_height))
        Plots.hline!(elem_ratio_norm_plt, [1.0], label=false,  left_margin=plt_lft_margin, linestyle=:dash)
        Plots.xlabel!(elem_ratio_norm_plt, L"\mathrm{Iterations}")   
        Plots.ylabel!(elem_ratio_norm_plt, L"\frac{\eta_{est}/\beta_{est}}{\eta_{gt}/\beta_{gt}}")

        # plots for height error
        elem_rel_height_error_plt = set_plot(fs, sz=(plt_width, plt_height))
        Plots.plot!(elem_rel_height_error_plt, [], left_margin=plt_lft_margin, label=false)
        Plots.xlabel!(elem_rel_height_error_plt, L"\mathrm{Time\;(s)}")
        Plots.ylabel!(elem_rel_height_error_plt, L"\mathrm{Relative\;Height\;Error}")
        Plots.xlims!(elem_rel_height_error_plt, 0, end_obs_win)

        elem_height_plt = set_plot(fs, sz=(plt_width, plt_height))
        Plots.plot!(elem_height_plt, [], left_margin=plt_lft_margin, label=false)
        Plots.xlabel!(elem_height_plt, L"\mathrm{Time\;(s)}")
        Plots.ylabel!(elem_height_plt, L"\mathrm{Height\;(mm)}")
        Plots.xlims!(elem_height_plt, 0, end_obs_win)
        
        for elem_size_folder_ in elem_size_folders
            if elem_size_folder_ == "post_analysis" || elem_size_folder_ == "Q2_2" || elem_size_folder_ == "Q2_4" || elem_size_folder_ == "Q2_8"
                continue
            end
            
            elem_size_folder = joinpath(filepath_dir, elem_size_folder_)
            printstyled("Processing element size folder: $(elem_size_folder)\n", color=:blue)  
            sim_time_folders = readdir(elem_size_folder)
            
            # figures for Simulation window vise camparison
            sim_window_η_plt = set_plot(fs, sz=(plt_width, plt_height))
            Plots.hline!(sim_window_η_plt, [η_gt[1]], label="Ground truth η")
            Plots.xlabel!(sim_window_η_plt,L"\mathrm{Iterations}")
            Plots.ylabel!(sim_window_η_plt,L"\eta\;\mathrm{(KPa\cdot s)}")
            
            sim_window_β_plt = set_plot(fs, sz=(plt_width, plt_height))
            Plots.hline!(sim_window_β_plt, [β_gt[1]], label="Ground truth β")
            Plots.xlabel!(sim_window_β_plt, L"\mathrm{Iterations}")
            Plots.ylabel!(sim_window_β_plt, L"\beta\;\mathrm{(mm^{-1})}")

            sim_window_ratio_plt = set_plot(fs, sz=(plt_width, plt_height))
            Plots.hline!(sim_window_ratio_plt, [η_gt[1]/β_gt[1]], label="Ground truth η/β",  left_margin=plt_lft_margin)
            Plots.xlabel!(sim_window_ratio_plt, L"\mathrm{Iterations}")   
            Plots.ylabel!(sim_window_ratio_plt, L"\eta/\beta\;(KPa\cdot s \cdot mm^{-1})")

            # plots for normalized values
            sim_window_η_norm_plt = set_plot(fs, sz=(plt_width, plt_height))
            Plots.hline!(sim_window_η_norm_plt, [1.0], label=false,  left_margin=plt_lft_margin, linestyle=:dash)
            Plots.xlabel!(sim_window_η_norm_plt,L"\mathrm{Iterations}")
            Plots.ylabel!(sim_window_η_norm_plt,L"\eta_{est}/\eta_{gt}")

            sim_window_β_norm_plt = set_plot(fs, sz=(plt_width, plt_height))
            Plots.hline!(sim_window_β_norm_plt, [1.0], label=false,  left_margin=plt_lft_margin, linestyle=:dash)
            Plots.xlabel!(sim_window_β_norm_plt, L"\mathrm{Iterations}")
            Plots.ylabel!(sim_window_β_norm_plt, L"\beta_{est}/\beta_{gt}")

            sim_window_ratio_norm_plt = set_plot(fs, sz=(plt_width, plt_height))
            Plots.hline!(sim_window_ratio_norm_plt, [1.0], label=false,  left_margin=plt_lft_margin, linestyle=:dash)
            Plots.xlabel!(sim_window_ratio_norm_plt, L"\mathrm{Iterations}")   
            Plots.ylabel!(sim_window_ratio_norm_plt, L"\frac{\eta_{est}/\beta_{est}}{\eta_{gt}/\beta_{gt}}")

            # plots for height error
            sim_window_rel_height_error_plt = set_plot(fs, sz=(plt_width, plt_height))
            Plots.plot!(sim_window_rel_height_error_plt, [], left_margin=plt_lft_margin, label=false)
            Plots.xlabel!(sim_window_rel_height_error_plt, L"\mathrm{Time\;(s)}")
            Plots.ylabel!(sim_window_rel_height_error_plt, L"\mathrm{Relative\;Height\;Error}")
            Plots.xlims!(sim_window_rel_height_error_plt, 0, end_obs_win)

            sim_window_height_plt = set_plot(fs, sz=(plt_width, plt_height))
            Plots.plot!(sim_window_height_plt, [], left_margin=plt_lft_margin, label=false)
            Plots.xlabel!(sim_window_height_plt, L"\mathrm{Time\;(s)}")
            Plots.ylabel!(sim_window_height_plt, L"\mathrm{Height\;(mm)}")
            Plots.xlims!(sim_window_height_plt, 0, end_obs_win)

            plt_slices = set_plot(fs, sz=(350, 750))
            Plots.vline!(plt_slices, [0.0], color=:blue, linestyle=:dash, label="Minimum",  lw=1)
            Plots.xlabel!(plt_slices, L"\mathrm{Distance\;from\;minimum\;(px)}")
            Plots.ylabel!(plt_slices, L"\mathrm{Cost}")
            Plots.ylims!(plt_slices, 0, 50)

            for sim_time_folder_ in sim_time_folders

                if sim_time_folder_ == "post_analysis_time" || sim_time_folder_ == "Results" || sim_time_folder_ == "simtime_2.0" || sim_time_folder_ == "simtime_30.0"
                    continue
                end
                
                height_error_plt = set_plot(fs, sz=(plt_width, plt_height))
                Plots.plot!(height_error_plt, [], label=false)
                Plots.xlabel!(height_error_plt, L"\mathrm{Time\;(s)}")
                Plots.ylabel!(height_error_plt, L"\mathrm{Height\;Error\;(mm)}")
                Plots.xlims!(height_error_plt, 0, end_obs_win)
                
                rel_height_error_plt = set_plot(fs, sz=(plt_width, plt_height))
                Plots.plot!(rel_height_error_plt, [], label=false)
                Plots.xlabel!(rel_height_error_plt, L"\mathrm{Time\;(s)}")
                Plots.ylabel!(rel_height_error_plt, L"\mathrm{Height\;Error\;(mm)}")
                Plots.xlims!(rel_height_error_plt, 0, end_obs_win)

                sim_time_folder = joinpath(elem_size_folder, sim_time_folder_, "noise_0.0")

                printstyled("Processing simulation time folder: $(sim_time_folder)\n", color=:cyan)

                window_dirs = readdir(sim_time_folder)
                for window_dir in window_dirs
                    if window_dir == "Results" || window_dir == "post_analysis_window" || window_dir == "single_window" || window_dir == "post_analysis_noise"
                        println("Skipping directory: $window_dir")
                        continue
                    end
                    win_exp_path = joinpath(filepath, elem_size_folder, sim_time_folder, window_dir)
                    
                    println("Processing window: $win_exp_path")
                    exp_params = read_json(joinpath(win_exp_path ,"Results","data","experiment_parameters.json"))
                    sim_time_exp = exp_params["sim_time_exp"]
                    data_type = exp_params["data_type"]
                    noise_level = exp_params["noise_level"]
                    ne = exp_params["ne_exp"]

                    if data_type  == "synthetic"
                        ObsDataList, splinexObs, splineyObs = read_csv(joinpath(filepath_gt,"data","img_data","contour_data"))  
                        @info "Data type $data_type Reading synthetic contour data of $(length(ObsDataList)) time steps"
                    elseif data_type == "simulated"
                        ObsDataList, splinexObs, splineyObs = read_csv(joinpath(filepath_gt,"data","sim_data","contour_data"))  
                        @info "Data type $data_type Reading simulated contour data from $(joinpath(filepath_gt,"data","sim_data","contour_data")) of $(length(ObsDataList)) time steps"
                    else
                        error("Unknown data type: $data_type")
                    end

                    obsBorderPts, nSplinex, nSpliney, pd = add_noise(ObsDataList, nFactor=0.0)

                    local η_gt::Vector{Float64} = float.(sim_params["η"])
                    local β_gt  = sim_params["β"]

                    est_ηpList = readdlm(joinpath(win_exp_path,"Results","data","est_η.csv"), ',', Float64)
                    est_βpList = readdlm(joinpath(win_exp_path,"Results","data","est_β.csv"), ',', Float64)
                    avg_ηList = readdlm(joinpath(win_exp_path,"Results","data","avg_η.csv"), ',', Float64)
                    est_h_list = readdlm(joinpath(win_exp_path,"Results","data","est_h.csv"), ',', Float64)
                    
                    data_ranges_ = get_time_windows(joinpath(win_exp_path,"Results","data","window_data","data_ranges.csv"))
                    t_windows = readdlm(joinpath(win_exp_path,"Results","data","window_data","t_windows.csv"),',',Float64)
                    time_windows = readdlm(joinpath(win_exp_path,"Results","data","window_data","time_windows.csv"),',',Float64)

                    symBorderPts, splinex, spliney = read_csv(joinpath(win_exp_path,"Results","data","sim_data","2D_border_points"))
                    
                    println("Time windows: $(time_windows)")
                    obs_time = sum(time_windows)
                    
                    if obs_time != sim_time
                        @warn "Observation time frame $obs_time is less than preset ground truth time frame $sim_time, updating time frame"
                        min_time = min(obs_time, sim_time)
                        sim_time = min_time
                        obs_time = min_time
                    end
                    
                    if obs_time < sim_time_exp
                        @warn "Observation time frame $obs_time is less than experimental simulation time frame $sim_time_exp, switching to observation time frame"
                        sim_time_exp = obs_time
                    end

                    data_point_len = round(Int, obs_time/t_steps)
                    η_gt = η_gt[1:data_point_len]
                    gt_h = gt_h_[1:(data_point_len+1)]
                    symBorderPts = symBorderPts[1:data_point_len+1, :]
                    obsBorderPts = obsBorderPts[1:data_point_len+1, :]
                    est_h = est_h_list[1:data_point_len+1, :]

                    height_error = abs.(est_h - gt_h)
                    rel_height_error = height_error ./ gt_h * 100.0

                    # normalize by ground-truth     
                    η_gt = η_gt                   
                    est_η_norm = est_ηpList./avg_ηList
                    est_β_norm = est_βpList./β_gt

                    ratio_est = est_ηpList ./ est_βpList
                    ratio_gt = η_gt / β_gt
                    normalized_ratio = (est_ηpList ./ est_βpList) * (β_gt / η_gt)
            
                    Plots.plot!(sim_window_η_plt, est_ηpList, label=string("Window - $(sim_time)s"), marker=1, legend=:outerbottom, legend_column=2)
                    Plots.plot!(sim_window_β_plt, est_βpList, label=string("Window - $(sim_time)s"), marker=1, legend=:outerbottom, legend_column=2)
                    Plots.plot!(sim_window_ratio_plt, ratio_est, label=string("Window - $(sim_time)s"), marker=1, legend=:outerbottom, legend_column=2)
                    
                    Plots.plot!(sim_window_rel_height_error_plt, time, rel_height_error, label=string("Window - $(sim_time)s"," - ne: ",ne), legend=:outerbottom, legend_column=2)
                    Plots.plot!(sim_window_height_plt, time, est_h, label=string("Window - $(sim_time)s"), legend=:outerbottom, legend_column=2)
                    Plots.plot!(sim_window_height_plt, time, gt_h, label=L"h_{gt}(t)", legend=:outerbottom, legend_column=2)
                    
                    Plots.plot!(sim_window_η_norm_plt, est_η_norm, label=string("Window - $(sim_time)s"), marker=1, legend=:outerbottom, legend_column=2)
                    Plots.plot!(sim_window_β_norm_plt, est_β_norm, label=string("Window - $(sim_time)s"), marker=1, legend=:outerbottom, legend_column=2)
                    Plots.plot!(sim_window_ratio_norm_plt, normalized_ratio, label=string("Window - $(sim_time)s"), marker=1, legend=:outerbottom, legend_column=2)

                    if sim_time_exp == 5.0
                        if ne == 6
                            t_prev = 0.1
                            for ti::Int in 1:(size(data_ranges_, 1)-1)
                                t = t_windows[ti]
                                data_range_ = data_ranges_[ti]
                                t_win = collect(range(start=t_prev, stop=t, step=t_steps))

                                Plots.vline!(η_norm_plot_5, [t], color=:gray, linestyle=:dash, label=false)
                                Plots.vline!(β_norm_plot_5, [t], color=:gray, linestyle=:dash, label=false)
                                Plots.vline!(ratio_plot_5, [t], color=:gray, linestyle=:dash, label=false)
                                Plots.vline!(ratio_norm_plot_5, [t], color=:gray, linestyle=:dash, label=false)

                                if ti == 1
                                    Plots.plot!(η_norm_plot_5, t_win, est_η_norm[data_range_], label=latexstring("\$\\beta_{gt}:$(β_gt[1])\\mathrm{\\frac{Pa\\cdot s}{m}}\$"), color=palette(:Set1)[(parse(Int, dir)-1) % length(palette(:Set1))+1])
                                    Plots.plot!(β_norm_plot_5, t_win, est_β_norm[data_range_], label=latexstring("\$\\beta_{gt}:$(β_gt[1])\\mathrm{\\frac{Pa\\cdot s}{m}}\$"), color=palette(:Set1)[(parse(Int, dir)-1) % length(palette(:Set1))+1])
                                    Plots.plot!(ratio_plot_5, t_win, ratio_est[data_range_], label=latexstring("\$\\beta_{gt}:$(β_gt[1])\\mathrm{\\frac{Pa\\cdot s}{m}}\$"), yscale=:log10, color=palette(:Set1)[(parse(Int, dir)-1) % length(palette(:Set1))+1])
                                    Plots.plot!(ratio_norm_plot_5, t_win, normalized_ratio[data_range_], label=latexstring("\$\\beta_{gt}:$(β_gt[1])\\mathrm{\\frac{Pa\\cdot s}{m}}\$"), color=palette(:Set1)[(parse(Int, dir)-1) % length(palette(:Set1))+1])
                                else
                                    Plots.plot!(η_norm_plot_5, t_win, est_η_norm[data_range_], label=false, color=palette(:Set1)[(parse(Int, dir)-1) % length(palette(:Set1))+1])
                                    Plots.plot!(β_norm_plot_5, t_win, est_β_norm[data_range_], label=false, color=palette(:Set1)[(parse(Int, dir)-1) % length(palette(:Set1))+1])
                                    Plots.plot!(ratio_plot_5, t_win, ratio_est[data_range_], label=false, yscale=:log10, color=palette(:Set1)[(parse(Int, dir)-1) % length(palette(:Set1))+1])
                                    Plots.plot!(ratio_norm_plot_5, t_win, normalized_ratio[data_range_], label=false, color=palette(:Set1)[(parse(Int, dir)-1) % length(palette(:Set1))+1])
                                end
                                t_prev = t+t_steps
                            end

                            Plots.plot!(h_plot_5, time, est_h, label=latexstring("\$\\beta_{gt}:$(β_gt[1])\\mathrm{\\frac{Pa\\cdot s}{m}}\$"), color=palette(:Set1)[(parse(Int, dir)-1) % length(palette(:Set1))+1])
                            Plots.plot!(h_plot_5, time, gt_h, label=false, linestyle=:dash, color=palette(:Set1)[(parse(Int, dir)-1) % length(palette(:Set1))+1])

                            t = collect(range(start=0, stop=sim_time, step=t_steps))
                            # Plots.plot!(ratio_plot_5, ratio_est, label=latexstring("\$\\beta_{gt}:$(β_gt[1])\\mathrm{\\frac{Pa\\cdot s}{m}}\$"), color=palette(:Set1)[(parse(Int, dir)-1) % length(palette(:Set1))+1])
                            # Plots.hline!(ratio_plot_5, [ratio_gt], linestyle=:dash, label=false, color=:black, yscale=:log10)
                            Plots.plot!(rel_height_error_glob_plot_5, t, rel_height_error, label=latexstring("\$\\beta_{gt}:$(β_gt[1])\\mathrm{\\frac{Pa\\cdot s}{m}}\$"), color=palette(:Set1)[(parse(Int, dir)-1) % length(palette(:Set1))+1])
                            Plots.plot!(h_norm_plot_5, t, est_h./gt_h, label=latexstring("\$\\beta_{gt}:$(β_gt[1])\\mathrm{\\frac{Pa\\cdot s}{m}}\$"), color=palette(:Set1)[(parse(Int, dir)-1) % length(palette(:Set1))+1])
                            if dir == "1"
                                for ti::Int in 1:(size(data_ranges_, 1)-1)
                                    t_win = t_windows[ti]
                                    data_range_ = data_ranges_[ti]

                                    Plots.vline!(h_norm_plot_5, [t_win], color=:gray, linestyle=:dash, label=false)
                                end
                            end
                            
                            # Plots.plot!(ratio_norm_plot_5, normalized_ratio, label=latexstring("\$\\beta_{gt}:$(β_gt[1])\\mathrm{\\frac{Pa\\cdot s}{m}}\$"), marker=1, legend=:outerbottom, legend_column=3)
                        end

                        # Plots.plot!(height_error_plt, time, height_error, label=string("Number of elements: ",ne), marker=1, legend=:outerbottom, legend_column=2)
                        # Plots.plot!(rel_height_error_plt, time, rel_height_error, label=string("Number of elements: ",ne), marker=1, legend=:outerbottom, legend_column=2)
                        
                        # Plots.plot!(elem_η_plt, est_η, label=string("Number of elements: ",ne), marker=1, legend=:outerbottom, legend_column=2)
                        # Plots.plot!(elem_β_plt, est_β, label=string("Number of elements: ",ne), marker=1, legend=:outerbottom, legend_column=2)
                        # Plots.plot!(elem_ratio_plt, ratio_est, label=string("Number of elements: ",ne), marker=1, legend=:outerbottom, legend_column=2)
                        
                        # Plots.plot!(elem_rel_height_error_plt, time, rel_height_error, label=string("Number of elements: ",ne), legend=:outerbottom, legend_column=2)
                        # Plots.plot!(elem_height_plt, time, est_h, label=string("Number of elements: ",ne), legend=:outerbottom, legend_column=2)
                        # Plots.plot!(elem_height_plt, time, gt_h, label=L"h_{gt}(t)", legend=:outerbottom, legend_column=2)
                    
                        # Plots.plot!(elem_η_norm_plt, est_η_norm, label=string("Number of elements: ",ne), marker=1, legend=:outerbottom, legend_column=2)
                        # Plots.plot!(elem_β_norm_plt, est_β_norm, label=string("Number of elements: ",ne), marker=1, legend=:outerbottom, legend_column=2)
                        # Plots.plot!(elem_ratio_norm_plt, normalized_ratio, label=string("Number of elements: ",ne), marker=1, legend=:outerbottom, legend_column=2)
                    
                    elseif sim_time == 10.0
                        if ne == 6
                            Plots.plot!(η_norm_plot_10, est_η_norm, label=latexstring("\$\\beta_{gt}:$(β_gt[1])\\mathrm{\\frac{Pa\\cdot s}{m}}\$"), marker=1, color=palette(:Set1)[(parse(Int, dir)-1) % length(palette(:Set1))+1])
                            Plots.plot!(β_norm_plot_10, est_β_norm, label=latexstring("\$\\beta_{gt}:$(β_gt[1])\\mathrm{\\frac{Pa\\cdot s}{m}}\$"), marker=1, color=palette(:Set1)[(parse(Int, dir)-1) % length(palette(:Set1))+1])
                            
                            Plots.plot!(ratio_plot_10, ratio_est, label=latexstring("\$\\beta_{gt}:$(β_gt[1])\\mathrm{\\frac{Pa\\cdot s}{m}}\$"), marker=1, yscale=:log10, color=palette(:Set1)[(parse(Int, dir)-1) % length(palette(:Set1))+1])
                            Plots.hline!(ratio_plot_10, [ratio_gt], linestyle=:dash, label=false, color=:black, yscale=:log10)

                            Plots.plot!(h_plot_10, time, est_h, label=latexstring("\$\\beta_{gt}:$(β_gt[1])\\mathrm{\\frac{Pa\\cdot s}{m}}\$"), color=palette(:Set1)[(parse(Int, dir)-1) % length(palette(:Set1))+1])
                            Plots.plot!(h_plot_10, time, gt_h, label=false, linestyle=:dash, color=palette(:Set1)[(parse(Int, dir)-1) % length(palette(:Set1))+1])
                            Plots.plot!(rel_height_error_glob_plot_10, time, rel_height_error, label=latexstring("\$\\beta_{gt}:$(β_gt[1])\\mathrm{\\frac{Pa\\cdot s}{m}}\$"), color=palette(:Set1)[(parse(Int, dir)-1) % length(palette(:Set1))+1])
                            Plots.plot!(ratio_norm_plot_10, normalized_ratio, label=latexstring("\$\\beta_{gt}:$(β_gt[1])\\mathrm{\\frac{Pa\\cdot s}{m}}\$"), marker=1, color=palette(:Set1)[(parse(Int, dir)-1) % length(palette(:Set1))+1])
                        end
                    end
                end
                plot_path_noise = joinpath(sim_time_folder,"post_analysis_noise","plots")
                set_file(plot_path_noise)
                @info "Saving plots to $plot_path_noise"
                Plots.savefig(height_error_plt, joinpath(plot_path_noise,"height_error.pdf"))
                Plots.savefig(rel_height_error_plt, joinpath(plot_path_noise,"relative_height_error.pdf"))
                @info "Saved plots to $plot_path_noise"
            end
            plot_path_sim_time = joinpath(elem_size_folder,"post_analysis_time","plots")
            set_file(plot_path_sim_time)
            
            Plots.savefig(sim_window_η_plt, joinpath(plot_path_sim_time,"η.pdf"))
            Plots.savefig(sim_window_β_plt, joinpath(plot_path_sim_time,"β.pdf"))
            Plots.savefig(sim_window_ratio_plt, joinpath(plot_path_sim_time,"η_β_ratio.pdf"))
            Plots.savefig(sim_window_rel_height_error_plt, joinpath(plot_path_sim_time,"relative_height_error.pdf"))
            Plots.savefig(sim_window_height_plt, joinpath(plot_path_sim_time,"height.pdf"))
            Plots.savefig(sim_window_η_norm_plt, joinpath(plot_path_sim_time,"η_normalized.pdf"))
            Plots.savefig(sim_window_β_norm_plt, joinpath(plot_path_sim_time,"β_normalized.pdf"))
            Plots.savefig(sim_window_ratio_norm_plt, joinpath(plot_path_sim_time,"η_β_ratio_normalized.pdf"))
            Plots.savefig(plt_slices, joinpath(plot_path_sim_time,"cost_slices_along_directions.pdf"))
            @info "Saved plots to $plot_path_sim_time"
        end
        
        plot_path_elems = joinpath(filepath_dir,"post_analysis","plots")
        set_file(plot_path_elems)

        Plots.savefig(elem_ratio_plt, joinpath(plot_path_elems,"η_β_ratio.pdf"))
        Plots.savefig(elem_η_plt, joinpath(plot_path_elems,"η.pdf"))
        Plots.savefig(elem_β_plt, joinpath(plot_path_elems,"β.pdf"))
        Plots.savefig(elem_ratio_norm_plt, joinpath(plot_path_elems,"η_β_ratio_normalized.pdf"))
        Plots.savefig(elem_η_norm_plt, joinpath(plot_path_elems,"η_normalized.pdf"))
        Plots.savefig(elem_β_norm_plt, joinpath(plot_path_elems,"β_normalized.pdf"))
        Plots.savefig(elem_rel_height_error_plt, joinpath(plot_path_elems,"relative_height_error.pdf"))
        Plots.savefig(elem_height_plt, joinpath(plot_path_elems,"height.pdf"))
        @info "Saved plots to $plot_path_elems"
    end
    plot_path_global = joinpath(filepath,"post_analysis_global","plots")
    set_file(plot_path_global)
    
    Plots.savefig(η_norm_plot_5, joinpath(plot_path_global,"η_normalized_5.pdf"))
    Plots.savefig(β_norm_plot_5, joinpath(plot_path_global,"β_normalized_5.pdf"))
    Plots.savefig(ratio_plot_5, joinpath(plot_path_global,"η_β_ratio_5.pdf"))
    Plots.savefig(rel_height_error_glob_plot_5, joinpath(plot_path_global,"relative_height_error_5.pdf"))
    Plots.savefig(h_plot_5, joinpath(plot_path_global,"height_comparison_5.pdf"))
    Plots.savefig(ratio_norm_plot_5, joinpath(plot_path_global,"η_β_ratio_normalized_5.pdf"))
    Plots.savefig(h_norm_plot_5, joinpath(plot_path_global,"height_normalized_5.pdf"))

    Plots.savefig(η_norm_plot_10, joinpath(plot_path_global,"η_normalized_10.pdf"))
    Plots.savefig(β_norm_plot_10, joinpath(plot_path_global,"β_normalized_10.pdf"))
    Plots.savefig(ratio_plot_10, joinpath(plot_path_global,"η_β_ratio_10.pdf"))
    Plots.savefig(rel_height_error_glob_plot_10, joinpath(plot_path_global,"relative_height_error_10.pdf"))
    Plots.savefig(ratio_norm_plot_10, joinpath(plot_path_global,"η_β_ratio_normalized_10.pdf"))
    Plots.savefig(h_norm_plot_10, joinpath(plot_path_global,"height_normalized_10.pdf"))
    Plots.savefig(h_plot_10, joinpath(plot_path_global,"height_comparison_10.pdf"))
    @info "Saved plots to $plot_path_global"
end

function post_analysis_real(filepath_gt_::String, filepath::String, avoid_list)
    dir_list = readdir(filepath)

    # for 1/2 linewidth
    # fs::Int = 12
    # plot_ht::Int = 350
    # plot_width ::Int = 477
    # plt_lft_marge = 1pt
    # plt_right_marge = 5pt
    # plt_top_margin = 1pt

    # for 1/3 linewidth
    # fs::Int = 12
    # plt_height::Int = 255
    # plt_width::Int = 318
    # plt_lft_margin = -5pt
    # plt_right_margin = -1pt
    # plt_top_margin = -5pt

    # figures to compare convergence with slip levels
    η_plot_5 = set_plot(fs, sz=(plt_width, plt_height), legend_column=3, right_margin=plt_right_margin, left_margin=plt_lft_margin)
    Plots.xlabel!(η_plot_5,L"\mathrm{Time\;(s)}")
    Plots.ylabel!(η_plot_5,latexstring("\$\\eta_{est}(t)\$ (KPa \$\\cdot\$ s)"))

    β_plot_5 = set_plot(fs, sz=(plt_width, plt_height), legend_column=3, right_margin=plt_right_margin, left_margin=plt_lft_margin)
    Plots.xlabel!(β_plot_5, L"\mathrm{Time\;(s)}")
    Plots.ylabel!(β_plot_5, L"\beta_{est}\;(mm^{-1})")

    η_plot_10 = set_plot(fs, sz=(plt_width, plt_height), legend_column=3, right_margin=plt_right_margin, left_margin=plt_lft_margin)
    Plots.xlabel!(η_plot_10,L"\mathrm{Time\;(s)}")
    Plots.ylabel!(η_plot_10,latexstring("\$\\eta_{est}(t)\$ (KPa \$\\cdot\$ s)"))

    β_plot_10 = set_plot(fs, sz=(plt_width, plt_height), legend_column=3, right_margin=plt_right_margin, left_margin=plt_lft_margin)
    Plots.xlabel!(β_plot_10, L"\mathrm{Time\;(s)}")
    Plots.ylabel!(β_plot_10, L"\beta_{est}\;(mm^{-1})")

    # height plots
    gt_h_plot = set_plot(fs, sz=(plt_width, plt_height), legend_column=3, right_margin=plt_right_margin, left_margin=plt_lft_margin)
    Plots.xlabel!(gt_h_plot,L"\mathrm{Time\;(s)}")
    Plots.ylabel!(gt_h_plot,L"h_{m}\;\mathrm{(mm)}")  
    Plots.xlims!(gt_h_plot, 0, end_obs_win)

    h_plot_5 = set_plot(fs, sz=(plt_width, plt_height), legend_column=3, right_margin=plt_right_margin, left_margin=plt_lft_margin)
    Plots.xlabel!(h_plot_5,L"\mathrm{Time\;(s)}")
    Plots.ylabel!(h_plot_5,L"h\;\mathrm{(mm)}")  
    Plots.xlims!(h_plot_5, 0, end_obs_win)

    h_plot_est_5 = set_plot(fs, sz=(plt_width, plt_height), legend_column=3, right_margin=plt_right_margin, left_margin=plt_lft_margin)
    Plots.xlabel!(h_plot_est_5,L"\mathrm{Time\;(s)}")
    Plots.ylabel!(h_plot_est_5,L"h_{est}\;\mathrm{(mm)}")  
    Plots.xlims!(h_plot_est_5, 0, end_obs_win)

    h_plot_10 = set_plot(fs, sz=(plt_width, plt_height), legend_column=3, right_margin=plt_right_margin, left_margin=plt_lft_margin)
    Plots.xlabel!(h_plot_10,L"\mathrm{Time\;(s)}")
    Plots.ylabel!(h_plot_10,L"h_{est}\;\mathrm{(mm)}") 
    Plots.xlims!(h_plot_10, 0, end_obs_win) 

    # normalised height plots
    h_norm_plot_5 = set_plot(fs, sz=(plt_width, plt_height), legend_column=3, right_margin=plt_right_margin, left_margin=plt_lft_margin)
    Plots.hline!(h_norm_plot_5, [1.0],  linestyle=:dash, label=false, color=:black)
    Plots.xlabel!(h_norm_plot_5,L"\mathrm{Time\;(s)}")
    Plots.ylabel!(h_norm_plot_5,L"h_{est}/h_{m}")  
    Plots.xlims!(h_norm_plot_5, 0, end_obs_win)
    Plots.ylims!(h_norm_plot_5, y_lims_h_norm)

    h_norm_plot_10 = set_plot(fs, sz=(plt_width, plt_height), legend_column=3, right_margin=plt_right_margin, left_margin=plt_lft_margin)
    Plots.hline!(h_norm_plot_10, [1.0],  linestyle=:dash, label=false, color=:black)
    Plots.xlabel!(h_norm_plot_10,L"\mathrm{Time\;(s)}")
    Plots.ylabel!(h_norm_plot_10,L"h_{est}/h_{m}") 
    Plots.xlims!(h_norm_plot_10, 0, end_obs_win)
    Plots.ylims!(h_norm_plot_10, y_lims_h_norm)

    # η\β ratio plot
    ratio_plot_5 = set_plot(fs, sz=(plt_width, plt_height), legend_column=3, right_margin=plt_right_margin, left_margin=plt_lft_margin)
    Plots.xlabel!(ratio_plot_5,L"\mathrm{Iterations}")
    Plots.ylabel!(ratio_plot_5,L"\eta_{est}/\beta_{est}")

    ratio_plot_10 = set_plot(fs, sz=(plt_width, plt_height), legend_column=3, right_margin=plt_right_margin, left_margin=plt_lft_margin)
    Plots.xlabel!(ratio_plot_10,L"\mathrm{Iterations}")
    Plots.ylabel!(ratio_plot_10,L"\eta_{est}/\beta_{est}")

    rel_height_error_glob_plot_5 = set_plot(fs, sz=(plt_width, plt_height), legend_column=3, right_margin=plt_right_margin, left_margin=plt_lft_margin)
    Plots.xlabel!(rel_height_error_glob_plot_5, L"\mathrm{Time\;(s)}")
    Plots.ylabel!(rel_height_error_glob_plot_5, latexstring("Relative Height Error (\$\\%\$)"))
    Plots.ylims!(rel_height_error_glob_plot_5, y_lims_rel_error)

    rel_height_error_glob_plot_10 = set_plot(fs, sz=(plt_width, plt_height), legend_column=3, right_margin=plt_right_margin, left_margin=plt_lft_margin)
    Plots.xlabel!(rel_height_error_glob_plot_10, L"\mathrm{Time\;(s)}")
    Plots.ylabel!(rel_height_error_glob_plot_10, latexstring("Relative Height Error (\$\\%\$)"))
    Plots.ylims!(rel_height_error_glob_plot_10, y_lims_rel_error)
    
    for dir in dir_list
        if dir in avoid_list || dir == "post_analysis_global"
            println("Skipping directory: ", dir)
            continue
        end
        filepath_dir = joinpath(filepath, dir)
        filepath_gt = joinpath(filepath_gt_, dir)

        printstyled("Processing directory: $(filepath_dir)\n", color=:green)
        sim_params = read_json(joinpath(filepath_gt,"data","sim_params.json"))
        sim_time = sim_params["simulation_time"]    
        t_steps = sim_params["time_steps"]
        gt_h_ = readdlm(joinpath(filepath_gt,"data","h.csv"), ',', Float64)
        time = collect(Float64, range(start=0, stop=sim_time, step=t_steps))

        elem_size_folders = readdir(filepath_dir)
        
        for elem_size_folder_ in elem_size_folders
            if elem_size_folder_ == "post_analysis" || elem_size_folder_ == "Q2_2" || elem_size_folder_ == "Q2_4" || elem_size_folder_ == "Q2_8"
                continue
            end
            
            elem_size_folder = joinpath(filepath_dir, elem_size_folder_)
            printstyled("Processing element size folder: $(elem_size_folder)\n", color=:blue)  
            sim_time_folders = readdir(elem_size_folder)
            
            # figures for Simulation window vise camparison
            sim_window_η_plt = set_plot(fs, sz=(plt_width, plt_height))
            Plots.xlabel!(sim_window_η_plt,L"\mathrm{Iterations}")
            Plots.ylabel!(sim_window_η_plt,L"\eta\;\mathrm{(KPa\cdot s)}")
            
            sim_window_β_plt = set_plot(fs, sz=(plt_width, plt_height))
            Plots.xlabel!(sim_window_β_plt, L"\mathrm{Iterations}")
            Plots.ylabel!(sim_window_β_plt, L"\beta\;\mathrm{(mm^{-1})}")

            sim_window_ratio_plt = set_plot(fs, sz=(plt_width, plt_height))
            Plots.xlabel!(sim_window_ratio_plt, L"\mathrm{Iterations}")   
            Plots.ylabel!(sim_window_ratio_plt, L"\eta/\beta\;(KPa\cdot s \cdot mm^{-1})")

            # plots for height error
            sim_window_rel_height_error_plt = set_plot(fs, sz=(plt_width, plt_height))
            Plots.plot!(sim_window_rel_height_error_plt, [], left_margin=plt_lft_margin, label=false)
            Plots.xlabel!(sim_window_rel_height_error_plt, L"\mathrm{Time\;(s)}")
            Plots.ylabel!(sim_window_rel_height_error_plt, L"\mathrm{Relative\;Height\;Error}")
            Plots.xlims!(sim_window_rel_height_error_plt, 0, end_obs_win)

            sim_window_height_plt = set_plot(fs, sz=(plt_width, plt_height))
            Plots.plot!(sim_window_height_plt, [], left_margin=plt_lft_margin, label=false)
            Plots.xlabel!(sim_window_height_plt, L"\mathrm{Time\;(s)}")
            Plots.ylabel!(sim_window_height_plt, L"\mathrm{Height\;(mm)}")
            Plots.xlims!(sim_window_height_plt, 0, end_obs_win)

            for sim_time_folder_ in sim_time_folders

                if sim_time_folder_ == "post_analysis_time" || sim_time_folder_ == "Results" || sim_time_folder_ == "simtime_2.0" #|| sim_time_folder_ == "simtime_10.0"
                    continue
                end
                
                height_error_plt = set_plot(fs, sz=(plt_width, plt_height))
                Plots.plot!(height_error_plt, [], label=false)
                Plots.xlabel!(height_error_plt, L"\mathrm{Time\;(s)}")
                Plots.ylabel!(height_error_plt, L"\mathrm{Height\;Error\;(mm)}")
                Plots.xlims!(height_error_plt, 0, end_obs_win)
                
                rel_height_error_plt = set_plot(fs, sz=(plt_width, plt_height))
                Plots.plot!(rel_height_error_plt, [], label=false)
                Plots.xlabel!(rel_height_error_plt, L"\mathrm{Time\;(s)}")
                Plots.ylabel!(rel_height_error_plt, L"\mathrm{Height\;Error\;(mm)}")
                Plots.xlims!(rel_height_error_plt, 0, end_obs_win)

                sim_time_folder = joinpath(elem_size_folder, sim_time_folder_, "noise_0.0")

                printstyled("Processing simulation time folder: $(sim_time_folder)\n", color=:cyan)

                window_dirs = readdir(sim_time_folder)
                for window_dir in window_dirs
                    if window_dir == "Results" || window_dir == "post_analysis_window" || window_dir == "single_window" || window_dir == "post_analysis_noise"
                        println("Skipping directory: $window_dir")
                        continue
                    end
                    win_exp_path = joinpath(filepath, elem_size_folder, sim_time_folder, window_dir)
                    
                    println("Processing window: $win_exp_path")
                    exp_params = read_json(joinpath(win_exp_path ,"Results","data","experiment_parameters.json"))
                    sim_time_exp = exp_params["sim_time_exp"]
                    data_type = exp_params["data_type"]
                    ne = exp_params["ne_exp"]

                    ObsDataList, splinexObs, splineyObs = read_csv(joinpath(filepath_gt,"data","img_data","contour_data"))  
                    obsBorderPts, nSplinex, nSpliney, pd = add_noise(ObsDataList, nFactor=0.0)

                    est_ηpList = readdlm(joinpath(win_exp_path,"Results","data","est_η.csv"), ',', Float64)
                    est_βpList = readdlm(joinpath(win_exp_path,"Results","data","est_β.csv"), ',', Float64)
                    est_h_list = readdlm(joinpath(win_exp_path,"Results","data","est_h.csv"), ',', Float64)
                    
                    data_ranges_ = get_time_windows(joinpath(win_exp_path,"Results","data","window_data","data_ranges.csv"))
                    t_windows = readdlm(joinpath(win_exp_path,"Results","data","window_data","t_windows.csv"),',',Float64)
                    time_windows = readdlm(joinpath(win_exp_path,"Results","data","window_data","time_windows.csv"),',',Float64)

                    symBorderPts, splinex, spliney = read_csv(joinpath(win_exp_path,"Results","data","sim_data","2D_border_points"))
                    
                    println("Time windows: $(time_windows)")
                    obs_time = sum(time_windows)
                    
                    if obs_time != sim_time
                        @warn "Observation time frame $obs_time is less than preset ground truth time frame $sim_time, updating time frame"
                        min_time = min(obs_time, sim_time)
                        sim_time = min_time
                        obs_time = min_time
                    end
                    
                    if obs_time < sim_time_exp
                        @warn "Observation time frame $obs_time is less than experimental simulation time frame $sim_time_exp, switching to observation time frame"
                        sim_time_exp = obs_time
                    end

                    data_point_len = round(Int, obs_time/t_steps)
                    gt_h = gt_h_[1:(data_point_len+1)]
                    symBorderPts = symBorderPts[1:data_point_len+1, :]
                    obsBorderPts = obsBorderPts[1:data_point_len+1, :]
                    est_h = est_h_list[1:data_point_len+1, :]
                    time = collect(range(start=0, stop=sim_time, step=t_steps))

                    height_error = abs.(est_h - gt_h)
                    rel_height_error = height_error ./ gt_h * 100.0

                    ratio_est = est_ηpList ./ est_βpList
            
                    Plots.plot!(sim_window_η_plt, est_ηpList, label=string("Window - $(sim_time)s"), marker=1, legend=:outerbottom, legend_column=2)
                    Plots.plot!(sim_window_β_plt, est_βpList, label=string("Window - $(sim_time)s"), marker=1, legend=:outerbottom, legend_column=2)
                    Plots.plot!(sim_window_ratio_plt, ratio_est, label=string("Window - $(sim_time)s"), marker=1, legend=:outerbottom, legend_column=2)
                    
                    Plots.plot!(sim_window_rel_height_error_plt, time, rel_height_error, label=string("Window - $(sim_time)s"," - ne: ",ne), legend=:outerbottom, legend_column=2)
                    Plots.plot!(sim_window_height_plt, time[1:length(est_h)], est_h, label=string("Window - $(sim_time)s"), legend=:outerbottom, legend_column=2)
                    Plots.plot!(sim_window_height_plt, time, gt_h, label=L"h_{gt}(t)", legend=:outerbottom, legend_column=2)
                    
                    if sim_time_exp == 0.5
                        if ne == 6
                            t_prev = 0.1
                            for ti::Int in 1:(size(data_ranges_, 1)-1)
                                t = t_windows[ti]
                                data_range_ = data_ranges_[ti]
                                t_win = collect(range(start=t_prev, stop=t, step=t_steps))

                                Plots.vline!(η_plot_5, [t], color=:gray, linestyle=:dash, label=false)
                                Plots.vline!(β_plot_5, [t], color=:gray, linestyle=:dash, label=false)
                                Plots.vline!(ratio_plot_5, [t], color=:gray, linestyle=:dash, label=false)

                                if ti == 1
                                    Plots.plot!(η_plot_5, t_win, est_ηpList[data_range_], label=string(L"\mathrm{Exp}:\;",dir), color=palette(:Set1)[(parse(Int, dir)-1) % length(palette(:Set1))+1])
                                    Plots.plot!(β_plot_5, t_win, est_βpList[data_range_], label=string(L"\mathrm{Exp}:\;",dir), color=palette(:Set1)[(parse(Int, dir)-1) % length(palette(:Set1))+1])
                                    Plots.plot!(ratio_plot_5, t_win, ratio_est[data_range_], label=string(L"\mathrm{Exp}:\;",dir), color=palette(:Set1)[(parse(Int, dir)-1) % length(palette(:Set1))+1])
                                else
                                    Plots.plot!(η_plot_5, t_win, est_ηpList[data_range_], label=false, color=palette(:Set1)[(parse(Int, dir)-1) % length(palette(:Set1))+1])
                                    Plots.plot!(β_plot_5, t_win, est_βpList[data_range_], label=false, color=palette(:Set1)[(parse(Int, dir)-1) % length(palette(:Set1))+1])
                                    Plots.plot!(ratio_plot_5, t_win, ratio_est[data_range_], label=false, color=palette(:Set1)[(parse(Int, dir)-1) % length(palette(:Set1))+1])
                                end
                                t_prev = t+t_steps
                            end
                            
                            Plots.plot!(rel_height_error_glob_plot_5, time, rel_height_error, label=string(L"\mathrm{Exp}:\;",dir), color=palette(:Set1)[(parse(Int, dir)-1) % length(palette(:Set1))+1])
                            Plots.plot!(h_norm_plot_5, time, est_h./gt_h, label=string(L"\mathrm{Exp}:\;",dir), color=palette(:Set1)[(parse(Int, dir)-1) % length(palette(:Set1))+1])
                            Plots.plot!(h_plot_est_5, time, est_h, label=string(L"\mathrm{Exp}:\;",dir), color=palette(:Set1)[(parse(Int, dir)-1) % length(palette(:Set1))+1])
                            Plots.plot!(gt_h_plot, time, gt_h, label=string(L"\mathrm{Exp}:\;",dir), color=palette(:Set1)[(parse(Int, dir)-1) % length(palette(:Set1))+1])
                            Plots.plot!(h_plot_5, time, gt_h, label=string(L"\mathrm{Exp}:\;",dir), color=palette(:Set1)[(parse(Int, dir)-1) % length(palette(:Set1))+1], linestyle=:dash)
                            Plots.plot!(h_plot_5, time, est_h, label=false, color=palette(:Set1)[(parse(Int, dir)-1) % length(palette(:Set1))+1])
                            if dir == "1"
                                for ti::Int in 1:(size(data_ranges_, 1))
                                    t_win = t_windows[ti]
                                    data_range_ = data_ranges_[ti]
                                    Plots.vline!(h_norm_plot_5, [t_win], color=:gray, linestyle=:dash, label=false)
                                    Plots.vline!(h_plot_est_5, [t_win], color=:gray, linestyle=:dash, label=false)
                                    Plots.vline!(h_plot_5, [t_win], color=:gray, linestyle=:dash, label=false)
                                    # Plots.vline!(gt_h_plot, [t_win], color=:gray, linestyle=:dash, label=false)
                                end
                            end
                        
                            # Plots.plot!(ratio_plot_5, ratio_est, label=latexstring("\$\\beta_{gt}:$(β_gt[1])\\mathrm{\\frac{Pa\\cdot s}{m}}\$"), legend=:outerbottom, legend_column=3, yscale=:log10)
                            # Plots.hline!(ratio_plot_5, [ratio_gt], linestyle=:dash, label=false, color=:black, yscale=:log10)
                            # Plots.plot!(ratio_norm_plot_5, normalized_ratio, label=latexstring("\$\\beta_{gt}:$(β_gt[1])\\mathrm{\\frac{Pa\\cdot s}{m}}\$"), marker=1, legend=:outerbottom, legend_column=3)
                        end

                        Plots.plot!(height_error_plt, time, height_error, label=string("Number of elements: ",ne), marker=1, legend=:outerbottom, legend_column=2)
                        Plots.plot!(rel_height_error_plt, time, rel_height_error, label=string("Number of elements: ",ne), marker=1, legend=:outerbottom, legend_column=2)
                        
                    elseif sim_time == 10.0
                        if ne == 6
                            Plots.plot!(η_norm_plot_10, est_η_norm, label=latexstring("\$\\beta_{gt}:$(β_gt[1])\\mathrm{\\frac{Pa\\cdot s}{m}}\$"), marker=1, legend=:outerbottom, legend_column=3)
                            Plots.plot!(β_norm_plot_10, est_β_norm, label=latexstring("\$\\beta_{gt}:$(β_gt[1])\\mathrm{\\frac{Pa\\cdot s}{m}}\$"), marker=1, legend=:outerbottom, legend_column=3)
                            
                            Plots.plot!(ratio_plot_10, ratio_est, label=latexstring("\$\\beta_{gt}:$(β_gt[1])\\mathrm{\\frac{Pa\\cdot s}{m}}\$"), marker=1, legend=:outerbottom, legend_column=3, yscale=:log10)
                            Plots.hline!(ratio_plot_10, [ratio_gt], linestyle=:dash, label=false, color=:black, yscale=:log10)

                            Plots.plot!(rel_height_error_glob_plot_10, time, rel_height_error, label=string("Window - $(sim_time)s"," - ne: ",ne), legend=:outerbottom, legend_column=2)
                            Plots.plot!(ratio_norm_plot_10, normalized_ratio, label=latexstring("\$\\beta_{gt}:$(β_gt[1])\\mathrm{\\frac{Pa\\cdot s}{m}}\$"), marker=1, legend=:outerbottom, legend_column=3)
                        end
                    end
                end
                # plot_path_noise = joinpath(sim_time_folder,"post_analysis_noise","plots")
                # set_file(plot_path_noise)
                # @info "Saving plots to $plot_path_noise"
                # Plots.savefig(height_error_plt, joinpath(plot_path_noise,"height_error.pdf"))
                # Plots.savefig(rel_height_error_plt, joinpath(plot_path_noise,"relative_height_error.pdf"))
                # @info "Saved plots to $plot_path_noise"
            end
            plot_path_sim_time = joinpath(elem_size_folder,"post_analysis_time","plots")
            set_file(plot_path_sim_time)
            
            Plots.savefig(sim_window_η_plt, joinpath(plot_path_sim_time,"η.pdf"))
            Plots.savefig(sim_window_β_plt, joinpath(plot_path_sim_time,"β.pdf"))
            Plots.savefig(sim_window_ratio_plt, joinpath(plot_path_sim_time,"η_β_ratio.pdf"))
            # Plots.savefig(sim_window_rel_height_error_plt, joinpath(plot_path_sim_time,"relative_height_error.pdf"))
            Plots.savefig(sim_window_height_plt, joinpath(plot_path_sim_time,"height.pdf"))
            @info "Saved plots to $plot_path_sim_time"
        end
    end
    plot_path_global = joinpath(filepath,"post_analysis_global","plots")
    set_file(plot_path_global)
    
    Plots.savefig(η_plot_5, joinpath(plot_path_global,"η_plot_5.pdf"))
    Plots.savefig(β_plot_5, joinpath(plot_path_global,"β_plot_5.pdf"))
    Plots.savefig(ratio_plot_5, joinpath(plot_path_global,"η_β_ratio_5.pdf"))
    Plots.savefig(h_plot_est_5, joinpath(plot_path_global,"height_5.pdf"))
    Plots.savefig(h_plot_5, joinpath(plot_path_global,"height_comparison_5.pdf"))
    Plots.savefig(rel_height_error_glob_plot_5, joinpath(plot_path_global,"relative_height_error_5.pdf"))
    Plots.savefig(h_norm_plot_5, joinpath(plot_path_global,"height_normalized_5.pdf"))
    Plots.savefig(gt_h_plot, joinpath(plot_path_global,"gt_height.pdf"))

    # Plots.savefig(η_norm_plot_10, joinpath(plot_path_global,"η_normalized_10.pdf"))
    # Plots.savefig(β_norm_plot_10, joinpath(plot_path_global,"β_normalized_10.pdf"))
    # Plots.savefig(ratio_plot_10, joinpath(plot_path_global,"η_β_ratio_10.pdf"))
    # Plots.savefig(rel_height_error_glob_plot_10, joinpath(plot_path_global,"relative_height_error_10.pdf"))
    # Plots.savefig(ratio_norm_plot_10, joinpath(plot_path_global,"η_β_ratio_normalized_10.pdf"))
    # Plots.savefig(h_norm_plot_10, joinpath(plot_path_global,"height_normalized_10.pdf"))
    @info "Saved plots to $plot_path_global"
end

function _get_borders(data_type::String, filepath_gt::String, exp_path::String)

    if data_type  == "synthetic"
        ObsDataList, splinexObs, splineyObs = read_csv(joinpath(filepath_gt,"data","img_data","contour_data"))  
        @info "Data type $data_type Reading synthetic contour data of $(length(ObsDataList)) time steps"
    elseif data_type == "simulated"
        ObsDataList, splinexObs, splineyObs = read_csv(joinpath(filepath_gt,"data","sim_data","contour_data"))  
        @info "Data type $data_type Reading simulated contour data from $(joinpath(filepath_gt,"data","sim_data","contour_data")) of $(length(ObsDataList)) time steps"
    else
        error("Unknown data type: $data_type")
    end

    obsBorderPts, nSplinex, nSpliney, pd = add_noise(ObsDataList, nFactor=0.0)
    symBorderPts, splinex, spliney = read_csv(joinpath(exp_path,"Results","data","sim_data","2D_border_points"))

    if length(obsBorderPts) != length(symBorderPts)
        n_time = min(length(obsBorderPts), length(symBorderPts))
        @warn "Mismatched number of time steps between observed border points ($(length(obsBorderPts))) and simulated border points ($(length(symBorderPts))). Truncating to $n_time time steps."
        obsBorderPts = obsBorderPts[1:n_time, :]
        symBorderPts = symBorderPts[1:n_time, :]

        nSplinex = nSplinex[1:n_time, :]
        nSpliney = nSpliney[1:n_time, :]

        splinex = splinex[1:n_time, :]
        spliney = spliney[1:n_time, :]
    end

    return obsBorderPts, symBorderPts, nSplinex, nSpliney, splinex, spliney
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
                                  gt_color::AbstractString = "heat",
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

function set_time_window(time_step_len::Float64, data::AbstractArray; method::String="linear", window_size::Float64=10.0, global_window_sz::Float64=30.0)
    windows::Vector{AbstractArray} = Vector{AbstractArray}()
    time_windows::Vector{Float64} = Vector{Float64}()
    data_ranges::Vector{AbstractArray} = Vector{AbstractArray}()
    t_windows::Vector{Float64} = Vector{Float64}()

    function get_t_window(window_size::Float64, step_len::Float64, iter::Int, method::String)::Float64
        t_window = 0.0
        if method == "linear"
            t_window = round(window_size*iter, digits=1)
        elseif method == "quadratic"
            t_window = round(window_size*iter^2, digits=1)
        elseif method == "exponential"
            if iter == 1
                t_window = round(window_size*exp(0.5*(iter-1)), digits=1)
            else 
                δt = window_size*round(exp(3*(iter-1)), digits=1) #*exp(0.5*(iter-1)) # round(window_size*exp(3*(iter-1)), digits=1)
                if δt < window_size
                    @warn "Computed time window increment $δt is less than minimum Window $window_size; using minimum."
                    δt = window_size
                end
                t_window = round(window_size*exp(3*(iter-1)), digits=1) #round(window_size*(exp(0.5*(iter-2)) + δt), digits=1)
            end
        end
        return t_window
    end

    iter::Int = 1
    start_point::Int = 1
    t_window_prev::Float64 = 0.0
    t_window_end::Float64 = get_t_window(window_size, time_step_len, iter, method)
    end_point::Int = round(Int,t_window_end*time_step_len)+1

    END_FLAG::Bool = false
    while true
            # If the computed end_point reaches or exceeds the data length, adjust.
            # Keep the original requested value for clearer messages.
            if end_point >= size(data, 1) || t_window_end >= global_window_sz
                if end_point == size(data, 1)
                    @info "Reached end of data at point $end_point."
                elseif end_point > size(data, 1)
                    requested_end = end_point
                    # adjust the time window end to the last available sample
                    t_window_end = round((size(data, 1)-1)/time_step_len, digits=1)
                    end_point = size(data, 1)
                    @warn "Requested end point $requested_end exceeds data size $(size(data,1)); adjusting to end of data (end_point=$end_point). Adjusted end time to $t_window_end seconds."
                elseif t_window_end >= global_window_sz
                    requested_end = t_window_end
                    # adjust the end point to match the global window size
                    t_window_end = global_window_sz
                    end_point = round(Int,t_window_end*time_step_len)+1
                    @warn "Requested time window end $requested_end seconds exceeds global window size $global_window_sz seconds; adjusting to global window size (end_point=$end_point)."
                end
                END_FLAG = true
            end

        data_range = start_point:end_point
        data_range_ = start_point:(end_point-1)
        
        println("Data frame : $data_range")
        println("time windows from : $t_window_prev to $t_window_end")
        println("time Window : $(t_window_end - t_window_prev) seconds")
        println("data length : $(size(data[data_range], 1))")
        println("data range size : $(length(data_range_))")
        println("----------")
        t_window_size = round(t_window_end - t_window_prev, digits=1)
        push!(time_windows, t_window_size)
        push!(windows, data[data_range])
        push!(data_ranges, data_range_)
        push!(t_windows, t_window_end)

        if END_FLAG == true
            break
        end
        iter = iter + 1
        start_point = end_point
        t_window_prev = t_window_end
        t_window_end = get_t_window(window_size, time_step_len, iter, method)
        end_point = round(Int,t_window_end*time_step_len)+1
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
                        LinearAlgebra.BLAS.set_num_threads(max_workers)
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
                                dest_dir = joinpath(params["filepath_res"], "Results", "logs")
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
    refine_list = [6] # [2, 3, 4, 5] # refinement levels, ne = ne_exp^refine
    control = "force" # "force" or "velocity"
    viscosity_type_list = ["bulk_viscosity","constant"]
    window = "multi_window"
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
            @info "Processing ground truth directory: $dir for $viscosity_type viscosity ..."
            filepath_gt = string(_filepath_gt,"/",dir)
            for ne in refine_list
                if ne == 6 && viscosity_type == "constant"
                    noise_level_list = [0.0]
                else
                    noise_level_list = [0.0]
                end
                for noise_level in noise_level_list 
                    if noise_level == 0.0 && viscosity_type == "constant" && ne != 6
                        sim_time_exp_list = [5.0, 2.0, 10.0] # simulation time in seconds
                    elseif viscosity_type == "constant" && ne == 6
                        # sim_time_exp_list = [5.0, 2.0, 10.0] # simulation time in seconds
                        sim_time_exp_list = [20.0, 30.0] # simulation time in seconds
                    elseif noise_level == 0.0 && viscosity_type == "bulk_viscosity"
                        sim_time_exp_list = [2.0, 5.0, 10.0] # simulation time in seconds
                    else
                        sim_time_exp_list = [2.0, 5.0, 10.0]  # simulation time in seconds
                    end
                    println("Simulation time experiments to run: $sim_time_exp_list")
                    for sim_time_exp::Float16 in sim_time_exp_list
                        # ne = ne_exp^ref
                        if viscosity_type == "constant"
                            window = ""
                        end
                        @info "Running optimization with ne = $ne and simulation time = $sim_time_exp with noise level = $noise_level"
                        for FunctionClass_x in FunctionClass_x_List
                            filepath_res = string("/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/experiments/sim_data/optimization/Stokes/$control/$viscosity_type/Q2_16/$dir/$(FunctionClass_x)_$(ne)/simtime_$(sim_time_exp)/noise_$(noise_level)/$window")
                            @info "Running optimization with FunctionClass_x = $FunctionClass_x with $ne elements"

                            exp_params = Dict("FunctionClass_x" => FunctionClass_x, "FunctionClass_u" => "Q2", "FunctionClass_p" => "Q1", "ne_exp" => ne, "sim_time_exp" => sim_time_exp, 
                            "filepath_res" => filepath_res, "filepath_gt"=>filepath_gt, "control" => control, "data_type"=>"simulated", "camera_matrix" => camera_matrix, "WRITE_GT"=> false,
                            "noise_level"=>noise_level, "mode"=>"multi_window")

                            push!(param_list, exp_params)
                        end
                    end
                end
            end
        end
        run_param_list(param_list; max_workers=15)
    end
end

function optimize_syn()

    FunctionClass_x_List = ["Q2"]
    # refine_list = [1, 2, 3] # refinement levels, ne = ne_exp^refine
    refine_list = [6] # [2, 3, 4, 5] # refinement levels, ne = ne_exp^refine
    control = "force" # "force" or "velocity"
    viscosity_type_list = ["bulk_viscosity"]
    window = "multi_window"
    camera_matrix::AbstractArray = [[2.39642674e+03, 0.0, 1.00429248e+03] [0.0, 2.40565353e+03, 7.57028161e+02] [0.0, 0.0, 1.0]]'
    filepath_res::String = ""
    param_list = Vector{Dict}(undef, 0)
    avoid_dirs = ["3_less_noise", "5", "s"]
    for viscosity_type in viscosity_type_list
        _filepath_gt = string("/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/ground_truth/sim_data/Stokes/$control/$viscosity_type/Q2_16")
        dir_list = readdir(_filepath_gt)
        for dir in dir_list
            if dir in avoid_dirs
                continue
                println("Skipping dir $dir")
            end
            @info "Processing ground truth directory: $dir for $viscosity_type viscosity with model type $model_type ..."
            filepath_gt = string(_filepath_gt,"/",dir)
            for ne in refine_list
                if ne == 6 && viscosity_type == "constant"
                    noise_level_list = [0.0]
                else
                    noise_level_list = [0.0]
                end
                for noise_level in noise_level_list 
                    if noise_level == 0.0 && viscosity_type == "constant" && ne != 6
                        sim_time_exp_list = [5.0, 2.0, 10.0] # simulation time in seconds
                    elseif viscosity_type == "constant" && ne == 6
                        # sim_time_exp_list = [5.0, 2.0, 10.0] # simulation time in seconds
                        sim_time_exp_list = [5.0] # simulation time in seconds
                    elseif noise_level == 0.0 && viscosity_type == "bulk_viscosity"
                        sim_time_exp_list = [5.0] # simulation time in seconds
                    else
                        sim_time_exp_list = [2.0, 5.0, 10.0]  # simulation time in seconds
                    end
                    println("Simulation time experiments to run: $sim_time_exp_list")
                    for sim_time_exp::Float16 in sim_time_exp_list
                        # ne = ne_exp^ref
                        if viscosity_type == "constant"
                            window = ""
                        end
                        @info "Running optimization with ne = $ne and simulation time = $sim_time_exp with noise level = $noise_level"
                        for FunctionClass_x in FunctionClass_x_List
                            filepath_res = string("/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/experiments/syn_data/optimization/Stokes/$control/$viscosity_type/Q2_16/$dir/$(FunctionClass_x)_$(ne)/simtime_$(sim_time_exp)/noise_$(noise_level)/$window")
                            @info "Running optimization with FunctionClass_x = $FunctionClass_x with $ne elements"

                            exp_params = Dict("FunctionClass_x" => FunctionClass_x, "FunctionClass_u" => "Q2", "FunctionClass_p" => "Q1", "ne_exp" => ne, "sim_time_exp" => sim_time_exp, 
                            "filepath_res" => filepath_res, "filepath_gt"=>filepath_gt, "control" => control, "data_type"=>"synthetic", "camera_matrix" => camera_matrix, "WRITE_GT"=> false,
                            "noise_level"=>noise_level, "mode"=>"multi_window")

                            push!(param_list, exp_params)
                        end
                    end
                end
            end
        end
        run_param_list(param_list; max_workers=15)
    end
end

function optimize_real()

    FunctionClass_x_List = ["Q2"]
    # refine_list = [1, 2, 3] # refinement levels, ne = ne_exp^refine
    refine_list = [6] # refinement levels, ne = ne_exp^refine
    control = "force" # "force" or "velocity"
    η_start = 1.0
    β_start = 1.0
    viscosity_type = "bulk_viscosity"
    
    r::Float64 = 25.0  # radius of the cylinder in mm
    h::Float64 = 38.5  # height of the cylinder in mm
    camera_matrix::AbstractArray = [[2.39642674e+03, 0.0, 1.00429248e+03] [0.0, 2.40565353e+03, 7.57028161e+02] [0.0, 0.0, 1.0]]'
    # sim_time_exp::Float64 = 5.0 # simulation time in seconds
    F_ext::Float64 = 9.812*1e3 # force applied to the cylinder in N
    model_type = "Stokes" # "carreau" or "Stokes"
    window = "multi_window" # "multi_window" or "single_window"
    filepath_res::String = ""
    param_list = Vector{Dict}(undef, 0)
    sim_time_exp_list = [0.5] # simulation time in seconds
    avoid_dirs = ["3_less_noise","s"]

    if model_type == "carreau" && viscosity_type == "bulk_viscosity"
        _filepath_gt = string("/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/ground_truth/sim_data/Carreau")
    else
        _filepath_gt = string("/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/ground_truth/physical_data")
    end

    dir_list = readdir(_filepath_gt)
    for dir in dir_list
        if dir in avoid_dirs
            continue
            println("Skipping dir $dir")
        end
        @info "Processing ground truth directory: $dir for $viscosity_type viscosity ..."
        filepath_gt = string(_filepath_gt,"/",dir)
        for ne in refine_list
            println("Simulation time experiments to run: $sim_time_exp_list")
            for sim_time_exp in sim_time_exp_list
                println("Simulation time: $sim_time_exp seconds")
                for ne in refine_list
                    for FunctionClass_x in FunctionClass_x_List
                        if model_type == "carreau" && viscosity_type == "bulk_viscosity"
                            filepath_res = string("/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/experiments/syn_data/optimization/Carreau/$dir/$(FunctionClass_x)_$(ne)/simtime_$(sim_time_exp)/noise_0.0/$window")
                            
                            exp_params = Dict("FunctionClass_x" => FunctionClass_x, "FunctionClass_u" => "Q2", "FunctionClass_p" => "Q1", "ne_exp" => ne, "sim_time_exp" => sim_time_exp, 
                            "η_start" => η_start, "β_start" => β_start, "filepath_res" => filepath_res, "filepath_gt"=>filepath_gt, "control" => control, "data_type"=>"synthetic", 
                            "camera_matrix" => camera_matrix, "WRITE_GT"=> false, "noise_level"=>0.0, "mode"=>window)
                         else
                            filepath_res = string("/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/experiments/physical_data/optimization/$dir/$(FunctionClass_x)_$(ne)/simtime_$(sim_time_exp)/noise_0.0/$window")
                            
                            exp_params = Dict("FunctionClass_x" => FunctionClass_x, "FunctionClass_u" => "Q2", "FunctionClass_p" => "Q1", "ne_exp" => ne, "sim_time_exp" => sim_time_exp, 
                            "η_start" => η_start, "β_start" => β_start, "filepath_res" => filepath_res, "filepath_gt"=>filepath_gt, "control" => control, "viscosity_type"=>viscosity_type, 
                            "data_type"=>"physical", "r" => r, "h" => h, "camera_matrix" => camera_matrix, "F_ext" => F_ext, "mode"=>window)
                        end
                        @info "Running optimization with FunctionClass_x = $FunctionClass_x with $ne elements"
                        
                        push!(param_list, exp_params)
                    end
                end
            end
        end
    end
    run_param_list(param_list; max_workers=5)
end

function plot_()

    control::String = "force" # "force" or "velocity"
    viscosity_type_list = ["bulk_viscosity"] #,"constant"]
    model_type::String = "Stokes" # "carreau" or "Stokes"
    avoid_dirs = ["3_less_noise", "s"]
    data_type_list = ["synthetic"] #,"simulated", "physical"]

    for data_type in data_type_list
        if data_type == "synthetic"
            if model_type == "carreau"
                base_path = string("/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/experiments/syn_data/optimization/Carreau")
            else
                base_path = "/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/experiments/syn_data/optimization/Stokes"
            end
        elseif data_type == "simulated"
            base_path = "/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/experiments/sim_data/optimization/Stokes"
        elseif data_type == "physical"
            base_path = "/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/experiments/physical_data/optimization"
        end
        for viscosity_type::String in viscosity_type_list
            if data_type == "physical"
                filepath_gt = string("/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/ground_truth/physical_data")
                filepath_res = base_path
            else
                filepath_gt = string("/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/ground_truth/sim_data/Stokes/$control/$viscosity_type/Q2_16")
                filepath_res = string("$base_path/$control/$viscosity_type/Q2_16")
            end
            if model_type == "carreau" && viscosity_type == "bulk_viscosity"
                filepath_gt = string("/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/ground_truth/sim_data/Carreau")
                filepath_res = base_path
            end
            
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
                post_analysis_const(filepath_gt, filepath_res, avoid_dirs)
            elseif viscosity_type == "bulk_viscosity" && model_type != "carreau" && data_type != "physical"
                post_analysis_bulk(filepath_gt, filepath_res, avoid_dirs)
            elseif data_type == "physical"
                post_analysis_real(filepath_gt, filepath_res, avoid_dirs)
            end
        end
    end
end

# main()
plot_()
# optimize_sim()
# optimize_syn()
# optimize_real()
