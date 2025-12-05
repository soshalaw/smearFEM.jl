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

function simulate_data(exp_params::Dict)
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
        
        gt_viscosity_type = params["viscosity_type"]

        η_gt = params["η"][1]
        β_gt = params["β"][1]

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
        
    else
        error("data_type should be either simulated or synthetic")
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

    if gt_viscosity_type == "constant"
        
        _range = 1:(round(Int,sim_time_exp/t_steps_exp)+1)
        @info "Considering from frame $(first(_range)) to frame $(last(_range)) in the observations"
        ObsDataList = ObsDataList[_range] # align the observation points with the simulation time

        conditions = Conditions(camera_matrix=camera_matrix, obj_pose=obj_pose, SIDES=SIDES, filepath=exp_path, ANIMATE=false)

        time = collect(Float64, range(start=0, stop=sim_time_exp, step=t_steps_exp))
        if noiseLevel == 0.0
                
            obsBorderPts, nSplinex, nSpliney, pd = add_noise(ObsDataList, nFactor=0.0)

            # simulate the gt model with the estimated parameters
            est_μ_list, gradList, borderPts2DList, fields, pos3D, pos2D, splinep, splineq = simulate(model_gt, scene_gt, conditions)

            animate_fields(filepath=string(exp_path,"/Results/plots"), p=splinep, q=splineq, pObs=nSplinex, qObs=nSpliney)

            # test the closest point function
            d, pairs = closest_point(borderPts2DList, obsBorderPts)

            write_csv(string(exp_path,"/Results/data/cost_iter"), costList)

            plot_error = set_plot(fs, sz=(1650, 1250))
            Plots.plot!(plot_error, time, d, label="Contour point error", dpi=400, lw=3)
            Plots.savefig(plot_error, string(exp_path,"/Results/plots/contour_error.pdf"))

        # @save string(exp_path,"/Results/data/sim_data/Cost_Matrices.jld2") ηList, βList, CostMat, ∂CostMat, ∂2CostMat
        else
            n_samples = 10
            error_list = Matrix{Float64}(undef, n_samples, round(Int,sim_time_exp/t_steps_exp)+1)
            ANIMATED = false

            set_file(string(exp_path,"/Results/plots"))
            for n::Int in 1:n_samples
                obsBorderPts, nSplinex, nSpliney, pd = add_noise(ObsDataList, nFactor=noiseLevel)

                # simulate the model with the estimated parameters
                est_μ_list, gradList, borderPts2DList, fields, pos3D, pos2D, splinep, splineq = simulate(model_gt, scene_gt, conditions)

                if ANIMATED == false
                    ANIMATED = true
                    animate_fields(filepath=string(exp_path,"/Results/plots"), p=splinep, q=splineq, pObs=nSplinex, qObs=nSpliney)
                    set_plot(fs, sz=(1650, 1250))
                    Plots.xlims!(-5,5)
                    plot!(x->pdf(pd, x), label="", dpi=400, lw=3)
                    savefig(string(exp_path,"/Results/plots/obs_pdf.pdf"))
                end
                error, pairs = closest_point(borderPts2DList, obsBorderPts)
                error_list[n, :] = error
            end
            
            write_csv(string(exp_path,"/Results/data/error_list"), error_list)
            plot_error = set_plot(fs, sz=(1650, 1250))
            StatsPlots.errorline!(plot_error, error_list', label="Contour point error with noise level $(noiseLevel)", dpi=400, lw=3)
            Plots.savefig(plot_error, string(exp_path,"/Results/plots/contour_error.pdf"))
        end

    elseif gt_viscosity_type == "bulk_viscosity"
        av_η::Float64 = 0.0
        gt_time_frame::Int = round(Int,sim_time_gt/t_steps_gt,)
        error_list = Vector{Float64}(undef,gt_time_frame)
        obsBorderPts, nSplinex, nSpliney, pd = add_noise(ObsDataList, nFactor=0.0)
        
        viscosity_type = "constant"
        
        conditions = Conditions(camera_matrix=camera_matrix, obj_pose=obj_pose)
        model, scene = def_problem(r, h, ne_exp, η_gt, ndim, FunctionClass_u, nDof_u, FunctionClass_p, nDof_p, FunctionClass_x, β_gt, F, control, viscosity_type, 
        sim_time_exp, t_steps_exp)
        conditions = Conditions(camera_matrix=camera_matrix, obj_pose=obj_pose)
        
        set_file(string(exp_path,"/Results/plots"))
        
        time_windows, windows, data_ranges_, t_windows = set_time_window(1/t_steps_exp, obsBorderPts, method="fixed", window_size=sim_time_exp)
        _, splinexObs_win, _, _ = set_time_window(1/t_steps_exp, splinexObs, method="fixed", window_size=sim_time_exp)
        _, splineyObs_win, _, _ = set_time_window(1/t_steps_exp, splineyObs, method="fixed", window_size=sim_time_exp)

        println("Time windows: $(time_windows)")

        for ti::Int in 1:length(windows)
            data_range_ = data_ranges_[ti]
            scene_gt.sim_time = time_windows[ti]
            scene_gt.cParam = F[data_range_]
            obsBorderPts_t = windows[ti] # align the observation points with the simulation time
            η_gt_ = model_gt.η[data_range_]
            av_η = mean(η_gt_)

            println("Data frame : $(data_range_)")
            println("Time frame : $(scene.sim_time)")
            @info "Time window $(t_windows[ti])"
            printstyled("Time window: $(ti), time frames: $(scene.sim_time)\n"; color = :blue)
            println("observation window size: $(size(obsBorderPts_t,1))")

            obsBorderPts_t = windows[ti]

            est_μ_list, gradList, borderPts2DList, fields, pos3D, pos2D, splinex, spliney = simulate(model_gt, scene_gt, conditions)

            error, pairs = closest_point(borderPts2DList, obsBorderPts)

            error_list[data_range_] = error
            
            update_model!(model_gt)
        end

        @info "Completed all time windows."

        plt_error = set_plot(fs, sz=(1650, 1250))
        for ti::Int in 1:size(data_ranges_, 1)
            t = t_windows[ti]
            Plots.vline!(plt_error, [t], color=:gray, lw=3, linestyle=:dash, label=false)
        end
        Plots.xlabel!(L"\mathrm{Time\;(s)}")
        Plots.ylabel!(L"\eta\mathrm{(t)\;(KPa\cdot s)}")
        Plots.savefig(plt_error, string(exp_path,"/Results/plots/η_gt.pdf"))
        t_prev = 0.1
        for ti::Int in 1:length(data_ranges_)
            t = t_windows[ti]
            data_range_ = data_ranges_[ti]
            t_win = collect(range(start=t_prev, stop=t, step=t_steps_gt))
            if ti == 1
                Plots.plot!(plt_error, t_win, error_list[data_range_], label="Estimated η(t)", lw=3, color=:orange)
            else
                Plots.plot!(plt_error, t_win, error_list[data_range_], lw=3, color=:orange, label=false)
            end
            t_prev = t+t_steps_gt
        end
        Plots.savefig(plt_error, string(exp_path,"/Results/plots/η.pdf"))
        
        viscosity_type = "bulk_viscosity"
        est_model, est_scene = def_problem(r, h, ne_exp, η_gt, ndim, FunctionClass_u, nDof_u, FunctionClass_p, nDof_p, FunctionClass_x, β_gt, F, control, viscosity_type, 
        sim_time_gt, t_steps_gt)
        est_μ_list, gradList, simBorderPts, splinex, spliney, pos2D = simulate(est_model, est_scene, conditions)
        
        error_glob, paris = closest_point(simBorderPts, obsBorderPts)
        animate_fields(filepath=string(exp_path,"/Results/plots"), p=splinex, q=spliney, pObs=splinexObs, qObs=splineyObs)

        write_csv(string(exp_path,"/Results/data/window_data/error_global"), error_glob)
        
    end
end

function compare()

    FunctionClass_x_List = ["Q2"]
    # refine_list = [1, 2, 3] # refinement levels, ne = ne_exp^refine
    refine_list = [4, 6, 8] # refinement levels, ne = ne_exp^refine
    noise_level_list = [0.0, 0.5, 1.0]
    control = "force" # "force" or "velocity"
    viscosity_type_list = ["constant"]
    sim_time_exp_list = [30.0]

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
            for noise_level in noise_level_list 
                for ne in refine_list
                    if ne == 4 && viscosity_type == "constant"
                        noise_level_list = [0.0]
                    else
                        noise_level_list = [0.0]
                    end
                    for sim_time_exp::Float16 in sim_time_exp_list
                        # ne = ne_exp^ref
                        @info "Running optimization with ne = $ne and simulation time = $sim_time_exp with noise level = $noise_level"
                        for FunctionClass_x in FunctionClass_x_List

                            start_time = Dates.now()
                            filepath_res = string("/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/experiments/sim_data/mesh_comparisons/Stokes/$control/$viscosity_type/Q2_16/$dir/$(FunctionClass_x)_$(ne)/simtime_$(sim_time_exp)/noise_$(noise_level)/multi_window")

                            exp_params = Dict("FunctionClass_x" => FunctionClass_x, "FunctionClass_u" => "Q2", "FunctionClass_p" => "Q1", "ne_exp" => ne, "sim_time_exp" => sim_time_exp, 
                            "filepath_res" => filepath_res, "filepath_gt"=>filepath_gt, "control" => control, "data_type"=>"simulated", "camera_matrix" => camera_matrix, "WRITE_GT"=> false,
                            "noise_level"=>noise_level, "mode"=>"multi_window")

                            simulate_data(exp_params)

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

compare()