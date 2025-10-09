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
    
    # simulation parameters for the ground truth
    r::Float64 = exp_params["r"]  # radius of the cylinder in mm
    h::Float64 = exp_params["h"]  # height of the cylinder in mm
    
    control::String = exp_params["control"]
    gt_viscosity_type::String = exp_params["viscosity_type"] # "constant" or "bulk_viscosity"
    filepath_gt::String = exp_params["filepath_gt"]
    filepath_res::String = exp_params["filepath_res"]

    obj_pose::AbstractArray = zeros(Float64, 4,4) # initial pose of the object
    camera_matrix::AbstractArray = exp_params["camera_matrix"]

    sim_time_exp::Float64 = exp_params["sim_time_exp"]
    if haskey(exp_params, "steps_exp")
        steps_exp = exp_params["steps_exp"]
    else
        steps_exp = sim_time_exp*30.0 # assuming 30 fps for the experiments
    end
    t_steps_exp::Float64 = sim_time_exp/steps_exp

    ne_exp::Int = exp_params["ne_exp"] # number of elements in the mesh for the experiment

    data_type ::String = exp_params["data_type"] # "simulated" or "physical" or "real"
    noiseLevel::Float64 = 0
    SIDES::Bool = false
    
    if data_type == "simulated" || data_type == "synthetic"

        WRITE_GT = exp_params["WRITE_GT"] 
        outlier_frames = Int[]

        if WRITE_GT == true # write the ground truth data
            @info "Writing ground truth gt data to with $ne_gt elements to $filepath_res"
            write_data(exp_params)
        end

        params = read_json(string(filepath_gt,"/data/sim_params.json"))
        
        r = params["r"]
        h = params["h"]

        η_gt = params["η"][1]
        β_gt = params["β"][1]
        gt_viscosity_type = params["viscosity_type"]
        F = params["cParam"]

        sim_time_gt = params["simulation_time"]
        t_steps_gt = params["time_steps"]

        camera_matrix = reshape(Array(params["camera_matrix"]), 3, 3)
        obj_pose = reshape(Array(params["obj_pose"])*1.0,3,1)
        control = params["control_type"]
        
        model_gt, scene_gt = def_problem(r, h, ne_exp, η_gt, ndim, FunctionClass_u, nDof_u, FunctionClass_p, nDof_p, FunctionClass_x, β_gt, F, control, viscosity_type, 
                                        sim_time_gt, t_steps_gt)

        ObsDataList, splinexObs, splineyObs = read_csv(string(filepath_gt,"/data/sim_data/contour_data"))  
        printstyled("Ground truth η: $(η_gt), ground truth β: $(β_gt)\n"; color = :green)

    elseif data_type == "physical"
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

    exp_path = string("$filepath_res/data_type/$FunctionClass_x","_","$ne_exp")
    set_file(exp_path)
    write_json(string(exp_path,"/Results/data/experiment_parameters"), exp_params)

    if t_steps_exp > t_steps_gt
        @info "time resolution of the ground truth $t_steps_gt is larger than the experimental $t_steps_exp, switching to ground truth resolution"
        t_steps_exp = t_steps_gt
    end

    # Read the gt data
    obsBorderPts, nSplinex, nSpliney, pd = add_noise(ObsDataList, nFactor=noiseLevel)
    
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

        # obsBorderPts = obsBorderPts[1:(round(Int,sim_time/t_steps)+1)] # align the observation points with the simulation time

        model, scene = def_problem(r, h, ne_exp, η_gt, ndim, FunctionClass_u, nDof_u, FunctionClass_p, nDof_p, FunctionClass_x, β_gt, F, control, gt_viscosity_type, 
                        sim_time, t_steps_exp)
        conditions = Conditions(camera_matrix=camera_matrix, obj_pose=obj_pose, SIDES=SIDES, filepath=exp_path, ANIMATE=false)

        stats = fit_model(model, scene, conditions, obsBorderPts, θ, outliers=outlier_frames)

        time = collect(Float64, range(start=0, stop=sim_time, step=t_steps_exp))

        iterList::Vector{Float64} = stats["iterList"]
        costList::Vector{Float64} = stats["cost_list"]
        ηpList::Vector{Float64} = stats["ηList"]
        βpList::Vector{Float64} = stats["βList"]

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

        animate_fields(filepath=string(exp_path,"/Results/plots"), p=splinep, q=splineq, pObs=splinexObs, qObs=splineyObs)

        est_h = get_height(est_μ_list, h)
        gt_h_ = readdlm(string(filepath_gt,"/data/h.csv"), ',', Float64)
        gt_h = gt_h_[1:(round(Int,sim_time/t_steps_exp)+1)]

        set_file(string(exp_path,"/Results/plots"))
        set_plot(fs)
        Plots.plot!(time, est_h, label="Estimated height", dpi=400, lw=3)
        Plots.plot!(time, gt_h, label="Ground truth height", dpi=400, lw=3)
        Plots.xticks!(0:1:round(Int,sim_time))
        Plots.xlabel!("Time (s)")
        Plots.ylabel!("Height (mm)")
        Plots.savefig(string(exp_path,"/Results/plots/h_est.pdf"))

        set_plot(fs)
        Plots.plot!(time, abs.(est_h-gt_h), label="Height estimation error", dpi=400, lw=3)
        Plots.xlabel!("Time (s)")
        Plots.ylabel!("Height Error (mm)")
        Plots.xticks!(0:1:round(Int,sim_time))
        Plots.savefig(string(exp_path,"/Results/plots/h_est_error.pdf"))

        # Plot the cost function with iterations
        set_plot(fs)
        Plots.plot!(iterList, costList, label="Cost", marker=2, dpi=400, yscale=:log10, xminorgrid = :false, lw=3)
        Plots.xlabel!("Iterations")
        Plots.ylabel!("Cost (px)")
        Plots.xticks!(minimum(iterList):2:maximum(iterList))
        Plots.savefig(string(exp_path,"/Results/plots/cost_steps.pdf"))

        # Plot the cost function with iterations
        set_plot(fs)
        Plots.plot!(iterList, costList, label="Cost", marker=2, dpi=400, yscale=:log10, xscale=:log10, lw=3)
        Plots.xlabel!("Iterations")
        Plots.ylabel!("Cost (px)")
        Plots.savefig(string(exp_path,"/Results/plots/cost_steps_log.pdf"))

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

        sampleNo = 5
        ηList = collect(range(η_start, stop=ηStop, length=sampleNo))
        βList = collect(range(β_start, stop=βStop, length=sampleNo))
        CostMat = zeros(size(ηList,1),size(βList,1))
        costη = zeros(size(ηList,1))
        costβ = zeros(size(βList,1))

        η_iter = 1:size(ηList,1)
        β_iter = 1:size(βList,1)

        for i::Int in η_iter
            η = ηList[i]
            for j::Int in β_iter
                β = βList[j]
                reset_model!(model)
                model.η = [η]
                scene.β = [β]
                μ_list, gradList, simBorderPts, fields, pos3D, pos2D, splinex, spliney = simulate(model, scene, conditions)
                
                # test the closest point function
                d_cp, pairs = closest_point(simBorderPts, obsBorderPts) 
                
                CostMat[i,j] = sum(d_cp)/length(d_cp)
            end
        end
        
        # Plot the cost function surface
        set_plot(fs)
        Plots.contour!(ηList, βList, CostMat, color=:turbo, fill=false, levels=100, xlabel="η", ylabel="β", dpi=400)
        Plots.plot!(ηpList, βpList, label="Estimations", ms=:4, m=:x, color=:royalblue, lw=3)
        Plots.plot!([η_gt], [β_gt], label="Ground truth", ms=:8, m=:star5, color=:indianred2, lw=3)
        Plots.xlabel!(L"\eta")
        Plots.ylabel!(L"\beta")
        Plots.savefig(string(exp_path,"/Results/plots/cost_surface_iter.pdf"))

        set_plot(fs)
        Plots.contourf!(ηList, βList, CostMat, color=:turbo, fill=false, levels=100, xlabel="η", ylabel="β", dpi=400)
        Plots.plot!([η_gt], [β_gt], label="Ground truth", ms=:8, m=:star5, color=:indianred2, lw=3)
        Plots.xlabel!(L"\eta")
        Plots.ylabel!(L"\beta")
        Plots.savefig(string(exp_path,"/Results/plots/cost_surface.pdf"))

        # Write the results to files
        contour_plot_params = Dict("η_list" => ηList, "β_list" => βList, "cost_mat" => CostMat)

        write_json(string(exp_path,"/Results/data/stats"), stats)
        write_json(string(exp_path,"/Results/data/contour_plot_params"), contour_plot_params)
        write_csv(string(exp_path,"/Results/data/η"), ηpList)
        write_csv(string(exp_path,"/Results/data/β"), βpList)
        write_csv(string(exp_path,"/Results/data/est_h"), est_h)
        write_csv(string(exp_path,"/Results/data/gt_h"), gt_h)
        write_csv(string(exp_path,"/Results/data/cost_iter"), costList)

        write_data(string(exp_path,"/Results/data/sim_data/2D_surface_points"), pos2D)
        write_data(string(exp_path,"/Results/data/sim_data/3D_points"), pos3D)
        write_data(string(exp_path,"/Results/data/sim_data/motion_fields "), fields)
        write_data(string(exp_path,"/Results/data/sim_data/2D_border_points"), borderPts2DList)

    elseif gt_viscosity_type == "bulk_viscosity"

        viscosity_type = "constant"

        if t_steps_exp > t_steps_gt 
            @info "time resolution of the ground truth $t_steps_gt is not enough switching to ground truth resolution"
            t_steps_exp = t_steps_gt
        end

        conditions = Conditions(camera_matrix=camera_matrix, obj_pose=obj_pose)
        model, scene = def_problem(r, h, ne_exp, η_gt, ndim, FunctionClass_u, nDof_u, FunctionClass_p, nDof_p, FunctionClass_x, β_gt, F, control, viscosity_type, 
                            sim_time_exp, t_steps_exp)

        gt_time_frame::Int = round(Int,sim_time_gt/t_steps_gt)
        sim_time_frame::Int = round(Int,sim_time_exp/t_steps_exp)
        window::Int = round(Int,gt_time_frame/sim_time_frame)

        est_ηpList = Vector{Float64}(undef,gt_time_frame)
        est_βpList = Vector{Float64}(undef,gt_time_frame)

        println("Number of time windows: $window")

        windows = collect(range(start=1, stop=sim_time_gt, length=window))
        titer::Int = 1
        set_file(string(exp_path,"/Results/plots"))
        # for ti::Int in 1:window
            ti = 1
            data_range = (round(Int,sim_time_frame*(titer-1))+1):(round(Int,sim_time_frame*(titer))+1)
            data_range_ = (round(Int,sim_time_frame*(titer-1))+1):(round(Int,sim_time_frame*(titer)))

            obsBorderPts_t = obsBorderPts[data_range] # align the observation points with the simulation time

            println("Time frame : $data_range")
            printstyled("Time window: $(ti)\n"; color = :green)

            if data_type != "physical"
                η_gt_ = model_gt.η[data_range_]
                av_η = mean(η_gt_)

                printstyled("Average ground truth η in the window: $(av_η), ground truth β: $(β_gt)\n"; color = :green)
            end

            stats = fit_model(model, scene, conditions, obsBorderPts_t, θ, outliers=outlier_frames)

            est_ηpList[data_range_] .= stats["η"]
            est_βpList[data_range_] .= stats["β"]
            titer = titer + 1

            θ[1] = stats["η"]
            θ[2] = stats["β"]

            update_model!(model)
      
            iterList = stats["iterList"]
            costList = stats["cost_list"]
    
            # Plot the cost function with iterations
            set_plot(fs)
            Plots.plot!(iterList, costList, label="Cost", marker=2, dpi=400, yscale=:log10, xminorgrid = :false, lw=3)
            Plots.xlabel!("Iterations")
            Plots.ylabel!("Cost (px)")
            Plots.xticks!(minimum(iterList):2:maximum(iterList))
            Plots.savefig(string(exp_path,"/Results/plots/cost_steps.pdf"))
    
            # Plot the cost function with iterations
            set_plot(fs)
            Plots.plot!(iterList, costList, label="Cost", marker=2, dpi=400, yscale=:log10, xscale=:log10, lw=3)
            Plots.xlabel!("Iterations")
            Plots.ylabel!("Cost (px)")
            Plots.savefig(string(exp_path,"/Results/plots/cost_steps_log.pdf"))
        # end

        plt_η = set_plot(fs)
        t_windows = collect(range(start=t_steps_gt, stop=sim_time_gt, step=t_steps_gt))
        if data_type != "physical"
            Plots.plot!(plt_η, t_windows, model_gt.η, label="Ground truth η(t)", dpi=400, lw=3)
        end
        Plots.plot!(plt_η, t_windows, est_ηpList, label="Estimated η(t)", lw=3)
        for t in windows
            Plots.vline!(plt_η, [t], color=:gray, lw=3, linestyle=:dash, label=false)
        end
        Plots.xlabel!("Time (s)")
        Plots.ylabel!("η(t)")
        Plots.savefig(plt_η, string(exp_path,"/Results/plots/η.pdf"))

        viscosity_type = "constant"
        est_model, est_scene = def_problem(r, h, ne_exp, η_gt, ndim, FunctionClass_u, nDof_u, FunctionClass_p, nDof_p, FunctionClass_x, β_gt, F, control, viscosity_type, 
                            sim_time_exp, t_steps_exp)
        est_model.η = est_ηpList
        est_μ_list, gradList, borderPts2DList, fields, pos3D, pos2D, splinex, spliney = simulate(est_model, est_scene, conditions)

        animate_fields(filepath=string(exp_path,"/Results/plots"), p=splinex, q=spliney, pObs=splinexObs[data_range], qObs=splineyObs[data_range])
        
        est_h_list = get_height(est_μ_list, h)

        if data_type != "physical"
            gt_h_list = get_height(gt_μ_list, h)

            plt_h = set_plot(fs)
            Plots.plot!(gt_h_list, label="Ground truth height", dpi=400, lw=3)
            Plots.plot!(est_h_list, label="Estimated height")
            Plots.xlabel!("Time (s)")
            Plots.ylabel!("Height (mm)")
            Plots.savefig(string(exp_path,"/Results/plots/h.pdf"))

            plt_error = set_plot(fs)
            Plots.plot!(abs.(est_h_list-gt_h_list), label="Height estimation error", dpi=400, lw=3)
            Plots.savefig(string(exp_path,"/Results/plots/h_est_error.pdf"))
            Plots.xlabel!("Time (s)")
            Plots.ylabel!("Height Error (px)")
            Plots.savefig(string(exp_path,"/Results/plots/error.pdf"))

            write_csv(string(exp_path,"/Results/data/η_gt"), model_gt.η)
            write_csv(string(exp_path,"/Results/data/gt_h"), gt_h_list)
            write_csv(string(exp_path,"/Results/data/β_gt"), β_gt)
        end

        write_csv(string(exp_path,"/Results/data/est_η"), est_ηpList)
        write_csv(string(exp_path,"/Results/data/est_β"), est_βpList)
        write_csv(string(exp_path,"/Results/data/est_h"), est_h_list)

        write_data(string(exp_path,"/Results/data/sim_data/2D_surface_points"), pos2D)
        write_data(string(exp_path,"/Results/data/sim_data/3D_points"), pos3D)
        write_data(string(exp_path,"/Results/data/sim_data/motion_fields "), fields)
        write_data(string(exp_path,"/Results/data/sim_data/2D_border_points"), borderPts2DList)
    end
end

function optimize_real()

    ne_exp::Int = 2 # number of elements in the mesh for the experiment 
    FunctionClass_x_List = ["Q2"]
    # refine_list = [1, 2, 3] # refinement levels, ne = ne_exp^refine
    refine_list = [2] # refinement levels, ne = ne_exp^refine
    control = "force" # "force" or "velocity"
    η_start = 5.0
    β_start = 1.0
    viscosity_type = "bulk_viscosity"

    r::Float64 = 25.0  # radius of the cylinder in mm
    h::Float64 = 40.0  # height of the cylinder in mm
    camera_matrix::AbstractArray = [[2.39642674e+03, 0.0, 1.00429248e+03] [0.0, 2.40565353e+03, 7.57028161e+02] [0.0, 0.0, 1.0]]'
    # sim_time_exp::Float64 = 5.0 # simulation time in seconds
    F_ext::Float64 = 1*9.812*1e3 # force applied to the cylinder in N

    sim_time_list = [1.0, 0.75, 0.5, 0.25, 0.1] # simulation time in seconds
    for sim_time_exp in sim_time_list
        println("Simulation time: $sim_time_exp seconds")
        file_id = 5
        filepath_gt = string("/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/ground_truth/physical_data/$file_id")
        filepath_res = string("/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/experiments/physical_data/integration_tests/single_window/$sim_time_exp")
        
        for ref in refine_list
            ne = ne_exp^ref
            @info "Running optimization with ne = $ne"
            for FunctionClass_x in FunctionClass_x_List
                @info "Running optimization with FunctionClass_x = $FunctionClass_x with $ne elements"

                exp_params = Dict("FunctionClass_x" => FunctionClass_x, "FunctionClass_u" => "Q2", "FunctionClass_p" => "Q1", "ne_exp" => ne, "sim_time_exp" => sim_time_exp, 
                "η_start" => η_start, "β_start" => β_start, "filepath_res" => filepath_res, "filepath_gt"=>filepath_gt, "control" => control, "viscosity_type"=>viscosity_type, 
                "data_type"=>"physical", "r" => r, "h" => h, "camera_matrix" => camera_matrix, "F_ext" => F_ext)

                optimize(exp_params)
            end
        end
        # post_analysis(filepath_gt, filepath_res)
        file_id = file_id + 1
    end
end

# main()
optimize_real()
