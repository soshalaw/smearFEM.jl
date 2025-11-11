using smearFEM
using LaTeXStrings
using DelimitedFiles
using Plots

global fs = 32  # font size for plots

function compare(exp_params::Dict)

    η_gt::Float64 = 0.0
    β_gt::Float64 = 0.0
    F_ext::Float64 = 0.0
    sim_time_gt::Float64 = 0.0
    t_steps_gt::Float64 = 0.0
    steps_exp::Float64 = 0.0
    outlier_frames::Vector{Int} = Int[]

    r::Float64 = 0.0
    h::Float64 = 0.0

    ndim::Int = 3
    nDof_p::Int = 1  # number of degree of freedom per node
    nDof_u::Int = ndim  # number of degree of freedom per node

    FunctionClass_u::String = exp_params["FunctionClass_u"]
    FunctionClass_p::String = exp_params["FunctionClass_p"]    
    FunctionClass_x::String = exp_params["FunctionClass_x"]
    
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
        
        params = read_json(string(filepath_gt,"/data/sim_params.json"))
        
        r = params["r"]
        h = params["h"]
        
        η_gt = params["η"][1]
        β_gt = params["β"][1]
        gt_viscosity_type = params["viscosity_type"]
        F = Array(params["cParam"])
        
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
        
    elseif data_type == "physical"
        error("Physical data type is not yet implemented")
    else
        error("data_type should be either simulated or physical")
    end

    exp_path = string(filepath_res)
    set_file(exp_path)
    write_json(string(exp_path,"/Results/data/experiment_parameters"), exp_params)

    if t_steps_exp > t_steps_gt
        @info "time resolution of the ground truth $t_steps_gt is larger than the experimental $t_steps_exp, switching to ground truth resolution"
        t_steps_exp = t_steps_gt
    end

    # Read the gt data
    obsBorderPts, nSplinex, nSpliney, pd = add_noise(ObsDataList, nFactor=noiseLevel)

    if gt_viscosity_type == "constant"

        printstyled("Ground truth η : $(η_gt), ground truth β: $(β_gt)\n"; color = :green)

        model, scene = def_problem(r, h, ne_exp, η_gt, ndim, FunctionClass_u, nDof_u, FunctionClass_p, nDof_p, FunctionClass_x, β_gt, F, control, gt_viscosity_type, 
                        sim_time_exp, t_steps_exp)
        conditions = Conditions(camera_matrix=camera_matrix, obj_pose=obj_pose, SIDES=SIDES, filepath=exp_path, ANIMATE=false)

        est_μ_list, gradList, borderPts2DList, fields, pos3D, pos2D, splinep, splineq = simulate(model, scene, conditions)

        error, _ = closest_point(borderPts2DList, obsBorderPts)

        # simulate the gt model with the estimated parameters

        animate_fields(filepath=string(exp_path,"/Results/plots"), p=splinep, q=splineq, pObs=splinexObs, qObs=splineyObs)

        est_h = get_height(est_μ_list, h)
        gt_h_ = readdlm(string(filepath_gt,"/data/h.csv"), ',', Float64)
        gt_h = gt_h_[1:(round(Int,sim_time_exp/t_steps_exp)+1)]

        set_file(string(exp_path,"/Results/plots"))

        set_plot(fs)
        time = collect(range(start=0.0, stop=sim_time_exp, step=t_steps_exp))
        Plots.plot!(time, error, label="Border points error", dpi=400, lw=3)
        Plots.xlabel!(L"\mathrm{Time\;(s)}")
        Plots.ylabel!("Border Points Error (px)")
        Plots.savefig(string(exp_path,"/Results/plots/border_points_error.pdf"))

        set_plot(fs)
        Plots.plot!(time, est_h, label="Estimated height", dpi=400, lw=3)
        Plots.plot!(time, gt_h, label="Ground truth height", dpi=400, lw=3)
        Plots.xlabel!(L"\mathrm{Time\;(s)}")
        Plots.ylabel!(L"\mathrm{Height\;(mm)}")
        Plots.savefig(string(exp_path,"/Results/plots/h_est.pdf"))

        set_plot(fs)
        Plots.plot!(time, abs.(est_h-gt_h), label="Height estimation error", dpi=400, lw=3)
        Plots.xlabel!(L"\mathrm{Time\;(s)}")
        Plots.ylabel!("Height Error (mm)")
        # Plots.xticks!(0:1:round(Int,sim_time_exp))
        Plots.savefig(string(exp_path,"/Results/plots/h_est_error.pdf"))

        write_csv(string(exp_path,"/Results/data/error"), error)
        write_csv(string(exp_path,"/Results/data/est_h"), est_h)
        write_csv(string(exp_path,"/Results/data/gt_h"), gt_h)

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
        for ti::Int in 1:window

            range = (round(Int,sim_time_frame*(titer-1))+1):(round(Int,sim_time_frame*(titer))+1)
            range_ = (round(Int,sim_time_frame*(titer-1))+1):(round(Int,sim_time_frame*(titer)))

            obsBorderPts_t = obsBorderPts[range] # align the observation points with the simulation time

            println("Time frame : $range")
            printstyled("Time window: $(ti)\n"; color = :green)

            if data_type != "physical"
                η_gt_ = model_gt.η[range_]
                av_η = mean(η_gt_)

                printstyled("Average ground truth η in the window: $(av_η), ground truth β: $(β_gt)\n"; color = :green)
            end

            stats = fit_model(model, scene, conditions, obsBorderPts_t, θ, outliers=outlier_frames)

            est_ηpList[range_] .= stats["η"]
            est_βpList[range_] .= stats["β"]
            titer = titer + 1

            θ[1] = stats["η"]
            θ[2] = stats["β"]

            update_model!(model)
        end
        set_file(string(exp_path,"/Results/plots"))

        plt_η = set_plot(fs)
        t_windows = collect(range(start=t_steps_gt, stop=sim_time_gt, step=t_steps_gt))
        if data_type != "physical"
            Plots.plot!(plt_η, t_windows, model_gt.η, label="Ground truth η(t)", dpi=400, lw=3)
        end
        Plots.plot!(plt_η, t_windows, est_ηpList, label="Estimated η(t)", lw=3)
        for t in windows
            Plots.vline!(plt_η, [t], color=:gray, lw=3, linestyle=:dash, label=false)
        end
        Plots.xlabel!(L"\mathrm{Time\;(s)}")
        Plots.ylabel!(L"\eta(\mathrm{t})\;\mathrm{(KPa\cdot s)}")
        Plots.savefig(plt_η, string(exp_path,"/Results/plots/η.pdf"))

        viscosity_type = "bulk_viscosity"
        est_model, est_scene = def_problem(r, h, ne_exp, η_gt, ndim, FunctionClass_u, nDof_u, FunctionClass_p, nDof_p, FunctionClass_x, β_gt, F, control, viscosity_type, 
                            sim_time_gt, t_steps_gt)
        est_model.η = est_ηpList
        est_μ_list, gradList, simBorderPts, splinex, spliney, pos2D = simulate(est_model, est_scene, conditions)

        animate_fields(filepath=string(exp_path,"/Results/plots"), p=splinex, q=spliney, pObs=splinexObs, qObs=splineyObs)
        
        est_h_list = get_height(est_μ_list, h)
        gt_h_list = get_height(gt_μ_list, h)

        plt_h = set_plot(fs)
        Plots.plot!(gt_h_list, label="Ground truth height", dpi=400, lw=3)
        Plots.plot!(est_h_list, label="Estimated height")
        Plots.xlabel!(L"\mathrm{Time\;(s)}")
        Plots.ylabel!(L"\mathrm{Height\;(mm)}")
        Plots.savefig(string(exp_path,"/Results/plots/h.pdf"))

        plt_error = set_plot(fs)
        Plots.plot!(abs.(est_h_list-gt_h_list), label="Height estimation error", dpi=400, lw=3)
        Plots.savefig(string(exp_path,"/Results/plots/h_est_error.pdf"))
        Plots.xlabel!(L"\mathrm{Time\;(s)}")
        Plots.ylabel!(L"\mathrm{Height\;Error\;(px)}")
        Plots.savefig(string(exp_path,"/Results/plots/error.pdf"))

        write_csv(string(exp_path,"/Results/data/est_η"), est_ηpList)
        write_csv(string(exp_path,"/Results/data/est_β"), est_βpList)
        write_csv(string(exp_path,"/Results/data/η_gt"), model_gt.η)
        write_csv(string(exp_path,"/Results/data/β_gt"), β_gt)
        write_csv(string(exp_path,"/Results/data/est_h"), est_h_list)
        write_csv(string(exp_path,"/Results/data/gt_h"), gt_h_list)

        write_data(string(exp_path,"/Results/data/sim_data/2D_surface_points"), pos2D)
        write_data(string(exp_path,"/Results/data/sim_data/3D_points"), pos3D)
        write_data(string(exp_path,"/Results/data/sim_data/motion_fields "), fields)
        write_data(string(exp_path,"/Results/data/sim_data/2D_border_points"), borderPts2DList)

    end
end


function test_model_accuracy()

    println("Running model accuracy tests...")
    FunctionClass_x_List = ["Q2"]
    refine_list = [2, 4, 6, 8, 12] #, 12, 16] # refinement levels, ne = ne_exp^refine
    noise_level_list = [0.0] # 0.5 1.0]
    control = "force" # "force" or "velocity"
    viscosity_type_list = ["constant"]
    data_type_list = ["synthetic", "simulated"] # "synthetic" or "simulated"
    
    camera_matrix::AbstractArray = [[2.39642674e+03, 0.0, 1.00429248e+03] [0.0, 2.40565353e+03, 7.57028161e+02] [0.0, 0.0, 1.0]]'
    sim_time_exp::Float64 = 20.0 # simulation time in seconds
    filepath_res::String = ""

    for data_type in data_type_list
        @info "Running tests for data type: $data_type"
        file_id = 1
        for viscosity_type in viscosity_type_list
            _filepath_gt = string("/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/ground_truth/sim_data/Stokes/$control/$viscosity_type/Q2_16/")
            dir_list = readdir(_filepath_gt)
            for dir in dir_list
                filepath_gt = string(_filepath_gt,"/",dir)
                for noise_level in noise_level_list
                    for ne in refine_list
                        filepath_res = string("/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/experiments/model_accuracy_tests/Stokes/$control/$viscosity_type/Q2_16/$dir/Q2_$ne/$data_type")
                        # ne = ne_exp^ref
                        @info "Running optimization with ne = $ne"
                        for FunctionClass_x in FunctionClass_x_List
                            @info "Running optimization with FunctionClass_x = $FunctionClass_x with $ne elements"

                            exp_params = Dict("FunctionClass_x" => FunctionClass_x, "FunctionClass_u" => "Q2", "FunctionClass_p" => "Q1", "ne_exp" => ne, "sim_time_exp" => sim_time_exp, 
                            "filepath_res" => filepath_res, "filepath_gt"=>filepath_gt, "control" => control, "data_type"=>data_type, "camera_matrix" => camera_matrix, "WRITE_GT"=> false,
                            "noise_level"=>noise_level)

                            compare(exp_params)
                        end
                    end
                file_id = file_id + 1
                end
            end
        end
    end
end

function analyse()

    FunctionClass_x_List = ["Q2"]
    refine_list = [2, 4, 6, 8, 12] #, 12, 16] # refinement levels, ne = ne_exp^refine
    noise_level_list = [0.0] # 0.5 1.0]
    control = "force" # "force" or "velocity"
    viscosity_type_list = ["constant"]
    data_type = ["synthetic", "simulated"] # "synthetic" or "simulated"
    
    camera_matrix::AbstractArray = [[2.39642674e+03, 0.0, 1.00429248e+03] [0.0, 2.40565353e+03, 7.57028161e+02] [0.0, 0.0, 1.0]]'
    sim_time_exp::Float64 = 20.0 # simulation time in seconds
    filepath_res::String = ""

   for viscosity_type in viscosity_type_list
        _filepath_gt = string("/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/ground_truth/sim_data/Stokes/$control/$viscosity_type/Q2_16/")
        dir_list = readdir(_filepath_gt)
        for dir in dir_list
            filepath_gt = string(_filepath_gt,"/",dir)
            for noise_level in noise_level_list
                for ne in refine_list
                    filepath_res = string("/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/experiments/model_accuracy_tests/Stokes/$control/$viscosity_type/Q2_16/$dir/Q2_$ne/")
                    # ne = ne_exp^ref
                    @info "Running optimization with ne = $ne"
                    for FunctionClass_x in FunctionClass_x_List
                        @info "Running optimization with FunctionClass_x = $FunctionClass_x with $ne elements"

                        post_analysis(filepath_res)
                    end
                end
            end
        end
    end

end

function post_analysis(filepath::String)

    run_filepath = readdir(filepath)
    println("Processing results in: ", filepath)

    
    fig1 = set_plot(fs)
    # Plots.ylims!(fig1, η_gt[1]*0.5, η_gt[1]*1.5)
    
    fig2 = set_plot(fs)
    # Plots.ylims!(fig2, abs(β_gt[1]-30), β_gt[1]*1.05)
    
    for run_folder_ in run_filepath
        if run_folder_ == "post_analysis"
            continue
        end
        run_folder = string(filepath,run_folder_)
        error = readdlm(string(run_folder,"/Results/data/error.csv"), ',', Float64)

        println("Processing folder: ", run_folder)
        params = read_json(string(run_folder,"/Results/data/experiment_parameters.json"))

        FunctionClass_x = params["FunctionClass_x"]
        ne = params["ne_exp"]

        # η = readdlm(string(run_folder,"/Results/data/η.csv"), ',', Float64)
        # β = readdlm(string(run_folder,"/Results/data/β.csv"), ',', Float64)

        Plots.plot!(fig1, error, label=string("Data type - ",run_folder_), marker=2, dpi=400, lw=3)
        xlabel!(fig1, L"\mathrm{Time\;(s)}")
        ylabel!(fig1, L"error\;(px)")

    end

    plot_path = string(filepath,"/post_analysis/plots")
    set_file(plot_path)

    @info "Saving plots to $plot_path"
    Plots.savefig(fig1, string(plot_path,"/error.pdf"))

end

test_model_accuracy()
analyse()