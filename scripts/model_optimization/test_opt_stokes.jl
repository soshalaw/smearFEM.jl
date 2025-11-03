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
    noiseLevel::Float64 = 0
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

    exp_path = string(filepath_res)
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
                        sim_time_exp, t_steps_exp)
        conditions = Conditions(camera_matrix=camera_matrix, obj_pose=obj_pose, SIDES=SIDES, filepath=exp_path, ANIMATE=false)

        stats = fit_model(model, scene, conditions, obsBorderPts, θ, outliers=outlier_frames)

        time = collect(Float64, range(start=0, stop=sim_time_exp, step=t_steps_exp))

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
        gt_h = gt_h_[1:(round(Int,sim_time_exp/t_steps_exp)+1)]

        set_file(string(exp_path,"/Results/plots"))
        set_plot(fs)
        Plots.plot!(time, est_h, label="Estimated height", dpi=400, lw=3)
        Plots.plot!(time, gt_h, label="Ground truth height", dpi=400, lw=3)
        Plots.xticks!(0:1:round(Int,sim_time_exp))
        Plots.xlabel!(L"\mathrm{Time\;(s)}")
        Plots.ylabel!(L"\mathrm{Height\;(mm)}")
        Plots.savefig(string(exp_path,"/Results/plots/h_est.pdf"))

        set_plot(fs)
        Plots.plot!(time, abs.(est_h-gt_h), label="Height estimation error", dpi=400, lw=3)
        Plots.xlabel!(L"\mathrm{Time\;(s)}")
        Plots.ylabel!("Height Error (mm)")
        Plots.xticks!(0:1:round(Int,sim_time_exp))
        Plots.savefig(string(exp_path,"/Results/plots/h_est_error.pdf"))

        # Plot the cost function with iterations
        set_plot(fs)
        Plots.plot!(iterList, costList, label="Cost", marker=2, dpi=400, yscale=:log10, xminorgrid = :false, lw=3)
        Plots.xlabel!(L"\mathrm{Iterations}")
        Plots.ylabel!(L"\mathrm{Cost\;(px)}")
        Plots.xticks!(minimum(iterList):2:maximum(iterList))
        Plots.savefig(string(exp_path,"/Results/plots/cost_steps.pdf"))

        # Plot the cost function with iterations
        set_plot(fs)
        Plots.plot!(iterList, costList, label="Cost", marker=2, dpi=400, yscale=:log10, xscale=:log10, lw=3)
        Plots.xlabel!(L"\mathrm{Iterations}")
        Plots.ylabel!(L"\mathrm{Cost\;(px)}")
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

        # for i::Int in η_iter
        #     η = ηList[i]
        #     for j::Int in β_iter
        #         β = βList[j]
        #         reset_model!(model)
        #         model.η = [η]
        #         scene.β = [β]
        #         μ_list, gradList, simBorderPts, fields, pos3D, pos2D, splinex, spliney = simulate(model, scene, conditions)
                
        #         # test the closest point function
        #         d_cp, pairs = closest_point(simBorderPts, obsBorderPts) 
                
        #         CostMat[i,j] = sum(d_cp)/length(d_cp)
        #     end
        # end
        
        # Plot the cost function surface
        set_plot(fs)
        Plots.contour!(ηList, βList, CostMat, color=:turbo, fill=false, levels=100, dpi=400)
        Plots.plot!(ηpList, βpList, label="Estimations", ms=:4, m=:x, color=:royalblue, lw=3)
        Plots.plot!([η_gt], [β_gt], label="Ground truth", ms=:8, m=:star5, color=:indianred2, lw=3)
        Plots.xlabel!(L"\eta\;\mathrm{(KPa\cdot s)}")
        Plots.ylabel!(L"\beta\;\mathrm{(L/mm^{-1})}")
        Plots.savefig(string(exp_path,"/Results/plots/cost_surface_iter.pdf"))

        set_plot(fs)
        Plots.contourf!(ηList, βList, CostMat, color=:turbo, fill=false, levels=100, dpi=400)
        Plots.plot!([η_gt], [β_gt], label="Ground truth", ms=:8, m=:star5, color=:indianred2, lw=3)
        Plots.xlabel!(L"\eta\;\mathrm{(KPa\cdot s)}")
        Plots.ylabel!(L"\beta\;\mathrm{(L/mm^{-1})}")
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

function compare_stats(filepath,filepath_gt)

    path_1 = []
    path_2 = []
    path_3 = []
    path_4 = []
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

        plot_η = set_plot(fs)
        plot_β = set_plot(fs)
        plot_cost = set_plot(fs)
        plot_error = set_plot(fs)

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

        plot_η = set_plot(fs)
        plot_β = set_plot(fs)
        plot_cost = set_plot(fs)
        plot_cost_log = set_plot(fs)
        plot_error = set_plot(fs)
        plot_h = set_plot(fs)

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
                xlabel!(plot_cost, LL"\mathrm{Iterations}")
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

        plot_η = set_plot(fs)
        plot_β = set_plot(fs)
        plot_cost = set_plot(fs)
        plot_cost_log = set_plot(fs)
        plot_error = set_plot(fs)
        plot_h = set_plot(fs)

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

        plot_η = set_plot(fs)
        plot_β = set_plot(fs)
        plot_cost = set_plot(fs)
        plot_cost_log = set_plot(fs)
        plot_error = set_plot(fs)
        plot_h = set_plot(fs)

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
end

function plot_(filepath, filepath_gt)
    
    run_filepath = readdir(string(filepath))

    sim_params = read_json(string(filepath_gt,"/data/sim_params.json")) 
    for exp_path_ in run_filepath
        if exp_path_ == "post_analysis"
            continue
        end

        exp_path = string(filepath, exp_path_)

        sim_time = 10.0 # simulation time in seconds 
        steps = 25.0 # number of time steps
        t_steps = sim_time/steps

        time = collect(Float64, range(start=0, stop=sim_time, step=t_steps))

        viscosity_type = sim_params["viscosity_type"]

        if viscosity_type == "constant"
            η_gt = sim_params["η"]
            β_gt = sim_params["β"]

            est_η = readdlm(string(exp_path,"/Results/data/η.csv"), ',', Float64)
            est_β = readdlm(string(exp_path,"/Results/data/β.csv"), ',', Float64)
            est_h = readdlm(string(exp_path,"/Results/data/est_h.csv"), ',', Float64)
            gt_h = readdlm(string(exp_path,"/Results/data/gt_h.csv"), ',', Float64)

            stats = read_json(string(exp_path,"/Results/data/stats.json")) 
            contour_plot_params = read_json(string(exp_path,"/Results/data/contour_plot_params.json")) 

            ηList = contour_plot_params["η_list"]
            βList = contour_plot_params["β_list"]
            CostMat = contour_plot_params["cost_mat"]

            costList = stats["cost_list"]
            iterList = stats["iterList"]

            set_plot(fs)
            Plots.plot!(time, est_h, label="Estimated height", dpi=400, lw=3)
            Plots.plot!(time, gt_h, label="Ground truth height", dpi=400, lw=3)
            Plots.xticks!(0:1:round(Int,sim_time))
            Plots.xlabel!(L"\mathrm{Time\;(s)}")
            Plots.ylabel!(L"\mathrm{Height\;(mm)}")
            Plots.savefig(string(exp_path,"/Results/plots/h_est.pdf"))

            set_plot(fs)
            Plots.plot!(time, abs.(est_h-gt_h), label="Height estimation error", dpi=400, lw=3)
            Plots.xlabel!(L"\mathrm{Time\;(s)}")
            Plots.ylabel!("Height Error (mm)")
            Plots.xticks!(0:1:round(Int,sim_time))
            Plots.savefig(string(exp_path,"/Results/plots/h_est_error.pdf"))

            set_plot(fs)
            Plots.plot!(est_η, label="Estimated η", dpi=400, lw=3)
            Plots.hline!([η_gt], label="Ground truth η", lw=3)
            Plots.xlabel!(L"\mathrm{Time\;(s)}")
            Plots.ylabel!(L"\eta\;\mathrm{(KPa\cdot s)}")
            Plots.savefig(string(exp_path,"/Results/plots/η.pdf"))

            set_plot(fs)
            Plots.plot!(est_β, label="Estimated η", dpi=400, lw=3)
            Plots.hline!([β_gt], label="Ground truth β", lw=3)
            Plots.xlabel!(L"\mathrm{Time\;(s)}")
            Plots.ylabel!(L"\beta\;\mathrm{mm^{-1}}")
            Plots.savefig(string(exp_path,"/Results/plots/β.pdf"))

            # Plot the cost function with iterations
            set_plot(fs)
            Plots.plot!(iterList, costList, label="Cost", marker=2, dpi=400, yscale=:log10, xminorgrid = :false, lw=3)
            Plots.xlabel!(L"\mathrm{Iterations}")
            Plots.ylabel!(L"\mathrm{Cost\;(px)}")
            Plots.xticks!(minimum(iterList):2:maximum(iterList))
            Plots.savefig(string(exp_path,"/Results/plots/cost_steps.pdf"))

            # Plot the cost function with iterations
            set_plot(fs)
            Plots.plot!(iterList, costList, label="Cost", marker=2, dpi=400, yscale=:log10, xscale=:log10, lw=3)
            Plots.xlabel!(L"\mathrm{Iterations}")
            Plots.ylabel!(L"\mathrm{Cost\;(px)}")
            Plots.savefig(string(exp_path,"/Results/plots/cost_steps_log.pdf"))

            # Plot the cost function surface
            set_plot(fs)
            Plots.contour!(ηList, βList, CostMat, color=:turbo, fill=false, levels=100, dpi=400)
            Plots.plot!(est_η, est_β, label="Estimations", ms=:4, m=:x, color=:royalblue, lw=3)
            Plots.plot!([η_gt], [β_gt], label="Ground truth", ms=:8, m=:star5, color=:indianred2, lw=3)
            Plots.xlabel!(L"\eta\;\mathrm{(KPa\cdot s)}")
            Plots.ylabel!(L"\beta\;\mathrm{mm^{-1}}")
            Plots.savefig(string(exp_path,"/Results/plots/cost_surface_iter.pdf"))

            set_plot(fs)
            Plots.contourf!(ηList, βList, CostMat, color=:turbo, fill=false, levels=100, dpi=400)
            Plots.plot!([η_gt], [β_gt], label="Ground truth", ms=:8, m=:star5, color=:indianred2, lw=3)
            Plots.xlabel!(L"\eta\;\mathrm{(KPa\cdot s)}")
            Plots.ylabel!(L"\beta\;\mathrm{mm^{-1}}")
            Plots.savefig(string(exp_path,"/Results/plots/cost_surface.pdf"))

        elseif viscosity_type == "bulk_viscosity"
            sim_time_gt = 40.0 # simulation time in seconds
            steps_gt = 100.0 # number of time steps
            t_steps_gt = sim_time_gt/steps_gt

            est_ηpList = readdlm(string(exp_path,"/Results/data/est_η.csv"), ',', Float64)
            est_βpList = readdlm(string(exp_path,"/Results/data/est_β.csv"), ',', Float64)
            η_gt = readdlm(string(exp_path,"/Results/data/η_gt.csv"), ',', Float64)
            β_gt = readdlm(string(exp_path,"/Results/data/β_gt.csv"), ',', Float64)
            est_h_list = readdlm(string(exp_path,"/Results/data/est_h.csv"), ',', Float64)
            gt_h_list = readdlm(string(exp_path,"/Results/data/gt_h.csv"), ',', Float64)

            gt_time_frame::Int = round(Int,sim_time_gt/t_steps_gt)
            sim_time_frame::Int = round(Int,sim_time/t_steps)
            window::Int = gt_time_frame/sim_time_frame

            println("Number of time windows: $window")

            windows = collect(range(start=t_steps_gt, stop=sim_time_gt, length=window+1))

            plt_η = set_plot(fs)
            t_windows = collect(range(start=t_steps_gt, stop=sim_time_gt, step=t_steps_gt))
            Plots.plot!(plt_η, t_windows, η_gt, label="Ground truth η(t)", dpi=400, lw=3)
            Plots.plot!(plt_η, t_windows, est_ηpList, label="Estimated η(t)", lw=3)
            for t in windows
                println(t)
                Plots.vline!(plt_η, [t], color=:gray, lw=3, linestyle=:dash, label=false)
            end
            Plots.xlabel!(L"\mathrm{Time\;(s)}")
            Plots.ylabel!(L"\eta(\mathrm{t})\;\mathrm{(KPa\cdot s)}")
            Plots.savefig(plt_η, string(exp_path,"/Results/plots/η.pdf"))

            plt_η_gt = set_plot(fs)
            t_windows = collect(range(start=t_steps_gt, stop=sim_time_gt, step=t_steps_gt))
            Plots.plot!(plt_η_gt, t_windows, η_gt, label="Ground truth η(t)", dpi=400, lw=3)
            # Plots.plot!(plt_η, t_windows, est_ηpList, label="Estimated η(t)", lw=3)
            for t in windows
                Plots.vline!(plt_η_gt, [t], color=:gray, lw=3, linestyle=:dash, label=false)
            end
            Plots.xlabel!(L"\mathrm{Time\;(s)}")
            Plots.ylabel!(L"\eta(\mathrm{t})\;\mathrm{(KPa\cdot s)}")
            Plots.savefig(plt_η_gt, string(exp_path,"/Results/plots/η_gt.pdf"))

            plt_h = set_plot(fs)
            Plots.plot!(gt_h_list, label="Ground truth height", dpi=400, lw=3)
            for t in windows
                Plots.vline!(plt_h, [t], color=:gray, lw=3, linestyle=:dash, label=false)
            end
            Plots.plot!(est_h_list, label="Estimated height")
            Plots.xlabel!(L"\mathrm{Time\;(s)}")
            Plots.ylabel!(L"\mathrm{Height\;(mm)}")
            Plots.savefig(string(exp_path,"/Results/plots/h.pdf"))

            plt_error = set_plot(fs)
            Plots.plot!(abs.(est_h_list-gt_h_list), label="Height estimation error", dpi=400, lw=3)
            for t in windows
                Plots.vline!(plt_error, [t], color=:gray, lw=3, linestyle=:dash, label=false)
            end
            Plots.savefig(string(exp_path,"/Results/plots/h_est_error.pdf"))
            Plots.xlabel!(L"\mathrm{Time\;(s)}")
            Plots.ylabel!(L"\mathrm{Height\;Error\;(px)}")
            Plots.savefig(string(exp_path,"/Results/plots/error.pdf"))

        end
    end

    return 
end

function post_analysis(filepath_gt::String, filepath::String)
    # ηList = readdlm(string(filepath_gt,"/data/η.csv"), ',', Float64)
    # βList = readdlm(string(filepath_gt,"/data/β.csv"), ',', Float64)

    params = read_json(string(filepath_gt,"/data/sim_params.json"))
    η_gt = params["η"]
    β_gt = params["β"]  

    println("Ground truth η: ", η_gt[1])
    η_gt_list = η_gt[1]*ones(30,1)
    β_gt_list = β_gt[1]*ones(30,1)

    run_filepath = readdir(string(filepath,"/runs/"))

    fig1 = Plots.plot(η_gt_list, label="Ground truth η", dpi=400)
    Plots.xlabel!(fig1,L"\mathrm{Iterations}")
    Plots.ylabel!(fig1,L"\eta\;\mathrm{(KPa\cdot s)}")
    # Plots.ylims!(fig1, η_gt[1]*0.5, η_gt[1]*1.5)


    fig2 = Plots.plot(β_gt_list, label="Ground truth β", dpi=400)
    Plots.xlabel!(fig2, L"\mathrm{Iterations}")
    Plots.ylabel!(fig2, L"\beta\;\mathrm{mm^{-1}}")
    # Plots.ylims!(fig2, abs(β_gt[1]-30), β_gt[1]*1.05)

    for run_folder_ in run_filepath
        if run_folder_ == "post_analysis"
            continue
        end
        run_folder = string(filepath,"/runs/",run_folder_)
        println("Processing folder: ", run_folder)
        params = read_json(string(run_folder,"/Results/data/sim_params.json"))

        FunctionClass_x = params["FunctionClass_x"]
        ne = params["ne_exp"]

        if ne != 2
            η = readdlm(string(run_folder,"/Results/data/η.csv"), ',', Float64)
            β = readdlm(string(run_folder,"/Results/data/β.csv"), ',', Float64)

            Plots.plot!(fig1, η, label=string("Basis - ",FunctionClass_x," - ne: ",ne))
            xlabel!(fig1, L"\mathrm{Iterations}")
            ylabel!(fig1, L"\eta\;\mathrm{(KPa\cdot s)}")

            Plots.plot!(fig2, β, label=string("Basis - ",FunctionClass_x," - ne: ",ne))
            xlabel!(fig2, L"\mathrm{Iterations}")
            ylabel!(fig2, L"\beta\;\mathrm{mm^{-1}}")

        end
    end

    plot_path = string(filepath,"/Results/plots")
    set_file(plot_path)

    @info "Saving plots to $plot_path"
    Plots.savefig(fig1, string(plot_path,"/η_comp_zoomed.pdf"))
    Plots.savefig(fig2, string(plot_path,"/β_comp_zoomed.pdf"))

end

function plot_results()
    ne_gt::Int = 8 # number of elements in the mesh for the ground truth
    ne_exp::Int = 2 # number of elements in the mesh for the experiment 
    # β_gt_list = [5, 10, 50, 100.0, 200.0, 500.0, 1000.0, 10000.0]
    # η_gt_list = [40.0]
    β_gt_list = [10.0, 50.0, 100.0, 1e3]
    # β_gt_list = [100.0]
    η_gt_list = [60.0]
    FunctionClass_x_List = ["S2", "Q2"]
    # refine_list = [1, 2, 3] # refinement levels, ne = ne_exp^refine
    refine_list = [2] # refinement levels, ne = ne_exp^refine
    control = "force" # "force" or "velocity"

    viscosity_type_list = ["constant"]
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
                        plot_(filepath, filepath_gt)
                        compare_stats(filepath,filepath_gt)
                    end
                    post_analysis(filepath_gt, filepath)
                    run_id = run_id + 1
                end
            end
        end
    end
end

function optimize_sim()
    ne_gt::Int = 4 # number of elements in the mesh for the ground truth
    ne_exp::Int = 2 # number of elements in the mesh for the experiment 
    # β_gt_list = [5, 10, 50, 100.0, 200.0, 500.0, 1000.0, 10000.0]
    # η_gt_list = [40.0]
    β_gt_list = [10.0, 50.0, 100.0, 1e3]
    # β_gt_list = [100.0]
    η_gt_list = [60.0]
    FunctionClass_x_List = ["Q2", "S2"]
    # refine_list = [1, 2, 3] # refinement levels, ne = ne_exp^refine
    refine_list = [2] # refinement levels, ne = ne_exp^refine
    control = "force" # "force" or "velocity"

    viscosity_type_list = ["bulk_viscosity"]
    FunctionClass_x_gt_list = ["S2", "Q2"] # Function space for the ground truth

    exp_size = size(FunctionClass_x_List,1)*size(refine_list,1)
    η_mat = zeros(30, exp_size)
    β_mat = zeros(30, exp_size)

    WRITE_GT = true
    for viscosity_type in viscosity_type_list
        for FunctionClass_x_gt in FunctionClass_x_gt_list
            run_id = 1
            for β_gt in β_gt_list
                for η_gt in η_gt_list

                    filepath_gt = string("/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/ground_truth/sim_data/Stokes/$control/$viscosity_type/$(FunctionClass_x_gt)_$(ne_gt)/$run_id")
                    filepath_res = string("/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/experiments/sim_data/optimization/Stokes/$control/$viscosity_type/$run_id")
                    for ref in refine_list
                        ne = ne_exp^ref
                        @info "Running optimization with ne = $ne"
                        for FunctionClass_x in FunctionClass_x_List
                            @info "Running optimization with FunctionClass_x = $FunctionClass_x with $ne elements"
    
                            exp_params = Dict("FunctionClass_x" => FunctionClass_x, "FunctionClass_u" => "Q2", "FunctionClass_p" => "Q1", "ne_gt" => ne_gt, "ne_exp" => ne, 
                                        "β_gt" => β_gt, "η_gt" => η_gt, "WRITE_GT" => WRITE_GT, "filepath" => filepath_res, "filepath_gt"=>filepath_gt, "control" => control, 
                                        "viscosity_type"=>viscosity_type, "FunctionClass_x_gt" => FunctionClass_x_gt)
    
                            optimize(exp_params)
                            WRITE_GT = false
                        end
                    end
                    WRITE_GT = true
                    post_analysis(filepath_gt, filepath_res)
                    run_id = run_id + 1
                end
            end
        end
    end
end

function optimize_real()

    ne_exp::Int = 2 # number of elements in the mesh for the experiment 
    FunctionClass_x_List = ["Q2", "S2"]
    # refine_list = [1, 2, 3] # refinement levels, ne = ne_exp^refine
    refine_list = [2] # refinement levels, ne = ne_exp^refine
    control = "force" # "force" or "velocity"
    η_start = 5.0
    β_start = 1.0
    viscosity_type_list = ["bulk_viscosity"]

    r::Float64 = 25.0  # radius of the cylinder in mm
    h::Float64 = 40.0  # height of the cylinder in mm
    camera_matrix::AbstractArray = [[2.39642674e+03, 0.0, 1.00429248e+03] [0.0, 2.40565353e+03, 7.57028161e+02] [0.0, 0.0, 1.0]]'
    sim_time_exp::Float64 = 5.0 # simulation time in seconds
    F_ext::Float64 = 1*9.812*1e3 # force applied to the cylinder in N

    for viscosity_type in viscosity_type_list
        file_id = 5
        filepath_gt = string("/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/ground_truth/physical_data/$file_id")
        filepath_res = string("/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/experiments/physical_data/optimization/Stokes/$control/$viscosity_type/$file_id")
        
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
        post_analysis(filepath_gt, filepath_res)
        file_id = file_id + 1
    end
end

function optimize_syn()

    ne_exp::Int = 2 # number of elements in the mesh for the experiment 
    FunctionClass_x_List = ["Q2"]
    # refine_list = [1, 2, 3] # refinement levels, ne = ne_exp^refine
    refine_list = [4, 6, 8] # refinement levels, ne = ne_exp^refine
    control = "force" # "force" or "velocity"
    η_start = 5.0
    β_start = 1.0
    viscosity_type_list = ["constant"]

    r::Float64 = 25.0  # radius of the cylinder in mm
    h::Float64 = 40.0  # height of the cylinder in mm
    camera_matrix::AbstractArray = [[2.39642674e+03, 0.0, 1.00429248e+03] [0.0, 2.40565353e+03, 7.57028161e+02] [0.0, 0.0, 1.0]]'
    sim_time_exp::Float64 = 20.0 # simulation time in seconds
    filepath_res::String = ""
    for viscosity_type in viscosity_type_list
        file_id = 2
        filepath_gt = string("/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/ground_truth/sim_data/Stokes/$control/$viscosity_type/Q2_16/$file_id")
        for ne in refine_list
            filepath_res = string("/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/experiments/syn_data/optimization/Stokes/$control/$viscosity_type/Q2_16/$file_id/Q2_$ne")
            # ne = ne_exp^ref
            @info "Running optimization with ne = $ne"
            for FunctionClass_x in FunctionClass_x_List
                @info "Running optimization with FunctionClass_x = $FunctionClass_x with $ne elements"

                exp_params = Dict("FunctionClass_x" => FunctionClass_x, "FunctionClass_u" => "Q2", "FunctionClass_p" => "Q1", "ne_exp" => ne, "sim_time_exp" => sim_time_exp, 
                "filepath_res" => filepath_res, "filepath_gt"=>filepath_gt, "control" => control, "data_type"=>"synthetic", "camera_matrix" => camera_matrix, "WRITE_GT"=> false)

                optimize(exp_params)
            end
        end
        # post_analysis(filepath_gt, filepath_res)
        file_id = file_id + 1
    end
end

function plot_syn()

    ne_exp::Int = 2 # number of elements in the mesh for the experiment 
    FunctionClass_x_List = ["Q2"]
    # refine_list = [1, 2, 3] # refinement levels, ne = ne_exp^refine
    refine_list = [4, 6, 8] # refinement levels, ne = ne_exp^refine
    control = "force" # "force" or "velocity"
    η_start = 5.0
    β_start = 1.0
    viscosity_type_list = ["constant"]

    r::Float64 = 25.0  # radius of the cylinder in mm
    h::Float64 = 40.0  # height of the cylinder in mm
    camera_matrix::AbstractArray = [[2.39642674e+03, 0.0, 1.00429248e+03] [0.0, 2.40565353e+03, 7.57028161e+02] [0.0, 0.0, 1.0]]'
    sim_time_exp::Float64 = 20.0 # simulation time in seconds
    filepath_res::String = ""
    for viscosity_type in viscosity_type_list
        file_id = 2
        filepath_gt = string("/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/ground_truth/sim_data/Stokes/$control/$viscosity_type/Q2_16/$file_id")
        for ne in refine_list
            filepath_res = string("/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/experiments/syn_data/optimization/Stokes/$control/$viscosity_type/Q2_16/$file_id/Q2_$ne/")
            # ne = ne_exp^ref
            @info "Running optimization with ne = $ne"
            for FunctionClass_x in FunctionClass_x_List
                @info "Running optimization with FunctionClass_x = $FunctionClass_x with $ne elements"

                exp_params = Dict("FunctionClass_x" => FunctionClass_x, "FunctionClass_u" => "Q2", "FunctionClass_p" => "Q1", "ne_exp" => ne, "sim_time_exp" => sim_time_exp, 
                "filepath_res" => filepath_res, "filepath_gt"=>filepath_gt, "control" => control, "data_type"=>"synthetic", "camera_matrix" => camera_matrix, "WRITE_GT"=> false)

                # optimize(exp_params)
                plot_(filepath_res, filepath_gt)
            end
        end
        # post_analysis(filepath_gt, filepath_res)
        file_id = file_id + 1
    end
end

# main()
# plot_syn()e
optimize_syn
