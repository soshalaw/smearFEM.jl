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

function test_opt_bulk()
    scale = 50
    r::Float64 = 0.5*scale  # radius of the cylinder in mm
    h::Float64 = 1*scale  # height of the cylinder in mm
    ndim::Int = 3
    FunctionClass_x::String = "Q2"
    FunctionClass_u::String = "Q2"
    nDof_u::Int = ndim  # number of degree of freedom per node
    FunctionClass_p::String = "Q1"
    nDof_p::Int = 1  # number of degree of freedom per node
    ne_gt::Int = 12 # number of elements in the mesh for the ground truth

    camera_matrix = [[8*2048/7.07, 0.0, 2048/2] [0.0, 8*1536/5.3, 1536/2] [0.0, 0.0, 1.0]]'
    camera_pose = scale*[0 -0.25 2]'   # camera position in mm

    F::Float64 = 50000.0
    # F::Float64 = 3.0 # force applied to the cylinder in N

    dev::Float64 = 0.1

    control::String = "force" # "force" or "velocity"
    viscosity_type::String = "bulk_viscosity" # "constant" or "bulk_viscosity"
    noiseLevel::Float64 = 0.0
    filepath::String = string("/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/sim_experiments/cost_function_test/optimization/Stokes/",control,"/test_bulk_viscosity")

    # simulation parameters for the ground truth
    sim_time_gt::Float64 = 90.0
    gt_steps::Float64 = 90.0
    t_steps_gt::Float64 = sim_time_gt/gt_steps
  
    β_gt::Float64 = 100.0
    η_gt::Float64 = 60.0

    # Write the ground truth
    printstyled("Ground truth η: $(η_gt), ground truth β: $(β_gt)\n"; color = :green)
    model_gt, scene_gt = def_problem(r, h, ne_gt, η_gt, ndim, FunctionClass_u, nDof_u, FunctionClass_p, nDof_p, β_gt, F, control, viscosity_type, sim_time_gt, t_steps_gt)
    write_sim_data(model_gt, scene_gt, camera_matrix, camera_pose, filepath)
    η_gt = model_gt.η[1]

    write_csv(string(filepath,"/Results/data/η_gt"), model_gt.η)

    # Read the gt data
    ObsDataList, splinexObs, splineyObs = read_csv(string(filepath,"/Results/contour_data"))  
    nScene, nSplinex, nSpliney, pd = add_noise(ObsDataList, nFactor=noiseLevel)
    ObsData = [nScene, nSplinex, nSpliney]
    obsBorderPts = ObsData[1]

    sim_time::Float64 = 15.0
    steps::Float64 = 15.0
    t_steps::Float64 = sim_time/steps

    ne_exp::Int = 4 # number of elements in the mesh for the experiment
    dev_η::Float64 = η_gt*dev
    dev_β::Float64 = β_gt*dev

    ηStart::Float64 = η_gt - dev_η
    βStart::Float64 = β_gt - dev_β

    viscosity_type = "constant"
    conditions = Conditions(camera_matrix=camera_matrix, camera_pose=camera_pose)
    model, scene = def_problem(r, h, ne_exp, η_gt, ndim, FunctionClass_u, nDof_u, FunctionClass_p, nDof_p, FunctionClass_x, β_gt, F, control, viscosity_type, 
                        sim_time, t_steps)
    θ::Vector{Float64} = [ηStart, βStart]
    
    gt_time_frame::Int = round(Int,sim_time_gt)
    sim_time_frame::Int = round(Int,sim_time)
    windows::Int = gt_time_frame/sim_time_frame

    est_ηpList = Vector{Float64}(undef,gt_time_frame)
    est_βpList = Vector{Float64}(undef,gt_time_frame)

    println("Number of time windows: ", windows)
    titer::Int = 1
    for ti::Int in 1:windows

        println((sim_time_frame*(titer-1)+1):(sim_time_frame*titer)+1)
        printstyled("Time window: $(ti)\n"; color = :green)
        η_gt_ = model_gt.η[round(Int,(sim_time*(titer-1))+1):round(Int,sim_time*(titer))]
        av_η = mean(η_gt_)
        printstyled("Average ground truth η in the window: $(av_η), ground truth β: $(β_gt)\n"; color = :green)

        stats = fit_model(model, scene, conditions, obsBorderPts[(sim_time_frame*(titer-1)+1):(sim_time_frame*titer)+1], θ)

        est_ηpList[(sim_time_frame*(titer-1)+1):(sim_time_frame*titer)] .= stats["η"]
        est_βpList[(sim_time_frame*(titer-1)+1):(sim_time_frame*titer)] .= stats["β"]
        titer = titer + 1

        θ[1] = stats["η"]
        θ[2] = stats["β"]

        update_model!(model)
    end
    set_file(string(filepath,"/Results/plots"))

    plt_η = set_plot(22, "Time (s)", "η")
    t_windows = collect(range(start=t_steps_gt, stop=sim_time_gt, step=t_steps_gt))
    Plots.plot!(model_gt.η, label="Ground truth η(t)", dpi=400)
    # Plots.plot!(t_windows, est_ηpList, label="Estimated η(t)")
    Plots.xlabel!("time (s)")
    Plots.ylabel!("η")
    Plots.savefig(string(filepath,"/Results/plots/η.pdf"))

    # run the simulation with the estimated parameters
    viscosity_type = "bulk_viscosity"
    est_model, est_scene = def_problem(r, h, ne_exp, η_gt, ndim, FunctionClass_u, nDof_u, FunctionClass_p, nDof_p, FunctionClass_x, β_gt, F, control, viscosity_type, 
                        sim_time_gt, t_steps_gt)
    est_model.η = est_ηpList
    est_μ_list, gradList, simBorderPts, splinex, spliney, pos2D = simulate(est_model, est_scene, conditions)

    gt_μ_list, gradList, simBorderPts, splinex, spliney, pos2D = simulate(model_gt, scene_gt, conditions)

    est_h_list = get_height(est_μ_list, h)
    gt_h_list = get_height(gt_μ_list, h)

    plt_h = set_plot(22, "Time (s)", "Height")
    Plots.plot!(gt_h_list, label="Ground truth height", dpi=400)
    Plots.plot!(est_h_list, label="Estimated height")
    Plots.savefig(string(filepath,"/Results/plots/h.pdf"))

    plt_error = set_plot(22, "Time (s)", "Error")
    Plots.plot!(abs.(est_h-gt_h), label="Height estimation error", dpi=400)
    Plots.savefig(string(filepath,"/Results/plots/h_est_error.pdf"))

    write_csv(string(filepath,"/Results/data/est_η"), est_ηpList)
    write_csv(string(filepath,"/Results/data/est_β"), est_βpList)
    write_csv(string(filepath,"/Results/data/η_gt"), model_gt.η)
    write_csv(string(filepath,"/Results/data/β_gt"), β_gt)
    write_csv(string(filepath,"/Results/data/est_h"), est_h_list)
    write_csv(string(filepath,"/Results/data/gt_h"), gt_h_list)
end

function test_opt_const(exp_params::Dict)
    WRITE_GT = exp_params["WRITE_GT"] 
    r::Float64 = 1  # radius of the cylinder in mm
    h::Float64 = 1  # height of the cylinder in
    filepath::String = ""
    β_gt::Float64 = 1.0
    η_gt::Float64 = 1.0
    sim_time_gt::Float64 = 1.0# simulation time in seconds
    steps_gt::Float64 = 1.0 # number of time steps
    t_steps_gt::Float64 = 1.0
    F::Vector{Float64} = ones(Float64, round(Int, (sim_time_gt/t_steps_gt))) # force applied to the cylinder in N

    ndim::Int = 3
    FunctionClass_x::String = ""
    FunctionClass_u::String = exp_params["FunctionClass_u"]
    nDof_u::Int = ndim  # number of degree of freedom per node
    FunctionClass_p::String = exp_params["FunctionClass_p"]
    nDof_p::Int = 1  # number of degree of freedom per node

    noiseLevel::Float64 = 0
    SIDES::Bool = false
    
    sim = true # true or false with simulated dummy data or real data
    if sim == true # with simulated dummy data
        # simulation parameters for the ground truth
        scale = 50
        r = 0.5*scale  # radius of the cylinder in mm
        h = 1*scale  # height of the cylinder in mm

        FunctionClass_x = exp_params["FunctionClass_x_gt"]
        control = exp_params["control"]
        viscosity_type = exp_params["viscosity_type"] # "constant" or "bulk_viscosity"
        filepath = exp_params["filepath"]
        filepath_gt = exp_params["filepath_gt"]

        sim_time_gt = 10.0 # simulation time in seconds
        steps_gt = 25.0 # number of time steps
        t_steps_gt = sim_time_gt/steps_gt

        β_gt = exp_params["β_gt"]
        η_gt = exp_params["η_gt"]
        ne_gt::Int = exp_params["ne_gt"] # number of elements in the mesh for the ground truth

        F_ext::Float64 = 1200.0*3*β_gt # force applied to the cylinder in N

        if viscosity_type == "bulk_viscosity"
            F_ext = 1500.0*3*β_gt # force applied to the cylinder in N
            sim_time_gt = 40.0 # simulation time in seconds
            steps_gt = 100.0 # number of time steps
            t_steps_gt = sim_time_gt/steps_gt
        end

        F = -F_ext*ones(Float64, round(Int, (sim_time_gt/t_steps_gt))) # force applied to the cylinder in N

        camera_matrix = [[8*2048/7.07, 0.0, 2048/2] [0.0, 8*1536/5.3, 1536/2] [0.0, 0.0, 1.0]]'
        camera_pose = scale*[0 -0.5 2.75]'   # camera position in mm

        # Write the ground truth
        printstyled("Ground truth η: $(η_gt), ground truth β: $(β_gt)\n"; color = :green)
        model_gt, scene_gt = def_problem(r, h, ne_gt, η_gt, ndim, FunctionClass_u, nDof_u, FunctionClass_p, nDof_p, FunctionClass_x, β_gt, F, control, viscosity_type, 
                        sim_time_gt, t_steps_gt)

        if WRITE_GT == true
            @info "Writing ground truth gt data to with $ne_gt elements to $filepath"
            write_sim_data(model_gt, scene_gt, camera_matrix, camera_pose, filepath_gt)
        end

        ObsDataList, splinexObs, splineyObs = read_csv(string(filepath_gt,"/data/sim_data/contour_data"))  
        params = read_json(string(filepath_gt,"/data/sim_params.json"))
    else
        filepath_gt = string("/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/sim_experiments/Synthtic_data/exp_1")
        params = read_json(string(filepath_gt,"/data/sim_params.json"))

        ObsDataList, splinexObs, splineyObs = read_csv(string(filepath,"/data/img_data/contour_data"))  
    end

    r = params["r"]
    # h = params["h"]
    h = 50.0
    η_gt = params["η"][1]
    β_gt = params["β"][1]
    camera_matrix = reshape(Array(params["camera_matrix"]), 3, 3)
    camera_pose = reshape(Array(params["camera_pose"])*1.0,3,1)
    control = params["control_type"]
    gt_viscosity_type = params["viscosity_type"]
    F = params["cParam"]

    sim_time_gt = params["simulation_time"]
    t_steps_gt = params["time_steps"]

    println(η_gt)
    printstyled("Ground truth η: $(η_gt), ground truth β: $(β_gt)\n"; color = :green)
    
    FunctionClass_x = exp_params["FunctionClass_x"]
    ne_exp::Int = exp_params["ne_exp"] # number of elements in the mesh for the experiment
    exp_path = string("$filepath/runs/$FunctionClass_x","_","$ne_exp")
    
    model_gt, scene_gt = def_problem(r, h, ne_exp, η_gt, ndim, FunctionClass_u, nDof_u, FunctionClass_p, nDof_p, FunctionClass_x, β_gt, F, control, viscosity_type, 
                    sim_time_gt, t_steps_gt)
    write_json(string(exp_path,"/Results/data/sim_params"), exp_params)

    # simulation parameters for the experiments
    sim_time::Float64 = 10.0 # simulation time in seconds 
    steps::Float64 = 25.0 # number of time steps
    t_steps::Float64 = sim_time/steps

    if t_steps > t_steps_gt 
        @info "time resolution of the ground truth is not enough switching to ground truth resolution"
        t_steps = t_steps_gt
    end

    # Read the gt data
    obsBorderPts, nSplinex, nSpliney, pd = add_noise(ObsDataList, nFactor=noiseLevel)
    
    dev::Float64 = 0.3

    dev_η::Float64 = dev*η_gt
    ηStart::Float64 = abs(η_gt - dev_η)

    dev_β::Float64 = dev*β_gt
    βStart::Float64 = abs(β_gt - dev_β)

    θ::Vector{Float64} = [ηStart, βStart]

    if gt_viscosity_type == "constant"

        # obsBorderPts = obsBorderPts[1:(round(Int,sim_time/t_steps)+1)] # align the observation points with the simulation time

        model, scene = def_problem(r, h, ne_exp, η_gt, ndim, FunctionClass_u, nDof_u, FunctionClass_p, nDof_p, FunctionClass_x, β_gt, F, control, gt_viscosity_type, 
                        sim_time, t_steps)
        conditions = Conditions(camera_matrix=camera_matrix, camera_pose=camera_pose, SIDES=SIDES, filepath=exp_path, ANIMATE=false)

        stats = fit_model(model, scene, conditions, obsBorderPts, θ)

        time = collect(Float64, range(start=0, stop=sim_time, step=t_steps))

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
        gt_h = gt_h_[1:(round(Int,sim_time/t_steps)+1)]

        set_file(string(exp_path,"/Results/plots"))
        set_plot(22)
        Plots.plot!(time, est_h, label="Estimated height", dpi=400)
        Plots.plot!(time, gt_h, label="Ground truth height", dpi=400)
        Plots.xticks!(0:1:round(Int,sim_time))
        Plots.savefig(string(exp_path,"/Results/plots/h_est.pdf"))

        set_plot(22)
        Plots.plot!(time, abs.(est_h-gt_h), label="Height estimation error", dpi=400)
        Plots.xlabel!("Time (s)")
        Plots.ylabel!("Error")
        Plots.xticks!(0:1:round(Int,sim_time))
        Plots.savefig(string(exp_path,"/Results/plots/h_est_error.pdf"))

        # Plot the cost function with iterations
        set_plot(22)
        Plots.plot!(iterList, costList, label="Cost", marker=2, dpi=400, yscale=:log10, xminorgrid = :false)
        Plots.xlabel!("Iterations")
        Plots.ylabel!("Cost")
        Plots.xticks!(minimum(iterList):2:maximum(iterList))
        Plots.savefig(string(exp_path,"/Results/plots/cost_steps.pdf"))

        # Plot the cost function with iterations
        set_plot(22)
        Plots.plot!(iterList, costList, label="Cost", marker=2, dpi=400, yscale=:log10, xscale=:log10)
        Plots.xlabel!("Iterations")
        Plots.ylabel!("Cost")
        Plots.savefig(string(exp_path,"/Results/plots/cost_steps_log.pdf"))

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

        sampleNo = 4
        ηList = collect(range(ηStart, stop=ηStop, length=sampleNo))
        βList = collect(range(βStart, stop=βStop, length=sampleNo))
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
                μ_list, gradList, simBorderPts, splinex, spliney, pos2D = simulate(model, scene, conditions)

                # test the closest point function
                d_cp, pairs = closest_point(simBorderPts, obsBorderPts) 
                
                CostMat[i,j] = sum(d_cp)/length(d_cp)
            end
        end
        
        # Plot the cost function surface
        set_plot(22)
        Plots.contour!(ηList, βList, CostMat, color=:turbo, fill=false, levels=100, xlabel="η", ylabel="β", dpi=400)
        Plots.plot!(ηpList, βpList, label="Estimations", ms=:4, m=:x, color=:royalblue)
        Plots.plot!([η_gt], [β_gt], label="Ground truth", ms=:8, m=:star5, color=:indianred2)
        Plots.xlabel!(L"\eta")
        Plots.ylabel!(L"\beta")
        Plots.savefig(string(exp_path,"/Results/plots/cost_surface_iter.pdf"))

        set_plot(22)
        Plots.contourf!(ηList, βList, CostMat, color=:turbo, fill=false, levels=100, xlabel="η", ylabel="β", dpi=400)
        Plots.plot!([η_gt], [β_gt], label="Ground truth", ms=:8, m=:star5, color=:indianred2)
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

        # simulation parameters for the experiments
        sim_time = 10.0 # simulation time in seconds 
        steps = 25.0 # number of time steps
        t_steps = sim_time/steps

        if t_steps > t_steps_gt 
            @info "time resolution of the ground truth is not enough switching to ground truth resolution"
            t_steps = t_steps_gt
        end

        conditions = Conditions(camera_matrix=camera_matrix, camera_pose=camera_pose)
        model, scene = def_problem(r, h, ne_exp, η_gt, ndim, FunctionClass_u, nDof_u, FunctionClass_p, nDof_p, FunctionClass_x, β_gt, F, control, viscosity_type, 
                            sim_time, t_steps)

        gt_time_frame::Int = round(Int,sim_time_gt/t_steps_gt)
        sim_time_frame::Int = round(Int,sim_time/t_steps)
        windows::Int = gt_time_frame/sim_time_frame

        est_ηpList = Vector{Float64}(undef,gt_time_frame)
        est_βpList = Vector{Float64}(undef,gt_time_frame)

        println("Number of time windows: $windows")

        titer::Int = 1
        for ti::Int in 1:windows

            range = (round(Int,sim_time_frame*(titer-1))+1):(round(Int,sim_time_frame*(titer))+1)
            range_ = (round(Int,sim_time_frame*(titer-1))+1):(round(Int,sim_time_frame*(titer)))

            obsBorderPts_t = obsBorderPts[range] # align the observation points with the simulation time

            println("Time frame : $range")
            printstyled("Time window: $(ti)\n"; color = :green)

            η_gt_ = model_gt.η[range_]
            av_η = mean(η_gt_)

            printstyled("Average ground truth η in the window: $(av_η), ground truth β: $(β_gt)\n"; color = :green)

            stats = fit_model(model, scene, conditions, obsBorderPts_t, θ)

            est_ηpList[range_] .= stats["η"]
            est_βpList[range_] .= stats["β"]
            titer = titer + 1

            θ[1] = stats["η"]
            θ[2] = stats["β"]

            update_model!(model)
        end
        set_file(string(exp_path,"/Results/plots"))

        plt_η = set_plot(22)
        t_windows = collect(range(start=t_steps_gt, stop=sim_time_gt, step=t_steps_gt))
        Plots.plot!(plt_η, t_windows, model_gt.η, label="Ground truth η(t)", dpi=400)
        Plots.plot!(plt_η, t_windows, est_ηpList, label="Estimated η(t)")
        Plots.xlabel!("Time (s)")
        Plots.ylabel!("η")
        Plots.savefig(plt_η, string(exp_path,"/Results/plots/η.pdf"))

        # run the simulation with the estimated parameters
        viscosity_type = "bulk_viscosity"
        est_model, est_scene = def_problem(r, h, ne_exp, η_gt, ndim, FunctionClass_u, nDof_u, FunctionClass_p, nDof_p, FunctionClass_x, β_gt, F, control, viscosity_type, 
                            sim_time_gt, t_steps_gt)
        est_model.η = est_ηpList
        est_μ_list, gradList, simBorderPts, splinex, spliney, pos2D = simulate(est_model, est_scene, conditions)

        gt_μ_list, gradList, simBorderPts, splinex, spliney, pos2D = simulate(model_gt, scene_gt, conditions)

        est_h_list = get_height(est_μ_list, h)
        gt_h_list = get_height(gt_μ_list, h)

        plt_h = set_plot(22)
        Plots.plot!(gt_h_list, label="Ground truth height", dpi=400)
        Plots.plot!(est_h_list, label="Estimated height")
        Plots.xlabel!("Time (s)")
        Plots.ylabel!("Height")
        Plots.savefig(string(exp_path,"/Results/plots/h.pdf"))

        plt_error = set_plot(22)
        Plots.plot!(abs.(est_h_list-gt_h_list), label="Height estimation error", dpi=400)
        Plots.savefig(string(exp_path,"/Results/plots/h_est_error.pdf"))
        Plots.xlabel!("Time (s)")
        Plots.ylabel!("Error")
        Plots.savefig(string(exp_path,"/Results/plots/error.pdf"))

        write_csv(string(exp_path,"/Results/data/est_η"), est_ηpList)
        write_csv(string(exp_path,"/Results/data/est_β"), est_βpList)
        write_csv(string(exp_path,"/Results/data/η_gt"), model_gt.η)
        write_csv(string(exp_path,"/Results/data/β_gt"), β_gt)
        write_csv(string(exp_path,"/Results/data/est_h"), est_h_list)
        write_csv(string(exp_path,"/Results/data/gt_h"), gt_h_list)
    end
end

function plot_()

    filepath = string("/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/sim_experiments/cost_function_test/optimization/Stokes/force/test_bulk_viscosity/Results/")
    
    est_η = readdlm(string(filepath,"data/est_η.csv"), ',', Float64)
    est_β = readdlm(string(filepath,"data/est_β.csv"), ',', Float64)
    η_gt = readdlm(string(filepath,"data/η_gt.csv"), ',', Float64)
    β_gt = readdlm(string(filepath,"data/β_gt.csv"), ',', Float64)
    est_h = readdlm(string(filepath,"data/est_h.csv"), ',', Float64)
    gt_h = readdlm(string(filepath,"data/gt_h.csv"), ',', Float64)

    set_file(string(filepath,"/plots"))
    Plots.plot(gt_h, label="Ground truth height", dpi=400)
    Plots.plot!(est_h, label="Estimated height")
    Plots.xlabel!("time (s)")
    Plots.ylabel!("h(t)")
    Plots.savefig(string(filepath,"/plots/h.pdf"))

    Plots.plot(abs.(est_h-gt_h), label="Height estimation error", dpi=400)
    Plots.xlabel!("Time (s)")
    Plots.ylabel!("Error")
    Plots.savefig(string(filepath,"/plots/h_est_error.pdf"))

    Plots.plot(est_η, label="Estimated η", dpi=400)
    Plots.plot!(η_gt, label="Ground truth η")
    Plots.xlabel!("Time (s)")
    Plots.ylabel!("η")
    Plots.savefig(string(filepath,"/plots/η.pdf"))

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
    Plots.xlabel!(fig1,"Iterations")
    Plots.ylabel!(fig1,"η")
    # Plots.ylims!(fig1, η_gt[1]*0.5, η_gt[1]*1.5)


    fig2 = Plots.plot(β_gt_list, label="Ground truth β", dpi=400)
    Plots.xlabel!(fig2, "Iterations")
    Plots.ylabel!(fig2, "β")
    # Plots.ylims!(fig2, abs(β_gt[1]-30), β_gt[1]*1.05)

    for run_folder_ in run_filepath
        println("Processing folder: ", run_folder_)
        run_folder = string(filepath,"/runs/",run_folder_)
        params = read_json(string(run_folder,"/Results/data/sim_params.json"))

        FunctionClass_x = params["FunctionClass_x"]
        ne = params["ne_exp"]

        if ne != 2
            η = readdlm(string(run_folder,"/Results/data/η.csv"), ',', Float64)
            β = readdlm(string(run_folder,"/Results/data/β.csv"), ',', Float64)

            Plots.plot!(fig1, η, label=string("Basis - ",FunctionClass_x," - ne: ",ne))
            Plots.plot!(fig2, β, label=string("Basis - ",FunctionClass_x," - ne: ",ne))
        end
    end

    plot_path = string(filepath,"/Results/plots")
    set_file(plot_path)
    Plots.savefig(fig1, string(plot_path,"/η_comp_zoomed.pdf"))
    Plots.savefig(fig2, string(plot_path,"/β_comp_zoomed.pdf"))

end

function main()
    ne_gt::Int = 4 # number of elements in the mesh for the ground truth
    ne_exp::Int = 2 # number of elements in the mesh for the experiment 
    # β_gt_list = [5, 10, 50, 100.0, 200.0, 500.0, 1000.0, 10000.0]
    # η_gt_list = [40.0]
    β_gt_list = [10.0, 50.0, 100.0, 1e3]
    # β_gt_list = [100.0]
    η_gt_list = [60.0]
    FunctionClass_x_List = ["Q2"]
    # refine_list = [1, 2, 3] # refinement levels, ne = ne_exp^refine
    refine_list = [2] # refinement levels, ne = ne_exp^refine
    control = "force" # "force" or "velocity"
    viscosity_type_list = ["constant"]
    FunctionClass_x_gt_list = ["Q2"] # Function space for the ground truth

    exp_size = size(FunctionClass_x_List,1)*size(refine_list,1)
    η_mat = zeros(30, exp_size)
    β_mat = zeros(30, exp_size)

    WRITE_GT = true
    for viscosity_type in viscosity_type_list
        for FunctionClass_x_gt in FunctionClass_x_gt_list
            run_id = 1
            for β_gt in β_gt_list
                for η_gt in η_gt_list
    
                    filepath = string("/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/sim_experiments/cost_function_test/optimization/Stokes/$control/model_validation/$viscosity_type/$ne_gt/$run_id")
                    filepath_gt = string("/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/sim_experiments/ground_truth/fem/Stokes/$control/$viscosity_type/$ne_gt/$run_id")
                    
                    for ref in refine_list
                        ne = ne_exp^ref
                        @info "Running optimization with ne = $ne"
                        for FunctionClass_x in FunctionClass_x_List
                            @info "Running optimization with FunctionClass_x = $FunctionClass_x with $ne elements"
    
                            exp_params = Dict("FunctionClass_x" => FunctionClass_x, "FunctionClass_u" => "Q2", "FunctionClass_p" => "Q1", "ne_gt" => ne_gt, "ne_exp" => ne, 
                                        "β_gt" => β_gt, "η_gt" => η_gt, "WRITE_GT" => WRITE_GT, "filepath" => filepath, "filepath_gt"=>filepath_gt, "control" => control, 
                                        "viscosity_type"=>viscosity_type, "FunctionClass_x_gt" => FunctionClass_x_gt)
    
                            test_opt_const(exp_params)
                            WRITE_GT = true
                        end
                    end
                    # post_analysis(filepath_gt, filepath)
                    run_id = run_id + 1
                end
            end
        end
    end
end



main()
