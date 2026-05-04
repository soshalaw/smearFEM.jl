using LinearAlgebra
using ProgressMeter
using SparseArrays
using Plots

using smearFEM
using StatsPlots
using Distributions
using Dates

function run_exp(exp_params::Dict)

    r::Float64 = 1  # radius of the cylinder in mm
    h::Float64 = 1  # height of the cylinder in
    filepath::String = ""
    β::Float64 = 1.0
    η::Float64 = 1.0
    sim_time_gt::Float64 = 1.0# simulation time in seconds
    steps_gt::Float64 = 1.0 # number of time steps
    t_steps_gt::Float64 = 1.0
    F::Vector{Float64} = ones(Float64, round(Int, (sim_time_gt/t_steps_gt))) # force applied to the cylinder in N

    ndim::Int = 3
    FunctionClass_u::String = "Q2"
    nDof_u::Int = ndim  # number of degree of freedom per node
    FunctionClass_p::String = "Q1"
    nDof_p::Int = 1  # number of degree of freedom per node

    noiseLevel::Float64 = 0
    SIDES::Bool = false
    scale = 50
    r = 0.5*scale  # radius of the cylinder in mm
    h = 1*scale  # height of the cylinder in mm

    control = exp_params["control"]
    viscosity_type = "constant" # "constant" or "bulk_viscosity"
    filepath = exp_params["filepath"]

    sim_time = 20.0 # simulation time in seconds
    steps = 20.0 # number of time steps
    t_steps = sim_time/steps

    time = collect(Float64, range(start=0, stop=sim_time, step=t_steps))

    β = exp_params["β_gt"]
    η = exp_params["η_gt"]
    ne::Int = exp_params["ne_gt"] # number of elements in the mesh for the ground truth
    WRITE_DATA = exp_params["WRITE_DATA"]

    if β == 1e-5
        F_ext = 1500
    else
        F_ext::Float64 = 1500.0*β # force applied to the cylinder in N
    end
    F = -F_ext*ones(Float64, round(Int, (sim_time/t_steps))) # force applied to the cylinder in N

    camera_matrix = [[8*2048/7.07, 0.0, 2048/2] [0.0, 8*1536/5.3, 1536/2] [0.0, 0.0, 1.0]]'
    camera_pose = scale*[0 -0.5 2.75]'   # camera position in mm

    # FEM data
    printstyled("η: $(η), β: $(β)\n"; color = :green)
    file_path_fem = string(filepath,"/fem")
    model_fem, scene_fem = def_problem(r, h, ne, η, ndim, FunctionClass_u, nDof_u, FunctionClass_p, nDof_p, "Q2", β, F, control, viscosity_type, 
                    sim_time, t_steps)
    if WRITE_DATA
        write_sim_data(model_fem, scene_fem, camera_matrix, camera_pose, file_path_fem)
    end

    ObsDataList_fem, splinexObs, splineyObs = read_csv(string(file_path_fem,"/data/sim_data/contour_data"))  
    obsBorderPts_fem, nSplinex_fem, nSpliney_fem, pd = add_noise(ObsDataList_fem, nFactor=noiseLevel)
    
    # IGA data
    file_path_iga = string(filepath,"/iga")
    model_iga, scene_iga = def_problem(r, h, ne, η, ndim, FunctionClass_u, nDof_u, FunctionClass_p, nDof_p, "S2", β, F, control, viscosity_type, 
                sim_time, t_steps)
    if WRITE_DATA
        write_sim_data(model_iga, scene_iga, camera_matrix, camera_pose, file_path_iga)
    end

    ObsDataList_iga, splinexObs, splineyObs = read_csv(string(file_path_iga,"/data/sim_data/contour_data"))  
    obsBorderPts_iga, nSplinex_iga, nSpliney_iga, pd = add_noise(ObsDataList_iga, nFactor=noiseLevel)

    animate_fields(filepath = string(filepath,"/Results/images"), p=nSplinex_fem, q=nSpliney_fem, pObs=nSplinex_iga, qObs=nSpliney_iga)
    
    d, pairs = closest_point(obsBorderPts_fem, obsBorderPts_iga)
    return d, time
end

function main()
    WRITE_DATA = true

    β_list = [1e-5 10 100 500 1000 1e5]
    exp_ne = [2, 4, 8]  # number of elements in radial and height directions,

    control = "force" # "force" or "velocity"
    η = 40
    error_list = zeros(length(β_list),length(exp_ne),21)
    error_sum_list = zeros(length(β_list),length(exp_ne))
    time = zeros(20,1)
    run_id = 1
    file_path = resolve_data_path("sim_experiments/convergence_analysis/stokes_convergence")
    
    for (i,β_gt) in enumerate(β_list)
        for (j,ne) in enumerate(exp_ne)
            filepath = string(file_path,"/",run_id)
            @info "Running optimization with ne = $ne"
            exp_params = Dict("ne_gt" => ne, "β_gt" => β_gt, "η_gt" => η, "filepath" => filepath, "control" => control, "WRITE_DATA"=> WRITE_DATA)
            d, time = run_exp(exp_params)
            error_list[i,j,:] = d
            error_sum_list[i,j] = sum(d)/length(d)
            run_id += 1
        end
    end

    plot_path = string(file_path,"/plot")
    set_file(plot_path)
    # set_plot(22)
    # Plots.plot!(time,error_list[1,1,:])
    # Plots.savefig(string(plot_path,"/contour_convergence.pdf"))
    set_plot(22)

    for (i,β) in enumerate(β_list)
        Plots.plot!(exp_ne, error_sum_list[i,:],label = "β = $β", lw=2)
    end
    Plots.xlabel!("Number of elements")
    Plots.ylabel!("Mean error (px/s)")
    Plots.savefig(string(plot_path,"/contour_convergence.pdf"))

    for (i,β) in enumerate(β_list)
        set_plot(22)
        array = error_list[i,:,:]
        iter = 1:size(array,1)
        for j in iter
            ne = exp_ne[j]
            Plots.plot!(time,array[j,:],label="ne = $ne", lw=2)
        end
        Plots.title!("β = $β")
        Plots.ylabel!("Contour Error (px)")
        Plots.xlabel!("Time (s)")
        Plots.savefig(string(plot_path,"/contour_convergence_with_ne_$i.pdf"))
    end

    for (i,ne) in enumerate(exp_ne)
        set_plot(22)
        array = error_list[:,i,:]
        iter = 1:size(array,1)
        for j in iter
            β = β_list[j]
            Plots.plot!(time,array[j,:], label="β = $β", lw=2)
        end
        Plots.title!("ne = $ne")
        Plots.ylabel!("Contour Error (px)")
        Plots.xlabel!("Time (s)")
        Plots.savefig(string(plot_path,"/contour_convergence_with_slip_$i.pdf"))
    end

    # for FunctionClass_x_gt in FunctionClass_x_gt_list
    #     run_id = 1
    #     file_path = resolve_data_path("sim_experiments/convergence_analysis/stokes_convergence")
    #     @info "Running optimization with FunctionClass_x_gt = $FunctionClass_x_gt with $ne_gt elements for the ground truth ..."
    #     for β_gt in β_gt_list
    #         for η_gt in η_gt_list
    #             filepath = string(file_path,"/",run_id)
    #             for ref in refine_list
    #                 ne = ne_exp^ref
    #                 @info "Running optimization with ne = $ne"
    #                 for FunctionClass_x in FunctionClass_x_List
    #                     @info "Running optimization with FunctionClass_x = $FunctionClass_x"
    #                     exp_params = Dict("FunctionClass_x" => FunctionClass_x, "FunctionClass_u" => "Q2", "FunctionClass_p" => "Q1", "ne_gt" => ne_gt, "ne_exp" => ne, 
    #                                 "β_gt" => β_gt, "η_gt" => η_gt, "filepath" => filepath, "control" => control, "FunctionClass_x_gt" => FunctionClass_x_gt)
    #                     d = run_exp(exp_params)
    #                     WRITE_GT = false
    #                 end
    #             end
    #             run_id += 1
    #         end
    #     end
    # end
end

main()