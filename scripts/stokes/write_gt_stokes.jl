using smearFEM
using Dates
using Plots, Plots.PlotMeasures
using LaTeXStrings
using StatsPlots

function write_data(exp_params::Dict)

    ndim::Int = 3
    FunctionClass_u::String = exp_params["FunctionClass_u"]
    FunctionClass_p::String = exp_params["FunctionClass_p"]    
    FunctionClass_x::String = exp_params["FunctionClass_x"]
    nDof_u::Int = ndim  # number of degree of freedom per node
    nDof_p::Int = 1  # number of degree of freedom per node

    # simulation parameters for the ground truth
    r::Float64 = exp_params["r"]  # radius of the cylinder in mm
    h::Float64 = exp_params["h"]  # height of the cylinder in mm

    control::String = exp_params["control"]
    viscosity_type::String = exp_params["viscosity_type"] # "constant" or "bulk_viscosity"
    filepath_gt::String = exp_params["filepath_gt"]
    obj_pose::AbstractArray = exp_params["obj_pose_gt"]
    camera_matrix::AbstractArray = exp_params["camera_matrix"]

    sim_time_gt::Float64 = exp_params["sim_time_gt"]
    steps_gt::Float64 = exp_params["steps_gt"]
    t_steps_gt::Float64 = sim_time_gt/steps_gt

    β_gt::Float64 = exp_params["β_gt"]
    η_gt::Float64 = exp_params["η_gt"]
    ne_gt::Int = exp_params["ne_gt"] # number of elements in the mesh for the ground truth

    F_ext::Float64 = exp_params["F_ext"]

    if viscosity_type == "bulk_viscosity"
        F_ext = 1500.0*3*β_gt # force applied to the cylinder in N
        sim_time_gt = 40.0 # simulation time in seconds
        steps_gt = 100.0 # number of time steps
        t_steps_gt = sim_time_gt/steps_gt
    end

    F = -F_ext*ones(Float64, round(Int, (sim_time_gt/t_steps_gt))) # force applied to the cylinder in N

    # Write the ground truth
    printstyled("Ground truth η: $(η_gt), ground truth β: $(β_gt)\n"; color = :green)
    model_gt, scene_gt = def_problem(r, h, ne_gt, η_gt, ndim, FunctionClass_u, nDof_u, FunctionClass_p, nDof_p, FunctionClass_x, β_gt, F, control, viscosity_type, 
                    sim_time_gt, t_steps_gt)

    @info "Writing ground truth gt data to with $ne_gt elements to $filepath_gt"
    write_sim_data(model_gt, scene_gt, camera_matrix, obj_pose, filepath_gt)

  end

  function main()
    # parameters for the optimization
    r::Float64 = 25.0  # radius of the cylinder in mm
    h::Float64 = 40.0  # height of the cylinder in mm
    ne_gt::Int = 16 # number of elements in the mesh for the ground truth

    β_gt_list = [10.0, 50.0, 100.0, 1e3]
    η_gt_list = [100.0]

    refine_list = [2] # refinement levels, ne = ne_exp^refine
    control = "force" # "force" or "velocity"

    viscosity_type_list = ["constant"] # "constant" or "bulk_viscosity"
    FunctionClass_x_gt_list = ["Q2"] # Function space for the ground truth

    F_ext_::Float64 = 0.15*9.812*1e3 # force applied to the cylinder in N
    sim_time_gt::Float64 = 20.0 # simulation time in seconds
    steps_gt::Float64 = 200.0 # number of time steps

    obj_pose = zeros(Float64, 4,4)
    obj_pose[1,1] = -1.0
    obj_pose[2,3] = -1.0
    obj_pose[3,2] = -1.0
    obj_pose[1:3,4] = [0.0, h/2, 150.0]
    camera_matrix::AbstractMatrix{Float64} = [[2.39642674e+03, 0.0, 1.00429248e+03] [0.0, 2.40565353e+03, 7.57028161e+02] [0.0, 0.0, 1.0]]'

    for viscosity_type in viscosity_type_list
        for FunctionClass_x in FunctionClass_x_gt_list
            run_id = 1
            for β_gt in β_gt_list
                for η_gt in η_gt_list
                    F_ext  = F_ext_*0.5*β_gt
                    filepath_gt = string("/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/sim_experiments/ground_truth/fem/Stokes/$control/$viscosity_type/$(FunctionClass_x)_$(ne_gt)/$run_id")

                    exp_params = Dict("FunctionClass_x" => FunctionClass_x, "FunctionClass_u" => "Q2", "FunctionClass_p" => "Q1", "ne_gt" => ne_gt, "β_gt" => β_gt, 
                                "η_gt" => η_gt, "filepath_gt"=>filepath_gt, "control" => control, "viscosity_type" => viscosity_type, "obj_pose_gt" => obj_pose, 
                                "F_ext" => F_ext, "sim_time_gt" => sim_time_gt, "steps_gt" => steps_gt, "r" => r, "h" => h, "camera_matrix" => camera_matrix)

                    write_data(exp_params)
                    run_id = run_id + 1
                end
            end
        end
    end
end

main()