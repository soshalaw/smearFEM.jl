using smearFEM
using Dates
using Plots, Plots.PlotMeasures
using LaTeXStrings
using StatsPlots
using LinearAlgebra

include("../ParallelExecution.jl")
using .ParallelExecution

function _get_F_ext(β_gt::Real)::Float64
    if β_gt <= 1.0
        return 9.813e3 * 0.85
    elseif β_gt <= 10.0
        return 9.813e3 * 0.93
    elseif β_gt <= 50.0
        return 9.813e3 * 0.955
    elseif β_gt <= 100.0
        return 9.813e3 * 0.97
    elseif β_gt <= 500.0
        return 9.813e3 * 0.99
    elseif β_gt <= 1e3
        return 9.813e3 * 1.0
    elseif β_gt <= 5e3
        return 9.813e3 * 1.01
    elseif β_gt <= 2e3
        return 9.813e3 * 0.995
    elseif β_gt == 1e4
        return 9.813e3 * 0.85 * 700
    elseif β_gt == 1e5
        return 9.813e3 * 0.85 * 7e3
    elseif β_gt == 1e10
        return 9.813e3 * 0.85 * 7e8
    else
        return 9.813e3 * 0.85
    end
    # if β_gt == 0.01
    #     return 9.813e3 * 0.638834  # F_ext=6268.88 N
    # elseif β_gt == 0.1
    #     return 9.813e3 * 0.637642  # F_ext=6257.18 N
    # elseif β_gt == 1.0
    #     return 9.813e3 * 0.62686  # F_ext=6151.38 N
    # elseif β_gt == 10.0
    #     return 9.813e3 * 0.687133  # F_ext=6742.83 N
    # elseif β_gt == 50.0
    #     return 9.813e3 * 0.664343  # F_ext=6519.19 N
    # elseif β_gt == 100.0
    #     return 9.813e3 * 0.671039  # F_ext=6584.91 N
    # elseif β_gt == 500.0
    #     return 9.813e3 * 0.685549  # F_ext=6727.29 N
    # elseif β_gt == 1000.0
    #     return 9.813e3 * 0.698172  # F_ext=6851.16 N
    # end
end

# Read final height from simulation output CSV
function get_final_height(filepath::String)::Float64
    height_file = joinpath(filepath, "data", "h.csv")
    
    if !isfile(height_file)
        @warn "Height file not found: $height_file"
        return NaN
    end
    
    try
        lines = readlines(height_file)
        if isempty(lines)
            return NaN
        end
        
        final_height = parse(Float64, strip(lines[end]))
        return final_height
    catch err
        @error "Failed to read final height from $height_file" exception=err
        return NaN
    end
end

# Iteratively adjust F_ext multipliers to achieve target height across all experiments
function calibrate_multipliers(param_list::Vector{Dict}; target_height::Float64=29.0, tolerance::Float64=1e-2, max_iterations::Int=3)
    base_force = 9.813e3
    multipliers = Dict(
        0.01 => 0.85,
        0.1 => 0.85,
        1.0 => 0.85,
        10.0 => 0.93,
        50.0 => 0.955,
        100.0 => 0.97,
        500.0 => 0.99,
        1e3 => 1.0
    )
    
    @info "Starting multiplier calibration (target height: $(target_height)mm ± $(tolerance)mm)"
    
    for iteration in 1:max_iterations
        @info "\n" * ("="^80)
        @info "CALIBRATION ITERATION $iteration / $max_iterations"
        @info ("="^80)
        
        for i in 1:length(param_list)
            β = param_list[i]["β_gt"]
            mult = get(multipliers, β, 0.85)
            param_list[i]["F_ext"] = base_force * mult
        end
        
        @info "Running experiments with current multipliers..."
        run_parallel_tasks(param_list, write_gt_data; max_workers=-1, memory_per_task_mb=512.0)
        
        @info "\nHeight measurements:"
        all_converged = true
        for i in 1:length(param_list)
            β = param_list[i]["β_gt"]
            h_final = get_final_height(param_list[i]["filepath_gt"])
            mult = multipliers[β]
            
            if isnan(h_final)
                @warn "  β=$β: Could not read height"
                all_converged = false
                continue
            end
            
            height_error = h_final - target_height
            @info "  β=$β: height=$(round(h_final, digits=4))mm (target=$target_height, error=$(round(height_error, digits=4))mm), multiplier=$mult"
            
            if abs(height_error) > tolerance
                all_converged = false
                adjustment_factor = target_height / h_final
                new_mult = mult * adjustment_factor * 0.9  # damping factor prevents oscillation
                multipliers[β] = new_mult
                @info "    → Adjusting multiplier to $(round(new_mult, digits=6))"
            else
                @info "    ✓ Converged"
            end
        end
        
        all_converged && break
    end
    
    return multipliers
end

# Get available system memory in MB (Linux: /proc/meminfo)


# Write error logs to file
function write_error_log(err::Exception, bt::Vector; params::Dict=Dict(), dest_dir::String=".")
    mkpath(dest_dir)
    log_file = joinpath(dest_dir, "error_log_$(now()).txt")
    open(log_file, "w") do f
        write(f, "Error Log\n")
        write(f, "=========\n\n")
        write(f, "Timestamp: $(now())\n")
        write(f, "Error: $(err)\n\n")
        write(f, "Parameters: $params\n\n")
        write(f, "Stacktrace:\n")
        write(f, string(bt))
    end
    @info "Error log written to $log_file"
end

# Handle errors in worker threads
function _handle_worker_error(err::Exception, idx::Int, params::Dict)
    bt = catch_backtrace()
    @error "write_gt_data failed for params index $idx" exception=(err, bt)
    
    try
        dest_dir = get(params, "filepath_res", ".") |> d -> joinpath(d, "Results", "logs")
        write_error_log(err, bt; params=params, dest_dir=dest_dir)
    catch ewrite
        @error "Failed to write error log" exception=ewrite
    end
end

function main(; use_parallel::Bool=true, calibrate::Bool=false, max_workers::Int=-1, memory_per_experiment_mb::Float64=512.0)
    # parameters for the optimization
    r::Float64 = 25.0  # radius of the cylinder in mm
    h::Float64 = 40.0  # height of the cylinder in mm
    ne_gt::Int = 16 # number of elements in the mesh for the ground truth

    β_gt_list = [0.01, 0.1, 1, 2, 4, 6, 8, 10, 20, 30, 40, 50, 100, 200, 300, 400,500, 1e3, 2e3, 5e3] # penalty parameters for the ground truth [2e3, 5e3, 1e4, 1e5, 1e10]
    η_gt_list = [1e2] # viscosity values for the ground truth in kg/(mm⋅s)

    control = "force" # "force" or "velocity"

    viscosity_type_list = ["constant"] # "constant" or "bulk_viscosity"
    FunctionClass_x_gt_list = ["Q2"] # Function space for the ground truth
    
    sim_time_gt::Float64 = 5.0 # simulation time in seconds
    steps_gt::Int =  50 # number of time steps[]

    obj_pose = zeros(Float64, 4, 4)
    obj_pose[1, 1] = -1.0
    obj_pose[2, 3] = -1.0
    obj_pose[3, 2] = -1.0
    obj_pose[1:3, 4] = [0.0, h/2, 150.0]
    camera_matrix::Matrix{Float64} = [2.39642674e+03  0.0  1.00429248e+03; 0.0  2.40565353e+03  7.57028161e+02; 0.0  0.0 1.0;]

    # Build list of all experiment parameters
    param_list = Dict[]
    for viscosity_type in viscosity_type_list
        for FunctionClass_x in FunctionClass_x_gt_list
            run_id = 1
            for β_gt in β_gt_list
                for η_gt in η_gt_list
                    F_ext = _get_F_ext(β_gt)
                    filepath_gt = string("/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/ground_truth/sim_data/Stokes/$control/$viscosity_type/$(FunctionClass_x)_$(ne_gt)/$run_id")

                    exp_params = Dict(
                        "FunctionClass_x" => FunctionClass_x,
                        "FunctionClass_u" => "Q2",
                        "FunctionClass_p" => "Q1",
                        "ne_gt" => ne_gt,
                        "β_gt" => β_gt,
                        "η_gt" => η_gt,
                        "filepath_gt" => filepath_gt,
                        "control" => control,
                        "viscosity_type" => viscosity_type,
                        "obj_pose_gt" => obj_pose,
                        "F_ext" => F_ext,
                        "sim_time_gt" => sim_time_gt,
                        "steps_gt" => steps_gt,
                        "r" => r,
                        "h" => h,
                        "camera_matrix" => camera_matrix,
                        "animate" => !use_parallel
                    )
                    push!(param_list, exp_params)
                    run_id = run_id + 1
                end
            end
        end
    end

    @info "Built parameter list with $(length(param_list)) experiment configurations"
    if !isempty(param_list)
        @info "First experiment: η=$(param_list[1]["η_gt"]), β=$(param_list[1]["β_gt"]), filepath=$(param_list[1]["filepath_gt"])"
    end

    # Step 1: Optionally calibrate forces to achieve ~29mm target height
    if calibrate
        @info "Starting force calibration to achieve target height of 29mm..."
        calibrated_multipliers = calibrate_multipliers(param_list; target_height=29.0, tolerance=1e-2, max_iterations=3)
        
        @info "\n" * ("="^80)
        @info "FINAL CALIBRATED MULTIPLIERS - Copy to _get_F_ext():"
        @info ("="^80)
        
        # Display in _get_F_ext format
        for β in sort(collect(keys(calibrated_multipliers)))
            mult = calibrated_multipliers[β]
            force = 9.813e3 * mult
            @info "  β=$β: multiplier=$(round(mult, digits=6)) → F_ext=$(round(force, digits=2)) N"
        end
        
        # Save to file
        open("calibrated_multipliers.txt", "w") do f
            write(f, "Calibrated multipliers for _get_F_ext function:\n")
            write(f, "Copy these into the function in write_gt_stokes.jl\n")
            write(f, ("="^70) * "\n\n")
            for β in sort(collect(keys(calibrated_multipliers)))
                mult = calibrated_multipliers[β]
                force = 9.813e3 * mult
                write(f, "elseif β_gt == $β\n")
                write(f, "    return 9.813e3 * $(round(mult, digits=6))  # F_ext=$(round(force, digits=2)) N\n")
            end
        end
        @info "\nCalibrated multipliers saved to calibrated_multipliers.txt"
    end
    
    # Step 2: Run all experiments with current/calibrated F_ext values
    if use_parallel
        @info "Running $(length(param_list)) experiments in parallel..."
        # Note: Gmsh is not thread-safe, mesh loading is serialized via lock in smearFEM
        run_parallel_tasks(param_list, write_gt_data; max_workers=max_workers, memory_per_task_mb=memory_per_experiment_mb)
    else
        @info "Running $(length(param_list)) experiments sequentially..."
        for (i, params) in enumerate(param_list)
            @info "Sequential execution: calling write_gt_data for index $i / $(length(param_list))"
            try
                write_gt_data(params)
                @info "Completed write_gt_data for index $i"
            catch err
                _handle_worker_error(err, i, params)
            end
        end
        println("All experiments completed.")
    end
end

# Usage: main(use_parallel=true, calibrate=false) or main(use_parallel=true, calibrate=true)
# Default: Run with current F_ext values (to calibrate: pass calibrate=true)
main(use_parallel=true, calibrate=false, memory_per_experiment_mb=1024.0)