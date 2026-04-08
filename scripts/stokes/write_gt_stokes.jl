using smearFEM
using Dates
using Plots, Plots.PlotMeasures
using LaTeXStrings
using StatsPlots
using LinearAlgebra

# Helper function to compute F_ext based on penalty parameter
function _get_F_ext(β_gt::Real)::Float64
    if β_gt <= 1.0
        return 9.813e3 * 0.85
    elseif β_gt == 10.0
        return 9.813e3 * 0.93
    elseif β_gt == 50.0
        return 9.813e3 * 0.955
    elseif β_gt == 100.0
        return 9.813e3 * 0.97
    elseif β_gt == 500.0
        return 9.813e3 * 0.99
    elseif β_gt == 1e3
        return 9.813e3 * 1.0
    elseif β_gt == 5e3
        return 9.813e3 * 1.01
    elseif β_gt == 2e3
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
end

function main(; use_parallel::Bool=true, max_workers::Int=8)
    # parameters for the optimization
    r::Float64 = 25.0  # radius of the cylinder in mm
    h::Float64 = 40.0  # height of the cylinder in mm
    ne_gt::Int = 16 # number of elements in the mesh for the ground truth

    β_gt_list = [0.01, 10.0, 50.0, 1e2, 500.0, 1e3] # penalty parameters for the ground truth [2e3, 5e3, 1e4, 1e5, 1e10]
    η_gt_list = [1e2] # viscosity values for the ground truth in kg/(mm⋅s)

    control = "force" # "force" or "velocity"

    viscosity_type_list = ["bulk_viscosity"] # "constant" or "bulk_viscosity"
    FunctionClass_x_gt_list = ["Q2"] # Function space for the ground truth
    
    sim_time_gt::Float64 = 30.0 # simulation time in seconds
    steps_gt::Int = 300 # number of time steps[]

    obj_pose = zeros(Float64, 4, 4)
    obj_pose[1, 1] = -1.0
    obj_pose[2, 3] = -1.0
    obj_pose[3, 2] = -1.0
    obj_pose[1:3, 4] = [0.0, h/2, 150.0]
    camera_matrix::Matrix{Float64} = [2.39642674e+03  0.0  1.00429248e+03; 0.0  2.40565353e+03  7.57028161e+02; 0.0  0.0 1.0;]

    # Build list of all experiment parameters
    params_list = Dict[]
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
                    push!(params_list, exp_params)
                    run_id = run_id + 1
                end
            end
        end
    end

    @info "Built parameter list with $(length(params_list)) experiment configurations"
    if !isempty(params_list)
        @info "First experiment: η=$(params_list[1]["η_gt"]), β=$(params_list[1]["β_gt"]), filepath=$(params_list[1]["filepath_gt"])"
    end

    # Run experiments in parallel or sequentially based on use_parallel flag
    if use_parallel
        run_param_list(params_list; max_workers=max_workers)
    else
        for (i, params) in enumerate(params_list)
            @info "Sequential execution: calling write_gt_data for index $i / $(length(params_list))"
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

# Helper: Display available cores and calculate batch information
function display_batch_info(n_experiments::Int, n_cores::Int=Threads.nthreads())
    n_batches = ceil(Int, n_experiments / n_cores)
    experiments_per_batch = n_experiments ÷ n_batches
    remaining = n_experiments % n_batches
    
    @info "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    @info "CPU CORES AVAILABLE: $n_cores cores"
    @info "TOTAL EXPERIMENTS: $n_experiments"
    @info "NUMBER OF BATCHES: $n_batches"
    if n_batches > 1
        @info "EXPERIMENTS PER BATCH: $experiments_per_batch (batch 1-$n_batches each have $(experiments_per_batch + (remaining > 0 ? 1 : 0)) experiments max)"
        @info "BATCH DISTRIBUTION: $remaining batch(es) with extra experiment(s)"
    end
    @info "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    
    return n_batches
end

# helper: run a vector of parameter Dicts with limited concurrent workers
function run_param_list(params_list::Vector{Dict}; max_workers::Int=-1)
    isempty(params_list) && return
    
    # Auto-detect available cores if max_workers not specified (default -1)
    n_available_cores = Threads.nthreads()
    if max_workers <= 0
        workers = n_available_cores
    else
        workers = max(1, min(n_available_cores, max_workers))
    end
    
    # Display batch information
    display_batch_info(length(params_list), workers)
    
    @info "Running $(length(params_list)) experiments with $workers workers"
    @warn "Note: Gmsh is not thread-safe. All Gmsh operations (mesh loading) are serialized via a lock.\n         If mesh generation is the bottleneck, consider using sequential execution: main(use_parallel=false)"
    
    # Atomic counter for tracking completed experiments
    completed = Threads.Atomic{Int}(0)
    
    # Set BLAS threads once before spawning tasks (not inside threads)
    LinearAlgebra.BLAS.set_num_threads(max(1, div(Threads.nthreads(), workers)))
    
    ch = Channel{Int}(length(params_list))
    @sync begin
        for i in 1:length(params_list)
            put!(ch, i)
        end
        close(ch)
        
        for w in 1:workers
            Threads.@spawn begin
                while true
                    try
                        idx = take!(ch)  # Blocks until item available or channel closed
                        
                        try
                            # Suppress output from write_gt_data
                            redirect_stdout(devnull) do
                                redirect_stderr(devnull) do
                                    write_gt_data(params_list[idx])
                                end
                            end
                            # Report progress without clutter
                            Threads.atomic_add!(completed, 1)
                            @info "Experiments completed: $(completed[]) / $(length(params_list))"
                        catch err
                            _handle_worker_error(err, idx, params_list[idx])
                        end
                    catch e
                        # EOFError or InvalidStateException when taking from closed/empty channel
                        if isa(e, EOFError) || isa(e, InvalidStateException)
                            break
                        else
                            @error "Worker $w encountered unexpected error" exception=e
                            break
                        end
                    end
                end
            end
        end
    end
    println("All experiments completed.")
end

# Helper function for consistent error handling in worker threads
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

main()