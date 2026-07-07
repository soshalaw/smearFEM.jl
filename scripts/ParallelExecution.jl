"""
    ParallelExecution

A reusable module for memory-aware parallel execution with thread-safe progress tracking.
Supports dynamic work-stealing queue-based execution with per-worker timing metrics.

Usage:
    using ParallelExecution
    
    # Define task function that accepts single parameter
    function my_task(params::Dict)
        # ... do work ...
    end
    
    # Run in parallel
    run_parallel_tasks(param_list, my_task; max_workers=-1, memory_per_task_mb=512.0)
"""

module ParallelExecution

using Dates
using LinearAlgebra

export get_available_memory_mb,
       allocate_workers,
       run_parallel_tasks,
       print_progress_spinner,
       display_batch_info

# Get available system memory in MB (Linux: /proc/meminfo)
function get_available_memory_mb()::Float64
    try
        meminfo = readlines("/proc/meminfo")
        for line in meminfo
            if startswith(line, "MemAvailable:")
                kb_str = split(line)[2]
                return parse(Float64, kb_str) / 1024.0
            end
        end
    catch e
        @debug "Failed to parse memory info: $(e.msg), using default 2048 MB"
    end
    return 2048.0
end

"""
    allocate_workers(n_tasks::Int, available_memory_mb::Float64, memory_per_task_mb::Float64=512.0, max_workers::Int=-1)

Allocate worker threads prioritizing 4 BLAS threads per worker when possible.
Returns (n_workers, blas_threads_per_worker, info_message)
"""
function allocate_workers(n_tasks::Int, available_memory_mb::Float64; 
                         memory_per_task_mb::Float64=512.0, max_workers::Int=-1)
    
    # Memory constraint: never plan to occupy more than 50% of available RAM
    safety_factor = 0.5
    max_by_memory = max(1, floor(Int, (available_memory_mb * safety_factor) / memory_per_task_mb))
    
    # PRIORITY: Prefer 4 threads per worker
    preferred_threads = 4
    workers = min(Threads.nthreads() ÷ preferred_threads, max_by_memory, n_tasks)
    
    # Fallback: Use all available threads if can't achieve 4 per worker
    if workers < 1
        workers = min(Threads.nthreads(), max_by_memory, n_tasks)
    end
    
    # Apply user override
    if max_workers > 0
        workers = min(workers, max_workers)
    end
    
    blas_threads = max(1, min(preferred_threads, Threads.nthreads() ÷ workers))
    
    info_msg = "Workers: $workers, BLAS threads/worker: $blas_threads, Available memory: $(round(Int, available_memory_mb))MB"
    
    return workers, blas_threads, info_msg
end

# Format task information for display
function format_task_info(idx::Int, params::Any)::String
    # Extract key parameters if they exist
    info_parts = ["Task $idx:"]
    
    # Try common parameter names
    param_keys = ["η_gt", "β_gt", "filepath_gt"]
    for key in param_keys
        if haskey(params, key)
            val = params[key]
            if isa(val, String) && length(val) > 60
                # Truncate long paths
                val = "..." * val[end-50:end]
            end
            push!(info_parts, "$key=$val")
        end
    end
    
    return join(info_parts, ", ")
end

# Print animated spinner with progress and optional live worker timing
function print_progress_spinner(completed::Int, total::Int, spinner_idx::Int; 
                               worker_timings::Dict{Int, Tuple{Int, Float64}}=Dict())
    spinners = ['◐', '◓', '◑', '◒']
    spinner = spinners[mod(spinner_idx, 4) + 1]
    pct = round(Int, (completed / total) * 100)
    
    # Build progress line
    progress_line = "$spinner Tasks: $completed / $total ($pct%)"
    
    # Add worker timing summary if available (append to same line)
    if !isempty(worker_timings)
        timing_strs = []
        for (w, (iters, total_time)) in sort(collect(worker_timings))
            if iters > 0
                avg_time = total_time / iters
                push!(timing_strs, "W$w: $(round(avg_time, digits=2))s")
            end
        end
        
        if !isempty(timing_strs)
            timing_summary = join(timing_strs, " ")
            print("\r$progress_line  |  $timing_summary")
        else
            print("\r$progress_line")
        end
    else
        print("\r$progress_line")
    end
    flush(stdout)
    sleep(0.01)  # Small delay to prevent overwriting
end

# Display batch information (minimal, informative output)
function display_batch_info(n_tasks::Int, n_workers::Int)
    @info "Executing $n_tasks tasks with $n_workers workers ($(Threads.nthreads()) threads available)"
end

"""
    run_parallel_tasks(task_list::Vector{Dict}, task_func::Function; 
                      max_workers::Int=-1, memory_per_task_mb::Float64=512.0)

Execute tasks in parallel using memory-aware worker allocation and work-stealing queue.

Args:
    task_list: Vector of Dict parameters, one per task
    task_func: Function that executes single task: task_func(param_dict)
    max_workers: Maximum workers to use (-1 = auto-detect)
    memory_per_task_mb: Estimated memory per task in MB

Returns:
    NamedTuple(completed::Int, worker_timings::Dict)
"""
function run_parallel_tasks(task_list::Vector, task_func::Function; 
                           max_workers::Int=-1, memory_per_task_mb::Float64=512.0)
    
    isempty(task_list) && return (completed=0, worker_timings=Dict())
    
    n_tasks = length(task_list)
    available_mb = get_available_memory_mb()
    workers, blas_threads, info_msg = allocate_workers(n_tasks, available_mb; 
                                                        memory_per_task_mb, max_workers)
    
    display_batch_info(n_tasks, workers)
    
    # Only warn if threads are critically low
    if blas_threads == 1 && workers > 1
        @warn "Only 1 BLAS thread per worker (undersubscribed). Consider reducing workers."
    end
    
    LinearAlgebra.BLAS.set_num_threads(blas_threads)
    
    completed = Threads.Atomic{Int}(0)
    progress_lock = ReentrantLock()
    worker_timings = Dict{Int, Tuple{Int, Float64}}()
    
    ch = Channel{Int}(n_tasks)
    @sync begin
        # Populate channel
        for i in 1:n_tasks
            put!(ch, i)
        end
        close(ch)
        
        # Spawn workers
        for w in 1:workers
            Threads.@spawn begin
                worker_start = time()
                worker_iters = 0
                timing_update_interval = max(1, workers)
                
                while true
                    try
                        idx = take!(ch)
                        
                        try
                            # Execute task with output redirection
                            redirect_stdout(devnull) do
                                redirect_stderr(devnull) do
                                    task_func(task_list[idx])
                                end
                            end
                            
                            Base.GC.gc(true)
                            worker_iters += 1
                            Threads.atomic_add!(completed, 1)
                            
                            # Log task completion
                            current_completed = completed[]
                            task_info = format_task_info(idx, task_list[idx])
                            if mod(current_completed, max(1, n_tasks ÷ 4)) == 0 || current_completed == n_tasks
                                @info "[$current_completed / $n_tasks] Completed: $task_info"
                            end
                            
                            lock(progress_lock) do
                                current = completed[]
                                # Update timings periodically
                                if mod(worker_iters, timing_update_interval) == 0
                                    elapsed = time() - worker_start
                                    worker_timings[w] = (worker_iters, elapsed)
                                end
                                print_progress_spinner(current, n_tasks, current; worker_timings=worker_timings)
                                flush(stdout)
                            end
                        catch err
                            bt = catch_backtrace()
                            @error "Task failed for index $idx" exception=(err, bt)
                        end
                    catch e
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
    
    # Print final summary (clear progress line and show completion)
    lock(progress_lock) do
        print("\r" * " "^100)  # Clear the progress line
        print("\r✓ All tasks completed\n")
    end
    
    return (completed=completed[], worker_timings=worker_timings)
end

end # module ParallelExecution
