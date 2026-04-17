#!/usr/bin/env julia
"""
Standalone test script to demonstrate ParallelExecution output
Run with: julia test_parallel_output.jl
"""

include("scripts/ParallelExecution.jl")
using .ParallelExecution

# Simple task function that takes time
function test_task(params::Dict)
    sleep(1.5)  # 1.5 seconds per task
end

# Create test tasks with realistic parameters
task_list = [
    Dict("η_gt" => 100.0, "β_gt" => 0.01, "filepath_gt" => "/data/ground_truth/sim_data/Stokes/force/constant/Q2_6/1"),
    Dict("η_gt" => 100.0, "β_gt" => 0.1, "filepath_gt" => "/data/ground_truth/sim_data/Stokes/force/constant/Q2_6/2"),
    Dict("η_gt" => 100.0, "β_gt" => 1.0, "filepath_gt" => "/data/ground_truth/sim_data/Stokes/force/constant/Q2_6/3"),
    Dict("η_gt" => 100.0, "β_gt" => 10.0, "filepath_gt" => "/data/ground_truth/sim_data/Stokes/force/constant/Q2_6/4"),
    Dict("η_gt" => 100.0, "β_gt" => 50.0, "filepath_gt" => "/data/ground_truth/sim_data/Stokes/force/constant/Q2_6/5"),
    Dict("η_gt" => 100.0, "β_gt" => 100.0, "filepath_gt" => "/data/ground_truth/sim_data/Stokes/force/constant/Q2_6/6"),
]

println("\n" * "="^80)
println("PARALLELEXECUTION OUTPUT TEST")
println("="^80)
println("\nThis test demonstrates all output messages from the ParallelExecution module.")
println("You should see:\n")
println("  1. Task info (parameters)")
println("  2. ✓ Task X completed (for each task)")
println("  3. Live progress bar with worker timing")
println("  4. Final completion message\n")
println("="^80 * "\n")

# Run with automatic worker allocation
run_parallel_tasks(task_list, test_task; max_workers=-1)

println("\n" * "="^80)
println("✓ TEST COMPLETE - All outputs demonstrated above")
println("="^80 * "\n")
