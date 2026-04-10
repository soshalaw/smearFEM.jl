"""
    run_all_benchmarks.jl

Comprehensive benchmarking suite runner for the Stokes solver.
Executes all benchmarking scripts in a logical order and collects results.

Run with: julia --project run_all_benchmarks.jl
Or with optimizations: julia -O3 --threads=auto --check-bounds=no --project run_all_benchmarks.jl
"""

using Printf
using Dates

# Color codes for output
const RESET = "\e[0m"
const BOLD = "\e[1m"
const GREEN = "\e[32m"
const BLUE = "\e[34m"
const YELLOW = "\e[33m"
const RED = "\e[31m"

println()
println("$(BOLD)$(BLUE)===== COMPREHENSIVE STOKES SOLVER BENCHMARKING SUITE =====$(RESET)")
println()

# System information
println("$(BOLD)SYSTEM INFORMATION$(RESET)")
println("  Julia version: $(VERSION)")
println("  Threads: $(Threads.nthreads())")
try
    using LinearAlgebra
    println("  BLAS threads: $(BLAS.get_num_threads())")
catch
    println("  BLAS threads: N/A")
end
println("  Timestamp: $(now())")
println("  Optimization flags: -O3 (if run with -O3), --check-bounds=no (if applicable)")
println()

# List of benchmarking scripts to run
scripts = [
    ("detailed_profiler.jl", "Component-level Performance Analysis"),
    ("profile_100hz.jl", "100Hz Optimization Validation"),
    ("iterative_solver_profiler.jl", "Iterative Solver Performance"),
    ("profileview_profiler.jl", "Flamegraph Visualization Analysis"),
]

results = []
failed = []

# Run each benchmarking script
for (i, (script, description)) in enumerate(scripts)
    println("---" * "-"^70)
    println("[$(i)/$(length(scripts))] Running: $description")
    println("       Script: $script")
    println("---" * "-"^70)
    
    start_time = time()
    try
        include(script)
        elapsed = time() - start_time
        push!(results, (script, description, elapsed, "SUCCESS"))
        println("\n$GREEN[OK]$RESET $script completed in $(round(elapsed, digits=2))s\n")
    catch e
        elapsed = time() - start_time
        push!(failed, (script, description, string(e)))
        println("\n$RED[FAIL]$RESET $script failed after $(round(elapsed, digits=2))s")
        println("  Error: $(string(e))\n")
    end
end

# Summary Report
println("\n")
println("$(BOLD)$(BLUE)===== BENCHMARKING SUMMARY =====$(RESET)")
println()

println("$(BOLD)RESULTS$(RESET)")
total_time = 0.0
for (script, description, elapsed, status) in results
    color = GREEN
    println("  $color[OK]$RESET $(lpad(script, 30)) .. $status ($(round(elapsed, digits=2))s)")
    total_time += elapsed
end

if !isempty(failed)
    println("\n$(BOLD)$(RED)FAILED BENCHMARKS$(RESET)")
    for (script, description, error) in failed
        println("  $RED[FAIL]$RESET $(lpad(script, 30)) .. FAILED")
    end
end

println("\n$(BOLD)PERFORMANCE SUMMARY$(RESET)")
println("  Total benchmarking time: $(round(total_time, digits=2))s")
println("  Successful: $(length(results))/$(length(scripts))")
println("  Failed: $(length(failed))/$(length(scripts))")

if length(failed) == 0
    println("\n$GREEN$(BOLD)[OK] All benchmarks completed successfully!$(RESET)")
else
    println("\n$YELLOW$(BOLD)[WARNING] Some benchmarks failed. See details above.$(RESET)")
end

println("\n$(BOLD)NEXT STEPS$(RESET)")
println("  1. Review individual benchmark output above")
println("  2. Check test/benchmarking/*.jl for detailed analysis")
println("  3. Use run_100hz.sh in performance/ for optimized execution")
println("  4. Compare results against baseline performance")
println()
