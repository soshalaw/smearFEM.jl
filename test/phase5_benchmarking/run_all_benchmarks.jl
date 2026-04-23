"""
    run_all_benchmarks.jl

Phase 5 Master Orchestrator

Runs all Phase 5 benchmarking suites (5.1-5.6) in sequence and
generates a comprehensive performance report.
"""

using smearFEM
using Dates

println("\n" * "="^80)
println("PHASE 5: COMPREHENSIVE PERFORMANCE BENCHMARKING SUITE")
println("="^80)
println("Start Time: $(now())\n")

# Define benchmark scripts
benchmarks = [
    ("Phase 5.1: Solver Performance Baseline", "profile_solver_performance.jl"),
    ("Phase 5.3: Assembly Breakdown Analysis", "profile_assembly_breakdown.jl"),
    ("Phase 5.4: GMRES Configuration Tuning", "profile_gmres_tuning.jl"),
    ("Phase 5.5: GPU vs CPU Comparison", "profile_gpu_vs_cpu.jl"),
    ("Phase 5.6: Memory Profiling", "profile_memory_usage.jl"),
]

results_summary = []

for (phase_name, script) in benchmarks
    println("\n" * "="^80)
    println("Running: $phase_name")
    println("="^80)
    
    try
        # Get script path relative to this file
        script_path = joinpath(dirname(@__FILE__), script)
        
        if isfile(script_path)
            println("\n[EXEC] $script\n")
            include(script_path)
            push!(results_summary, (phase_name, "✓ PASS"))
            println("\n[SUCCESS] $phase_name completed successfully\n")
        else
            println("\n[SKIP] Script not found: $script_path\n")
            push!(results_summary, (phase_name, "⊘ SKIP"))
        end
    catch e
        println("\n[ERROR] $phase_name failed: $e\n")
        push!(results_summary, (phase_name, "✗ FAIL"))
    end
end

# Final Report
println("\n" * "="^80)
println("PHASE 5 BENCHMARKING SUITE - FINAL REPORT")
println("="^80)

println("\nExecution Summary:")
for (phase, status) in results_summary
    println("  $status  $phase")
end

println("\nEnd Time: $(now())")

# Performance Summary
println("\n" * "="^80)
println("KEY PERFORMANCE METRICS")
println("="^80)

println("\n✓ Phase 5 Benchmarking Complete!")
println("\nNext Steps:")
println("  1. Review assembly breakdown for optimization targets")
println("  2. Apply GMRES tuning recommendations to solver_integration.jl")
println("  3. Monitor GPU acceleration once hardware available")
println("  4. Integrate real mesh systems from squeeze_stokes.jl")
println("  5. Proceed to Phase 6: Deployment & Production Tuning")
