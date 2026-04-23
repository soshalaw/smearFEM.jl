"""
    profile_solver_performance.jl

Phase 5: Solver Performance Profiling (Simplified)

Benchmarks CPU solver performance on matrices of varying sizes,
measuring convergence and timing breakdown.
"""

using smearFEM
using Statistics
using LinearAlgebra
using SparseArrays
using Random

println("\n" * "="^80)
println("PHASE 5: SOLVER PERFORMANCE PROFILING")
println("="^80)

# Configuration
DOF_SIZES = [1000, 5000, 10000]
N_ITERATIONS = 30
WARMUP_ITERS = 5

results = []

for dof in DOF_SIZES
    println("\n[TEST] Profiling $dof DOF system...")
    
    # Setup
    config = cpu_fallback_config()  # CPU only
    
    # Create test system (random SPD matrix for Stokes-like problem)
    Random.seed!(42)
    A_diag = rand(dof) .+ dof  # Ensure positive definite
    A = spdiagm(0 => A_diag)
    
    # Add off-diagonal elements (sparse pattern)
    for i in 1:min(20, dof-1)
        for j in i+1:min(i+5, dof)
            val = 0.1 * rand()
            A[i, j] = val
            A[j, i] = val
        end
    end
    
    b = rand(dof)
    
    println("[MEMORY] Matrix: $(size(A)[1])×$(size(A)[2]), $(nnz(A)) nonzeros")
    
    # Warmup
    println("[WARMUP] Running $WARMUP_ITERS warmup iterations...")
    for _ in 1:WARMUP_ITERS
        x = solve_system(A, b, config)
    end
    
    # Benchmark
    println("[BENCHMARK] Running $N_ITERATIONS iterations...")
    solve_times = []
    
    for iter in 1:N_ITERATIONS
        t_start = time()
        x = solve_system(A, b, config)
        t_end = time()
        
        solve_ms = (t_end - t_start) * 1000
        push!(solve_times, solve_ms)
        
        if iter % 10 == 0
            println("[ITER $iter] $(round(solve_ms, digits=3))ms")
        end
    end
    
    # Statistics
    mean_time = mean(solve_times)
    std_time = std(solve_times)
    min_time = minimum(solve_times)
    max_time = maximum(solve_times)
    p95_time = quantile(solve_times, 0.95)
    
    println("\n[RESULT] $dof DOF:")
    println("  Mean:   $(round(mean_time, digits=3))ms")
    println("  Std:    $(round(std_time, digits=3))ms")
    println("  Min:    $(round(min_time, digits=3))ms")
    println("  Max:    $(round(max_time, digits=3))ms")
    println("  P95:    $(round(p95_time, digits=3))ms")
    
    # Store results
    for (i, t) in enumerate(solve_times)
        push!(results, (dof=dof, iteration=i, solver_ms=t))
    end
end

# Summary
println("\n" * "="^80)
println("BENCHMARKING COMPLETE")
println("="^80)

println("\nPerformance Summary by DOF:")
for dof in DOF_SIZES
    subset_times = [r.solver_ms for r in results if r.dof == dof]
    mean_time = mean(subset_times)
    println("  $dof DOF: $(round(mean_time, digits=3))ms avg")
end

println("\n✓ Solver performance profiling complete!")
