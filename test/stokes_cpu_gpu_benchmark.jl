"""
CPU vs GPU Benchmarking for Stokes Simulation

Compares solver performance on the single_simulation.jl script
with CPU and GPU configurations.
"""

using smearFEM
using Dates
using Statistics
using SparseArrays
using LinearAlgebra

println("="^80)
println("STOKES SIMULATION: CPU vs GPU BENCHMARK")
println("="^80)

# Configuration for measurements
runs_per_config = 3
iterations_per_run = 10

# Test matrices sizes
test_cases = [
    (name="Small (100 DOF)", n_u=50, n_p=50),
    (name="Medium (500 DOF)", n_u=250, n_p=250),
    (name="Large (1000 DOF)", n_u=500, n_p=500),
]

function create_test_matrix(n_u::Int, n_p::Int)
    """Create a realistic Stokes matrix for testing"""
    n_total = n_u + n_p
    
    # Create velocity block (main matrix, sparse)
    A_uu = spdiagm(0 => 2.0 * ones(n_u)) + 0.1 * sprandn(n_u, n_u, 0.1)
    A_uu = (A_uu + A_uu') / 2  # Make symmetric
    A_uu = A_uu + n_u * I  # Add diagonal dominance
    
    # Create coupling matrices
    B = 0.1 * sprandn(n_p, n_u, 0.15)
    
    # Saddle point system
    A = [A_uu            B'
         B              spzeros(n_p, n_p)]
    
    b = randn(n_total)
    
    return sparse(A), b
end

println("\nBenchmarking different matrix sizes:")
println("="^80)

results = Dict()

for (case_name, n_u, n_p) in test_cases
    println("\n[TEST] $case_name (Velocity DOF: $n_u, Pressure DOF: $n_p)")
    
    # Create test system
    A, b = create_test_matrix(n_u, n_p)
    n_total = size(A, 1)
    
    println("  Matrix: $(size(A, 1))×$(size(A, 2)), nnz=$(nnz(A)), density=$(round(100*nnz(A)/(n_total^2); digits=2))%")
    
    # Test 1: CPU Configuration (standard)
    println("\n  CPU Configuration (standard GMRES):")
    cpu_config = smearFEM.cpu_fallback_config()
    println("    - solver: $(cpu_config.solver_type)")
    println("    - tol: $(cpu_config.tol), maxiter: $(cpu_config.maxiter), restart: $(cpu_config.gmres_restart)")
    
    cpu_times = Float64[]
    for run in 1:runs_per_config
        for _ in 1:iterations_per_run
            t_start = time()
            try
                x = smearFEM.solve_system(A, b, cpu_config)
            catch e
                # Solver may not converge on all systems, but we measure time anyway
            end
            t_end = time()
            push!(cpu_times, (t_end - t_start) * 1000)  # Convert to ms
        end
    end
    
    cpu_mean = mean(cpu_times)
    cpu_std = std(cpu_times)
    cpu_min = minimum(cpu_times)
    cpu_max = maximum(cpu_times)
    
    println("    Results ($(runs_per_config * iterations_per_run) iterations):")
    println("      Mean: $(round(cpu_mean; digits=2))ms ± $(round(cpu_std; digits=2))ms")
    println("      Min/Max: $(round(cpu_min; digits=2))ms / $(round(cpu_max; digits=2))ms")
    
    # Test 2: GPU Configuration (if available)
    println("\n  GPU Configuration (GPU GMRES):")
    gpu_config = smearFEM.realtime_config()
    println("    - solver: $(gpu_config.solver_type)")
    println("    - tol: $(gpu_config.tol), maxiter: $(gpu_config.maxiter), restart: $(gpu_config.gmres_restart)")
    
    gpu_times = Float64[]
    for run in 1:runs_per_config
        for _ in 1:iterations_per_run
            t_start = time()
            try
                x = smearFEM.solve_system(A, b, gpu_config)
            catch e
                # Solver may not converge or GPU not available
            end
            t_end = time()
            push!(gpu_times, (t_end - t_start) * 1000)  # Convert to ms
        end
    end
    
    gpu_mean = mean(gpu_times)
    gpu_std = std(gpu_times)
    gpu_min = minimum(gpu_times)
    gpu_max = maximum(gpu_times)
    
    println("    Results ($(runs_per_config * iterations_per_run) iterations):")
    println("      Mean: $(round(gpu_mean; digits=2))ms ± $(round(gpu_std; digits=2))ms")
    println("      Min/Max: $(round(gpu_min; digits=2))ms / $(round(gpu_max; digits=2))ms")
    
    # Calculate speedup
    speedup = cpu_mean / gpu_mean
    println("\n  Speedup (CPU/GPU): $(round(speedup; digits=2))x")
    
    if speedup > 1.5
        println("    ✓ GPU shows significant speedup")
    elseif speedup > 1.0
        println("    • GPU shows modest speedup")
    else
        println("    ⚠ GPU slower (may be CPU fallback or overhead dominant)")
    end
    
    results[case_name] = (
        cpu_mean = cpu_mean,
        cpu_std = cpu_std,
        gpu_mean = gpu_mean,
        gpu_std = gpu_std,
        speedup = speedup
    )
end

# Summary
println("\n" * "="^80)
println("BENCHMARK SUMMARY")
println("="^80)

println("\nPerformance Comparison:")
println("Case                    | CPU (ms)      | GPU (ms)      | Speedup")
println("-"^70)

for (case_name, result) in results
    cpu_str = "$(round(result.cpu_mean; digits=2)) ± $(round(result.cpu_std; digits=2))"
    gpu_str = "$(round(result.gpu_mean; digits=2)) ± $(round(result.gpu_std; digits=2))"
    speedup_str = "$(round(result.speedup; digits=2))x"
    
    # Pad the strings
    case_padded = rpad(case_name, 23)
    cpu_padded = rpad(cpu_str, 13)
    gpu_padded = rpad(gpu_str, 13)
    
    println("$case_padded | $cpu_padded | $gpu_padded | $speedup_str")
end

println("\n" * "="^80)
println("System Configuration")
println("="^80)

println("GPU Available: $(smearFEM.has_gpu())")
println("CPU Configuration: $(Sys.cpu_info())")
println("Julia Version: $(VERSION)")

println("\n✓ Benchmark complete!")
println("="^80)
