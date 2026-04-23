"""
    profile_gpu_vs_cpu.jl

Phase 5.5: GPU vs CPU Performance Comparison

Compares solver performance on GPU vs CPU across different DOF sizes,
measuring speedup factors and PCIe transfer overhead.
"""

using smearFEM
using Statistics
using LinearAlgebra
using SparseArrays
using Random

println("\n" * "="^80)
println("PHASE 5.5: GPU VS CPU PERFORMANCE COMPARISON")
println("="^80)

# Configuration
DOF_SIZES = [1000, 5000, 10000]
N_RUNS = 5

has_gpu_device = has_gpu()
println("\nGPU Available: $has_gpu_device\n")

results = []

for dof in DOF_SIZES
    println("[TEST] Profiling $dof DOF - CPU vs GPU")
    
    # Create test system
    Random.seed!(42)
    A_diag = rand(dof) .+ dof
    A = spdiagm(0 => A_diag)
    
    for i in 1:min(20, dof-1)
        for j in i+1:min(i+5, dof)
            val = 0.1 * rand()
            A[i, j] = val
            A[j, i] = val
        end
    end
    
    b = rand(dof)
    
    # CPU benchmarking
    cpu_config = cpu_fallback_config()
    cpu_times = []
    for _ in 1:N_RUNS
        t1 = time()
        x_cpu = solve_system(A, b, cpu_config)
        t2 = time()
        push!(cpu_times, (t2 - t1) * 1000)
    end
    
    cpu_mean = mean(cpu_times)
    cpu_std = std(cpu_times)
    
    println("  CPU: $(round(cpu_mean, digits=3))ms ± $(round(cpu_std, digits=3))ms")
    
    # GPU benchmarking (if available)
    gpu_mean = NaN
    gpu_std = NaN
    speedup = NaN
    
    if has_gpu_device
        try
            gpu_config = realtime_config()
            gpu_times = []
            for _ in 1:N_RUNS
                t1 = time()
                x_gpu = solve_system(A, b, gpu_config)
                t2 = time()
                push!(gpu_times, (t2 - t1) * 1000)
            end
            
            gpu_mean = mean(gpu_times)
            gpu_std = std(gpu_times)
            speedup = cpu_mean / gpu_mean
            
            println("  GPU: $(round(gpu_mean, digits=3))ms ± $(round(gpu_std, digits=3))ms")
            println("  Speedup: $(round(speedup, digits=2))x")
        catch e
            println("  GPU: ERROR - $e")
        end
    else
        println("  GPU: Not available (CPU fallback only)")
    end
    
    push!(results, (
        dof = dof,
        cpu_ms = cpu_mean,
        gpu_ms = gpu_mean,
        speedup = speedup
    ))
end

# Summary
println("\n" * "="^80)
println("GPU VS CPU PERFORMANCE SUMMARY")
println("="^80)

println("\nDOF      | CPU (ms) | GPU (ms) | Speedup")
println("-" ^ 45)
for r in results
    if isnan(r.speedup)
        println("$(lpad(r.dof, 8)) | $(lpad(round(r.cpu_ms, digits=2), 8)) | N/A      | N/A")
    else
        println("$(lpad(r.dof, 8)) | $(lpad(round(r.cpu_ms, digits=2), 8)) | $(lpad(round(r.gpu_ms, digits=2), 8)) | $(lpad(round(r.speedup, digits=1), 7))x")
    end
end

if has_gpu_device
    gpu_times = [r.speedup for r in results if !isnan(r.speedup)]
    if !isempty(gpu_times)
        println("\nAverage GPU speedup: $(round(mean(gpu_times), digits=2))x")
        println("Max speedup: $(round(maximum(gpu_times), digits=2))x at $(results[argmax(gpu_times)].dof) DOF")
    end
else
    println("\n[NOTE] GPU acceleration not available on this system")
end

println("\n✓ GPU vs CPU performance comparison complete!")
