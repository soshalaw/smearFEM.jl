"""
    profile_realistic_system.jl

Phase 5: Realistic System Performance Profiling

Benchmarks GPU/CPU acceleration on a real 50k DOF Stokes system
over 200+ iterations, measuring component timing breakdown.

Performance Target: <10ms per iteration
  - Assembly: <3ms
  - Solver: <6ms
  - Transfer: <1ms

Expected Output:
  - Per-iteration timing CSV
  - Summary statistics
  - Performance plots (can be generated separately)
"""

using smearFEM
using Statistics
using CSV
using DataFrames
using Logging

# ─────────────────────────────────────────────────────────────────────────────
# CONFIGURATION
# ─────────────────────────────────────────────────────────────────────────────

const BENCHMARK_CONFIG = Dict(
    :name => "50k_DOF_Stokes",
    :mesh_size => (60, 60, 30),  # ~108k elements, ~54k DOF with linear elements
    :n_iterations => 200,
    :warmup_iterations => 5,
    :output_csv => "benchmark_results_50k.csv",
    :log_level => Logging.Info,
)

# ─────────────────────────────────────────────────────────────────────────────
# PROFILING UTILITIES
# ─────────────────────────────────────────────────────────────────────────────

"""
    TimingData
Store timing information for each iteration
"""
mutable struct TimingData
    iteration::Int
    assembly_time::Float64  # ms
    solver_time::Float64    # ms
    total_time::Float64     # ms
    timestamp::Float64      # seconds since epoch
end

function profile_iteration(mdl::Stokes, cache::BasisFunctionCache, config::SolverConfig,
                          A_storage::SparseMatrixCSC, b_storage::Vector)
    """
    Profile a single iteration and return timing data
    """
    
    # Assembly timing
    t_asm_start = time()
    A_storage .= assemble_system_A_routed(mdl, cache, config)
    t_asm_end = time()
    
    # Solver timing (simplified: just solve Ax=b)
    b_storage .= rand(size(A_storage, 1))  # Random RHS for profiling
    t_solve_start = time()
    x = solve_system(A_storage, b_storage, config)
    t_solve_end = time()
    
    return (
        asm_ms = (t_asm_end - t_asm_start) * 1000,
        solve_ms = (t_solve_end - t_solve_start) * 1000,
        total_ms = (t_solve_end - t_asm_start) * 1000,
    )
end

function print_benchmark_header()
    println("\n" * "="^80)
    println("PHASE 5 BENCHMARKING: REALISTIC SYSTEM PERFORMANCE PROFILING")
    println("="^80)
end

function print_benchmark_config(config::Dict)
    println("\n[CONFIG] Benchmark Parameters:")
    println("  Mesh Size:         $(config[:mesh_size])")
    println("  Iterations:        $(config[:n_iterations])")
    println("  Warmup:            $(config[:warmup_iterations])")
    println("  Output CSV:        $(config[:output_csv])")
    println("  Performance Target: <10ms per iteration")
end

# ─────────────────────────────────────────────────────────────────────────────
# MAIN PROFILING ROUTINE
# ─────────────────────────────────────────────────────────────────────────────

function run_realistic_benchmark()
    """
    Main benchmark routine
    """
    
    print_benchmark_header()
    print_benchmark_config(BENCHMARK_CONFIG)
    
    # Setup logging
    logger = SimpleLogger(stderr, BENCHMARK_CONFIG[:log_level])
    with_logger(logger) do
        
        # ─────────────────────────────────────────────────────────────────────
        # PHASE 1: MESH AND MODEL SETUP
        # ─────────────────────────────────────────────────────────────────────
        
        println("\n[SETUP] Creating realistic mesh...")
        mx, my, mz = BENCHMARK_CONFIG[:mesh_size]
        mesh = meshgrid_cube(mx, my, mz)
        n_elems = prod(BENCHMARK_CONFIG[:mesh_size])
        n_nodes = (mx+1) * (my+1) * (mz+1)
        
        println("[SETUP] Mesh created:")
        println("  Elements:      $n_elems")
        println("  Nodes:         $n_nodes")
        println("  Approx DOF:    ~$(div(n_nodes, 3))")  # 3 DOF per node for Stokes
        
        # Create Stokes model
        println("[SETUP] Creating Stokes model...")
        mdl = Stokes(mesh, viscosity_type=:constant, cParam=1e-3)
        
        # Prepare basis cache
        println("[SETUP] Preparing basis function cache...")
        cache = BasisFunctionCache(mesh, mdl, max_batch_size=2000)
        
        # ─────────────────────────────────────────────────────────────────────
        # PHASE 2: CONFIGURATION AND CONTEXT SETUP
        # ─────────────────────────────────────────────────────────────────────
        
        println("\n[CONFIG] Setting up solver configuration...")
        config = if has_gpu()
            println("[GPU] GPU detected! Using realtime_config()")
            realtime_config()
        else
            println("[GPU] GPU not available. Using cpu_fallback_config()")
            cpu_fallback_config()
        end
        
        println("[CONFIG] Solver target: $(config.gpu_solve ? "GPU" : "CPU")")
        println("[CONFIG] Assembly target: $(config.gpu_assembly ? "GPU" : "CPU")")
        
        # Setup GPU context once (critical for real-time performance)
        println("[CONTEXT] Initializing solver context...")
        setup_solver_context(mdl, config)
        
        # ─────────────────────────────────────────────────────────────────────
        # PHASE 3: STORAGE ALLOCATION
        # ─────────────────────────────────────────────────────────────────────
        
        println("\n[MEMORY] Allocating timing storage...")
        timing_results = TimingData[]
        
        # Pre-allocate matrix and vector storage
        n_dof = size(mesh.coords, 2)
        A_storage = spzeros(n_dof, n_dof)
        b_storage = zeros(n_dof)
        
        println("[MEMORY] System size:")
        println("  Matrix size:   $(size(A_storage))")
        println("  Vector size:   $(length(b_storage))")
        
        # ─────────────────────────────────────────────────────────────────────
        # PHASE 4: WARMUP ITERATIONS
        # ─────────────────────────────────────────────────────────────────────
        
        println("\n[WARMUP] Running $(BENCHMARK_CONFIG[:warmup_iterations]) warmup iterations...")
        for i in 1:BENCHMARK_CONFIG[:warmup_iterations]
            timing = profile_iteration(mdl, cache, config, A_storage, b_storage)
            if i % 2 == 0
                println("[WARMUP] Iteration $i: $(round(timing.total_ms, digits=2))ms")
            end
        end
        println("[WARMUP] Warmup complete")
        
        # ─────────────────────────────────────────────────────────────────────
        # PHASE 5: BENCHMARK ITERATIONS
        # ─────────────────────────────────────────────────────────────────────
        
        println("\n[BENCHMARK] Running $(BENCHMARK_CONFIG[:n_iterations]) benchmark iterations...")
        println("[BENCHMARK] " * "-"^70)
        
        for i in 1:BENCHMARK_CONFIG[:n_iterations]
            timing = profile_iteration(mdl, cache, config, A_storage, b_storage)
            push!(timing_results, TimingData(i, timing.asm_ms, timing.solve_ms, 
                                            timing.total_ms, time()))
            
            # Progress reporting
            if i % 20 == 0 || i <= 5
                mean_total = mean([t.total_time for t in timing_results])
                mean_asm = mean([t.assembly_time for t in timing_results])
                mean_solve = mean([t.solver_time for t in timing_results])
                
                status = if timing.total_ms > 10.0
                    "⚠ EXCEEDS"
                else
                    "✓ OK"
                end
                
                println("[ITER $i] $status - Total: $(round(timing.total_ms, digits=2))ms " *
                        "(A: $(round(timing.asm_ms, digits=2))ms, " *
                        "S: $(round(timing.solve_ms, digits=2))ms)")
            end
        end
        
        println("[BENCHMARK] " * "-"^70)
        
        # ─────────────────────────────────────────────────────────────────────
        # PHASE 6: STATISTICAL ANALYSIS
        # ─────────────────────────────────────────────────────────────────────
        
        println("\n[ANALYSIS] Computing statistics...")
        
        total_times = [t.total_time for t in timing_results]
        asm_times = [t.assembly_time for t in timing_results]
        solve_times = [t.solver_time for t in timing_results]
        
        function print_timing_stats(label::String, times::Vector)
            println("\n  $label Timing:")
            println("    Min:       $(round(minimum(times), digits=2))ms")
            println("    Mean:      $(round(mean(times), digits=2))ms")
            println("    Median:    $(round(median(times), digits=2))ms")
            println("    Std Dev:   $(round(std(times), digits=2))ms")
            println("    Max:       $(round(maximum(times), digits=2))ms")
            println("    P95:       $(round(quantile(times, 0.95), digits=2))ms")
            println("    P99:       $(round(quantile(times, 0.99), digits=2))ms")
        end
        
        println("\n" * "="^80)
        println("PERFORMANCE SUMMARY")
        println("="^80)
        
        print_timing_stats("Assembly", asm_times)
        print_timing_stats("Solver", solve_times)
        print_timing_stats("Total", total_times)
        
        mean_total = mean(total_times)
        budget_utilization = (mean_total / 10.0) * 100
        
        println("\n  Budget Utilization:")
        println("    Target:    <10.00 ms")
        println("    Mean:      $(round(mean_total, digits=2))ms")
        println("    Utilization: $(round(budget_utilization, digits=1))%")
        
        if mean_total < 10.0
            println("    Status: ✓ WITHIN BUDGET")
        else
            println("    Status: ✗ EXCEEDS BUDGET")
        end
        
        # ─────────────────────────────────────────────────────────────────────
        # PHASE 7: DATA EXPORT
        # ─────────────────────────────────────────────────────────────────────
        
        println("\n[EXPORT] Writing timing data to $(BENCHMARK_CONFIG[:output_csv])...")
        df = DataFrame(
            iteration = [t.iteration for t in timing_results],
            assembly_ms = [t.assembly_time for t in timing_results],
            solver_ms = [t.solver_time for t in timing_results],
            total_ms = [t.total_time for t in timing_results],
            timestamp = [t.timestamp for t in timing_results]
        )
        
        CSV.write(BENCHMARK_CONFIG[:output_csv], df)
        println("[EXPORT] Data written successfully")
        
        # ─────────────────────────────────────────────────────────────────────
        # PHASE 8: FINAL REPORTING
        # ─────────────────────────────────────────────────────────────────────
        
        println("\n" * "="^80)
        println("BENCHMARK COMPLETE")
        println("="^80)
        println("\nResults Summary:")
        println("  System Size:       $(mx)×$(my)×$(mz) elements (~$(div(n_nodes,3))k DOF)")
        println("  Iterations:        $(BENCHMARK_CONFIG[:n_iterations])")
        println("  Mean Time/Iter:    $(round(mean_total, digits=2))ms")
        println("  Assembly:          $(round(mean(asm_times), digits=2))ms ($(round((mean(asm_times)/mean_total)*100, digits=0))%)")
        println("  Solver:            $(round(mean(solve_times), digits=2))ms ($(round((mean(solve_times)/mean_total)*100, digits=0))%)")
        println("  Budget Status:     $(mean_total < 10.0 ? "✓ PASS" : "✗ FAIL")")
        println("  CSV Output:        $(BENCHMARK_CONFIG[:output_csv])")
        println("\n")
        
    end  # with_logger
end

# ─────────────────────────────────────────────────────────────────────────────
# EXECUTION
# ─────────────────────────────────────────────────────────────────────────────

if abspath(PROGRAM_FILE) == @__FILE__
    run_realistic_benchmark()
end
