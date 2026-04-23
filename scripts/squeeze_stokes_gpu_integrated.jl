"""
    squeeze_stokes_gpu_integrated.jl

Phase 4b Example: GPU-Integrated Assembly and Solving

Demonstrates how to use the new routed assembly functions (assemble_system_A_routed,
assemble_system_B_routed, apply_boundary_conditions_routed) in a time-stepping loop
with GPU/CPU configuration.

This example shows the minimal modifications needed to convert existing squeeze_stokes.jl
code to use the GPU acceleration layer.

Performance Target: <10ms per iteration (3ms assembly + 6ms solver + 1ms transfer)
"""

using smearFEM
using Logging

# ─────────────────────────────────────────────────────────────────────────────
# 1. SETUP PHASE: Initialize model and configuration
# ─────────────────────────────────────────────────────────────────────────────

println("[PHASE 4B] GPU-Integrated Assembly Example")
println("=" ^ 70)

# Create squeeze flow scenario on coarse mesh (for example)
println("\n[SETUP] Creating squeeze flow scenario...")
mesh = meshgrid_cube(10, 10, 10)  # 10x10x10 elements = 1331 DOF for linear elements
scene = SqueezeFlow(mesh, β=1.0, q_d=0.01, control=:velocity)

println("[SETUP] Creating Stokes model...")
mdl = Stokes(mesh, viscosity_type=:constant, cParam=1e-3)

# Prepare basis function cache for time-stepping reuse
println("[SETUP] Preparing basis function cache...")
cache = BasisFunctionCache(mesh, mdl, max_batch_size=1000)

# ─────────────────────────────────────────────────────────────────────────────
# 2. CONFIGURATION PHASE: Choose GPU or CPU
# ─────────────────────────────────────────────────────────────────────────────

# Configuration preset based on GPU availability
if has_gpu()
    println("\n[CONFIG] GPU detected! Using realtime_config()...")
    config = realtime_config()
else
    println("\n[CONFIG] GPU not available. Using cpu_fallback_config()...")
    config = cpu_fallback_config()
end

# Print configuration status
println("[CONFIG] Assembly target: $(config.gpu_assembly ? "GPU" : "CPU")")
println("[CONFIG] Solver target: $(config.gpu_solve ? "GPU" : "CPU")")
println("[CONFIG] Keep on GPU: $(config.keep_on_gpu)")
println("[CONFIG] Cache preconditioner: $(config.cache_precond)")

# Pre-allocate GPU context once (critical for budget adherence)
println("[CONTEXT] Setting up solver context...")
setup_solver_context(mdl, config)

# ─────────────────────────────────────────────────────────────────────────────
# 3. TIME-STEPPING PHASE: Loop with routed assembly
# ─────────────────────────────────────────────────────────────────────────────

println("\n[TIME-STEP] Starting simulation loop...")
println("-" ^ 70)

# Time parameters
n_steps = 5  # Short run for demonstration
dt = 0.001
time_steps = range(0, step=dt, length=n_steps)

# Storage (reused across time steps)
n_dof = size(mesh.coords, 2)
_A_bar = spzeros(n_dof, n_dof)
B = spzeros(n_dof, n_dof)
b = zeros(n_dof)
u = zeros(n_dof)

timing_log = []

# Main time-stepping loop with routed assembly
for (step, t) in enumerate(time_steps[1:end-1])
    step_time = @timed begin
        # ─────────────────────────────────────────────────────────────────
        # ROUTED ASSEMBLY (Phase 4b)
        # ─────────────────────────────────────────────────────────────────
        
        asm_A_time = @timed _A_bar .= assemble_system_A_routed(mdl, cache, config)
        asm_B_time = @timed B .= assemble_system_B_routed(mdl, cache, config)
        bc_time = @timed b .= apply_boundary_conditions_routed(mdl, cache, config)
        
        # Combine matrices (simplified - see squeeze_stokes.jl for full details)
        η = mdl.η[end]
        A_combined = _A_bar + (1/η) * B  # Simplified; real code has more terms
        
        # ─────────────────────────────────────────────────────────────────
        # UNIFIED SOLVER (Phase 4a)
        # ─────────────────────────────────────────────────────────────────
        solve_time = @timed u .= solve_system(A_combined, b, config)
        
        # Store timing for this step
        push!(timing_log, (
            step = step,
            t = t,
            asm_A = asm_A_time.time * 1000,      # ms
            asm_B = asm_B_time.time * 1000,      # ms
            bc = bc_time.time * 1000,            # ms
            solve = solve_time.time * 1000,      # ms
        ))
    end
    
    # Print per-step timing
    if step % 1 == 0
        asm_total = timing_log[step].asm_A + timing_log[step].asm_B + timing_log[step].bc
        solve = timing_log[step].solve
        total = asm_total + solve
        
        println("[STEP $step] t=$(round(t, digits=4)) | " *
                "Assembly: $(round(asm_total, digits=2))ms | " *
                "Solve: $(round(solve, digits=2))ms | " *
                "Total: $(round(total, digits=2))ms")
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# 4. REPORTING PHASE: Print performance summary
# ─────────────────────────────────────────────────────────────────────────────

println("\n" * "=" ^ 70)
println("[SUMMARY] Phase 4b GPU Integration Performance")
println("=" ^ 70)

if !isempty(timing_log)
    asm_times = [log.asm_A + log.asm_B + log.bc for log in timing_log]
    solve_times = [log.solve for log in timing_log]
    total_times = asm_times .+ solve_times
    
    println("\nAssembly Timing:")
    println("  Min:    $(round(minimum(asm_times), digits=2)) ms")
    println("  Mean:   $(round(mean(asm_times), digits=2)) ms")
    println("  Max:    $(round(maximum(asm_times), digits=2)) ms")
    
    println("\nSolver Timing:")
    println("  Min:    $(round(minimum(solve_times), digits=2)) ms")
    println("  Mean:   $(round(mean(solve_times), digits=2)) ms")
    println("  Max:    $(round(maximum(solve_times), digits=2)) ms")
    
    println("\nTotal Per-Iteration:")
    println("  Min:    $(round(minimum(total_times), digits=2)) ms")
    println("  Mean:   $(round(mean(total_times), digits=2)) ms")
    println("  Max:    $(round(maximum(total_times), digits=2)) ms")
    println("  Budget: <10.00 ms ✓" * (mean(total_times) < 10 ? "" : " ✗ BUDGET EXCEEDED"))
end

# Print assembly router status
println("\n[CONFIG] Final Assembly Configuration:")
print_assembly_status(config)

println("\n[COMPLETE] Phase 4b GPU Integration Example Finished!")
println("=" ^ 70)

# ─────────────────────────────────────────────────────────────────────────────
# USAGE INSTRUCTIONS FOR INTEGRATION
# ─────────────────────────────────────────────────────────────────────────────

"""
To integrate GPU acceleration into existing squeeze_stokes.jl:

1. At the top of squeeze_stokes.jl, add:
   using smearFEM
   config = has_gpu() ? realtime_config() : cpu_fallback_config()
   setup_solver_context(mdl, config)

2. Replace assembly calls in the time-stepping loop:
   
   OLD:
   _A_bar .= assemble_system_A(mdl, cache)
   B .= assemble_system_B(mdl, cache)
   b .= apply_boundary_conditions(mdl, cache)
   
   NEW:
   _A_bar .= assemble_system_A_routed(mdl, cache, config)
   B .= assemble_system_B_routed(mdl, cache, config)
   b .= apply_boundary_conditions_routed(mdl, cache, config)

3. Replace solve call:
   
   OLD:
   u = A_combined \\ b  # LU or other solver
   
   NEW:
   u = solve_system(A_combined, b, config)

4. For monitoring, add timing around the loop:
   time_step_with_config!(mdl, A, b, config)
   print_real_time_report()

The routed functions automatically select GPU or CPU based on config,
with fallback chains ensuring convergence at all levels.
"""
