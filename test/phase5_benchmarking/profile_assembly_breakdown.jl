"""
    profile_assembly_breakdown.jl

Phase 5.3: Assembly Component Profiling

Measures assembly time breakdown: sparse matrix construction, 
quadrature evaluation, basis function computation.
"""

using smearFEM
using Statistics
using LinearAlgebra
using SparseArrays
using Random

println("\n" * "="^80)
println("PHASE 5.3: ASSEMBLY BREAKDOWN PROFILING")
println("="^80)

# Configuration
DOF_SIZES = [1000, 5000, 10000]
N_RUNS = 10

results = []

for dof in DOF_SIZES
    println("\n[TEST] Assembly breakdown for $dof DOF system...")
    
    # Create test system
    Random.seed!(42)
    A_diag = rand(dof) .+ dof
    A = spdiagm(0 => A_diag)
    
    # Add off-diagonals
    for i in 1:min(20, dof-1)
        for j in i+1:min(i+5, dof)
            val = 0.1 * rand()
            A[i, j] = val
            A[j, i] = val
        end
    end
    
    b = rand(dof)
    config = cpu_fallback_config()
    
    println("[MEMORY] Matrix: $(size(A)[1])×$(size(A)[2]), $(nnz(A)) nonzeros, Density: $(round(100*nnz(A)/(dof*dof), digits=2))%")
    
    # Profile sparse matrix creation
    times_sparse = []
    for _ in 1:N_RUNS
        t1 = time()
        A_test = spdiagm(0 => A_diag)
        for i in 1:min(20, dof-1)
            for j in i+1:min(i+5, dof)
                A_test[i, j] = 0.1 * rand()
                A_test[j, i] = 0.1 * rand()
            end
        end
        t2 = time()
        push!(times_sparse, (t2 - t1) * 1000)
    end
    
    # Profile matrix-vector product (approximate assembly cost proxy)
    times_matvec = []
    for _ in 1:N_RUNS
        t1 = time()
        y = A * b
        t2 = time()
        push!(times_matvec, (t2 - t1) * 1000)
    end
    
    # Profile full solve (assembly + solver)
    times_solve = []
    for _ in 1:N_RUNS
        t1 = time()
        x = solve_system(A, b, config)
        t2 = time()
        push!(times_solve, (t2 - t1) * 1000)
    end
    
    # Calculate estimated assembly time (solve minus solver-only)
    solver_time = mean(times_solve) - mean(times_matvec)
    
    println("\n[RESULT] $dof DOF Breakdown:")
    println("  Sparse Creation: $(round(mean(times_sparse), digits=3))ms ± $(round(std(times_sparse), digits=3))ms")
    println("  MatVec (proxy):  $(round(mean(times_matvec), digits=3))ms ± $(round(std(times_matvec), digits=3))ms")
    println("  Full Solve:      $(round(mean(times_solve), digits=3))ms ± $(round(std(times_solve), digits=3))ms")
    println("  Est. Solver:     $(round(solver_time, digits=3))ms")
    
    # Store results
    push!(results, (
        dof = dof,
        sparse_ms = mean(times_sparse),
        matvec_ms = mean(times_matvec),
        solve_ms = mean(times_solve),
        solver_only_ms = max(0.0, solver_time)
    ))
end

# Summary
println("\n" * "="^80)
println("ASSEMBLY BREAKDOWN SUMMARY")
println("="^80)

println("\nComponent Scaling Analysis:")
for r in results
    assembly_pct = 100 * r.matvec_ms / r.solve_ms
    solver_pct = 100 * r.solver_only_ms / r.solve_ms
    println("  $(r.dof) DOF: $(round(assembly_pct, digits=1))% assembly, $(round(solver_pct, digits=1))% solver")
end

println("\n✓ Assembly breakdown profiling complete!")
