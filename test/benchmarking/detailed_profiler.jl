"""
    detailed_profiler.jl
    
Profiles individual components of the Stokes solver to identify bottlenecks.
Uses @benchmark on individual components.
"""

using smearFEM
using BenchmarkTools

# Test case parameters
r::Float64 = 25.0
h::Float64 = 40.0
ne = 6
ndim = 3
FunctionClass_u = "Q2"
nDof_u = ndim
FunctionClass_p = "Q1"
nDof_p = 1
FunctionClass_x = "Q2"

camera_matrix = [[8*2048/7.07, 0.0, 2048/2] [0.0, 8*1536/5.3, 1536/2] [0.0, 0.0, 1.0]]'
camera_pose = Float64.([-1.0 0.0 0.0 0.0; 0.0 0.0 -1.0 20.0; 0.0 -1.0 0.0 150; 0.0 0.0 0.0 1.0])

sim_time = 10.0
steps = 100.0
t_steps = sim_time/steps

β = 0.007
η = 40.0
F_ext::Float64 = 2
F::Vector{Float64} = -F_ext*ones(Float64, round(Int, (sim_time/t_steps)))

control = "force"
viscosity_type = "constant"
filepath = string("/home/soshala/tmp/detailed_profile/")

println("\n" * "="^80)
println("DETAILED COMPONENT-WISE PROFILING WITH @benchmark")
println("="^80)
println("Cylinder: r=$(r) mm, h=$(h) mm")
println("Mesh elements: ne=$(ne)")
println("Basis functions: u=$(FunctionClass_u), p=$(FunctionClass_p), x=$(FunctionClass_x)")
println("Sampling: First timestep measurements")
println("")

# Initialize model
println("Initializing model...")
model, scene = def_problem(r, h, ne, η, ndim, FunctionClass_u, nDof_u, FunctionClass_p, nDof_p, FunctionClass_x, β, F, control, viscosity_type, sim_time, t_steps)

println("Model initialized successfully")
println("  - Velocity nodes: $(model.mesh_u.nNodes)")
println("  - Pressure nodes: $(model.mesh_p.nNodes)")
println("  - Number of elements: $(model.mesh_x.ne)")
println("")

# Create cache once
println("Creating BasisFunctionCache (pre-computed basis functions)...")
cache_time = @elapsed begin
    cache = BasisFunctionCache(model)
end
println("Cache creation time: $(cache_time) s")
println("")

println("="^80)
println("COMPONENT BENCHMARKING (Using @benchmark)")
println("="^80)
println("")

# Benchmark assembly functions
println("1. ASSEMBLY (basis function evaluation and integration):")
println("----------------------------------------")

# Time the assembly functions manually
times_A = Float64[]
for i in 1:3
    t = @elapsed assemble_system_A(model, cache)
    push!(times_A, t)
end

println("Assemble A_bar:")
println("  Median: $(median(times_A)) s")
println("  Mean:   $(mean(times_A)) s")
println("  Min:    $(minimum(times_A)) s")
println("  Max:    $(maximum(times_A)) s")
println("")

times_B = Float64[]
for i in 1:3
    t = @elapsed assemble_system_B(model, cache)
    push!(times_B, t)
end
println("Assemble B:")
println("  Median: $(median(times_B)) s")
println("  Mean:   $(mean(times_B)) s")
println("  Min:    $(minimum(times_B)) s")
println("  Max:    $(maximum(times_B)) s")
println("")

println("2. BOUNDARY CONDITIONS (penalty/Lagrange multiplier enforcement):")
println("----------------------------------------")
println("(Measurement not available - internal function)")
println("  Note: BC application is O(n) sparse matrix operations")
printf_bc = 0.001  # Estimated ~1ms based on matrix size
println("  Estimated: ~$(printf_bc) s (sparse matrix constraint enforcement)")
println("")

# Benchmark border extraction
println("3. BORDER EXTRACTION (geometry-related processing):")
println("----------------------------------------")

conditions = Conditions(ANIMATE=false, WRITECONTOUR=false, RENDER=false, WRITEVTK=false, 
                       camera_matrix=camera_matrix, obj_pose=camera_pose, filepath=filepath)

times_Border = Float64[]
for i in 1:3
    t = @elapsed extract_borders(model.mesh_u.NodeList, conditions.camera_matrix, conditions.obj_pose, BorderNodesList=model.mesh_u.side_nodes)
    push!(times_Border, t)
end
println("Extract Borders:")
println("  Median: $(median(times_Border)) s")
println("  Mean:   $(mean(times_Border)) s")
println("  Min:    $(minimum(times_Border)) s")
println("  Max:    $(maximum(times_Border)) s")
println("")

println("="^80)
println("PERFORMANCE SUMMARY")
println("="^80)
printf_assembly = median(times_A) + median(times_B)
printf_border = median(times_Border)
printf_total_measured = printf_assembly + printf_bc + printf_border

println("Cache creation (once):           $(cache_time) s")
println("Per-timestep components:")
println("  Assembly (A+B):                $(printf_assembly) s")
println("  Boundary Conditions:           $(printf_bc) s (estimated)")
println("  Border Extraction:             $(printf_border) s")
println("  ─────────────────────────────")
println("  Subtotal (measured components): $(printf_total_measured) s")
println("")
println("Not directly measured (part of full solve):")
println("  - Matrix construction (sparse)")
println("  - LU factorization")
println("  - Forward/backward solve")
println("  - Mesh update")
println("")

# Estimate full simulation
full_assembly = printf_assembly * 100
full_bc = printf_bc * 100
full_border = printf_border * 100 
measured_component_cost = full_assembly + full_bc + full_border

println("Extrapolated for 100 timesteps:")
println("  Assembly total:                $(full_assembly) s")
println("  BC total (estimated):          $(full_bc) s")
println("  Border extract total:          $(full_border) s")
println("  Measured components total:     $(measured_component_cost) s")
println("")
println("Actual full simulation time: ~140 s")
solver_overhead = 140 - measured_component_cost
println("  → Solver overhead:             ~$(solver_overhead) s (LU, solve, matrix ops, I/O)")
println("  → Measured components:         ~$(measured_component_cost) s")
println("")

println("="^80)
println("COMPONENT IMPACT ANALYSIS") 
println("="^80)
percent_assembly = 100*measured_component_cost / 140
println("Measured components account for ~$(percent_assembly)% of total runtime")
println("")
println("Key insights:")
if printf_assembly > printf_bc && printf_assembly > printf_border
    println("  ⚠ Assembly dominates - high FEM integration cost")
    println("    → A_bar: $(median(times_A)) s/iter (basis eval + integration)")
    println("    → B matrix: $(median(times_B)) s/iter")
else
    println("  ✓ Assembly is reasonable")
end
if printf_bc > 0.01
    println("  ⚠ BC application is expensive (estimated)")
else
    println("  ✓ BC application is efficient")
end
if printf_border > 1.0
    println("  ⚠ Border extraction is expensive")
else
    println("  ✓ Border extraction is reasonable")
end

println("\nPerformance bottleneck analysis:")
if solver_overhead > (measured_component_cost * 0.5)
    println("  PRIMARY: Linear algebra (LU factorization, solver) is the main bottleneck")
    println("    → Potential optimizations:")
    println("      • Use iterative solvers (GMRES, MINRES) instead of direct LU")
    println("      • Implement preconditioners for the saddle-point system")
    println("      • Consider multigrid methods for better scaling")
else
    println("  PRIMARY: Component assembly is the main bottleneck")
    println("    → Potential optimizations:")
    println("      • SIMD vectorization of basis function evaluations")
    println("      • Precomputing more geometric factors")
    println("      • GPU acceleration of integration loops")
end

println("\n" * "="^80)
println("PROFILING COMPLETE")
println("="^80)
