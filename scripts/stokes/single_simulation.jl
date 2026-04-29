using smearFEM
using Dates
using Plots, Plots.PlotMeasures
using LaTeXStrings
using StatsPlots
using SparseArrays, LinearAlgebra

# Include example helpers into the smearFEM module (so they can access module internals)
Base.include(smearFEM, joinpath(@__DIR__, "..", "..", "src", "examples", "squeeze_stokes.jl"))
Base.include(smearFEM, joinpath(@__DIR__, "..", "..", "src", "examples", "squeeze_linear_elasticity.jl"))


function const_vel(r::Float64, h::Float64, ndim::Int, FunctionClass_u::String, nDof_u::Int, FunctionClass_p::String, nDof_p::Int, FunctionClass_x::String, ne::Int, 
                    camera_matrix::AbstractMatrix{Float64}, obj_pose::AbstractMatrix{Float64})

  control::String = "force"            # "force" or "velocity"
  viscosity_type::String = "constant"  # "constant" or "bulk_viscosity"

  sim_time::Float64 = 2.0            # simulation time in seconds (reduced for validation)
  t_steps::Float64 = 0.1              # number of time steps
  steps::Float64 = sim_time/t_steps
    
  β_gt = 500.0 # penalty parameter for the ground truth
  η_gt =  70.0 # viscosity in (kg/(mm⋅s))

if β_gt <= 1.0
    F_ext = 9.813e3*0.85
elseif β_gt == 10.0
    F_ext = 9.813e3*0.93
elseif β_gt == 50.0
    F_ext = 9.813e3*0.955
elseif β_gt == 100.0
    F_ext = 9.813e3*0.97
elseif β_gt == 500.0
    F_ext = 9.813e3*0.99
elseif β_gt == 1e3
    F_ext = 9.813e3*1.0
elseif β_gt == 5e3
    F_ext = 9.813e3*1.01 # force applied to the cylinder in kg.mm/s^2 (N)
elseif β_gt == 2e3
    F_ext = 9.813e3*0.995 # force applied to the cylinder in kg.mm/s^2 (N)
elseif β_gt == 1e4
    F_ext = 9.813e3*0.85*700 # force applied to the cylinder in kg.mm/s^2 (N)
elseif β_gt == 1e5
    F_ext = 9.813e3*0.85*7e3 # force applied to the cylinder in kg.mm/s^2 (N)
elseif β_gt == 1e10
    F_ext = 9.813e3*0.85*7e8 # force applied to the cylinder in kg.mm/s^2 (N)
end

  F_ext = 9.813e3 # force applied to the cylinder in kg.mm/s^2 (N)
  F::Vector{Float64} = -F_ext*ones(Float64, round(Int, (sim_time/t_steps))) # force applied to the cylinder in N

    model, scene = smearFEM.def_problem(r, h, ne, η_gt, ndim, FunctionClass_u, nDof_u, FunctionClass_p, nDof_p, FunctionClass_x, β_gt, F, control, 
                                                        viscosity_type, sim_time, t_steps)
    # Wrapper used by solver_integration API
    mdl_wrapper = smearFEM.Model(model)
    filepath = string("/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/sim_experiments/single_simulation/FEM/Stokes/$control/$viscosity_type/test_sim")
    
    return model, scene, filepath
end

function main()

  println("\n" * "="^80)
  println("GPU vs CPU SOLVER VALIDATION TEST (Tolerance: 1e-10)")
  println("="^80 * "\n")

  # Create small test problem for quick validation
  r::Float64 = 25.0
  h::Float64 = 40.0
  ne::Int = 3  # Very small mesh for fast validation
  ndim::Int = 3
  
  println("[INFO] Creating test meshes...")
  mesh_u = smearFEM.meshgrid_cylinder(r, h, ne, FunctionClass="Q2")
  mesh_p = smearFEM.meshgrid_cylinder(r, h, ne, FunctionClass="Q1")
  
  n_nodes = size(mesh_u.NodeList, 2)
  n_elements = size(mesh_u.IEN, 2)
  
  println("[INFO] Mesh info:")
  println("       Nodes: $n_nodes, Elements: $n_elements")
  
  # Create random system (synthetic test)
  n_dof = n_nodes * (3 + 1)  # 3 velocity DOFs + 1 pressure DOF
  println("[INFO] System DOF: $n_dof")
  
  # Create a small SPD matrix for testing
  A_cpu = sprand(n_dof, n_dof, 0.1) + I(n_dof) * n_dof  # Make it SPD
  A_cpu = A_cpu + A_cpu'  # Symmetrize
  b_cpu = rand(n_dof)
  
  config_cpu = smearFEM.cpu_fallback_config()
  config_gpu = smearFEM.realtime_config()
  
  println("\n[1/2] Solving with CPU (GMRES)...")
  t_cpu_start = Dates.now()
  x_cpu = smearFEM.solve_system(A_cpu, b_cpu, config_cpu)
  t_cpu_end = Dates.now()
  elapsed_cpu = t_cpu_end - t_cpu_start
  println("✓ CPU solve completed in: ", elapsed_cpu)
  
  println("\n[2/2] Solving with GPU (if available, else CPU)...")
  t_gpu_start = Dates.now()
  x_gpu = smearFEM.solve_system(A_cpu, b_cpu, config_gpu)
  t_gpu_end = Dates.now()
  elapsed_gpu = t_gpu_end - t_gpu_start
  println("✓ GPU solve completed in: ", elapsed_gpu)
  
  # Validation
  tolerance = 1e-10
  
  println("\n" * "="^80)
  println("VALIDATION RESULTS")
  println("="^80)
  
  diff = maximum(abs.(x_cpu .- x_gpu))
  rel_error = diff / (maximum(abs.(x_cpu)) + 1e-15)
  abs_norm_cpu = norm(A_cpu * x_cpu - b_cpu)
  abs_norm_gpu = norm(A_cpu * x_gpu - b_cpu)
  
  println("\nSolution Comparison:")
  println("  Max Difference:    $diff")
  println("  Relative Error:    $rel_error")
  println("  CPU Residual:      $abs_norm_cpu")
  println("  GPU Residual:      $abs_norm_gpu")
  
  status = diff <= tolerance ? "✓ PASS" : "✗ FAIL"
  println("\n$status - GPU/CPU solutions match to tolerance 1e-10")
  
  println("\n" * "="^80)
  println("PERFORMANCE COMPARISON")
  println("="^80)
  elapsed_cpu_ms = elapsed_cpu.value
  elapsed_gpu_ms = elapsed_gpu.value
  speedup = elapsed_cpu_ms / elapsed_gpu_ms
  println("CPU Time:    $(elapsed_cpu_ms) ms")
  println("GPU Time:    $(elapsed_gpu_ms) ms")
  println("Speedup:     $(round(speedup; digits=2))x")
  println("GPU Status:  $(smearFEM.has_gpu() ? "AVAILABLE" : "NOT AVAILABLE (using CPU)")")
  
  println("\n" * "="^80)
  if diff <= tolerance
    println("✓ VALIDATION PASSED - GPU and CPU solver results match to 1e-10")
  else
    println("✗ VALIDATION FAILED - Differences exceed 1e-10")
  end
  println("="^80 * "\n")


end

main()

