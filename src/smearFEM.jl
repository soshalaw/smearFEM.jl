module smearFEM

using LinearAlgebra

# Performance optimization: enable multi-threaded BLAS at module load time
BLAS.set_num_threads(Threads.nthreads())

abstract type AbstractMeshgrid end

# abstract type model end
export AbstractMeshgrid, Meshgrid, Meshgrid1D, Meshgrid2D, Meshgrid3D, MeshgridCylinder # Meshes.jl
export Model, LinearElasticity, Stokes # models.jl
export SqueezeFlow # scenarios.jl
export EnvConditions, Conditions # types.jl

export meshgrid_line, meshgrid_square, meshgrid_cube, inflate_cylinder, meshgrid_ring, meshgrid_cylinder # Meshes.jl
export gaussian_quadrature, basis_function, BasisFunctionCache, get_basis_volume_functions, get_surface_basis_functions # fem.jl
export fit_curve, extract_borders, filter_points, rearrange, add_noise, project_to, back_project, ∇π, get_height, plot_covariance, eval_on_cylinder, get_lagrange_proj, get_lagrange_pts, get_nurbs_2_lagrange_proj, detect_outlier_observations, get_pose # PostProcess.jl
export closest_point, height_sample, match_points, fit_model # smearOptimize.jl
export reset_model!, update_model! # models.jl

export simulate, write_sim_data, test, simulate_single_tstep, compare, readData, initialize_mesh, write_gt_data, plot_rad_norm_vel_vs_slip, plot_rad_norm_vel_vs_visc
export simulate_single_tstep_stokes, simulate_stokes, test_stokes, write_sim_data_stokes, stokes_single_step_force
export set_boundary_conditions, simulate, set_file, assemble_system, set_slip_conditions, get_cMat, get_volume
export set_boundary_conditions_dense, assemble_system_dense, set_slip_conditions_dense
export assemble_system_A, assemble_system_B, apply_boundary_conditions_stokes, set_boundary_cond_stokes, compare_stokes, def_problem, reset_mesh
export set_boundary_cond_flow_cube, set_boundary_cond_flow_cyl

export read_csv, write_vtk, write_scene, write_csv, write_json, write_data, read_h5, read_json, read_perception_data, read_hdf5_generic, get_time_windows, write_stokes_scene # io.jl
export PlotGrid, plot_mesh, animate_fields, plot_matches, plot_matches_h, set_plot, set_subplot # plotting.jl
export plot_noise_covariance, plot_height_vs_slip, plot_field_at_height, arrow0!, get_norm, plot_data, plot_covariance! # analysis_plots.jl
export get_mesh_data # gmsh_utils.jl

export mat_nan_inf_check, write_time_log, dataframe_2_vec

# GPU acceleration utilities (Phase 1)
export has_gpu, alloc_gpu_storage, keep_on_gpu_strategy, setup_gpu_kernel_workspace
export query_gpu_memory, log_gpu_status, reset_gpu_memory
export GPUContext, reset_gpu_context, is_gpu_available, need_precond_recompute

# GPU assembly and solver (Phase 2)
export assemble_stokes_gpu!, assemble_elasticity_gpu!
export solve_stokes_gpu!, solve_elasticity_gpu!
export SolverConfig, realtime_config, cpu_fallback_config
export prepare_basis_cache

# GPU acceleration utilities
export has_gpu, alloc_gpu_storage, keep_on_gpu_strategy, setup_gpu_kernel_workspace, query_gpu_memory, log_gpu_status, reset_gpu_memory
export GPUContext

include("fem/models.jl")
include("fem/fem.jl")
include("fem/Meshes.jl")
include("fem/PostProcess.jl")
include("fem/types.jl")
include("fem/scenarios.jl")

include("io/plotting.jl")
include("io/io.jl")
include("io/analysis_plots.jl")
include("io/gmsh_utils.jl")

include("optimization/smearOptimize.jl")

# GPU acceleration utilities (must be before examples since they may use GPU features)
include("utils/gpu_utils.jl")
include("utils/gpu_memory.jl")

# GPU assembly and solver kernels (Phase 2 GPU acceleration)
include("fem/gpu_assembly.jl")
include("fem/gpu_solver.jl")

# Example files are not auto-included to avoid importing unnecessary dependencies like LinearSolve
# Users can include them explicitly when needed
# include("examples/squeeze_stokes.jl")
# include("examples/fluid_flow_stokes.jl")
# include("examples/squeeze_linear_elasticity.jl")
# include("examples/run_example.jl")

include("utils.jl")

# Automatically route package logging to stderr so progress output on stdout
# (sent via the package progress manager) isn't clobbered by Julia log lines.
# This runs when the module is loaded. If you prefer manual control, remove
# or comment out this call.
try
	setup_stderr_logger()
catch e
	# swallow any error to avoid breaking module load; it's safe to ignore
end

# Provide a safe alias to Base.isatty so external code that references
# `smearFEM.isatty` (incorrectly) will still work. If Base.isatty is not
# available in this Julia build, fall back to a safe no-op that returns false.
if isdefined(Base, :isatty)
	const isatty = Base.isatty
else
	const isatty = (io->false)
end

# GPU status reporting at module load
try
	log_gpu_status()
catch e
	# Silently ignore if GPU status reporting fails
end

end # module smearFem
