module smearFEM

using LinearAlgebra

# Performance optimization: enable multi-threaded BLAS at module load time
BLAS.set_num_threads(Threads.nthreads())

abstract type AbstractMeshgrid end

# abstract type model end
export AbstractMeshgrid, MeshgridLine, MeshgridDisk, MeshgridSquare, MeshgridCuboid, MeshgridCylinder # Meshes.jl
export AbstractGeometry, Cylinder, Cuboid, Disk, Square, Segment, ndim # geometries.jl
export Model, LinearElasticity, Stokes # models.jl
export SqueezeFlow # scenarios.jl
export EnvConditions, Conditions # types.jl

export meshgrid_cylinder, meshgrid_cuboid, meshgrid_square, meshgrid_disk, meshgrid_line # Meshes.jl
export reset_mesh!, update_initial_state! # Meshes.jl
export gaussian_quadrature, basis_function, get_quadrature, BasisFunctionCache, get_basis_volume_functions, get_surface_basis_functions # fem.jl
export fit_curve, extract_borders, filter_points, rearrange, add_noise, project_to, back_project, ∇π, get_height, plot_covariance, eval_on_cylinder, get_lagrange_proj, get_lagrange_pts, get_nurbs_2_lagrange_proj, detect_outlier_observations, get_pose # PostProcess.jl
export closest_point, match_points, fit_model # smearOptimize.jl
export reset_model!, update_model! # models.jl

export simulate, write_sim_data, readData, initialize_mesh, write_gt_data, plot_rad_norm_vel_vs_slip, plot_rad_norm_vel_vs_visc
export simulate_single_tstep_stokes, stokes_single_step_force
export assemble_system_A, assemble_system_B, def_problem, set_model

export read_csv, write_vtk, write_scene, write_csv, write_json, write_data, write_2d_data, read_h5, read_json, read_perception_data, get_time_windows, write_stokes_scene, set_file # io.jl
export plot_mesh, animate_fields, plot_matches, plot_matches_h, set_plot, set_subplot # plotting.jl
export plot_noise_covariance, plot_height_vs_slip, plot_field_at_height, arrow0!, get_norm, plot_data, plot_covariance! # analysis_plots.jl

export mat_nan_inf_check, write_time_log, dataframe_2_vec, get_cMat # utils.jl

export get_data_dir, get_mesh_dir, get_scratch_dir, resolve_data_path, resolve_mesh_path, create_output_dir, show_config # config.jl

include("fem/models.jl")
include("fem/fem.jl")
include("fem/meshes/Meshes.jl")
include("fem/PostProcess.jl")
include("fem/types.jl")
include("fem/scenarios.jl")

include("io/plotting.jl")
include("io/io.jl")
include("io/analysis_plots.jl")
include("io/gmsh_utils.jl")

include("optimization/smearOptimize.jl")

include("fem/geometries.jl")
include("examples/stokes_solver.jl")
include("examples/stokes_setup.jl")
include("examples/run_example.jl")

include("config.jl")
include("utils.jl")

# Automatically route package logging to stderr so progress output on stdout
# (sent via the package progress manager) isn't clobbered by Julia log lines.
# This runs when the module is loaded. If you prefer manual control, remove
# or comment out this call.
try
	setup_stderr_logger()
catch e
	@warn "Failed to configure logger: $(e.msg)"
end

# Provide a safe alias to Base.isatty so external code that references
# `smearFEM.isatty` (incorrectly) will still work. If Base.isatty is not
# available in this Julia build, fall back to a safe no-op that returns false.
if isdefined(Base, :isatty)
	const isatty = Base.isatty
else
	const isatty = (io->false)
end

end # module smearFem
