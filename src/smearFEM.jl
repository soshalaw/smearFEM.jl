module smearFEM

abstract type AbstractMeshgrid end

# abstract type model end
export AbstractMeshgrid, Meshgrid, Meshgrid1D, Meshgrid2D, Meshgrid3D, MeshgridCylinder # Meshes.jl
export Model, LinearElasticity, Stokes # models.jl
export SqueezeFlow # scenarios.jl
export EnvConditions, Conditions # types.jl

export meshgrid_line, meshgrid_square, meshgrid_cube, inflate_cylinder, meshgrid_ring, meshgrid_cylinder # Meshes.jl
export gaussian_quadrature, basis_function # fem.jl
export fit_curve, extract_borders, filter_points, rearrange, add_noise, project_to, back_project, ∇π, get_height, plot_covariance, eval_on_cylinder, get_lagrange_proj, get_lagrange_pts, get_nurbs_2_lagrange_proj, detect_outlier_observations, get_pose # PostProcess.jl
export closest_point, height_sample, match_points, fit_model # smearOptimize.jl
export reset_model!, update_model! # models.jl

export simulate, write_sim_data, test, simulate_single_tstep, compare, readData, initialize_mesh, write_gt_data
export simulate_single_tstep_stokes, simulate_stokes, test_stokes, write_sim_data_stokes
export set_boundary_conditions, simulate, set_file, assemble_system, set_slip_conditions, get_cMat, get_volume
export set_boundary_conditions_dense, assemble_system_dense, set_slip_conditions_dense
export assemble_system_A, assemble_system_B, apply_boundary_conditions_stokes, set_boundary_cond_stokes, compare_stokes, def_problem, reset_mesh
export set_boundary_cond_flow_cube, set_boundary_cond_flow_cyl

export read_csv, write_vtk, write_scene, write_csv, write_json, write_data, read_h5, read_json, read_perception_data, read_hdf5_generic # io.jl
export PlotGrid, plot_mesh, animate_fields, plot_matches, plot_matches_h, set_plot, set_subplot # plotting.jl
export plot_noise_covariance, plot_height_vs_slip, plot_field_at_height, arrow0!, get_norm, plot_data, plot_covariance! # analysis_plots.jl

export mat_nan_inf_check, write_time_log

include("fem/models.jl")
include("fem/fem.jl")
include("fem/Meshes.jl")
include("fem/PostProcess.jl")
include("fem/types.jl")
include("fem/scenarios.jl")

include("visualization/plotting.jl")
include("visualization/io.jl")
include("visualization/analysis_plots.jl")

include("optimization/smearOptimize.jl")

include("examples/squeeze_stokes.jl")
include("examples/fluid_flow_stokes.jl")
include("examples/squeeze_linear_elasticity.jl")
include("examples/run_example.jl")

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

end # module smearFem
