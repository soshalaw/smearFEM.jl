module smearFEM

using LinearAlgebra

# Performance optimization: enable multi-threaded BLAS at module load time
BLAS.set_num_threads(Threads.nthreads())

"""
    AbstractMeshgrid

Supertype for all mesh containers. Every subtype carries `volume_element_shape::Symbol` and
`basis_order::Int` as two independent fields, never a combined "Q2"-style string.
"""
abstract type AbstractMeshgrid end

# abstract type model end
export AbstractMeshgrid, MeshgridLine, MeshgridDisk, MeshgridSquare, MeshgridCuboid, MeshgridCylinder # Meshes.jl
export Model, LinearElasticity, lame_parameters # models.jl
export EnvConditions, Conditions # types.jl

export meshgrid_cylinder, meshgrid_cuboid, meshgrid_square, meshgrid_disk, meshgrid_line # Meshes.jl
export reset_mesh!, update_initial_state! # Meshes.jl
export gaussian_quadrature, basis_function, get_quadrature, BasisFunctionCache, get_volume_basis_functions, get_surface_basis_functions # fem.jl
export fit_curve, upsample_contour, extract_borders, filter_points, rearrange, add_noise, project_to, back_project, ∇π, get_height, eval_on_cylinder, get_lagrange_proj, get_lagrange_pts, get_nurbs_2_lagrange_proj, detect_outlier_observations, get_pose, project_to_camera_frame # PostProcess.jl
export reset_model!, update_model! # models.jl

export simulate # solver/elastic_squeeze.jl
export volume_update! # solver/volume_update.jl

export read_csv, sorted_csv_files, write_vtk, write_scene, write_csv, write_json, write_data, write_2d_data, read_h5, read_json, read_perception_data, get_time_windows, set_file # io.jl
export plot_mesh, animate_fields, set_plot, set_subplot, set_plot_from_config # plotting.jl

export mat_nan_inf_check, write_time_log, dataframe_2_vec, get_c_mat, get_camera_matrix # utils.jl

export get_data_dir, get_mesh_dir, get_scratch_dir, resolve_data_path, resolve_mesh_path, create_output_dir, show_config # config.jl

export hausdorff_distance, chamfer_distance, hausdorff_distance_kdtree, directed_hausdorff_kdtree, chamfer_distance_kdtree, chamfer_sq_distance_kdtree, closest_point_distance_kdtree, compare_pt_clouds # analysis/pointcloud_metrics.jl

include("fem/models.jl")
include("fem/fem.jl")
include("fem/meshes/Meshes.jl")
include("fem/PostProcess.jl")
include("fem/types.jl")

include("io/plotting.jl")
include("io/io.jl")
include("io/gmsh_utils.jl")

include("solver/sparse_lu.jl")
include("solver/elastic_system.jl")
include("solver/elastic_squeeze.jl")
include("solver/volume_update.jl")

include("config.jl")
include("utils.jl")

include("analysis/pointcloud_metrics.jl")

# Automatically route package logging to stderr so progress output on stdout
# (sent via the package progress manager) isn't clobbered by Julia log lines.
# This runs when the module is loaded. If you prefer manual control, remove
# or comment out this call.
try
	setup_stderr_logger()
catch e
	@warn "Failed to configure logger: $(e.msg)"
end


end # module smearFEM
