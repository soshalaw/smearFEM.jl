module smearFEM

abstract type AbstractMeshgrid end

# abstract type model end
export AbstractMeshgrid, Meshgrid, Meshgrid1D, Meshgrid2D, Meshgrid3D, MeshgridCylinder # Meshes.jl
export Model, LinearElasticity, Stokes # models.jl
export SqueezeFlow # scenarios.jl
export EnvConditions, Conditions # types.jl

export meshgrid_line, meshgrid_square, meshgrid_cube, inflate_cylinder, meshgrid_ring, meshgrid_cylinder # Meshes.jl
export gaussian_quadrature, basis_function # fem.jl
export fit_curve, extract_borders, filter_points, rearrange, add_noise, project_to, back_project, ∇π, get_height, plot_covariance # PostProcess.jl
export read_csv, write_vtk, write_scene, write_csv, write_json, write_contour_data, read_h5, read_json # io.jl
export PlotGrid, plot_mesh, animate_fields, plot_matches, plot_matches_h # plotting.jl
export closest_point, height_sample, match_points, fit_model # smearOptimize.jl
export reset_model!, update_model! # models.jl
export reset_model!, reset_model! # Meshes.jl
export simulate, write_sim_data, test, simulate_single_tstep, compare, readData
export simulate_single_tstep_stokes, simulate_stokes, test_stokes, write_sim_data_stokes
export set_boundary_conditions, simulate, set_file, initialize_mesh_test, assemble_system, set_slip_conditions, get_cMat, get_volume
export set_boundary_conditions_dense, assemble_system_dense, set_slip_conditions_dense
export assemble_system_A, assemble_system_B, apply_boundary_conditions_stokes, set_boundary_cond_stokes, compare_stokes, def_problem, reset_mesh
export set_boundary_cond_flow_cube, set_boundary_cond_flow_cyl

export mat_nan_inf_check

include("fem/models.jl")
include("fem/fem.jl")
include("fem/Meshes.jl")
include("fem/PostProcess.jl")
include("fem/types.jl")
include("fem/scenarios.jl")

include("visualization/plotting.jl")
include("visualization/io.jl")

include("optimization/smearOptimize.jl")

include("examples/squeeze_stokes.jl")
include("examples/fluid_flow_stokes.jl")
include("examples/squeeze_linear_elasticity.jl")
include("examples/run_example.jl")

include("utils.jl")

end # module smearFem
