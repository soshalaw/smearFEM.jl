module smearFEM

export gaussian_quadrature, basis_function #fem.jl
export fit_curve, extract_borders, filter_points, rearrange #PostProcess.jl
export read_csv, write_vtk, write_scene, write_csv, write_json, write_contour_data, read_h5 #io.jl
export PlotGrid, plot_mesh, animate_fields, plot_matches, plot_matches_h #plotting.jl
export closest_point, height_sample, match_points
export meshgrid_line, meshgrid_square, meshgrid_cube, inflate_cylinder, meshgrid_ring #Meshes.jl
export def_model, linearElasticity, stokes #models.jl

export GradDescent, update

export simulate, write_sim_data, test, simulate_single_tstep, compare, readData
export simulate_single_tstep_stokes
export setboundaryCond, simulate, set_file, initialize_mesh_test, assemble_system, apply_boundary_conditions, get_cMat, get_volume
export setboundaryCond_dense, assemble_system_dense, apply_boundary_conditions_dense
export assemble_system_A, assemble_system_B, apply_boundary_conditions_stokes, set_boundary_cond_stokes
export set_boundary_cond_flow_cube, set_boundary_cond_flow_cyl

include("fem/models.jl")
include("fem/fem.jl")
include("fem/Meshes.jl")
include("fem/PostProcess.jl")
include("visualization/plotting.jl")
include("visualization/io.jl")
include("optimization/smearOptimize.jl")
include("optimization/GradDescent.jl")

include("examples/squeeze_stokes.jl")
include("examples/fluid_flow_stokes.jl")
include("examples/squeeze_linear_elasticity.jl")
include("examples/run_example.jl")

end # module smearFem
