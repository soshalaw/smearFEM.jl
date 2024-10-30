module smearFEM

export gaussian_quadrature, basis_function
export fit_curve, extract_borders, filter_points, rearrange
export readCSV, write_vtk, write_scene, writeCSV
export PlotGrid, plot_mesh, animate_fields, plot_matches, plot_matches_h
export closest_point, height_sample, match_points
export GradDescent, update
export setboundaryCond, simulate, set_file, initialize_mesh_test, assemble_system, apply_boundary_conditions, get_cMat
export meshgrid_line, meshgrid_square, meshgrid_cube, inflate_cylinder, meshgrid_ring
export simulate, write_sim_data, test, simulate_single_tstep
export assemble_system_A, assemble_system_B, apply_boundary_conditions_stokes, set_boundary_cond_stokes
export def_model, linearElasticity, stokes

include("models.jl")
include("fem.jl")
include("Meshes.jl")
include("PostProcess.jl")
include("plotting.jl")
include("io.jl")
include("smearOptimize.jl")
include("GradDescent.jl")
include("squeeze_stokes.jl")
include("squeeze_linear_elasticity.jl")

end # module smearFem
