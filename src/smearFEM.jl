module smearFEM

export gaussian_quadrature, basis_function, greet_fem, def_model
export fit_curve, extract_borders, greet_pp, inflate_sphere, average_pts, rearrange
export readCSV, write_vtk, write_scene, writeCSV
export PlotGrid, plot_mesh, animate_fields, plot_matches, plot_matches_h
export closest_point, height_sample, match_points
export GradDescent, update
export meshgrid, setboundaryCond, simulate, set_file, initialize_mesh_test, assemble_system, apply_boundary_conditions, get_cMat

export simulate, write_sim_data, test

include("fem.jl")
include("PostProcess.jl")
include("plotting.jl")
include("io.jl")
include("smearOptimize.jl")
include("GradDescent.jl")

include("squeezeFlow.jl")

end # module smearFem
