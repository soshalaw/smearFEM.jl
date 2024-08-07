module smearFEM

export assemble_system, gaussian_quadrature, basis_function, greet_fem
export write_scene, animate_fields, fit_curve, extract_borders, greet_pp, inflate_sphere

include("fem.jl")
include("PostProcess.jl")

end # module smearFem

