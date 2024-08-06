module smearFEM

export gaussian_quadrature, basis_function, assemble_system, greet_fem
export greet_pp

include("fem.jl")
include("PostProcess.jl")

end
