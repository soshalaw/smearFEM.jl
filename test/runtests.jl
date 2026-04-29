using smearFEM
using Test
using Aqua

@testset "smearFEM.jl" begin

    include("linear_elasticity_test.jl")
    include("stokes_test.jl")
    include("qa.jl")
    include("fem_test.jl")
    # include("IGA_tests/extraction_test.jl")
    # stokes_optimization_test requires LinearSolve and external data files; skip for now
    # include("stokes_optimization_test.jl")
    include("shape_opt_test.jl")
    
end
