using smearFEM
using Test
using Aqua

@testset "smearFEM.jl" begin

    include("stokes_test.jl")
    include("qa.jl")
    include("fem_test.jl")
    include("cost_functions_test.jl")
    include("stokes_optimization_test.jl")
    include("mesh_templates_test.jl")
    
end
