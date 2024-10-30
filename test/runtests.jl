using smearFEM
using Test
using Aqua

@testset "smearFEM.jl" begin

    include("linear_elasticity_test.jl")
    include("qa.jl")
    include("fem_test.jl")
    
end
