using smearFEM
using Test
using Aqua

@testset "smearFEM.jl" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(smearFEM)
    end
    # Write your tests here.
end
