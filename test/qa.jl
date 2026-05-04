using smearFEM
using Test
using Aqua

@testset "Code quality (Aqua.jl)" begin
    Aqua.test_all(smearFEM;
    ambiguities=false,
    deps_compat = false,)
end