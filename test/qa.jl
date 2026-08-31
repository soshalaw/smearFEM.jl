using smearFEM
using Test
using Aqua

@testset "Code quality (Aqua.jl)" begin
    Aqua.test_all(smearFEM;
    ambiguities=true,
    deps_compat = false,)
end