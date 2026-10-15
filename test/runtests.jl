using smearFEM
using Test
using Aqua

@testset "smearFEM.jl" begin
    include("qa.jl")
    include("fem_test.jl")
    include("mesh_templates_test.jl")
    # also the guard against the simplex quadrature weights drifting out of step with the
    # master element `_basis_tet`/`_basis_tri` assume, which silently stiffens every Tet solve
    include("linear_elasticity_test.jl")
end
