using smearFEM
using Test
using Aqua

@testset "smearFEM.jl" begin

    include("stokes_test.jl")
    include("qa.jl")
    include("fem_test.jl")
    include("cost_functions_test.jl")
    include("stokes_optimization_test.jl")
    include("tikhonov_test.jl")
    include("mesh_templates_test.jl")
    # Slow (~8.5 min): six full squeeze-flow solves across a beta sweep. It is the only guard against the simplex
    # quadrature weights drifting out of step with the master element `_basis_tet`/`_basis_tri`
    # assume — a mismatch that silently over-stiffens every Tet/Tri solve.
    include("verify_meshes_test.jl")

end
