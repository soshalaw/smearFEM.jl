using smearFEM
using Test
using Aqua

# "fast", "slow" or "all" (default). CI runs the groups as separate jobs to stay inside its
# per-job time limit; coverage roughly doubles the slow group's runtime.
const TEST_GROUP = get(ENV, "SMEAR_TEST_GROUP", "all")

@testset "smearFEM.jl" begin

    if TEST_GROUP in ("all", "fast")
        include("stokes_test.jl")
        include("qa.jl")
        include("fem_test.jl")
        include("cost_functions_test.jl")
        include("stokes_optimization_test.jl")
        include("tikhonov_test.jl")
        include("mesh_templates_test.jl")
    end
    if TEST_GROUP in ("all", "slow")
        # Slow (~8.5 min): six full squeeze-flow solves across a beta sweep. It is the only guard against the simplex
        # quadrature weights drifting out of step with the master element `_basis_tet`/`_basis_tri`
        # assume — a mismatch that silently over-stiffens every Tet/Tri solve.
        include("verify_meshes_test.jl")
    end

end
