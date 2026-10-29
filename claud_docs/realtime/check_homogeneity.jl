# Checks η ∂sol/∂η + β ∂sol/∂β = -[u; 0] on the force system (u = velocity unknowns q_f and μ).
# usage: julia --project=<repo> check_homogeneity.jl
using smearFEM, LinearAlgebra, SparseArrays
const S = smearFEM
mesh_path = joinpath(dirname(dirname(@__DIR__)), "mesh_files")

F = -9518.61 * ones(10)
for (η, β) in ((100.0, 100.0), (37.0, 2500.0))
    model, scene = def_problem(Cylinder(25.0, 40.0), 6, η, :Hex, 2, 3, :Hex, 1, 1, :Hex, 2,
                               β, F, "force", "constant", 1.0, 0.1; mesh_path=mesh_path)
    cache = S.BasisFunctionCache(model)
    sys = S._squeeze_system(model, scene)
    S._assemble_squeeze_system!(sys, model, cache, 1.0, 0.0); K10 = copy(sys.K)
    S._assemble_squeeze_system!(sys, model, cache, 0.0, 0.0); Ā = K10 - sys.K   # viscous part
    S._assemble_squeeze_system!(sys, model, cache, η, β)
    n = size(sys.K, 1)
    r = zeros(n); r[end] = F[1]
    Fk = lu(sys.K)
    sol = Fk \ r
    sη = Fk \ (-(Ā * sol))
    sβ = Fk \ (-(sys.dKdβ * sol))
    w = -copy(sol); w[sys.nf+1:sys.nf+sys.np] .= 0
    println("η=$η β=$β: ‖η sη + β sβ - w‖/‖w‖ = ", norm(η * sη + β * sβ - w) / norm(w),
            ",  ‖sη‖ = ", norm(sη), ", ‖sβ‖ = ", norm(sβ))
end
