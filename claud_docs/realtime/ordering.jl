# UMFPACK strategy/ordering vs fill, factor and solve time on the benchmark force system.
# usage: julia -O3 --project=<repo> ordering.jl
using smearFEM, LinearAlgebra, SparseArrays
const S = smearFEM
const U = SparseArrays.UMFPACK
mesh_path = joinpath(dirname(dirname(@__DIR__)), "mesh_files")

F = -9518.61 * ones(10)
model, scene = def_problem(Cylinder(25.0, 40.0), 6, 100.0, :Hex, 2, 3, :Hex, 1, 1, :Hex, 2,
                           100.0, F, "force", "constant", 1.0, 0.1; mesh_path=mesh_path)
cache = S.BasisFunctionCache(model)
sys = S._squeeze_system(model, scene)
S._assemble_squeeze_system!(sys, model, cache, 100.0, 100.0)
K = sys.K
x = rand(size(K, 1))
ms(t) = round(1e3 * t, digits=1)

for (sname, strat) in (("auto", 0), ("unsym", 1), ("sym", 3)),
    (oname, ord) in (("amd", 1), ("metis", 3), ("best", 4))
    ctrl = U.get_umfpack_control(Float64, Int)
    ctrl[SparseArrays.LibSuiteSparse.UMFPACK_STRATEGY + 1] = strat
    ctrl[U.JL_UMFPACK_ORDERING] = ord
    try
        t_first = @elapsed Fk = lu(K; control=ctrl)
        t_re = minimum(@elapsed(lu!(Fk, K)) for _ in 1:3)
        t_s = minimum(@elapsed(Fk \ x) for _ in 1:10)
        res = norm(K * (Fk \ x) - x) / norm(x)
        nlu = nnz(Fk.L) + nnz(Fk.U)
        println(rpad("$sname/$oname", 12), " nnz(L+U) = ", lpad(nlu, 8), "  first lu ", lpad(ms(t_first), 7),
                " ms  lu! ", lpad(ms(t_re), 7), " ms  solve ", lpad(ms(t_s), 5), " ms  res ", res)
    catch err
        println(rpad("$sname/$oname", 12), " failed: ", err)
    end
end
