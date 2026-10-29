# Throughput of two UMFPACK triangular solves run concurrently (one numeric factorization,
# two workspaces via copy(F)) vs one after the other.
# usage: julia -O3 --threads=auto --project=<repo> concurrent_solves.jl
using smearFEM, LinearAlgebra, SparseArrays
const S = smearFEM
mesh_path = joinpath(dirname(dirname(@__DIR__)), "mesh_files")
F = -9518.61 * ones(10)
model, scene = def_problem(Cylinder(25.0, 40.0), 6, 100.0, :Hex, 2, 3, :Hex, 1, 1, :Hex, 2,
                           100.0, F, "force", "constant", 1.0, 0.1; mesh_path=mesh_path)
cache = S.BasisFunctionCache(model)
sys = S._squeeze_system(model, scene)
S._assemble_squeeze_system!(sys, model, cache, 100.0, 100.0)
F1 = lu(sys.K); F2 = copy(F1)
n = size(sys.K, 1); b1, b2 = rand(n), rand(n); x1, x2 = zeros(n), zeros(n)
seq() = (ldiv!(x1, F1, b1); ldiv!(x2, F1, b2))
function par()
    t = Threads.@spawn ldiv!(x2, F2, b2)
    ldiv!(x1, F1, b1)
    wait(t)
end
for f in (seq, par); f(); end
ts = minimum(@elapsed(seq()) for _ in 1:30)
tp = minimum(@elapsed(par()) for _ in 1:30)
println("two solves: sequential $(round(1e3ts, digits=2)) ms, concurrent $(round(1e3tp, digits=2)) ms")
par(); y1 = copy(x1); y2 = copy(x2); seq()
println("concurrent results identical to sequential: ", y1 == x1 && y2 == x2)
