# Stress test: concurrent UMFPACK lu! and solves on SEPARATE factorization objects from
# different threads must not crash and must give the same numbers as serial use.
# usage: julia -O3 --threads=auto --project=<repo> umfpack_concurrency.jl [rounds]
using smearFEM, LinearAlgebra, SparseArrays
const S = smearFEM
rounds = length(ARGS) > 0 ? parse(Int, ARGS[1]) : 30
mesh_path = joinpath(dirname(dirname(@__DIR__)), "mesh_files")
F = -9518.61 * ones(10)
model, scene = def_problem(Cylinder(25.0, 40.0), 6, 100.0, :Hex, 2, 3, :Hex, 1, 1, :Hex, 2,
                           100.0, F, "force", "constant", 1.0, 0.1; mesh_path=mesh_path)
cache = S.BasisFunctionCache(model)
sys = S._squeeze_system(model, scene)
S._assemble_squeeze_system!(sys, model, cache, 100.0, 100.0)
K1 = copy(sys.K)
model.mesh_x.NodeList = model.mesh_x.NodeList .* [1.01, 1.01, 0.98]
S._assemble_squeeze_system!(sys, model, cache, 120.0, 80.0)
K2 = copy(sys.K)
n = size(K1, 1); b = rand(n)

F1 = lu(K1); x_ref = F1 \ b                      # serial references
F2 = lu(K2); lu!(F2, K2); F2_ref = F2 \ b
mismatch = 0
for r in 1:rounds
    t = Threads.@spawn begin                     # background: refactor the second object
        for _ in 1:2
            lu!(F2, (r % 2 == 0) ? K2 : K2)
        end
        F2 \ b
    end
    xs = [F1 \ b for _ in 1:40]                  # foreground: solves on the first
    y = fetch(t)
    global mismatch += count(x -> x != x_ref, xs) + (y != F2_ref)
end
println("rounds $rounds, threads $(Threads.nthreads()): mismatches = $mismatch (expect 0)")

# several workers, each with its own objects, doing full lu (symbolic + numeric), lu! and solves
Ks = [K1, K2, K1 + 0.5 * K2, 2.0 * K2]
refs = [lu(K) \ b for K in Ks]
bad = Threads.Atomic{Int}(0)
@sync for (i, K) in enumerate(Ks)
    Threads.@spawn for _ in 1:rounds ÷ 4
        Fi = lu(K)
        lu!(Fi, K)
        (Fi \ b) == refs[i] || Threads.atomic_add!(bad, 1)
    end
end
println("4 concurrent workers × $(rounds ÷ 4) rounds: mismatches = $(bad[]) (expect 0)")
