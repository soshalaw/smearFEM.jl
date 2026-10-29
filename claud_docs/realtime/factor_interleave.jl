# Pardiso factorization and solve times when interleaved with the Julia-threaded assembly, as
# inside `simulate`, vs back to back. Run with different KMP_BLOCKTIME / MKL_NUM_THREADS.
# usage: julia -O3 --threads=auto --project=<repo> factor_interleave.jl
using smearFEM, LinearAlgebra, Statistics
const S = smearFEM
mesh_path = joinpath(dirname(dirname(@__DIR__)), "mesh_files")
F = -9518.61 * ones(10)
model, scene = def_problem(Cylinder(25.0, 40.0), 6, 100.0, :Hex, 2, 3, :Hex, 1, 1, :Hex, 2,
                           100.0, F, "force", "constant", 1.0, 0.1; mesh_path=mesh_path)
cache = S.BasisFunctionCache(model)
sys = S._squeeze_system(model, scene)
S._assemble_squeeze_system!(sys, model, cache, 100.0, 100.0)
S._refactor!(sys)
n = size(sys.K, 1); b = rand(n); x = zeros(n)
ms(v) = round(1e3 * median(v), digits=2)
fac() = S._factorize!(sys.lu, sys.K)
sol() = S._solve!(x, sys.lu, b)
asm() = S._assemble_squeeze_system!(sys, model, cache, 100.0, 100.0)
for _ in 1:3; fac(); sol(); asm(); end
tf_b2b = [@elapsed(fac()) for _ in 1:15]
ts_b2b = [@elapsed(sol()) for _ in 1:30]
tf_int = [(asm(); @elapsed(fac())) for _ in 1:15]
ts_int = [(asm(); @elapsed(sol())) for _ in 1:30]
ta_int = [(sol(); @elapsed(asm())) for _ in 1:30]
tf_idle = [(sleep(0.3); @elapsed(fac())) for _ in 1:5]
println("KMP_BLOCKTIME=$(get(ENV, "KMP_BLOCKTIME", "default")) MKL threads=$(S.Pardiso.get_nprocs_mkl()): ",
        "factor b2b $(ms(tf_b2b)) / after assembly $(ms(tf_int)) / after 0.3 s idle $(ms(tf_idle)) ms; ",
        "solve b2b $(ms(ts_b2b)) / after assembly $(ms(ts_int)) ms; assembly after solve $(ms(ta_int)) ms")
