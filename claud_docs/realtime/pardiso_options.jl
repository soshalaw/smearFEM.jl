# Pardiso numeric factorization / solve time and accuracy for a few iparm settings on the
# benchmark system, interleaved with assembly as inside `simulate`.
# usage: julia -O3 --threads=auto --project=<repo> pardiso_options.jl
using smearFEM, LinearAlgebra, SparseArrays, Statistics, Pardiso
const S = smearFEM
mesh_path = joinpath(dirname(dirname(@__DIR__)), "mesh_files")
F = -9518.61 * ones(10)
model, scene = def_problem(Cylinder(25.0, 40.0), 6, 100.0, :Hex, 2, 3, :Hex, 1, 1, :Hex, 2,
                           100.0, F, "force", "constant", 1.0, 0.1; mesh_path=mesh_path)
cache = S.BasisFunctionCache(model)
sys = S._squeeze_system(model, scene)
asm() = S._assemble_squeeze_system!(sys, model, cache, 100.0, 100.0)
asm()
K = sys.K; n = size(K, 1); b = rand(n)
ms(v) = round(1e3 * median(v), digits=2)
for (name, opts, T) in (("ours (2=2, 8=0, 11=1, 13=1, 34=n)", Dict(), Float64),
                         ("+ 24=1 two-level factorization", Dict(24 => 1), Float64),
                         ("13=0 no weighted matching", Dict(13 => 0), Float64),
                         ("11=0, 13=0 no scaling/matching", Dict(11 => 0, 13 => 0), Float64),
                         ("28=1 single precision", Dict(28 => 1), Float32),
                         ("28=1 + 24=1", Dict(28 => 1, 24 => 1), Float32))
    f = S._Pardiso(K)
    for (k, v) in opts; set_iparm!(f.ps, k, v); end
    A = SparseMatrixCSC{T,Int32}(f.A.m, f.A.n, f.A.colptr, f.A.rowval, T.(f.A.nzval))
    x, bb = zeros(T, n), zeros(T, n)
    fill_A!() = (for (i, k) in enumerate(f.from); A.nzval[i] = k == 0 ? 0 : T(K.nzval[k]); end)
    fill_A!()
    set_phase!(f.ps, Pardiso.ANALYSIS); pardiso(f.ps, x, A, bb)
    fac() = (set_phase!(f.ps, Pardiso.NUM_FACT); pardiso(f.ps, x, A, bb))
    sol() = (copyto!(bb, b); set_phase!(f.ps, Pardiso.SOLVE_ITERATIVE_REFINE); pardiso(f.ps, x, A, bb))
    fac(); sol()
    tf = [(asm(); fill_A!(); @elapsed(fac())) for _ in 1:10]
    ts = [(asm(); @elapsed(sol())) for _ in 1:10]
    res = norm(K * Float64.(x) - b) / norm(b)
    println(rpad(name, 36), "numeric $(ms(tf)) ms, solve after assembly $(ms(ts)) ms, rel. residual $(round(res, sigdigits=2)), perturbed pivots $(get_iparm(f.ps, 14))")
end
