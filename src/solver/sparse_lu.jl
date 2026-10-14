using LinearAlgebra
using SparseArrays
using Pardiso

# The sparse direct solver behind `_lagged_solve!`, for the symmetric saddle-point matrix of a
# `SqueezeSystem`. MKL Pardiso where MKL is available (x86-64): on the benchmark system
# (n = 5653) it factorizes in ~9 ms on 8 threads and solves in ~1.3 ms, where UMFPACK — the
# fallback everywhere else — takes ~150 ms and ~3.6 ms. Both reach relative residuals ~1e-14.

const _UmfpackLU = SparseArrays.UMFPACK.UmfpackLU{Float64,Int}

mutable struct _Umfpack
    F::Union{Nothing,_UmfpackLU}
    F2::Union{Nothing,_UmfpackLU}       # copy(F): own workspace, so two solves can run at once
end

mutable struct _Pardiso
    ps::MKLPardisoSolver
    A::SparseMatrixCSC{Float64,Int32}   # lower triangle of K plus every diagonal entry, which
                                        # Pardiso requires; Int32 is Pardiso.jl's index type here,
                                        # anything else is converted on every call
    from::Vector{Int}                   # K.nzval index of each A.nzval entry, 0 = added diagonal
    analysed::Bool
    x::Vector{Float64}                  # every solve goes through these two fixed buffers: MKL's
    b::Vector{Float64}                  # results depend on the alignment of its arguments
end

_sparse_lu(K::SparseMatrixCSC{Float64,Int}) = Pardiso.mkl_is_available() ? _Pardiso(K) : _Umfpack(nothing, nothing)

function _Pardiso(K::SparseMatrixCSC{Float64,Int})
    n = size(K, 1)
    colptr = ones(Int32, n + 1); rowval = Int32[]; from = Int[]
    for c in 1:n
        rows = [(K.rowval[k], k) for k in K.colptr[c]:(K.colptr[c+1] - 1) if K.rowval[k] >= c]
        any(r == c for (r, _) in rows) || push!(rows, (c, 0))
        sort!(rows)
        append!(rowval, first.(rows)); append!(from, last.(rows))
        colptr[c+1] = length(rowval) + 1
    end
    ps = MKLPardisoSolver()
    set_matrixtype!(ps, Pardiso.REAL_SYM_INDEF)
    pardisoinit(ps)
    fix_iparm!(ps, :N)
    set_iparm!(ps, 2, 2)      # METIS ordering
    set_iparm!(ps, 8, 0)      # no refinement inside Pardiso: `_lagged_solve!` refines
    set_iparm!(ps, 11, 1)     # scaling and weighted matching, recommended for saddle points
    set_iparm!(ps, 13, 1)
    set_iparm!(ps, 34, Pardiso.get_nprocs_mkl())   # fixed thread count: bit-identical repeats
    f = _Pardiso(ps, SparseMatrixCSC(n, n, colptr, rowval, zeros(length(rowval))), from, false, zeros(n), zeros(n))
    # Pardiso's memory lives outside Julia's heap
    finalizer(f) do f
        f.analysed || return
        set_phase!(f.ps, Pardiso.RELEASE_ALL)
        try pardiso(f.ps) catch end
    end
    return f
end

function _factorize!(f::_Pardiso, K::SparseMatrixCSC{Float64,Int})
    nz, Kz, from = f.A.nzval, K.nzval, f.from
    @inbounds for i in eachindex(from)
        k = from[i]
        nz[i] = k == 0 ? 0.0 : Kz[k]
    end
    if !f.analysed           # the pattern never changes: one symbolic analysis per system
        set_phase!(f.ps, Pardiso.ANALYSIS)
        pardiso(f.ps, f.x, f.A, f.b)
        f.analysed = true
    end
    set_phase!(f.ps, Pardiso.NUM_FACT)
    pardiso(f.ps, f.x, f.A, f.b)
    return f
end

function _factorize!(f::_Umfpack, K::SparseMatrixCSC{Float64,Int})
    lock(SPARSE_LU_LOCK) do
        f.F === nothing ? (f.F = lu(K)) : lu!(f.F, K)
    end
    f.F2 = copy(f.F)         # lu! gives F a new numeric factorization; a copy keeps the old one
    return f
end

function _solve!(x::AbstractVector{Float64}, f::_Pardiso, b::AbstractVector{Float64})
    copyto!(f.b, b)
    set_phase!(f.ps, Pardiso.SOLVE_ITERATIVE_REFINE)
    pardiso(f.ps, f.x, f.A, f.b)
    return copyto!(x, f.x)
end

_solve!(x::AbstractVector{Float64}, f::_Umfpack, b::AbstractVector{Float64}) =
    lock(() -> ldiv!(x, f.F, b), SPARSE_LU_LOCK)

# x = F \ b and y = F \ c. UMFPACK runs the two at once, each on its own workspace (1.7× the
# throughput of one after the other); Pardiso already spreads one solve over every core.
function _solve_pair!(x, y, f::_Pardiso, b, c)
    _solve!(x, f, b)
    _solve!(y, f, c)
    return nothing
end

function _solve_pair!(x, y, f::_Umfpack, b, c)
    lock(SPARSE_LU_LOCK) do
        t = Threads.@spawn ldiv!(y, f.F2, c)
        ldiv!(x, f.F, b)
        wait(t)
    end
    return nothing
end

# Pardiso refactorizes in about three refinement rounds' time, so once a round shows the
# factorization has gone stale — contracting the residual by less than 20×, which takes ~25 steps
# of mesh deformation on the benchmark — it is replaced at the start of the next step, before any
# round is spent on it. UMFPACK takes ~40 rounds' worth, so it never does that; only a solve that
# needs more than five rounds refactorizes, there and then (Pardiso: if refining stops working).
_stale(::_Pardiso, ρ::Float64) = ρ > 0.05
_stale(::_Umfpack, ρ::Float64) = false
_refactor_now(::_Pardiso, ρ::Float64, rounds::Int) = ρ > 0.5
_refactor_now(::_Umfpack, ρ::Float64, rounds::Int) = rounds > 5
