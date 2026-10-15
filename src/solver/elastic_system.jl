using LinearAlgebra
using SparseArrays

# Fixed-pattern assembly of the linear elastic stiffness
#
#     K = μ A + λ L + β b,   A = ∫ 2ε(v):ε(u),   L = ∫ div v div u,   b = tangential slip mass,
#
# on a chosen set of unknowns: `dof[i, a]` is the unknown of displacement component i of node a,
# 0 when the component is eliminated; several nodal components may share one unknown (a rigid
# plate). The mesh topology and the unknowns are fixed per system, so the sparsity pattern is
# too: each K.nzval entry sums its element contributions through precomputed gather lists, and
# the factorization's symbolic analysis is reused. L and b are kept on their own as ∂K/∂λ and
# ∂K/∂β (b on its own, much sparser pattern); ∂K/∂μ is never needed, since K is homogeneous of
# degree one in (λ, μ, β).

"""
    ElasticSystem

Linear elastic stiffness matrix on a fixed sparsity pattern, with the factorization the solves
reuse. Build it with `_elastic_system`; refill it at the current geometry with `_assemble!`.
"""
mutable struct ElasticSystem
    K::SparseMatrixCSC{Float64,Int}
    nz_L::Vector{Float64}                 # ∂K/∂λ on K's pattern, when assembled with sensitivities
    dKdβ::SparseMatrixCSC{Float64,Int}    # ∂K/∂β = b
    slip_nz::Vector{Int}                  # K.nzval index of each dKdβ.nzval entry
    ptr::Vector{Int32}                    # nzval entry p sums element-matrix entries
    src::Vector{Int32}                    # src[ptr[p]:ptr[p+1]-1] (Int32: gathering is bandwidth bound)
    slip_faces::Vector{Matrix{Int}}       # IEN of each face group with Navier slip
    slip_pos::Vector{Matrix{Int}}         # (2nS², faces) tangential slip entry -> dKdβ.nzval index
    Ke::Matrix{Float64}                   # (nUd², ne) element matrices μA + λL, filled in parallel
    Ke_L::Matrix{Float64}                 # (nUd², ne) element matrices L
    dof::Matrix{Int}                      # (3, nNodes) unknown of each nodal component, 0 = eliminated
    cache::BasisFunctionCache
    lu::Union{Nothing,_Pardiso,_Umfpack}  # factorization the solves use (sparse_lu.jl)
    stale::Bool                           # refactorize before the next solve
    buf::Matrix{Float64}                  # (n, 2) work vectors of `_refined_solve!`
    scratch::Vector{Tuple{Matrix{Float64},Array{Float64,3},Array{Float64,3}}}  # per assembly chunk
end

# Invert an element-entry -> nzval position map `pos` (0 = eliminated) into per-entry lists:
# entry p gathers the element-matrix values at linear indices src[ptr[p]:ptr[p+1]-1]. They come
# out in element order, so a parallel per-entry sum adds in the same order as a serial scatter.
function _gather_lists(pos::Matrix{Int}, nnz::Int)
    ptr = zeros(Int32, nnz + 1)
    for p in pos
        p != 0 && (ptr[p + 1] += 1)
    end
    ptr[1] = 1
    cumsum!(ptr, ptr)
    next = ptr[1:nnz]
    src = zeros(Int32, ptr[end] - 1)
    for (k, p) in enumerate(pos)
        p == 0 && continue
        src[next[p]] = k
        next[p] += 1
    end
    return ptr, src
end

# Position of (r, c) in a CSC matrix's nzval; 0 if either index is eliminated.
function _nz_index(K::SparseMatrixCSC, r::Int, c::Int)
    (r == 0 || c == 0) && return 0
    rng = K.colptr[c]:(K.colptr[c+1]-1)
    k = searchsortedfirst(view(K.rowval, rng), r)
    return first(rng) + k - 1
end

"""
    _elastic_system(mesh::AbstractMeshgrid, dof::Matrix{Int}; slip_faces=Matrix{Int}[]) -> ElasticSystem

Build the sparsity pattern and element-to-nzval maps of the stiffness on the unknowns `dof`,
with Navier slip on the tangential (x, y) components of every face group in `slip_faces`.
Values are filled by `_assemble!`.
"""
function _elastic_system(mesh::AbstractMeshgrid, dof::Matrix{Int}; slip_faces::Vector{Matrix{Int}}=Matrix{Int}[])
    IEN::Matrix{Int} = mesh.IEN
    nN, ne = size(IEN)
    nUd = 3nN
    n = maximum(dof)

    ru = reshape([dof[i, IEN[a, e]] for i in 1:3, a in 1:nN, e in 1:ne], nUd, ne)
    Is = Int[]; Js = Int[]
    sizehint!(Is, nUd^2 * ne); sizehint!(Js, nUd^2 * ne)
    for e in 1:ne, c in view(ru, :, e), r in view(ru, :, e)
        r != 0 && c != 0 && (push!(Is, r); push!(Js, c))
    end
    # every diagonal: Pardiso needs them, and a prescribed unknown becomes an identity row
    append!(Is, 1:n); append!(Js, 1:n)
    K = sparse(Is, Js, ones(length(Is)), n, n)
    fill!(K.nzval, 0.0)

    pos = zeros(Int, nUd^2, ne)
    for e in 1:ne, j in 1:nUd, i in 1:nUd
        pos[(j-1)*nUd + i, e] = _nz_index(K, ru[i, e], ru[j, e])
    end

    # slip acts on the two tangential components of a horizontal face
    function slip_positions(IEN_s)
        nS = size(IEN_s, 1)
        p = zeros(Int, 2nS^2, size(IEN_s, 2))
        for f in axes(IEN_s, 2), i in 1:2, b in 1:nS, a in 1:nS
            p[((i-1)*nS + (b-1))*nS + a, f] = _nz_index(K, dof[i, IEN_s[a, f]], dof[i, IEN_s[b, f]])
        end
        p
    end

    # b on its own pattern: a product with it is then far cheaper than with K's
    slip_pos = [slip_positions(f) for f in slip_faces]
    slip_nz = sort(unique(filter(!=(0), reduce(vcat, vec.(slip_pos); init=Int[]))))
    slip_index = zeros(Int, length(K.nzval))
    slip_index[slip_nz] = eachindex(slip_nz)
    dKdβ = sparse(K.rowval[slip_nz], [searchsortedlast(K.colptr, p) for p in slip_nz],
                  ones(length(slip_nz)), n, n)
    fill!(dKdβ.nzval, 0.0)

    nzv = length(K.nzval)
    return ElasticSystem(K, zeros(nzv), dKdβ, slip_nz, _gather_lists(pos, nzv)...,
                         slip_faces, [map(p -> p == 0 ? 0 : slip_index[p], sp) for sp in slip_pos],
                         zeros(nUd^2, ne), zeros(nUd^2, ne), dof, BasisFunctionCache(mesh),
                         nothing, false, zeros(n, 2),
                         Tuple{Matrix{Float64},Array{Float64,3},Array{Float64,3}}[])
end

# Inverse (row-major entries) and determinant of J = coords * ΔN, the 3×3 Jacobian of the
# isoparametric map at one Gauss point.
function _inverse_jacobian(coords::Matrix{Float64}, ΔN::Matrix{Float64})
    a = b = c = d = e = f = g = h = k = 0.0
    @inbounds for n in axes(ΔN, 1)
        x, y, z = coords[1, n], coords[2, n], coords[3, n]
        d1, d2, d3 = ΔN[n, 1], ΔN[n, 2], ΔN[n, 3]
        a += x*d1; b += x*d2; c += x*d3
        d += y*d1; e += y*d2; f += y*d3
        g += z*d1; h += z*d2; k += z*d3
    end
    A, B, C = e*k - f*h, -(d*k - f*g), d*h - e*g
    detJ = a*A + b*B + c*C
    inv = (A/detJ, -(b*k - c*h)/detJ, (b*f - c*e)/detJ,
           B/detJ, (a*k - c*g)/detJ, -(a*f - c*d)/detJ,
           C/detJ, -(a*h - b*g)/detJ, (a*e - b*d)/detJ)
    return inv, detJ
end

# Element matrices μA + λL and L of element e, written per node pair so no B matrix is formed:
#     A[(a,i),(b,j)] = Σ w (δij ∇Na·∇Nb + ∂jNa ∂iNb),   L[(a,i),(b,j)] = Σ w ∂iNa ∂jNb.
# The physical gradients of all Gauss points come first, stored Gauss point fastest, so each 3×3
# block is summed over the Gauss points in registers.
function _element_matrices!(KK, KL, coords::Matrix{Float64}, dN::Array{Float64,3}, wdN::Array{Float64,3},
                            e::Int, X::Matrix{Float64}, IEN::Matrix{Int},
                            ΔN_gp::Vector{Matrix{Float64}}, wpoints::Vector{Float64}, λ::Float64, μ::Float64)
    ngp, _, nN = size(dN)
    nUd = 3nN
    @inbounds for a in 1:nN, i in 1:3
        coords[i, a] = X[i, IEN[a, e]]
    end
    @inbounds for gp in 1:ngp
        (i11, i12, i13, i21, i22, i23, i31, i32, i33), detJ = _inverse_jacobian(coords, ΔN_gp[gp])
        w = wpoints[gp] * abs(detJ)
        ΔN = ΔN_gp[gp]
        for n in 1:nN
            x, y, z = ΔN[n, 1], ΔN[n, 2], ΔN[n, 3]
            g1 = x*i11 + y*i21 + z*i31
            g2 = x*i12 + y*i22 + z*i32
            g3 = x*i13 + y*i23 + z*i33
            dN[gp, 1, n], dN[gp, 2, n], dN[gp, 3, n] = g1, g2, g3
            wdN[gp, 1, n], wdN[gp, 2, n], wdN[gp, 3, n] = w*g1, w*g2, w*g3
        end
    end
    @inbounds for b in 1:nN, a in 1:nN
        p11 = p12 = p13 = p21 = p22 = p23 = p31 = p32 = p33 = 0.0   # pij = Σ w ∂iNa ∂jNb
        @simd for gp in 1:ngp
            a1, a2, a3 = wdN[gp, 1, a], wdN[gp, 2, a], wdN[gp, 3, a]
            b1, b2, b3 = dN[gp, 1, b], dN[gp, 2, b], dN[gp, 3, b]
            p11 += a1*b1; p12 += a1*b2; p13 += a1*b3
            p21 += a2*b1; p22 += a2*b2; p23 += a2*b3
            p31 += a3*b1; p32 += a3*b2; p33 += a3*b3
        end
        μg = μ * (p11 + p22 + p33)
        r = 3(a-1)
        c1, c2, c3 = (3(b-1))*nUd, (3(b-1) + 1)*nUd, (3(b-1) + 2)*nUd
        KK[c1 + r + 1] = μg + (μ + λ)*p11; KK[c2 + r + 1] = μ*p21 + λ*p12;  KK[c3 + r + 1] = μ*p31 + λ*p13
        KK[c1 + r + 2] = μ*p12 + λ*p21;  KK[c2 + r + 2] = μg + (μ + λ)*p22; KK[c3 + r + 2] = μ*p32 + λ*p23
        KK[c1 + r + 3] = μ*p13 + λ*p31;  KK[c2 + r + 3] = μ*p23 + λ*p32;  KK[c3 + r + 3] = μg + (μ + λ)*p33
        KL[c1 + r + 1] = p11;     KL[c2 + r + 1] = p12;     KL[c3 + r + 1] = p13
        KL[c1 + r + 2] = p21;     KL[c2 + r + 2] = p22;     KL[c3 + r + 2] = p23
        KL[c1 + r + 3] = p31;     KL[c2 + r + 3] = p32;     KL[c3 + r + 3] = p33
    end
    return nothing
end

"""
    _assemble!(sys::ElasticSystem, X::Matrix{Float64}, IEN::Matrix{Int}, λ, μ, β; sensitivities=false)

Fill `sys.K = μA + λL + βb` and `sys.dKdβ = b` at the geometry `X`, and with `sensitivities`
also `sys.nz_L = L`.
"""
function _assemble!(sys::ElasticSystem, X::Matrix{Float64}, IEN::Matrix{Int}, λ::Real, μ::Real, β::Real;
                    sensitivities::Bool=false)
    ΔN_gp, wpoints = sys.cache.ΔN, sys.cache.w
    nN, ne = size(IEN)
    ngp = length(wpoints)

    chunks = collect(Iterators.partition(1:ne, cld(ne, Threads.nthreads())))
    while length(sys.scratch) < length(chunks)
        push!(sys.scratch, (zeros(3, nN), zeros(ngp, 3, nN), zeros(ngp, 3, nN)))
    end
    @sync for (chunk, (coords, dN, wdN)) in zip(chunks, sys.scratch)
        Threads.@spawn for e in chunk
            _element_matrices!(view(sys.Ke, :, e), view(sys.Ke_L, :, e), coords, dN, wdN, e,
                               X, IEN, ΔN_gp, wpoints, Float64(λ), Float64(μ))
        end
    end

    nzv = length(sys.K.nzval)
    @sync for chunk in Iterators.partition(1:nzv, cld(nzv, Threads.nthreads()))
        Threads.@spawn begin
            _gather_entries!(sys.K.nzval, sys.Ke, sys.ptr, sys.src, chunk)
            sensitivities && _gather_entries!(sys.nz_L, sys.Ke_L, sys.ptr, sys.src, chunk)
        end
    end

    nz_b = sys.dKdβ.nzval
    fill!(nz_b, 0.0)
    for (IEN_s, pos) in zip(sys.slip_faces, sys.slip_pos)
        _slip_faces!(nz_b, pos, IEN_s, X, sys.cache.N_s, sys.cache.ΔN_s, sys.cache.w_s)
    end
    nzval = sys.K.nzval
    @inbounds for (i, p) in enumerate(sys.slip_nz)
        nzval[p] += β * nz_b[i]
    end
    return sys
end

function _gather_entries!(nz::Vector{Float64}, Ke::Matrix{Float64}, ptr::Vector{Int32}, src::Vector{Int32},
                          ps::UnitRange{Int})
    @inbounds for p in ps
        v = 0.0
        for k in ptr[p]:(ptr[p+1] - 1)
            v += Ke[src[k]]
        end
        nz[p] = v
    end
    return nothing
end

# Tangential slip mass ∫ Na Nb of a group of horizontal faces, scattered into nz.
function _slip_faces!(nz::Vector{Float64}, pos::Matrix{Int}, IEN_s::Matrix{Int}, X::Matrix{Float64},
                      N_gp::Vector{Vector{Float64}}, ΔN_gp::Vector{Matrix{Float64}}, wpoints::Vector{Float64})
    nS = size(IEN_s, 1)
    @inbounds for f in axes(IEN_s, 2), gp in eachindex(wpoints)
        ΔN = ΔN_gp[gp]
        t1x = t1y = t1z = t2x = t2y = t2z = 0.0
        for a in 1:nS
            n = IEN_s[a, f]
            t1x += X[1,n]*ΔN[a,1]; t1y += X[2,n]*ΔN[a,1]; t1z += X[3,n]*ΔN[a,1]
            t2x += X[1,n]*ΔN[a,2]; t2y += X[2,n]*ΔN[a,2]; t2z += X[3,n]*ΔN[a,2]
        end
        w = wpoints[gp] * sqrt((t1y*t2z - t1z*t2y)^2 + (t1z*t2x - t1x*t2z)^2 + (t1x*t2y - t1y*t2x)^2)
        N = N_gp[gp]
        for i in 1:2, b in 1:nS, a in 1:nS
            p = pos[((i-1)*nS + (b-1))*nS + a, f]
            p != 0 && (nz[p] += w * N[a] * N[b])
        end
    end
    return nothing
end

# y = α M x for a symmetric matrix with K's pattern and values `nz`, computed as α M'x: M'x
# splits by column across threads without write conflicts. Anything that makes K unsymmetric
# must change this.
function _symv!(y::AbstractVector{Float64}, K::SparseMatrixCSC{Float64,Int}, nz::Vector{Float64},
                x::AbstractVector{Float64}, α::Real=1.0)
    n = size(K, 2)
    @sync for cols in Iterators.partition(1:n, cld(n, Threads.nthreads()))
        Threads.@spawn begin
            cp, rv = K.colptr, K.rowval
            @inbounds for j in cols
                s = 0.0
                for k in cp[j]:(cp[j+1] - 1)
                    s += nz[k] * x[rv[k]]
                end
                y[j] = α * s
            end
        end
    end
    return y
end

function _residual!(res, sys::ElasticSystem, rhs, x)
    _symv!(res, sys.K, sys.K.nzval, x)
    @. res = rhs - res
    return res
end

# Numeric factorization of the current K, reusing the symbolic analysis when there is one.
function _refactor!(sys::ElasticSystem)
    sys.lu === nothing && (sys.lu = _sparse_lu(sys.K))
    t = @elapsed _factorize!(sys.lu, sys.K)
    @debug "factorization: $(round(1e3 * t, digits=1)) ms"
    _singular(sys.lu) && @warn "the stiffness is singular: a rigid-body mode is not removed by the boundary conditions, so the solution carries an arbitrary amount of it"
    sys.stale = false
    return sys
end

"""
    _refined_solve!(sys::ElasticSystem, x, rhs; tol, max_rounds=10) -> rounds

Solve `sys.K x = rhs` to relative residual `tol` by iterative refinement against the current
factorization, which may belong to an earlier geometry, starting from the value in `x`: the
mesh moves little per step, so a stale factorization still converges in a few cheap triangular
solves. Refactorizes when the factorization has been found stale (`_stale`, `_refactor_now`).
Returns the number of refinement rounds.
"""
function _refined_solve!(sys::ElasticSystem, x::Vector{Float64}, rhs::Vector{Float64};
                         tol::Float64, max_rounds::Int=10)
    nb = norm(rhs)
    nb == 0 && return (fill!(x, 0.0); 0)
    res, dx = view(sys.buf, :, 1), view(sys.buf, :, 2)
    if sys.lu === nothing || sys.stale
        _refactor!(sys)
        _solve!(x, sys.lu, rhs)
    end
    _residual!(res, sys, rhs, x)
    r = norm(res)
    rounds, ρ, refactored = 0, 0.0, false
    while r > tol * nb
        rounds > 2max_rounds && error("refined solve did not converge: relative residual $(r / nb)")
        if !refactored && rounds > 0 && (rounds == max_rounds || _refactor_now(sys.lu, ρ, rounds))
            _refactor!(sys)
            refactored = true
            _solve!(x, sys.lu, rhs)
        else
            _solve!(dx, sys.lu, res)
            x .+= dx
        end
        _residual!(res, sys, rhs, x)
        ρ, r = norm(res) / r, norm(res)
        _stale(sys.lu, ρ) && (sys.stale = true)
        rounds += 1
    end
    @debug "refined solve: $rounds rounds"
    return rounds
end

# Nodal displacements (3 × nNodes) of a solution x: eliminated components are zero.
function _nodal!(u::Matrix{Float64}, x::AbstractVector{Float64}, dof::Matrix{Int})
    @inbounds for i in eachindex(u, dof)
        k = dof[i]
        u[i] = k == 0 ? 0.0 : x[k]
    end
    return u
end

# Starting guess from the last `k` (≤ 4) solutions in `hist`, newest first, by polynomial
# extrapolation over the steps.
function _extrapolate!(x::Vector{Float64}, hist::Vector{Vector{Float64}}, k::Int)
    h1, h2, h3, h4 = hist
    if k == 0
        fill!(x, 0.0)
    elseif k == 1
        copyto!(x, h1)
    elseif k == 2
        @. x = 2h1 - h2
    elseif k == 3
        @. x = 3h1 - 3h2 + h3
    else
        @. x = 4h1 - 6h2 + 4h3 - h4
    end
    return x
end

# Make `x` the newest entry of `hist`, dropping the oldest, without allocating.
function _push_history!(hist::Vector{Vector{Float64}}, x::Vector{Float64})
    hist[1], hist[2], hist[3], hist[4] = hist[4], hist[1], hist[2], hist[3]
    copyto!(hist[1], x)
    return hist
end
