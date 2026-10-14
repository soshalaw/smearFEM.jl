using LinearAlgebra
using SparseArrays

# Fixed-pattern assembly and lagged-factorization solve of the squeeze-flow Stokes system
#
#     M = [T'AT  T'B ; B'T  0],   A = η Ā + β b,   T = [C_uc  q_top]
#
# which is the matrix `simulate` used to build with sparse products and block concatenation.
# Unknowns keep that ordering: [q_f; p; μ], μ being the top-plate velocity — an unknown under
# force control, prescribed under velocity control (`_prescribe_plate_velocity!`). The
# boundary conditions are the same in both modes. The mesh topology never changes, so the sparsity
# pattern is fixed: each `K.nzval` entry gathers its element contributions through
# precomputed lists. b is also kept on its own, as dM/dβ; dM/dη is never needed (see the
# homogeneity note in `simulate`).

"""
    SqueezeSystem

Saddle-point system of the squeeze flow, with a fixed sparsity pattern and a factorization
that is reused across time steps. Build it with `_squeeze_system`.
"""
mutable struct SqueezeSystem
    K::SparseMatrixCSC{Float64,Int}       # M at the current geometry
    dKdβ::SparseMatrixCSC{Float64,Int}    # ∂M/∂β = b, on its own (much sparser) pattern
    nz_slip::Vector{Float64}              # dKdβ.nzval
    slip_pos::Vector{Int}                 # K.nzval index of each nz_slip entry
    visc_ptr::Vector{Int32}               # nzval entry p sums Ke_visc[visc_src[visc_ptr[p]:visc_ptr[p+1]-1]]
    visc_src::Vector{Int32}               # (Int32: assembly is memory-bandwidth bound)
    div_ptr::Vector{Int32}                # same for Ke_div (B and its B' mirror)
    div_src::Vector{Int32}
    pos_slip_top::Matrix{Int}             # (2*nS², ne_surf) tangential slip entry -> nz_slip index
    pos_slip_btm::Matrix{Int}
    Ke_visc::Matrix{Float64}              # (81², ne) element matrices, filled in parallel
    Ke_div::Matrix{Float64}               # (nP*81, ne)
    lu::Union{Nothing,_Pardiso,_Umfpack}  # factorization the solves use (sparse_lu.jl)
    stale::Bool                           # refactorize before the next solve
    steps_since_factor::Int               # lagged solves since it was factorized
    η::Float64                            # viscosity K was last assembled at
    η_factor::Float64                     # viscosity F was factorized at
    nf::Int                               # free velocity unknowns
    np::Int                               # pressure unknowns
    node_dofs::Matrix{Int}                # (nDof, nNodes) unknown of each nodal velocity, 0 = eliminated
    plate_row::Vector{Int}                # nzval indices of row n (μ's equation) outside column n
    buf::Matrix{Float64}                  # (n, 7) work vectors of `_lagged_solve!`
    scratch::Vector{Tuple{Matrix{Float64},Array{Float64,3},Array{Float64,3}}}  # per assembly chunk
end

# Invert an element-entry -> nzval position map `pos` (0 = eliminated) into per-entry lists:
# entry p gathers the element-matrix values at linear indices src[ptr[p]:ptr[p+1]-1]. They
# come out in element order, so a parallel per-entry sum adds in the same order as a serial
# scatter over elements. A trailing dimension of `pos` maps the same element values twice.
# `src_of` (shaped like `pos`) redirects an entry to another element-matrix value.
function _gather_lists(pos::AbstractArray{Int}, nnz::Int; src_of::Union{Nothing,AbstractArray{Int}}=nothing)
    nsrc = size(pos, 1) * size(pos, 2)
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
        src[next[p]] = src_of === nothing ? mod1(k, nsrc) : src_of[k]
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
    _squeeze_system(mdl::Stokes, scene::SqueezeFlow) -> SqueezeSystem

Build the sparsity pattern and element-to-nzval maps of the squeeze-flow system (both control
modes). Values are filled by `_assemble_squeeze_system!`.
"""
function _squeeze_system(mdl::Stokes, scene::SqueezeFlow)
    # the mesh fields are abstractly typed; concrete locals keep the loops below fast
    IEN::Matrix{Int} = mdl.mesh_u.IEN
    ID::Matrix{Int} = mdl.mesh_u.ID
    IEN_top::Matrix{Int} = mdl.mesh_u.IEN_top
    IEN_bottom::Matrix{Int} = mdl.mesh_u.IEN_bottom
    IEN_p::Matrix{Int} = mdl.mesh_p.IEN
    nDof = mdl.nDof_u
    nU, ne = size(IEN)
    nP = size(IEN_p, 1)
    nS = size(IEN_top, 1)

    C_uc = sparse(scene.C_uc)
    q_top = sparse(scene.q_d[:top])[:]
    nf = size(C_uc, 2)
    np = mdl.mesh_p.nNodes
    iμ = nf + np + 1

    # reduced index of every global velocity dof: free, top-plate μ, or 0 (eliminated)
    red = zeros(Int, size(C_uc, 1))
    rows, cols, _ = findnz(C_uc)
    red[rows] .= cols
    red[findnz(q_top)[1]] .= iμ

    ru = [red[ID[i, IEN[a, e]]] for i in 1:nDof, a in 1:nU, e in 1:ne]   # (nDof, nU, ne)
    ru = reshape(ru, nDof * nU, ne)
    rp = IEN_p .+ nf

    Is = Int[]; Js = Int[]
    for e in 1:ne
        for c in ru[:, e], r in ru[:, e]
            r != 0 && c != 0 && (push!(Is, r); push!(Js, c))
        end
        for k in rp[:, e], r in ru[:, e]
            r != 0 && (push!(Is, r, k); push!(Js, k, r))
        end
    end
    n = iμ
    K = sparse(Is, Js, ones(length(Is)), n, n)

    nUd = nDof * nU
    pos_visc = zeros(Int, nUd^2, ne)
    pos_div = zeros(Int, nP * nUd, ne, 2)
    for e in 1:ne
        for j in 1:nUd, i in 1:nUd
            pos_visc[(j-1)*nUd + i, e] = _nz_index(K, ru[i, e], ru[j, e])
        end
        for l in 1:nUd, k in 1:nP
            pos_div[(l-1)*nP + k, e, 1] = _nz_index(K, ru[l, e], rp[k, e])
            pos_div[(l-1)*nP + k, e, 2] = _nz_index(K, rp[k, e], ru[l, e])
        end
    end

    # slip acts on the two tangential components only
    function slip_positions(IEN_s)
        pos = zeros(Int, 2 * nS^2, size(IEN_s, 2))
        for f in axes(IEN_s, 2), i in 1:2, b in 1:nS, a in 1:nS
            r = red[ID[i, IEN_s[a, f]]]
            c = red[ID[i, IEN_s[b, f]]]
            pos[((i-1)*nS + (b-1))*nS + a, f] = _nz_index(K, r, c)
        end
        pos
    end

    # the element kernel only computes node pairs a ≤ b of the symmetric viscous block; the
    # lower node blocks read their mirror
    src_visc = zeros(Int, nUd^2, ne)
    for e in 1:ne, c in 1:nUd, r in 1:nUd
        upper = div(r - 1, nDof) <= div(c - 1, nDof)
        src_visc[(c-1)*nUd + r, e] = (e - 1) * nUd^2 + (upper ? (c-1)*nUd + r : (r-1)*nUd + c)
    end

    # b on its own pattern: a product with it is then ~100x cheaper than with K's pattern
    pos_top, pos_btm = slip_positions(IEN_top), slip_positions(IEN_bottom)
    slip_pos = sort(unique(filter(!=(0), vcat(vec(pos_top), vec(pos_btm)))))
    slip_index = zeros(Int, length(K.nzval))
    slip_index[slip_pos] = eachindex(slip_pos)
    to_slip(pos) = map(p -> p == 0 ? 0 : slip_index[p], pos)
    dKdβ = sparse(K.rowval[slip_pos], [searchsortedlast(K.colptr, p) for p in slip_pos],
                  zeros(length(slip_pos)), n, n)
    nnz(dKdβ) == length(slip_pos) || error("slip pattern lost entries")

    nzv = length(K.nzval)
    return SqueezeSystem(K, dKdβ, dKdβ.nzval, slip_pos,
                       _gather_lists(pos_visc, nzv; src_of=src_visc)..., _gather_lists(pos_div, nzv)...,
                       to_slip(pos_top), to_slip(pos_btm),
                       zeros(nUd^2, ne), zeros(nP * nUd, ne),
                       nothing, false, 0, 1.0, 1.0, nf, np,
                       [red[ID[i, a]] for i in 1:nDof, a in axes(ID, 2)],
                       [k for c in 1:(n - 1) for k in K.colptr[c]:(K.colptr[c+1] - 1) if K.rowval[k] == n],
                       zeros(n, 7), Tuple{Matrix{Float64},Array{Float64,3},Array{Float64,3}}[])
end

# The pattern, position maps and symbolic factorization depend only on the mesh topology and
# the boundary conditions — both fixed per model by `def_problem` — so each model keeps one
# SqueezeSystem across `simulate` calls (a fit calls `simulate` on the same model many times).
const _SQUEEZE_SYSTEMS = WeakKeyDict{Stokes,SqueezeSystem}()

_cached_squeeze_system(mdl::Stokes, scene::SqueezeFlow) =
    get!(() -> _squeeze_system(mdl, scene), _SQUEEZE_SYSTEMS, mdl)

# Velocity control: the plate velocity μ is given, so its row and column of K become the
# identity and the column moves to the right-hand side, written to `rhs` for a unit plate
# velocity. The pattern is unchanged, so both control modes share the factorization machinery.
function _prescribe_plate_velocity!(sys::SqueezeSystem, rhs::Vector{Float64})
    K = sys.K
    n = size(K, 1)
    fill!(rhs, 0.0)
    for k in K.colptr[n]:(K.colptr[n+1] - 1)
        r = K.rowval[k]
        r != n && (rhs[r] = -K.nzval[k])
        K.nzval[k] = r == n ? 1.0 : 0.0
    end
    K.nzval[sys.plate_row] .= 0.0
    rhs[n] = 1.0
    return rhs
end

# Numeric refactorization at the current K, reusing the symbolic analysis when there is one.
function _refactor!(sys::SqueezeSystem)
    sys.lu === nothing && (sys.lu = _sparse_lu(sys.K))
    t = @elapsed _factorize!(sys.lu, sys.K)
    @debug "refactorization after $(sys.steps_since_factor) lagged solves, $(round(1e3 * t, digits=1)) ms"
    sys.stale = false
    sys.steps_since_factor = 0
    sys.η_factor = sys.η
    return sys
end

# dx = F \ res, and dy = F \ res_y when given, for the factorization F made at viscosity
# η_factor while K is at η. K(η, β) = D K(1, β/η) D with D = diag(√η velocity, 1/√η pressure),
# so rescaling around the solve leaves only the β/η and geometry changes as the factorization's
# error — bulk viscosity changes η every step. `tmp` and `tmp_y` are scratch.
function _corrections!(sys::SqueezeSystem, dx, res, tmp, dy=nothing, res_y=nothing, tmp_y=nothing)
    s = sqrt(sys.η_factor / sys.η)
    if s != 1
        res = _scale_unknowns!(copyto!(tmp, res), sys, s, 1/s)
        dy === nothing || (res_y = _scale_unknowns!(copyto!(tmp_y, res_y), sys, s, 1/s))
    end
    dy === nothing ? _solve!(dx, sys.lu, res) : _solve_pair!(dx, dy, sys.lu, res, res_y)
    if s != 1
        _scale_unknowns!(dx, sys, s, 1/s)
        dy === nothing || _scale_unknowns!(dy, sys, s, 1/s)
    end
    return nothing
end

# Inverse (row-major entries) and determinant of J = coords * ΔNx, the 3×3 Jacobian of the
# isoparametric map at one Gauss point.
function _inverse_jacobian(coords::Matrix{Float64}, ΔNx::Matrix{Float64})
    a = b = c = d = e = f = g = h = k = 0.0
    @inbounds for n in axes(ΔNx, 1)
        x, y, z = coords[1, n], coords[2, n], coords[3, n]
        d1, d2, d3 = ΔNx[n, 1], ΔNx[n, 2], ΔNx[n, 3]
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

# Element viscous and divergence matrices of element e. The viscous block is 2∫ε(u):ε(v) — the
# same form as assemble_system_A's 2wB C B' with the "standard" C — written per node pair so no
# B matrix is formed. It is symmetric, so only node pairs a ≤ b are computed;
# `_squeeze_system`'s gather lists mirror the rest. The physical gradients of all Gauss points
# come first, stored Gauss point fastest, so each 3×3 block is summed over the Gauss points in
# registers instead of being read and written once per Gauss point.
# Arguments are concretely typed on purpose: BasisFunctionCache stores NTuple{7, Any}.
function _element_matrices!(Kv, Kd, coords::Matrix{Float64}, dN::Array{Float64,3}, wdN::Array{Float64,3},
                            e::Int, X::Matrix{Float64}, IEN_x::Matrix{Int},
                            ΔN_u_gp::Vector{Matrix{Float64}}, Np::Matrix{Float64},
                            ΔN_x_gp::Vector{Matrix{Float64}}, wpoints::Vector{Float64})
    ngp, _, nU = size(dN)
    nP = size(Np, 2)
    nUd = 3nU
    @inbounds for a in axes(IEN_x, 1), i in 1:3
        coords[i, a] = X[i, IEN_x[a, e]]
    end
    @inbounds for gp in 1:ngp
        (i11, i12, i13, i21, i22, i23, i31, i32, i33), detJ = _inverse_jacobian(coords, ΔN_x_gp[gp])
        w = wpoints[gp] * abs(detJ)
        ΔN = ΔN_u_gp[gp]
        for n in 1:nU
            x, y, z = ΔN[n, 1], ΔN[n, 2], ΔN[n, 3]
            g1 = x*i11 + y*i21 + z*i31
            g2 = x*i12 + y*i22 + z*i32
            g3 = x*i13 + y*i23 + z*i33
            dN[gp, 1, n], dN[gp, 2, n], dN[gp, 3, n] = g1, g2, g3
            wdN[gp, 1, n], wdN[gp, 2, n], wdN[gp, 3, n] = w*g1, w*g2, w*g3
        end
    end
    # K[(a,i),(b,j)] = Σ w (δij ∇Na·∇Nb + ∂jNa ∂iNb)
    @inbounds for b in 1:nU, a in 1:b
        g = k11 = k12 = k13 = k21 = k22 = k23 = k31 = k32 = k33 = 0.0
        @simd for gp in 1:ngp
            a1, a2, a3 = wdN[gp, 1, a], wdN[gp, 2, a], wdN[gp, 3, a]
            b1, b2, b3 = dN[gp, 1, b], dN[gp, 2, b], dN[gp, 3, b]
            g += a1*b1 + a2*b2 + a3*b3
            k11 += a1*b1; k12 += a2*b1; k13 += a3*b1
            k21 += a1*b2; k22 += a2*b2; k23 += a3*b2
            k31 += a1*b3; k32 += a2*b3; k33 += a3*b3
        end
        r = 3(a-1)
        c1, c2, c3 = (3(b-1))*nUd, (3(b-1) + 1)*nUd, (3(b-1) + 2)*nUd
        Kv[c1 + r + 1] = k11 + g; Kv[c2 + r + 1] = k12;     Kv[c3 + r + 1] = k13
        Kv[c1 + r + 2] = k21;     Kv[c2 + r + 2] = k22 + g; Kv[c3 + r + 2] = k23
        Kv[c1 + r + 3] = k31;     Kv[c2 + r + 3] = k32;     Kv[c3 + r + 3] = k33 + g
    end
    # Kd[(a,i), k] = -Σ w ∂iNa Npk
    @inbounds for a in 1:nU, i in 1:3
        base = (3(a-1) + i - 1) * nP
        for k in 1:nP
            acc = 0.0
            @simd for gp in 1:ngp
                acc += wdN[gp, i, a] * Np[gp, k]
            end
            Kd[base + k] = -acc
        end
    end
    return nothing
end

# Tangential slip mass of the top and bottom faces, scattered into nz_slip (= dKdβ.nzval).
function _assemble_slip!(sys::SqueezeSystem, mdl::Stokes, cache::BasisFunctionCache)
    N_s_gp, ΔN_s_gp, _, _, _, _, wpoints = cache.bf_surface
    X::Matrix{Float64} = mdl.mesh_x.NodeList
    fill!(sys.nz_slip, 0.0)
    _slip_faces!(sys.nz_slip, sys.pos_slip_top, mdl.mesh_x.IEN_top, X, N_s_gp, ΔN_s_gp, wpoints)
    _slip_faces!(sys.nz_slip, sys.pos_slip_btm, mdl.mesh_x.IEN_bottom, X, N_s_gp, ΔN_s_gp, wpoints)
    return nothing
end

function _slip_faces!(nz::Vector{Float64}, pos::Matrix{Int}, IEN_s::Matrix{Int}, X::Matrix{Float64},
                      N_s_gp::Vector{Vector{Float64}}, ΔN_s_gp::Vector{Matrix{Float64}}, wpoints::Vector{Float64})
    nS = length(N_s_gp[1])
    @inbounds for f in axes(IEN_s, 2), gp in eachindex(wpoints)
        ΔN = ΔN_s_gp[gp]
        t1x = t1y = t1z = t2x = t2y = t2z = 0.0
        for a in 1:nS
            n = IEN_s[a, f]
            t1x += X[1,n]*ΔN[a,1]; t1y += X[2,n]*ΔN[a,1]; t1z += X[3,n]*ΔN[a,1]
            t2x += X[1,n]*ΔN[a,2]; t2y += X[2,n]*ΔN[a,2]; t2z += X[3,n]*ΔN[a,2]
        end
        w = wpoints[gp] * sqrt((t1y*t2z - t1z*t2y)^2 + (t1z*t2x - t1x*t2z)^2 + (t1x*t2y - t1y*t2x)^2)
        N = N_s_gp[gp]
        for i in 1:2, b in 1:nS, a in 1:nS
            p = pos[((i-1)*nS + (b-1))*nS + a, f]
            p != 0 && (nz[p] += w * N[a] * N[b])
        end
    end
    return nothing
end

"""
    _assemble_squeeze_system!(sys::SqueezeSystem, mdl::Stokes, cache::BasisFunctionCache, η, β)

Assemble M = [T'AT T'B; B'T 0] at the current geometry `mdl.mesh_x.NodeList` into `sys.K`.
"""
function _assemble_squeeze_system!(sys::SqueezeSystem, mdl::Stokes, cache::BasisFunctionCache, η::Real, β::Real)
    _, ΔN_u_gp, N_p_gp, _, _, ΔN_x_gp, wpoints = cache.bf_volume
    X::Matrix{Float64} = mdl.mesh_x.NodeList
    IEN_x::Matrix{Int} = mdl.mesh_x.IEN
    nX, ne = size(IEN_x)
    nU = size(ΔN_u_gp[1], 1)
    ngp = length(wpoints)
    Np = [N_p_gp[gp][k] for gp in 1:ngp, k in eachindex(N_p_gp[1])]
    Kv_all, Kd_all = sys.Ke_visc, sys.Ke_div

    chunks = collect(Iterators.partition(1:ne, cld(ne, Threads.nthreads())))
    while length(sys.scratch) < length(chunks)
        push!(sys.scratch, (zeros(3, nX), zeros(ngp, 3, nU), zeros(ngp, 3, nU)))
    end
    @sync for (chunk, (coords, dN, wdN)) in zip(chunks, sys.scratch)
        Threads.@spawn for e in chunk
            _element_matrices!(view(Kv_all, :, e), view(Kd_all, :, e), coords, dN, wdN, e,
                               X, IEN_x, ΔN_u_gp, Np, ΔN_x_gp, wpoints)
        end
    end

    _assemble_slip!(sys, mdl, cache)

    sys.η = η
    nzv = length(sys.K.nzval)
    @sync for chunk in Iterators.partition(1:nzv, cld(nzv, Threads.nthreads()))
        Threads.@spawn _gather_entries!(sys, chunk, η)
    end
    nzval, nz_slip = sys.K.nzval, sys.nz_slip
    @inbounds for (i, p) in enumerate(sys.slip_pos)
        nzval[p] += β * nz_slip[i]
    end
    return sys
end

# Sum the element contributions of nzval entries `ps` into K (slip is added afterwards).
function _gather_entries!(sys::SqueezeSystem, ps::UnitRange{Int}, η::Real)
    Kv, Kd = sys.Ke_visc, sys.Ke_div
    vp, vs, dp, ds = sys.visc_ptr, sys.visc_src, sys.div_ptr, sys.div_src
    nzval = sys.K.nzval
    @inbounds for p in ps
        v = 0.0
        for k in vp[p]:(vp[p+1] - 1)
            v += Kv[vs[k]]
        end
        d = 0.0
        for k in dp[p]:(dp[p+1] - 1)
            d += Kd[ds[k]]
        end
        nzval[p] = η * v + d
    end
    return nothing
end

# Starting guess from the last `k` (≤ 4) solutions in `hist`, newest first, by polynomial
# extrapolation in time. From exact solutions on the benchmark the cubic guess has relative
# residual ~5e-11 (constant viscosity; ~1e-8 bulk), quadratic ~7e-9, linear ~2e-6, the last
# solution alone ~7e-4. The stored solutions carry their own tolerance-level residuals, which
# the cubic amplifies up to 15×, so in practice a solve still takes about one correction.
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

# Nodal velocities (nDof × nNodes) of a solution x = [q_f; p; μ]: the same values as
# C_uc*q_f + μ*q_top, which only select entries; eliminated components are zero.
function _nodal_velocity!(v::Matrix{Float64}, x::Vector{Float64}, sys::SqueezeSystem)
    @inbounds for i in eachindex(v, sys.node_dofs)
        k = sys.node_dofs[i]
        v[i] = k == 0 ? 0.0 : x[k]
    end
    return v
end
_nodal_velocity(x::Vector{Float64}, sys::SqueezeSystem) = _nodal_velocity!(similar(sys.node_dofs, Float64), x, sys)

# Multiply the velocity unknowns [q_f; μ] of x by su and the pressure unknowns by sp.
function _scale_unknowns!(x::AbstractVector{Float64}, sys::SqueezeSystem, su::Real, sp::Real)
    nf, np = sys.nf, sys.np
    @views x[1:nf] .*= su
    @views x[(nf + 1):(nf + np)] .*= sp
    x[end] *= su
    return x
end

# res = rhs - K x, computed as rhs - K'x: K is symmetric (to rounding), and K'x splits by
# column across threads without write conflicts (4x faster than the serial K*x here).
# Anything that makes K unsymmetric must change this.
function _residual!(res::AbstractVector{Float64}, K::SparseMatrixCSC{Float64,Int},
                   rhs::AbstractVector{Float64}, x::AbstractVector{Float64})
    n = size(K, 2)
    @sync for cols in Iterators.partition(1:n, cld(n, Threads.nthreads()))
        Threads.@spawn _residual_cols!(res, K, rhs, x, cols)
    end
    return res
end

function _residual_cols!(res, K, rhs, x, cols::UnitRange{Int})
    cp, rv, nz = K.colptr, K.rowval, K.nzval
    @inbounds for j in cols
        s = 0.0
        for k in cp[j]:(cp[j+1] - 1)
            s += nz[k] * x[rv[k]]
        end
        res[j] = rhs[j] - s
    end
    return nothing
end

# Relative-residual tolerances of the lagged solves. The state accumulates into the geometry
# over the horizon; the sensitivities only feed the Gauss-Newton step.
const STATE_TOL = 1e-9
const SENS_TOL = 1e-8

"""
    _lagged_solve!(sys::SqueezeSystem, z, rhs_z, w=nothing, rhs_w! =nothing; max_rounds=10)

Solve `sys.K z = rhs_z` to relative residual STATE_TOL and, when `w` is given, `sys.K w = r`
to SENS_TOL, where `rhs_w!(r, z)` builds r from z — it is rebuilt whenever z changes. Both
are found by iterative refinement against the factorization of an earlier geometry, starting
from the values in `z` and `w`: the mesh moves little per step, so a stale factorization
converges in a few cheap triangular solves. When the factorization has been found stale
(`_stale`, `_refactor_now`), the system is refactorized and solved directly instead. Returns
the number of refinement rounds.
"""
function _lagged_solve!(sys::SqueezeSystem, z::Vector{Float64}, rhs_z::Vector{Float64},
                        w::Union{Nothing,Vector{Float64}}=nothing, rhs_w! =nothing;
                        max_rounds::Int=10)
    sys.lu === nothing && _refactor!(sys)
    res_z, res_w, dz, dw, tmp_z, tmp_w, rhs_w = eachcol(sys.buf)
    if sys.stale
        _refactor!(sys)
        _solve!(z, sys.lu, rhs_z)
        if w !== nothing
            rhs_w!(rhs_w, z)
            _solve!(w, sys.lu, rhs_w)
        end
    end
    sys.steps_since_factor += 1
    _residual!(res_z, sys.K, rhs_z, z)
    if w !== nothing
        rhs_w!(rhs_w, z)
        _residual!(res_w, sys.K, rhs_w, w)
    end
    @debug "lagged solve: initial residual $(norm(res_z) / norm(rhs_z)), $(sys.steps_since_factor) steps since factorization"
    rounds, ρ, refactored = 0, 0.0, false
    while true
        r_z, r_w = norm(res_z), w === nothing ? 0.0 : norm(res_w)
        todo_z = r_z > STATE_TOL * norm(rhs_z)
        todo_w = w !== nothing && r_w > SENS_TOL * norm(rhs_w)
        (todo_z || todo_w) || break
        rounds > 2max_rounds && error("lagged solve did not converge: relative residuals $(r_z / norm(rhs_z)), $(r_w / norm(rhs_w))")
        rounds > 0 && _stale(sys.lu, ρ) && (sys.stale = true)
        if !refactored && rounds > 0 && (rounds == max_rounds || _refactor_now(sys.lu, ρ, rounds))
            _refactor!(sys)
            refactored = true
            _solve!(z, sys.lu, rhs_z)
            _residual!(res_z, sys.K, rhs_z, z)
            if w !== nothing
                rhs_w!(rhs_w, z)
                _solve!(w, sys.lu, rhs_w)
                _residual!(res_w, sys.K, rhs_w, w)
            end
            continue
        end
        if todo_z && todo_w
            _corrections!(sys, dz, res_z, tmp_z, dw, res_w, tmp_w)
        elseif todo_z
            _corrections!(sys, dz, res_z, tmp_z)
        else
            _corrections!(sys, dw, res_w, tmp_w)
        end
        if todo_z
            z .+= dz
            _residual!(res_z, sys.K, rhs_z, z)
            w !== nothing && rhs_w!(rhs_w, z)
            ρ = norm(res_z) / r_z
        end
        if w !== nothing
            todo_w && (w .+= dw)
            _residual!(res_w, sys.K, rhs_w, w)
            todo_z || (ρ = norm(res_w) / r_w)
        end
        rounds += 1
    end
    @debug "lagged solve: $rounds corrections"
    return rounds
end
