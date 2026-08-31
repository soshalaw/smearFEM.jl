# Contour cost functions: how a simulated contour is scored against an observed one, and the
# derivatives the Gauss-Newton optimizers need from it. Add a cost by writing a `_residuals`
# method; normalization, per-frame accumulation and outlier skipping are shared.
#
# Correspondences come from a nearest-neighbour search and are held fixed while
# differentiating (the standard ICP linearization). `knn` is piecewise constant in θ, so it
# must not be differentiated through.

using LinearAlgebra
using ArgCheck
using NearestNeighbors

"""
    ContourCost

Abstract supertype for contour cost functions. A subtype is a stateless tag selecting a
`_residuals` method.

Every cost here is squared (`Δx² + Δy²`, never `√`) — Gauss-Newton needs `uᵀu` for `Jᵀu`
and `JᵀJ` to be its gradient and Hessian. So none of them equal the same-named metrics in
`analysis/pointcloud_metrics.jl`, which are rooted for interpretability: compare a cost
only against another cost.
"""
abstract type ContourCost end

"""
    ClosestPointCost()

One-directional cost: every *simulated* point is matched to its nearest observed point and
contributes the squared offset. Observed points with no simulated point near them are not
penalized, so a contour covering only part of the observation still scores well.

The historical behavior, and the default everywhere, so earlier results reproduce exactly.
"""
struct ClosestPointCost <: ContourCost end

"""
    ChamferCost()

Symmetric squared-Chamfer cost: the [`ClosestPointCost`](@ref) term plus its reverse, where
every *observed* point is matched to its nearest simulated point. Penalizes an incomplete
simulated contour, which the one-directional cost cannot see.

The reverse term reuses simulated points — several observed points may share a nearest
simulated neighbour, so those rows of `∂u/∂θ` repeat. Correct, but it reweights `JᵀJ`, so
Hessian conditioning differs from `ClosestPointCost` on the same data.
"""
struct ChamferCost <: ContourCost end

"""
    match_points(p_sim::AbstractMatrix{Float64}, p_obs::AbstractMatrix{Float64})

Match each point of `p_sim` to its nearest point in `p_obs` using KDTree spatial indexing.
One-directional; call it twice with the arguments swapped for a symmetric correspondence.

# Arguments:
- `p_sim::AbstractMatrix{Float64}` : Simulation points matrix of size [2 × n_sim]
- `p_obs::AbstractMatrix{Float64}` : Observation points matrix of size [2 × n_obs]

# Returns:
- `pairs::Matrix{Int64}` : Point pairs array of size [n_sim × 2] where pairs[i,:] = [i, j] matches sim point i to obs point j
"""
function match_points(p_sim::AbstractMatrix{Float64}, p_obs::AbstractMatrix{Float64})::Matrix{Int64}
    simSize = size(p_sim, 2)
    pairs = zeros(Int64, simSize, 2)
    tree = KDTree(p_obs; leafsize=10)
    for sim_counter in 1:simSize
        idx, _ = knn(tree, p_sim[:, sim_counter], 1)
        pairs[sim_counter, :] = [sim_counter, idx[1]]
    end
    return pairs
end

"""
    _as_2xN(x) -> Matrix{Float64}

Coerce one frame of contour points to a 2 × n_points matrix, unwrapping a single-element
`Vector{Matrix}` wrapper and transposing an n_points × 2 layout.
"""
function _as_2xN(x)
    raw = x isa AbstractVector ? x[1] : x   # unwrap Vector{Matrix} wrapper if present
    m = Matrix{Float64}(collect(raw))
    return size(m, 1) == 2 ? m : Matrix{Float64}(m')
end

"""
    _as_dudθ(x) -> Array{Float64,3}

Coerce one frame of contour-point sensitivities to a 2 × n_points × n_params array.
"""
_as_dudθ(x) = Array{Float64}(x isa AbstractVector ? x[1] : x)

"""
    _residuals(cost, sim_t, obs_t, du_tdθ) -> (u, J, u_ref, ndiv, pairs)

Build one frame's least-squares pieces for `cost`: the residual `u`, its Jacobian `J`, the
matched observation coordinates `u_ref` (for the scale-normalized cost), the normalization
denominator `ndiv`, and the raw correspondence(s) `pairs`.

`du_tdθ` is 2 × n_sim × n_params, or `nothing` when only the cost is wanted, in which case
`J` is `nothing` too.

`ndiv` counts residual *components*, not pairs — `length(pairs) == 2·n_pairs`. That
doubling is inherited and kept deliberately: changing it would move every recorded cost
value and every `λ` calibrated against one.
"""
function _residuals(::ClosestPointCost, sim_t::AbstractMatrix{Float64}, obs_t::AbstractMatrix{Float64}, du_tdθ)
    pairs = match_points(sim_t, obs_t)
    s, o = pairs[:, 1], pairs[:, 2]

    u     = [(sim_t[1, s] - obs_t[1, o]); (sim_t[2, s] - obs_t[2, o])]
    u_ref = [obs_t[1, o]; obs_t[2, o]]
    J     = isnothing(du_tdθ) ? nothing : [du_tdθ[1, s, :]; du_tdθ[2, s, :]]

    return u, J, u_ref, length(pairs), pairs
end

function _residuals(::ChamferCost, sim_t::AbstractMatrix{Float64}, obs_t::AbstractMatrix{Float64}, du_tdθ)
    fwd = match_points(sim_t, obs_t)   # fwd[:,1] sim index, fwd[:,2] obs index
    rev = match_points(obs_t, sim_t)   # rev[:,1] obs index, rev[:,2] sim index
    fs, fo = fwd[:, 1], fwd[:, 2]
    ro, rs = rev[:, 1], rev[:, 2]

    u = [(sim_t[1, fs] - obs_t[1, fo]); (sim_t[2, fs] - obs_t[2, fo]);
         (sim_t[1, rs] - obs_t[1, ro]); (sim_t[2, rs] - obs_t[2, ro])]
    u_ref = [obs_t[1, fo]; obs_t[2, fo]; obs_t[1, ro]; obs_t[2, ro]]
    J = isnothing(du_tdθ) ? nothing :
        [du_tdθ[1, fs, :]; du_tdθ[2, fs, :]; du_tdθ[1, rs, :]; du_tdθ[2, rs, :]]

    return u, J, u_ref, length(fwd) + length(rev), (fwd, rev)
end

"""
    contour_cost(sim_frames, obs_frames; outliers=[], cost=ClosestPointCost())

Score simulated contours against observed ones, without derivatives. Used by the line
searches to evaluate a trial step.

# Arguments:
- `sim_frames::AbstractArray` : Vector of simulation point coordinates [2 × n_points for each frame]
- `obs_frames::AbstractArray` : Vector of observation point coordinates [2 × n_points for each frame]

# Keyword Arguments:
- `outliers::AbstractArray` : Frame indices to skip
- `cost::ContourCost` : which cost to evaluate (default: [`ClosestPointCost`](@ref))

# Returns:
- `cost_list::Vector{Float64}` : Mean squared error for each frame
- `[pairsList, norm_cost_list]` : Point correspondences and scale-normalized cost per frame
"""
function contour_cost(sim_frames::AbstractArray, obs_frames::AbstractArray;
                       outliers::AbstractArray=[], cost::ContourCost=ClosestPointCost())
    cost_list = Float64[]
    norm_cost_list = Float64[]
    pairsList = []

    @argcheck length(sim_frames) == length(obs_frames) "Size of the simulation and observation scenes should be the same"
    for (frame_idx, (obs_t, _sim_t)) in enumerate(zip(obs_frames, sim_frames))

        if frame_idx in outliers
            @info "Skipping frame $frame_idx as it is marked as an outlier."
            continue
        end

        sim_t = _as_2xN(_sim_t)
        obs_t = _as_2xN(obs_t)

        u, _, u_ref, ndiv, pairs = _residuals(cost, sim_t, obs_t, nothing)

        tcost = u' * u
        push!(cost_list, tcost / (2 * ndiv))
        push!(norm_cost_list, tcost / (u_ref' * u_ref))
        push!(pairsList, pairs)
    end
    return cost_list, [pairsList, norm_cost_list]
end

"""
    contour_cost(sim_frames, obs_frames, dudθ; outliers=[], cost=ClosestPointCost())

Score simulated contours against observed ones, with the derivatives the Gauss-Newton
optimizers need: per frame, `uᵀu`, `Jᵀu` and the Gauss-Newton Hessian `JᵀJ`.

# Arguments:
- `sim_frames::AbstractArray` : Vector of simulation point coordinates
- `obs_frames::AbstractArray` : Vector of observation point coordinates
- `dudθ::AbstractArray` : Gradient of simulation points w.r.t. parameters [2 × n_points × n_params per frame]

# Keyword Arguments:
- `outliers::AbstractArray` : Frame indices to skip
- `cost::ContourCost` : which cost to differentiate (default: [`ClosestPointCost`](@ref))

# Returns:
- `cost_list::Vector{Float64}` : Mean squared error for each frame
- `dcost_list::Vector{Matrix}` : First derivative (gradient) w.r.t. parameters for each frame
- `dcost2List::Vector{Matrix}` : Second derivative (Gauss-Newton Hessian) for each frame
- `pairsList` : Point correspondences for each frame
"""
function contour_cost(sim_frames::AbstractArray, obs_frames::AbstractArray, dudθ::AbstractArray;
                       outliers::AbstractArray=[], cost::ContourCost=ClosestPointCost())::Tuple{Vector{Float64}, Vector, Vector, Vector}
    cost_list::Vector{Float64} = Float64[]
    dcost_list::Vector = []
    dcost2List::Vector = []
    pairsList::Vector = []

    @argcheck length(sim_frames) == length(obs_frames) "Size of the simulation and observation scenes should be the same"
    @argcheck length(sim_frames) == length(dudθ) "Size of the simulation and observation scenes should be the same"

    for (frame_idx, (obs_t, _sim_t, _du_tdθ)) in enumerate(zip(obs_frames, sim_frames, dudθ))

        if frame_idx in outliers
            @info "Skipping frame $frame_idx as it is marked as an outlier."
            continue
        end

        sim_t  = _as_2xN(_sim_t)
        obs_t  = _as_2xN(obs_t)
        du_tdθ = _as_dudθ(_du_tdθ)

        @argcheck size(sim_t,2) == size(du_tdθ,2) "Number of the border points and the gradient points should be the same"

        mat_nan_inf_check(du_tdθ[:,:,1])
        mat_nan_inf_check(du_tdθ[:,:,2])

        u, Jmat, _, ndiv, pairs = _residuals(cost, sim_t, obs_t, du_tdθ)

        tcost   = u' * u
        dtcost  = Jmat' * u
        dt2cost = Jmat' * Jmat

        mCost   = tcost/(2*ndiv)  # 1/2m(Σ(√(xi-x_obs)^2+(yi-y_obs)^2))^2 (mean error)
        dmcost  = dtcost/ndiv     # 1/m(Σ(xi-x_obs)∂x/∂θ_i +(yi-y_obs)∂y/∂θ_i))
        dm2cost = dt2cost/ndiv    # 1/m(Σ(∂x2/∂2θ_i + ∂y2/∂2θ_i))

        push!(cost_list, mCost)
        push!(dcost_list, dmcost)
        push!(dcost2List, dm2cost)
        push!(pairsList, pairs)
    end
    return cost_list, dcost_list, dcost2List, pairsList
end

"""
    closest_point(args...; kwargs...)

Deprecated alias for [`contour_cost`](@ref), which was renamed because it evaluates
whichever [`ContourCost`](@ref) it is given. Kept for scripts outside this repository;
delete once nothing calls it.
"""
function closest_point(args...; kwargs...)
    Base.depwarn("`closest_point` has been renamed to `contour_cost`, since it evaluates " *
                 "whichever ContourCost it is passed, not only ClosestPointCost.", :closest_point)
    return contour_cost(args...; kwargs...)
end
