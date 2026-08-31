using LinearAlgebra
using Plots
using ArgCheck

# Contour cost functions (`match_points`, `contour_cost`, the `ContourCost` types and their
# residuals) live in `optimization/cost_functions.jl`, included before this file.


"""
    init_cylinder()

Initialize and run a complete optimization cycle for cylinder geometry fitting.

Creates a ground truth cylinder, initializes an optimization model, and iteratively refines
cylinder radius and height using Newton's method until convergence or max iterations.

# Returns:
None. Prints optimization progress and final results to console.
"""
function init_cylinder()::Nothing
    scale::Int = 100
    ne::Int = 4
    camera_matrix::Matrix{Float64} = get_camera_matrix()
    camera_pose::Vector{Float64} = scale * [0 -0.25 2]'

    _box = meshgrid_cuboid(1.0, 1.0, 1.0; mesh_type=:structured, ne=ne, element_shape=:Hex, basis_order=2)
    NodeList = _box.NodeList
    nNodes = 2*ne + 1
    BorderNodes = [_box.side_nodes, _box.bottom_nodes, _box.top_nodes]

    r_gt = 0.25 * scale
    h_gt = 0.5 * scale
    NodeListCyl_gt = _inflate_cylinder(NodeList, -0.5, 0.5, -0.5, 0.5, r_gt, h_gt)
    side_nodes = BorderNodes[1]
    obsBorderPts, _ = extract_borders(NodeListCyl_gt, camera_matrix, camera_pose, BorderNodesList=side_nodes, nNodes)

    # Optimizer
    r = 1 * scale * ones(ne)
    h = 1 * scale
    NodeListCyl, ∇NodeListCyl = _inflate_cylinder(NodeList, -0.5, 0.5, -0.5, 0.5, r, h, GRAD=true)
    simBorderPts, ∇BorderPts2D = extract_borders(NodeListCyl, camera_matrix, camera_pose, nNodes, BorderNodesList=side_nodes, GRAD=true, dqdθ=∇NodeListCyl, )

    d, ∂d, ∂2d, _ = contour_cost([simBorderPts], [obsBorderPts], [∇BorderPts2D])
    totdinit = sum(d) / length(d)
    θ = vcat(r, h)
    iter = 1
    printstyled("Ground truth: r=$(round(r_gt, sigdigits=4)), h=$(round(h_gt, sigdigits=4))\n", color=:green)
    printstyled("Initial error=$(round(totdinit, sigdigits=4))\n", color=:yellow)
    while true
        printstyled("\n================ Iteration $iter ================", color=:blue)
        @debug "Current parameters: r=$(round(θ[1], sigdigits=4)), h=$(round(θ[2], sigdigits=4)), cost=$(round(totdinit, sigdigits=4))"

        t∂2d = zeros(size(∂2d[1]))
        t∂d = zeros(size(∂d[1]))
        for i in 1:length(d)
            t∂2d += ∂2d[i]
            t∂d  += ∂d[i]
        end

        p = t∂2d\t∂d
        @debug "Newton step: [$(round(p[1], sigdigits=4)), $(round(p[2], sigdigits=4))]"
        
        α = 1
        θ = θ - α*p
        val_check(θ)
        
        r = θ[1]
        h = θ[2]

        NodeListCyl, ∇NodeListCyl = _inflate_cylinder(NodeList, -0.5, 0.5, -0.5, 0.5, r, h, GRAD=true)
        simBorderPts, ∇BorderPts2D = extract_borders(NodeListCyl, camera_matrix, camera_pose, nNodes, BorderNodesList=side_nodes, GRAD=true, dqdθ=∇NodeListCyl, )

        d, ∂d, ∂2d, _ = contour_cost([simBorderPts],[obsBorderPts],[∇BorderPts2D])

        totd = sum(d)/length(d)
        c_grad = abs(totdinit - totd)
        c_grad_rel = c_grad / (abs(totdinit) + 1e-12)
        totdinit = totd

        iter = iter + 1
        
        @debug "Result: r = $(round(θ[1], sigdigits=4)), h = $(round(θ[2], sigdigits=4)), cost = $(round(totd, sigdigits=4))"
        @debug "Delta: Δcost = $(round(c_grad, sigdigits=3)) (rel: $(round(c_grad_rel, sigdigits=3)))"

        if c_grad < 0.05
            printstyled("[CONVERGED] Cost change = $(round(c_grad, sigdigits=3))\n", color=:green)
            break
        elseif iter ≥ 100
            printstyled("[WARNING] MAX ITERATIONS (100) REACHED\n", color=:yellow)
            break
        end
    end
    printstyled("\n========== Optimization Complete =========\n", color=:blue)
end

"""
    armijo_line_search(model, scene, conditions, obsBorderPts, θ_prev, p_damped, ∂d_prev, cost_prev; c=1e-4, max_backtracks=10, outliers=Int[], penalty, ∇penalty)

Armijo line search with sufficient descent condition for optimization.

Implements Armijo backtracking to ensure acceptable step sizes that satisfy:
cost(θ - α·p) ≤ cost_prev - c·α·∇cost·p

# Arguments:
- `model::Stokes` : Finite element model
- `scene::SqueezeFlow` : Squeeze flow problem setup
- `conditions::Conditions` : Simulation conditions and settings
- `obsBorderPts::Vector{AbstractArray}` : Observed border points
- `θ_prev::Vector{Float64}` : Previous parameter values
- `p_damped::Vector{Float64}` : Damped Newton step direction
- `∂d_prev::Vector` : Previous gradient of cost
- `cost_prev::Float64` : Previous cost value
- `c::Float64` : Sufficient decrease parameter (default 1e-4, range [1e-4, 0.1])
- `max_backtracks::Int` : Maximum backtracking iterations (default 10)
- `outliers::Vector{Int}` : Frame indices to skip
- `penalty::Function` : Regularization term `R(θ)` added to the data misfit (default: zero)
- `∇penalty::Function` : Gradient `∇R(θ)` of the regularization term (default: zero)

# Returns:
- `θ_trial::Vector{Float64}` : New parameter values
- `d_trial::Vector{Float64}` : Cost values at new point
- `∂d::Vector` : Gradient of cost at new point
- `∂2d::Vector` : Hessian of cost at new point
- `cost_trial::Float64` : Cost value at trial point (data misfit + `penalty`)
- `accepted::Bool` : Whether step was accepted
- `simBorderPts` : Simulated border points at the returned parameters by Armijo condition
- `simBorderPts` : Simulated border points at the returned parameters

The Armijo test uses the same objective that produced `p_damped`: pass the same
`penalty`/`∇penalty` here that were used to build the regularized gradient and Hessian.
"""
function armijo_line_search(model::Stokes, scene::SqueezeFlow, conditions::Conditions, obsBorderPts::Vector{AbstractArray},
                            θ_prev::Vector{Float64}, p_damped::Vector{Float64}, ∂d_prev::Vector, cost_prev::Float64;
                            c::Float64=1e-4, max_backtracks::Int=10, outliers::Vector{Int}=Int[],
                            penalty::Function = _ -> 0.0, ∇penalty::Function = θ -> zero(θ),
                            cost::ContourCost=ClosestPointCost())
    # Armijo line search with sufficient descent condition.
    # Accepts step if: cost(θ - α·p) ≤ cost_prev - c·α·∇f·p
    # f is the regularized objective: mean data misfit + penalty(θ).
    len_d = length(∂d_prev)
    α = 1.0
    t_grad = sum(∂d_prev) / len_d + ∇penalty(θ_prev)   # ∇f of the same objective that produced p
    grad_prod = dot(t_grad, p_damped)     # directional derivative ∇fᵀ·p
    
    # Initialize variables for scope after loop
    θ_trial = θ_prev - α * p_damped
    d_trial = similar(∂d_prev[1])
    ∂d = similar(∂d_prev)
    ∂2d = zeros(0)  # Will be overwritten
    cost_trial = cost_prev
    simBorderPts_accepted = nothing
    
    for backtrack_iter in 1:max_backtracks
        θ_trial = θ_prev - α * p_damped
        val_check(θ_trial)
        
        model.η = [θ_trial[1]]
        scene.β = [θ_trial[2]]
        μ_list, gradList, simBorderPts, _, _, _, _, _, _, _ = simulate(model, scene, conditions)
        
        # Use cost-only version during line search (fast)
        simBorderPts_accepted = simBorderPts

        d_trial_costs, _ = contour_cost(simBorderPts, obsBorderPts, outliers=outliers, cost=cost)
        cost_trial = sum(d_trial_costs) / len_d + penalty(θ_trial)
        
        # Armijo condition: cost_trial ≤ cost_prev - c·α·∇f·p
        Δcost = cost_prev - cost_trial
        sufficient_decrease = c * α * grad_prod
        
        status = Δcost ≥ sufficient_decrease ? "[ACCEPT]" : "[REJECT]"
        @debug "      Trial α    | Cost         | Δcost        | Required Decrease | Status"
        @debug "      $("-"^73)"
        @debug "      $(round(α, sigdigits=5)) | $(round(cost_trial, sigdigits=3)) | $(round(Δcost, sigdigits=3)) | $(round(sufficient_decrease, sigdigits=3)) | $status"
        
        if cost_trial ≤ cost_prev - sufficient_decrease
            # Compute full gradients/Hessians only for accepted step
            d_trial, ∂d, ∂2d, pairs = contour_cost(simBorderPts, obsBorderPts, gradList, outliers=outliers, cost=cost)
            @debug "      $("-"^73)"
            printstyled("      [ACCEPT] Step accepted: α = $(round(α, sigdigits=5)), cost = $(round(cost_trial, sigdigits=4))\n", color=:green)
            return θ_trial, d_trial, ∂d, ∂2d, cost_trial, true, simBorderPts_accepted
        end
        
        α *= 0.5  # Halve step size and retry
    end
    
    # Fallback: compute gradients/Hessians for final trial (even if rejected)
    model.η = [θ_trial[1]]
    scene.β = [θ_trial[2]]
    μ_list, gradList, simBorderPts, _, _, _, _, _, _, _ = simulate(model, scene, conditions)
    d_trial, ∂d, ∂2d, pairs = contour_cost(simBorderPts, obsBorderPts, gradList, outliers=outliers, cost=cost)
    simBorderPts_accepted = simBorderPts
    
    @debug "      $("-"^73)"
    printstyled("      [WARNING] Max backtracks reached. Accepting last trial (α=$(round(α, sigdigits=5)), cost=$(round(cost_trial, sigdigits=4)))\n", color=:yellow)
    return θ_trial, d_trial, ∂d, ∂2d, cost_trial, false, simBorderPts_accepted
end

"""
    backtrack_line_search(model, scene, conditions, obsBorderPts, θ_prev, p_damped, cost_prev; outliers=Int[], penalty)

Simple backtracking line search that tries full step, then half steps.

A basic line search method for reference and comparison. Attempts α=1.0, then α=0.5 if cost didn't improve.
No principled descent guarantee like Armijo condition.

# Arguments:
- `model::Stokes` : Finite element model
- `scene::SqueezeFlow` : Squeeze flow problem setup
- `conditions::Conditions` : Simulation conditions and settings  
- `obsBorderPts::Vector{AbstractArray}` : Observed border points
- `θ_prev::Vector{Float64}` : Previous parameter values
- `p_damped::Vector{Float64}` : Damped Newton step direction
- `cost_prev::Float64` : Previous cost value
- `outliers::Vector{Int}` : Frame indices to skip
- `penalty::Function` : Regularization term `R(θ)` added to the data misfit (default: zero)

# Returns:
- `θ_trial::Vector{Float64}` : New parameter values
- `d_trial::Vector{Float64}` : Cost values at new point
- `∂d::Vector` : Gradient of cost at new point
- `∂2d::Vector` : Hessian of cost at new point
- `cost_trial::Float64` : Cost value at trial point
- `accepted::Bool` : Whether step was accepted
- `simBorderPts` : Simulated border points at the returned parameters
"""
function backtrack_line_search(model::Stokes, scene::SqueezeFlow, conditions::Conditions, obsBorderPts::Vector{AbstractArray},
                               θ_prev::Vector{Float64}, p_damped::Vector{Float64}, cost_prev::Float64; outliers::Vector{Int}=Int[],
                               penalty::Function = _ -> 0.0, cost::ContourCost=ClosestPointCost())
    # Simple backtracking line search (original method - kept for reference/comparison).
    # Takes full step, then half step if needed. No principled descent guarantee.
    len_d_calc = 0
    
    @debug "      Strategy: Try α=1.0 (full), then α=0.5 (half)"
    
    # Try full step
    @debug "      Trying α = 1.0 (full step)..."
    θ_trial = θ_prev - p_damped
    val_check(θ_trial)
    model.η = [θ_trial[1]]
    scene.β = [θ_trial[2]]
    μ_list, gradList, simBorderPts, _, _, _, _, _, _, _ = simulate(model, scene, conditions)
    
    # Use cost-only version first (fast)
    d_trial_costs, _ = contour_cost(simBorderPts, obsBorderPts, outliers=outliers, cost=cost)
    len_d_calc = length(d_trial_costs)
    cost_trial = sum(d_trial_costs) / len_d_calc + penalty(θ_trial)
    Δcost_full = cost_prev - cost_trial
    @debug "        cost = $(round(cost_trial, sigdigits=4)), Δcost = $(round(Δcost_full, sigdigits=4))"
    
    if cost_trial < cost_prev
        # Compute full gradients/Hessians only for accepted step
        d_trial, ∂d, ∂2d, pairs = contour_cost(simBorderPts, obsBorderPts, gradList, outliers=outliers, cost=cost)
        printstyled("        [ACCEPT] Cost decreased, accepting\n", color=:green)
        return θ_trial, d_trial, ∂d, ∂2d, cost_trial, true, simBorderPts
    end
    
    # Try half step
    @debug "      Trying α = 0.5 (half step)..."
    θ_trial = θ_prev - 0.5 * p_damped
    val_check(θ_trial)
    model.η = [θ_trial[1]]
    scene.β = [θ_trial[2]]
    μ_list, gradList, simBorderPts, _, _, _, _, _, _, _ = simulate(model, scene, conditions)
    
    # Use cost-only version first (fast)
    d_trial_costs, _ = contour_cost(simBorderPts, obsBorderPts, outliers=outliers, cost=cost)
    cost_trial_half = sum(d_trial_costs) / len_d_calc + penalty(θ_trial)
    Δcost_half = cost_prev - cost_trial_half
    @debug "        cost = $(round(cost_trial_half, sigdigits=4)), Δcost = $(round(Δcost_half, sigdigits=4))"
    
    # Compute full gradients/Hessians for final result (accepted or rejected)
    d_trial, ∂d, ∂2d, pairs = contour_cost(simBorderPts, obsBorderPts, gradList, outliers=outliers, cost=cost)
    
    if cost_trial_half < cost_prev
        printstyled("        [ACCEPT] Cost decreased, accepting\n", color=:green)
        return θ_trial, d_trial, ∂d, ∂2d, cost_trial_half, true, simBorderPts
    else
        printstyled("        [WARNING] Neither step decreased cost. Accepting half step anyway.\n", color=:yellow)
        return θ_trial, d_trial, ∂d, ∂2d, cost_trial_half, false, simBorderPts
    end
end


"""
    _fit_model_GN(model, scene, conditions, obsBorderPts, θ; outliers, line_search_method)
        -> (ηpList, βpList, cost_list, iterList)

Run the Gauss-Newton parameter identification loop for a Stokes model.

Internal implementation called by `fit_model` when `method=:GN`. Updates `model.η`
and `scene.β` in-place and returns per-iteration histories.

# Arguments
- `model::Stokes`: mutable Stokes FEM model (reset at each iteration).
- `scene::SqueezeFlow`: squeeze-flow scenario holding boundary conditions and timing.
- `conditions::Conditions`: observation conditions (camera matrix and object pose).
- `obsBorderPts::Vector{AbstractArray}`: observed 2D contour points per time step.
- `θ::Vector{Float64}`: initial parameter vector `[η, β]`.

# Keyword Arguments
- `outliers::Vector{Int}`: frame indices to exclude from the cost (default: `[]`).
- `line_search_method::Symbol`: `:armijo` or `:backtrack` (default: `:armijo`).

# Returns
- `ηpList::Vector{Float64}`: viscosity history.
- `βpList::Vector{Float64}`: friction coefficient history.
- `cost_list::Vector{Float64}`: cost per iteration.
- `iterList::Vector{Float64}`: iteration indices.
"""
function _fit_model_GN(model::Stokes, scene::SqueezeFlow, conditions::Conditions, obsBorderPts::Vector{AbstractArray}, θ::Vector{Float64};
                   outliers::Vector{Int}=Int[], line_search_method::Symbol=:armijo,
                   store_border_pts::Bool=false, cost::ContourCost=ClosestPointCost())
    reset_model!(model)
    model.η = [θ[1]]
    scene.β = [θ[2]]

    ηpList::Vector{Float64} = Float64[]
    βpList::Vector{Float64} = Float64[]
    cost_list::Vector{Float64} = Float64[]
    iterList::Vector{Float64} = Float64[]
    # Simulated contours at each accepted iterate; index-aligned with ηpList/cost_list.
    # Opt-in: one full contour set per iteration is large (frames × 2 × border points).
    simBorderPtsList::Vector{Any} = Any[]
    
    printstyled("Initializing simulation with η: $(round(θ[1], sigdigits=4)), β: $(round(θ[2], sigdigits=4))\n", color=:cyan)
    μ_list, gradList, simBorderPts, _, _, _, _, _, _, _ = simulate(model, scene, conditions)

    @debug begin
        path = joinpath(get_scratch_dir(), "optimization_animation")
        animate_fields(filepath=path, sim_border_nodes_2d=simBorderPts, obs_border_nodes_2d=obsBorderPts)
        "Saved optimization animation to $path"
    end

    d, ∂d, ∂2d, _ = contour_cost(simBorderPts, obsBorderPts, gradList, outliers=outliers, cost=cost)
    totdinit::Float64 = sum(d)/length(d)

    push!(ηpList,θ[1])
    push!(βpList,θ[2])
    push!(cost_list,totdinit)
    push!(iterList,1)
    store_border_pts && push!(simBorderPtsList, simBorderPts)

    iter::Int = 1
    c_grad::Float64 = 1.0
    len_d::Int = length(d)
    t∂2d::Matrix{Float64} = zeros(size(∂2d[1]))
    t∂d::Vector{Float64} = zeros(size(∂d[1]))
    printstyled("Initial error = $totdinit\n", color=:yellow)
    while true
        printstyled("\n==================== Iteration $iter ===================\n", color=:blue)

        reset_model!(model)
        t∂2d .= 0.0
        t∂d .= 0.0

        @debug "Current parameters: η = $(round(θ[1], sigdigits=4)), β = $(round(θ[2], sigdigits=4)), cost = $(round(totdinit, sigdigits=4))"
        for i::Int in 1:len_d
            t∂2d = t∂2d + ∂2d[i]
            t∂d = t∂d + ∂d[i]
        end

        # Hessian diagnostics and damping
        condition_number::Float64 = cond(t∂2d)
        H_η::Float64 = abs(t∂2d[1, 1])
        H_β::Float64 = abs(t∂2d[2, 2])
        H_ratio::Float64 = H_η / (H_β + 1e-12)
        
        
        if condition_number > 1e2 || H_ratio > 1e2
            printstyled("  [WARNING] Hessian: κ=$(round(condition_number, sigdigits=3)), H_β=$H_β, H_η=$H_η, ratio=$(round(H_ratio, sigdigits=3))\n", color=:yellow)
            reg_param = 1e-6 * norm(t∂2d, 2)
            t∂2d_reg = t∂2d + reg_param * I
            p = t∂2d_reg \ t∂d
            printstyled("  [DAMPED] Using regularized Hessian for step\n", color=:yellow)
        else
            p = t∂2d \ t∂d
        end

        # Compute Newton step
        @debug "Newton step (damped): [$(round(p[1], sigdigits=4)), $(round(p[2], sigdigits=4))]"

        θ_prev::Vector{Float64} = copy(θ)
        cost_prev::Float64 = totdinit
        
        if line_search_method == :armijo
            θ, d, ∂d, ∂2d, totd, _, simBorderPts = armijo_line_search(model, scene, conditions, obsBorderPts, θ_prev, p, ∂d, cost_prev, outliers=outliers, cost=cost)
        else
            θ, d, ∂d, ∂2d, totd, _, simBorderPts = backtrack_line_search(model, scene, conditions, obsBorderPts, θ_prev, p, cost_prev, outliers=outliers, cost=cost)
        end
        
        c_grad = abs(cost_prev - totd)
        c_grad_rel = c_grad / (abs(cost_prev) + 1e-12)
        totdinit = totd
        
        Δη_rel = abs(θ[1] - θ_prev[1]) / (abs(θ_prev[1]) + 1e-12)
        Δβ_rel = abs(θ[2] - θ_prev[2]) / (abs(θ_prev[2]) + 1e-12)

        iter = iter + 1
        push!(ηpList, θ[1])
        push!(βpList, θ[2])
        push!(cost_list, totd)
        push!(iterList, iter)
        store_border_pts && push!(simBorderPtsList, simBorderPts)
        @debug "Result: η = $(round(θ[1], sigdigits=4)), β = $(round(θ[2], sigdigits=4)), cost = $(round(totd, sigdigits=4))"
        @debug "Deltas: Δη/η = $(round(Δη_rel, sigdigits=3)), Δβ/β = $(round(Δβ_rel, sigdigits=3)), Δcost = $(round(c_grad, sigdigits=3)) (rel: $(round(c_grad_rel, sigdigits=3)))"

        if c_grad_rel < 1e-3 && c_grad < 1e-3
            printstyled("[CONVERGED] Relative cost change = $(round(c_grad_rel, sigdigits=3))\n 
                        η = $(round(θ[1], sigdigits=4)), β = $(round(θ[2], sigdigits=4))", color=:green)
            break
        end
        
        if iter ≥ 100
            printstyled("[WARNING] MAX ITERATIONS (100) REACHED\n", color=:yellow)
            break
        end
    end
    
    printstyled("\n========= Optimization Complete =========\n", color=:blue)

    stats = Dict(
        "η" => θ[1],
        "β" => θ[2],
        "ηList" => ηpList,
        "βList" => βpList,
        "cost_list" => cost_list,
        "iterList" => iterList,
        "simBorderPtsList" => simBorderPtsList,
    )
    return stats
end

"""
    _fit_model_GN_tikhonov(model, scene, conditions, obsBorderPts, θ; outliers, λ_scale, θ_p, Γ, line_search_method)

Gauss-Newton parameter identification with Tikhonov regularization toward a prior `θ_p`.

Minimizes `mean_data_misfit(θ) + R(θ)` with
`R(θ) = ½·(θ - θ_p)ᵀ ΓᵀΛΓ (θ - θ_p)`, where `Λ = Diagonal(λ)`.

`λ` holds **one entry per parameter**, so η and β can be regularized at different
strengths — useful when one of them is far better constrained by the data than the other.
A scalar `λ_scale` is broadcast to every parameter, which recovers the usual single-weight
Tikhonov penalty. It is fixed before the first iteration and held there for the whole fit:
recomputing it per iteration would make the objective move under the line search, and the
convergence test would stop meaning anything.

`λ` is an **absolute** weight, not a multiple of the data misfit, so it does not transfer
between experiments with different units, noise levels or window lengths on its own — the
first-iteration prior/data precision diagnostic is the thing to read when choosing it.

The default is L2 shrinkage toward the origin — `θ_p = 0` with `Γ = diag(1 ./ |θ₀|)`, the
initial parameters — so `R(θ) = ½·[λ₁·(η/η₀)² + λ₂·(β/β₀)²]`. η and β are assumed
correlated, which leaves the data misfit near-flat along one direction; the penalty
resolves that degeneracy by selecting the smallest-magnitude point along the valley
instead of letting the fit drift.

Normalizing by the initial guess is what makes the two shrink comparably: an identity `Γ`
measures magnitude in raw units, so with η and β carrying different units and differing by
orders of magnitude it would crush the numerically larger one and leave the other
essentially free. It also makes `λ` dimensionless and directly comparable between the two
parameters, since `‖Γθ₀‖² = length(θ)`. Pass `Γ = Matrix{Float64}(I, 2, 2)` for literal
identity weighting, or a non-zero `θ_p` to regularize toward a prior estimate (e.g. the
previous time window) rather than toward zero.

# Arguments
- `model::Stokes` : mutable Stokes FEM model (reset each iteration).
- `scene::SqueezeFlow` : squeeze-flow scenario holding boundary conditions and timing.
- `conditions::Conditions` : observation conditions (camera matrix and object pose).
- `obsBorderPts::Vector{AbstractArray}` : observed 2D contour points per time step.
- `θ::Vector{Float64}` : initial parameter vector `[η, β]`.

# Keyword Arguments
- `outliers::Vector{Int}` : frame indices excluded from the cost (default: `Int[]`).
- `λ_scale::Union{Real,AbstractVector{<:Real}}` : prior strength, either one entry per
  parameter or a scalar broadcast to all of them (default: `1.0`). Entries must be finite
  and non-negative; all-zero disables regularization and recovers plain Gauss-Newton, and a
  single zero entry leaves that one parameter unregularized.
- `θ_p::Vector{Float64}` : prior mean (default: `zeros(length(θ))` — shrinkage toward the
  origin, i.e. prefer the smallest parameters that still fit the data).
- `Γ::Matrix{Float64}` : penalty shaping matrix (default: `diag(1 ./ |θ₀|)`, i.e. L2 on the
  parameters relative to their initial values). Must be finite and square.
- `line_search_method::Symbol` : `:armijo` or `:backtrack` (default: `:armijo`).

# Returns
- `stats::Dict` with keys `"η"`, `"β"`, `"ηList"`, `"βList"`, `"cost_list"`, `"iterList"`,
  plus `"H"` (converged regularized Gauss-Newton Hessian) and `"λ"` (the frozen
  per-parameter weight vector used).
  `"cost_list"` holds the regularized objective, so values are not comparable across fits
  with different `λ`.
"""
function _fit_model_GN_tikhonov(model::Stokes, scene::SqueezeFlow, conditions::Conditions, obsBorderPts::Vector{AbstractArray}, θ::Vector{Float64};
                              outliers::Vector{Int}=Int[], λ_scale::Union{Real,AbstractVector{<:Real}}=1.0,
                              θ_p::Vector{Float64}=zeros(length(θ)),
                              Γ::Matrix{Float64}=Matrix{Float64}(I, length(θ), length(θ)) ./ abs.(θ),
                              line_search_method::Symbol=:armijo, store_border_pts::Bool=false,
                              cost::ContourCost=ClosestPointCost())
    @argcheck length(θ_p) == length(θ) "Prior mean θ_p must match the length of θ"
    @argcheck all(isfinite, θ_p) "Prior mean θ_p must be finite"
    @argcheck size(Γ) == (length(θ), length(θ)) "Penalty matrix Γ must be square of side length(θ)"
    @argcheck all(isfinite, Γ) "Penalty matrix Γ must be finite (the default scaling needs all(θ .!= 0))"
    @argcheck λ_scale isa Real || length(λ_scale) == length(θ) "λ_scale must be a scalar or one entry per parameter"
    @argcheck all(≥(0), λ_scale) "λ_scale entries must be non-negative"
    @argcheck all(isfinite, λ_scale) "λ_scale entries must be finite"

    reset_model!(model)
    model.η = [θ[1]]
    scene.β = [θ[2]]

    ηpList::Vector{Float64} = Float64[]
    βpList::Vector{Float64} = Float64[]
    cost_list::Vector{Float64} = Float64[]
    iterList::Vector{Float64} = Float64[]
    # Simulated contours at each accepted iterate; index-aligned with ηpList/cost_list.
    # Opt-in: one full contour set per iteration is large (frames × 2 × border points).
    simBorderPtsList::Vector{Any} = Any[]
    
    printstyled("Initializing simulation with η: $(round(θ[1], sigdigits=4)), β: $(round(θ[2], sigdigits=4))\n", color=:cyan)
    μ_list, gradList, simBorderPts, _, _, _, _, _, _, _ = simulate(model, scene, conditions)

    @debug begin
        path = joinpath(get_scratch_dir(), "optimization_animation")
        animate_fields(filepath=path, sim_border_nodes_2d=simBorderPts, obs_border_nodes_2d=obsBorderPts)
        "Saved optimization animation to $path"
    end

    d, ∂d, ∂2d, _ = contour_cost(simBorderPts, obsBorderPts, gradList, outliers=outliers, cost=cost)
    data_cost_init::Float64 = sum(d)/length(d)

    # One λ per parameter, so η and β can be regularized at different strengths. A scalar
    # λ_scale means "the same strength for all of them". Held fixed for the whole fit:
    # recomputing it per iteration would make the objective move under the line search and
    # the convergence test would stop meaning anything.
    λ::Vector{Float64} = λ_scale isa Real ? fill(float(λ_scale), length(θ)) :
                                            Vector{Float64}(λ_scale)

    # Tikhonov penalty R(θ) = ½·(θ - θ_p)ᵀ ΓᵀΛΓ (θ - θ_p) with Λ = Diagonal(λ), i.e.
    # ½·Σᵢ λᵢ·[Γ(θ - θ_p)]ᵢ². Defined once so that the cost, the Newton step and the line
    # search all refer to the same objective. Λ sits between Γᵀ and Γ rather than
    # multiplying through, so a per-parameter λ weights the penalty in the basis Γ maps
    # into — with the default diagonal Γ that is exactly λᵢ·((θᵢ - θ_pᵢ)/θ₀ᵢ)².
    ΓᵀΛΓ::Matrix{Float64} = Γ' * Diagonal(λ) * Γ
    R(v::Vector{Float64})  = 0.5 * dot(v - θ_p, ΓᵀΛΓ * (v - θ_p))
    ∇R(v::Vector{Float64}) = ΓᵀΛΓ * (v - θ_p)
    ∇²R::Matrix{Float64}   = ΓᵀΛΓ

    printstyled("Prior: θ_p = $(round.(θ_p, sigdigits=4)), λ = $(round.(λ, sigdigits=4)) " *
                "(initial data misfit $(round(data_cost_init, sigdigits=4)))\n", color=:cyan)

    totdinit::Float64 = data_cost_init + R(θ)

    push!(ηpList,θ[1])
    push!(βpList,θ[2])
    push!(cost_list,totdinit)
    push!(iterList,1)
    store_border_pts && push!(simBorderPtsList, simBorderPts)

    iter::Int = 1
    c_grad::Float64 = 1.0
    len_d::Int = length(d)
    t∂2d::Matrix{Float64} = zeros(size(∂2d[1]))
    t∂d::Vector{Float64} = zeros(size(∂d[1]))
    t∂2d_reg::Matrix{Float64} = zeros(size(∂2d[1]))
    printstyled("Initial error = $totdinit\n", color=:yellow)
    while true
        printstyled("\n==================== Iteration $iter ===================\n", color=:blue)

        reset_model!(model)
        t∂2d .= 0.0
        t∂d .= 0.0

        @debug "Current parameters: η = $(round(θ[1], sigdigits=4)), β = $(round(θ[2], sigdigits=4)), cost = $(round(totdinit, sigdigits=4))"
        for i::Int in 1:len_d
            t∂2d = t∂2d + ∂2d[i]
            t∂d = t∂d + ∂d[i]
        end

        # Hessian diagnostics and damping
        condition_number::Float64 = cond(t∂2d)
        H_η::Float64 = abs(t∂2d[1, 1])
        H_β::Float64 = abs(t∂2d[2, 2])
        H_ratio::Float64 = H_η / (H_β + 1e-12)
          
        if condition_number > 1e2 || H_ratio > 1e2
            printstyled("  [WARNING] Hessian: κ=$(round(condition_number, sigdigits=3)), H_β=$H_β, H_η=$H_η, ratio=$(round(H_ratio, sigdigits=3))\n", color=:yellow)
        end

        t∂d_reg = t∂d/len_d + ∇R(θ)
        t∂2d_reg = t∂2d/len_d + ∇²R

        if iter == 1 && any(>(0), λ)
            # How hard the prior pulls relative to what the data actually constrains.
            # A flat direction (tiny data curvature) is dominated by even a small penalty,
            # which the cost-vs-penalty comparison alone does not reveal.
            prior_ratio = diag(∇²R) ./ max.(abs.(diag(t∂2d/len_d)), 1e-300)
            dom = maximum(prior_ratio) > 1 ? "  [prior dominates]" : ""
            printstyled("  Prior/data precision: η = $(round(prior_ratio[1], sigdigits=3)), " *
                        "β = $(round(prior_ratio[2], sigdigits=3))$dom\n", color=:cyan)
        end

        p = t∂2d_reg \ t∂d_reg

        # Compute Newton step
        @debug "Newton step (damped): [$(round(p[1], sigdigits=4)), $(round(p[2], sigdigits=4))]"

        θ_prev::Vector{Float64} = copy(θ)
        cost_prev::Float64 = totdinit
        
        if line_search_method == :armijo
            θ, d, ∂d, ∂2d, totd, _, simBorderPts = armijo_line_search(model, scene, conditions, obsBorderPts, θ_prev, p, ∂d, cost_prev,
                                                        outliers=outliers, penalty=R, ∇penalty=∇R, cost=cost)
        else
            θ, d, ∂d, ∂2d, totd, _, simBorderPts = backtrack_line_search(model, scene, conditions, obsBorderPts, θ_prev, p, cost_prev,
                                                           outliers=outliers, penalty=R, cost=cost)
        end
        
        # totd already includes R(θ): the line search minimizes the regularized objective.
        c_grad = abs(cost_prev - totd)
        c_grad_rel = c_grad / (abs(cost_prev) + 1e-12)
        totdinit = totd
        
        Δη_rel = abs(θ[1] - θ_prev[1]) / (abs(θ_prev[1]) + 1e-12)
        Δβ_rel = abs(θ[2] - θ_prev[2]) / (abs(θ_prev[2]) + 1e-12)

        iter = iter + 1
        push!(ηpList, θ[1])
        push!(βpList, θ[2])
        push!(cost_list, totd)
        push!(iterList, iter)
        store_border_pts && push!(simBorderPtsList, simBorderPts)
        @debug "Result: η = $(round(θ[1], sigdigits=4)), β = $(round(θ[2], sigdigits=4)), cost = $(round(totd, sigdigits=4))"
        @debug "Deltas: Δη/η = $(round(Δη_rel, sigdigits=3)), Δβ/β = $(round(Δβ_rel, sigdigits=3)), Δcost = $(round(c_grad, sigdigits=3)) (rel: $(round(c_grad_rel, sigdigits=3)))"

        if c_grad_rel < 1e-3 && c_grad < 1e-3
            printstyled("[CONVERGED] Relative cost change = $(round(c_grad_rel, sigdigits=3))\n 
                        η = $(round(θ[1], sigdigits=4)), β = $(round(θ[2], sigdigits=4))", color=:green)
            break
        end
        
        if iter ≥ 100
            printstyled("[WARNING] MAX ITERATIONS (100) REACHED\n", color=:yellow)
            break
        end
    end
    
    printstyled("\n========= Optimization Complete =========\n", color=:blue)

    stats = Dict(
        "η" => θ[1],
        "β" => θ[2],
        "ηList" => ηpList,
        "βList" => βpList,
        "cost_list" => cost_list,
        "iterList" => iterList,
        "simBorderPtsList" => simBorderPtsList,
        "H" => t∂2d_reg,
        "λ" => λ,
    )
    return stats
end

"""
    fit_model_LM(model::Stokes, scene::SqueezeFlow, conditions::Conditions,
                 obsBorderPts::Vector{AbstractArray}, θ::Vector{Float64};
                 outliers::Vector{Int}=Int[], λ::Float64=1e-3)

Optimize model parameters using Levenberg-Marquardt algorithm with adaptive damping.

# Arguments
- `model::Stokes` : Stokes flow model to optimize
- `scene::SqueezeFlow` : Squeeze flow scene configuration
- `conditions::Conditions` : Boundary and initial conditions
- `obsBorderPts::Vector{AbstractArray}` : Observed border point coordinates
- `θ::Vector{Float64}` : Initial parameters [η, β]
- `outliers::Vector{Int}` : Frame indices to skip (default: [])
- `λ::Float64` : Initial damping parameter (default: 1e-3)

# Returns
- `stats::Dict` : Optimization results with keys "η", "β", "ηList", "βList", "cost_list", "iterList"

# Details
Implements adaptive trust-region Levenberg-Marquardt with gain ratio ρ for step acceptance.
Decreases λ on accepted steps (→Newton), increases λ on rejected steps (→gradient descent).
Stops on convergence (relative cost change < 1e-4) or max iterations (100).
"""
function _fit_model_LM(model::Stokes, scene::SqueezeFlow, conditions::Conditions, obsBorderPts::Vector{AbstractArray}, θ::Vector{Float64};
                      outliers::Vector{Int}=Int[], λ::Float64=1e-3, cost::ContourCost=ClosestPointCost())
                      
    reset_model!(model)
    model.η = [θ[1]]
    scene.β = [θ[2]]

    ηpList::Vector{Float64} = Float64[]
    βpList::Vector{Float64} = Float64[]
    cost_list::Vector{Float64} = Float64[]
    iterList::Vector{Float64} = Float64[]
    
    printstyled("Initializing simulation with η: $(round(θ[1], sigdigits=4)), β: $(round(θ[2], sigdigits=4))\n", color=:cyan)
    μ_list, gradList, simBorderPts, _, _, _, _, _, _, _ = simulate(model, scene, conditions)
    d, ∂d, ∂2d, _ = contour_cost(simBorderPts, obsBorderPts, gradList, outliers=outliers, cost=cost)
    cost_prev::Float64 = sum(d)/length(d)

    push!(ηpList, θ[1])
    push!(βpList, θ[2])
    push!(cost_list, cost_prev)
    push!(iterList, 1)

    len_d::Int = length(d)
    iter::Int = 1
    t∂2d::Matrix{Float64} = zeros(size(∂2d[1]))
    t∂d::Vector{Float64} = zeros(size(∂d[1]))
    printstyled("Initial error = $cost_prev\n", color=:yellow)
    while true
        printstyled("\n==================== Iteration $iter ===================\n", color=:blue)

        reset_model!(model)
        @debug "Current parameters: η = $(round(θ[1], sigdigits=4)), β = $(round(θ[2], sigdigits=4)), cost = $(round(cost_prev, sigdigits=4))"
        
        t∂2d .= 0.0
        t∂d .= 0.0
        for i::Int in 1:len_d
            t∂2d = t∂2d + ∂2d[i]
            t∂d = t∂d + ∂d[i]
        end
        if iter == 1
            λ = 1e-3 * maximum(diag(t∂2d))
        end

        p::Vector{Float64} = (t∂2d + λ * Diagonal(diag(t∂2d))) \ t∂d
        θ_prev::Vector{Float64} = copy(θ)
        
        # Execute and evaluate step
        θ = θ - p
        θ = val_check(θ)
        
        reset_model!(model)
        model.η = [θ[1]]
        scene.β = [θ[2]]
        _, gradList, simBorderPts, _, _, _, _, _, _, _ = simulate(model, scene, conditions)
        d, ∂d, ∂2d, _ = contour_cost(simBorderPts, obsBorderPts, gradList, outliers=outliers, cost=cost)
        totd::Float64 = sum(d) / len_d
        
        c_grad::Float64 = abs(cost_prev - totd)
        c_grad_rel::Float64 = c_grad / (abs(cost_prev) + 1e-12)
        
        # Gain ratio: actual vs predicted cost reduction
        ρ::Float64 = (cost_prev - totd) / (0.5 * p' * (λ * (diag(t∂2d) .* p) + t∂d) + 1e-12)
        
        if ρ > 0
            λ *= max(1 / 3, 1 - (2 * ρ - 1)^3)
            λ = max(λ, 1e-7)
            printstyled("  [ACCEPT] ρ=$(round(ρ, sigdigits=3)), λ→$(round(λ, sigdigits=3))\n", color=:green)
            cost_prev = totd
        else
            λ *= 2
            λ = min(λ, 1e7)
            θ = θ_prev
            totd = cost_prev
            printstyled("  [REJECT] ρ=$(round(ρ, sigdigits=3)), λ→$(round(λ, sigdigits=3))\n", color=:yellow)
        end

        Δη_rel::Float64 = abs(θ[1] - θ_prev[1]) / (abs(θ_prev[1]) + 1e-12)
        Δβ_rel::Float64 = abs(θ[2] - θ_prev[2]) / (abs(θ_prev[2]) + 1e-12)

        if ρ > 0
            iter = iter + 1
            push!(ηpList, θ[1])
            push!(βpList, θ[2])
            push!(cost_list, totd)
            push!(iterList, iter)
            @debug "Result: η = $(round(θ[1], sigdigits=4)), β = $(round(θ[2], sigdigits=4)), cost = $(round(totd, sigdigits=4))"
            @debug "Deltas: Δη/η = $(round(Δη_rel, sigdigits=3)), Δβ/β = $(round(Δβ_rel, sigdigits=3)), Δcost = $(round(c_grad, sigdigits=3)) (rel: $(round(c_grad_rel, sigdigits=3)))"
            
            if c_grad_rel < 1e-4
                printstyled("[CONVERGED] Relative cost change = $(round(c_grad_rel, sigdigits=3))\n", color=:green)
                break
            end
            if iter ≥ 100
                printstyled("[WARNING] MAX ITERATIONS (100) REACHED\n", color=:yellow)
                break
            end
        else
            @debug "Retrying with λ = $(round(λ, sigdigits=3))..."
            continue
        end
    end
    
    printstyled("\n========= Optimization Complete =========\n", color=:blue)

    stats = Dict(
        "η" => θ[1],
        "β" => θ[2],
        "ηList" => ηpList, 
        "βList" => βpList,
        "cost_list" => cost_list,
        "iterList" => iterList
    )
    return stats
end

"""
    fit_model(model::Stokes, scene::SqueezeFlow, conditions::Conditions,
              obsBorderPts::Vector{AbstractArray}, θ::Vector{Float64};
              outliers::Vector{Int}=Int[], method::Symbol=:gn, line_search_method::String="armijo")

Unified interface for parameter optimization with multiple algorithm selection.

Dispatches to appropriate optimization method based on `method` parameter.

# Arguments
- `model::Stokes` : Stokes flow model to optimize
- `scene::SqueezeFlow` : Squeeze flow scene configuration
- `conditions::Conditions` : Boundary and initial conditions
- `obsBorderPts::Vector{AbstractArray}` : Observed border point coordinates
- `θ::Vector{Float64}` : Initial parameters [η, β]
- `outliers::Vector{Int}` : Frame indices to skip (default: [])
- `method::Symbol` : Optimization algorithm (default: `:gn`)
  - `:gn` - Gauss-Newton with damping and line search
  - `:gn_tikhonov` - Gauss-Newton regularized toward a prior `θ_p`
  - `:lm` - Levenberg-Marquardt with adaptive damping
- `kwargs...` : forwarded verbatim to the chosen implementation, which owns their defaults.
  Passing one the method does not accept raises a `MethodError` rather than being ignored,
  so a stray `λ_scale` alongside `method=:lm` fails loudly. The accepted keywords are:
  - all methods — `outliers::Vector{Int}` (frame indices to skip).
  - `:gn`, `:gn_tikhonov` — `line_search_method::Symbol` (`:armijo`, the default, applies
    the sufficient-descent condition; `:backtrack` tries the full step then halves it) and
    `store_border_pts::Bool` (collect the simulated contours at every accepted iterate into
    `stats["simBorderPtsList"]`; `false` by default because each entry holds a full
    frames × 2 × border-points contour set).
  - `:gn_tikhonov` — `λ_scale`, `θ_p` and `Γ`; see [`_fit_model_GN_tikhonov`](@ref) for
    what they mean and what they default to.
  - `:lm` — `λ::Float64`, the initial Levenberg-Marquardt damping.

# Returns
- `stats::Dict` : Optimization results with keys "η", "β", "ηList", "βList", "cost_list", "iterList"

# Throws
- `ArgumentError` : If `method` is not `:gn`, `:gn_tikhonov` or `:lm`

# Example
```julia
# Using Gauss-Newton (default)
result_gn = fit_model(model, scene, conditions, obs_pts, θ_init)

# Using Levenberg-Marquardt
result_lm = fit_model(model, scene, conditions, obs_pts, θ_init; method=:lm)

# Regularized toward a prior (e.g. the previous time window's estimate)
result_tk = fit_model(model, scene, conditions, obs_pts, θ_init;
                      method=:gn_tikhonov, θ_p=θ_prev_window, λ_scale=0.1)

# Using Gauss-Newton with backtracking
result_bt = fit_model(model, scene, conditions, obs_pts, θ_init; 
                      method=:gn, line_search_method=:backtrack)
```
"""
function fit_model(model::Stokes, scene::SqueezeFlow, conditions::Conditions, obsBorderPts::Vector{AbstractArray}, θ::Vector{Float64};
                   outliers::Vector{Int}=Int[], method::Symbol=:gn, kwargs...)

    # Validate method parameter
    @assert method in [:gn, :lm, :gn_tikhonov] "Invalid optimization method: $method. Must be :gn (Gauss-Newton), :gn_tikhonov (regularized Gauss-Newton) or :lm (Levenberg-Marquardt)"

    # Validate line_search_method
    if method in (:gn, :gn_tikhonov)
        line_search_method = get(kwargs, :line_search_method, :armijo)
        @assert line_search_method in [:armijo, :backtrack] "Invalid line search method: $line_search_method. Must be :armijo or :backtrack for Gauss-Newton"
    end

    # Everything except `method` and `outliers` is forwarded untouched to the chosen
    # implementation, which owns the defaults. Two consequences worth knowing:
    #
    #   * an argument left out here is not defaulted here either, so `θ_p`/`Γ` fall back to
    #     `_fit_model_GN_tikhonov`'s own defaults — no `nothing` sentinel needed; and
    #   * a keyword the chosen method does not accept is a MethodError rather than a silent
    #     no-op, so `method=:lm, λ_scale=…` fails loudly instead of quietly ignoring λ_scale.
    if method == :lm
        return _fit_model_LM(model, scene, conditions, obsBorderPts, θ; outliers=outliers, kwargs...)
    elseif method == :gn
        return _fit_model_GN(model, scene, conditions, obsBorderPts, θ; outliers=outliers, kwargs...)
    elseif method == :gn_tikhonov
        return _fit_model_GN_tikhonov(model, scene, conditions, obsBorderPts, θ; outliers=outliers, kwargs...)
    end
end

"""
    val_check(v::Vector{Float64})

Enforce physical bounds on optimization parameters to prevent unbounded growth.

Clamps parameter values to physically reasonable ranges and corrects negative values.

# Arguments:
- `v::Vector{Float64}` : Parameter vector [η, β] where:
  - η is viscosity (Pa·s): clamped to [1e-3, 1e5]
  - β is slip penalty parameter (L/mm): clamped to [1e-3, 1e8]

# Returns:
- `v::Vector{Float64}` : Bounded parameter vector

# Bounds:
- η (viscosity): typically 1e-3 to 1e5 Pa·s
- β (penalty parameter): can be very large in no-slip cases (1e-3 to 1e8 L/mm)
"""
function val_check(v::Vector{Float64})::Vector{Float64}
    # Enforce physical bounds on parameters
    η_min::Float64, η_max::Float64 = 1e-12, 1e12
    β_min::Float64, β_max::Float64 = 1e-12, 1e12
    
    sz::Int = size(v, 1)
    for i::Int in 1:sz
        if i == 1  # η parameter
            if v[i] < 0
                v[i] = abs(v[i])
                printstyled("  [WARNING] η: Negative value, taking absolute\n", color=:yellow)
            elseif v[i] < η_min
                v[i] = η_min
                printstyled("  [WARNING] η: Below minimum ($η_min), clamped\n", color=:yellow)
            elseif v[i] > η_max
                v[i] = η_max
                printstyled("  [WARNING] η: Above maximum ($η_max), clamped\n", color=:yellow)
            end
        elseif i == 2  # β parameter (penalty for slip condition)
            if v[i] < 0
                v[i] = abs(v[i])
                printstyled("  [WARNING] β: Negative value, taking absolute\n", color=:yellow)
            elseif v[i] < β_min
                v[i] = β_min
                printstyled("  [WARNING] β: Below minimum ($β_min), clamped\n", color=:yellow)
            elseif v[i] > β_max
                v[i] = β_max
                printstyled("  [WARNING] β: Above maximum ($β_max, no-slip limit), clamped\n", color=:yellow)
            end
        end
    end
    return v
end