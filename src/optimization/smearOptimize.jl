using LinearAlgebra
using Plots
using ArgCheck
using NearestNeighbors

"""
    match_points(p_sim::AbstractMatrix{Float64}, p_obs::AbstractMatrix{Float64})

Match simulation points to observation points using KDTree spatial indexing.

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
    closest_point(sim_frames::AbstractArray, obs_frames::AbstractArray; outliers::AbstractArray=[])

Compute point-wise distance between simulation and observation frames without gradient information.

For each frame pair, matches simulation points to observation points and computes the mean squared error.

# Arguments:
- `sim_frames::AbstractArray` : Vector of simulation point coordinates [2 × n_points for each frame]
- `obs_frames::AbstractArray` : Vector of observation point coordinates [2 × n_points for each frame]
- `outliers::AbstractArray` : Frame indices to skip

# Returns:
- `cost_list::Vector{Float64}` : Mean squared error for each frame
- `[pairsList, norm_cost_list]` : Tuple containing point pairs and normalized costs for each frame
"""
function _as_2xN(x)
    raw = x isa AbstractVector ? x[1] : x   # unwrap Vector{Matrix} wrapper if present
    m = Matrix{Float64}(collect(raw))
    return size(m, 1) == 2 ? m : Matrix{Float64}(m')
end

_as_dudθ(x) = Array{Float64}(x isa AbstractVector ? x[1] : x)

function closest_point(sim_frames::AbstractArray, obs_frames::AbstractArray; outliers::AbstractArray=[])
    # Define the cost function
    cost_list = Float64[]
    norm_cost_list = Float64[]
    pairsList = []

    # frame_counter = 1
    @argcheck length(sim_frames) == length(obs_frames) "Size of the simulation and observation scenes should be the same"
    for (frame_idx, (obs_t, _sim_t)) in enumerate(zip(obs_frames, sim_frames)) # iterate over the scenes

        if frame_idx in outliers
            @info "Skipping frame $frame_idx as it is marked as an outlier."
            continue
        end

        sim_t = _as_2xN(_sim_t)
        obs_t = _as_2xN(obs_t)
        
        pairs = match_points(sim_t, obs_t) # match the points using the first border
        
        pSim, qSim = sim_t[1, :], sim_t[2, :]
        pObs, qObs = obs_t[1, :], obs_t[2, :]
        
        u = [(pSim[pairs[:, 1]] - pObs[pairs[:, 2]]); (qSim[pairs[:, 1]] - qObs[pairs[:, 2]])]
        cost = u' * u

        u_ = [(pObs[pairs[:, 2]]); (qObs[pairs[:, 2]])]
        denom = u_' * u_
        
        mean_cost = cost / (2 * length(pairs))
        norm_cost = cost / denom
        push!(cost_list, mean_cost)
        push!(pairsList, pairs)
        push!(norm_cost_list, norm_cost)
    end
    return cost_list, [pairsList, norm_cost_list]
end

"""
    closest_point(sim_frames::AbstractArray, obs_frames::AbstractArray, dudθ::AbstractArray; outliers::AbstractArray=[])

Compute point-wise distance with first and second derivatives for gradient-based optimization.

For each frame pair, matches simulation points to observation points and computes the mean squared error
along with its first and second derivatives with respect to optimization parameters.

# Arguments:
- `sim_frames::AbstractArray` : Vector of simulation point coordinates
- `obs_frames::AbstractArray` : Vector of observation point coordinates
- `dudθ::AbstractArray` : Gradient of simulation points w.r.t. parameters [2 × n_points × n_params per frame]
- `outliers::AbstractArray` : Frame indices to skip

# Returns:
- `cost_list::Vector{Float64}` : Mean squared error for each frame
- `dcost_list::Vector{Matrix}` : First derivative (gradient) w.r.t. parameters for each frame
- `dcost2List::Vector{Matrix}` : Second derivative (Hessian) w.r.t. parameters for each frame
- `pairsList` : Point correspondence pairs for each frame
"""
function closest_point(sim_frames::AbstractArray, obs_frames::AbstractArray, dudθ::AbstractArray; outliers::AbstractArray=[])::Tuple{Vector{Float64}, Vector, Vector, Vector}
    cost_list::Vector{Float64} = Float64[]
    dcost_list::Vector = []
    dcost2List::Vector = []
    pairsList::Vector = []
    
    @argcheck length(sim_frames) == length(obs_frames) "Size of the simulation and observation scenes should be the same"
    @argcheck length(sim_frames) == length(dudθ) "Size of the simulation and observation scenes should be the same"

    for (frame_idx, (obs_t, _sim_t, _du_tdθ)) in enumerate(zip(obs_frames, sim_frames, dudθ)) # iterate over the scenes

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

        nθ::Int = size(du_tdθ, 3)
        tcost::Float64 = 0.0
        dtcost::Vector{Float64} = zeros(Float64, nθ)
        dt2cost::Matrix{Float64} = zeros(Float64, nθ, nθ)
        
        pairs::Matrix{Int64} = match_points(sim_t, obs_t)
        
        pSim::Vector{Float64} = sim_t[1, :]
        qSim::Vector{Float64} = sim_t[2, :]
        dpSim::Matrix{Float64} = du_tdθ[1, :, :]
        dqSim::Matrix{Float64} = du_tdθ[2, :, :]
        pObs::Vector{Float64} = obs_t[1, :]
        qObs::Vector{Float64} = obs_t[2, :]

        u = [(pSim[pairs[:,1]] - pObs[pairs[:,2]]); (qSim[pairs[:,1]] - qObs[pairs[:,2]])]
        Jmat = [dpSim[pairs[:,1],:]; dqSim[pairs[:,1],:]]

        tcost = u'*u
        dtcost = Jmat'*u
        dt2cost = Jmat'*Jmat

        mCost = tcost/(2*length(pairs)) # 1/2m(Σ(√(xi-x_obs)^2+(yi-y_obs)^2))^2 (mean error)
        dmcost = dtcost/length(pairs)   # 1/m(Σ(xi-x_obs)∂x/∂θ_i +(yi-y_obs)∂y/∂θ_i))
        dm2cost = dt2cost/length(pairs) # 1/m(Σ(∂x2/∂2θ_i + ∂y2/∂2θ_i))

        push!(cost_list, mCost)
        push!(dcost_list,dmcost)
        push!(dcost2List,dm2cost)
        push!(pairsList,pairs)
    end
    return cost_list, dcost_list, dcost2List, pairsList
end

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
    camera_matrix::Matrix{Float64} = [[8 * 2048 / 7.07, 0.0, 2048 / 2] [0.0, 8 * 1536 / 5.3, 1536 / 2] [0.0, 0.0, 1.0]]'
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

    d, ∂d, ∂2d, _ = closest_point([simBorderPts], [obsBorderPts], [∇BorderPts2D])
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

        d, ∂d, ∂2d, _ = closest_point([simBorderPts],[obsBorderPts],[∇BorderPts2D])

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
    armijo_line_search(model, scene, conditions, obsBorderPts, θ_prev, p_damped, ∂d_prev, cost_prev; c=1e-4, max_backtracks=10, outliers=Int[])

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

# Returns:
- `θ_trial::Vector{Float64}` : New parameter values
- `d_trial::Vector{Float64}` : Cost values at new point
- `∂d::Vector` : Gradient of cost at new point
- `∂2d::Vector` : Hessian of cost at new point
- `cost_trial::Float64` : Cost value at trial point
- `accepted::Bool` : Whether step was accepted by Armijo condition
"""
function armijo_line_search(model::Stokes, scene::SqueezeFlow, conditions::Conditions, obsBorderPts::Vector{AbstractArray},
                            θ_prev::Vector{Float64}, p_damped::Vector{Float64}, ∂d_prev::Vector, cost_prev::Float64;
                            c::Float64=1e-4, max_backtracks::Int=10, outliers::Vector{Int}=Int[])
    # Armijo line search with sufficient descent condition.
    # Accepts step if: cost(θ - α·p) ≤ cost_prev - c·α·∇f·p
    len_d = length(∂d_prev)
    α = 1.0
    t_grad = sum(∂d_prev)                 # total gradient ∇f = Σ ∂d[j]
    grad_prod = dot(t_grad, p_damped)     # directional derivative ∇fᵀ·p
    
    # Initialize variables for scope after loop
    θ_trial = θ_prev - α * p_damped
    d_trial = similar(∂d_prev[1])
    ∂d = similar(∂d_prev)
    ∂2d = zeros(0)  # Will be overwritten
    cost_trial = cost_prev
    simBorderPts_accepted = nothing
    gradList_accepted = nothing
    
    for backtrack_iter in 1:max_backtracks
        θ_trial = θ_prev - α * p_damped
        val_check(θ_trial)
        
        model.η = [θ_trial[1]]
        scene.β = [θ_trial[2]]
        μ_list, gradList, simBorderPts, _, _, _, _, _, _, _ = simulate(model, scene, conditions)
        
        # Use cost-only version during line search (fast)
        d_trial_costs, _ = closest_point(simBorderPts, obsBorderPts, outliers=outliers)
        cost_trial = sum(d_trial_costs) / len_d
        
        # Armijo condition: cost_trial ≤ cost_prev - c·α·∇f·p
        Δcost = cost_prev - cost_trial
        sufficient_decrease = c * α * grad_prod
        
        status = Δcost ≥ sufficient_decrease ? "[ACCEPT]" : "[REJECT]"
        @debug "      Trial α    | Cost         | Δcost        | Required Decrease | Status"
        @debug "      $("-"^73)"
        @debug "      $(round(α, sigdigits=5)) | $(round(cost_trial, sigdigits=3)) | $(round(Δcost, sigdigits=3)) | $(round(sufficient_decrease, sigdigits=3)) | $status"
        
        if cost_trial ≤ cost_prev - sufficient_decrease
            # Compute full gradients/Hessians only for accepted step
            d_trial, ∂d, ∂2d, pairs = closest_point(simBorderPts, obsBorderPts, gradList, outliers=outliers)
            @debug "      $("-"^73)"
            printstyled("      [ACCEPT] Step accepted: α = $(round(α, sigdigits=5)), cost = $(round(cost_trial, sigdigits=4))\n", color=:green)
            return θ_trial, d_trial, ∂d, ∂2d, cost_trial, true
        end
        
        α *= 0.5  # Halve step size and retry
    end
    
    # Fallback: compute gradients/Hessians for final trial (even if rejected)
    model.η = [θ_trial[1]]
    scene.β = [θ_trial[2]]
    μ_list, gradList, simBorderPts, _, _, _, _, _, _, _ = simulate(model, scene, conditions)
    d_trial, ∂d, ∂2d, pairs = closest_point(simBorderPts, obsBorderPts, gradList, outliers=outliers)
    
    @debug "      $("-"^73)"
    printstyled("      [WARNING] Max backtracks reached. Accepting last trial (α=$(round(α, sigdigits=5)), cost=$(round(cost_trial, sigdigits=4)))\n", color=:yellow)
    return θ_trial, d_trial, ∂d, ∂2d, cost_trial, false
end

"""
    backtrack_line_search(model, scene, conditions, obsBorderPts, θ_prev, p_damped, cost_prev; outliers=Int[])

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

# Returns:
- `θ_trial::Vector{Float64}` : New parameter values
- `d_trial::Vector{Float64}` : Cost values at new point
- `∂d::Vector` : Gradient of cost at new point
- `∂2d::Vector` : Hessian of cost at new point
- `cost_trial::Float64` : Cost value at trial point
- `accepted::Bool` : Whether step was accepted
"""
function backtrack_line_search(model::Stokes, scene::SqueezeFlow, conditions::Conditions, obsBorderPts::Vector{AbstractArray},
                               θ_prev::Vector{Float64}, p_damped::Vector{Float64}, cost_prev::Float64; outliers::Vector{Int}=Int[])
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
    d_trial_costs, _ = closest_point(simBorderPts, obsBorderPts, outliers=outliers)
    len_d_calc = length(d_trial_costs)
    cost_trial = sum(d_trial_costs) / len_d_calc
    Δcost_full = cost_prev - cost_trial
    @debug "        cost = $(round(cost_trial, sigdigits=4)), Δcost = $(round(Δcost_full, sigdigits=4))"
    
    if cost_trial < cost_prev
        # Compute full gradients/Hessians only for accepted step
        d_trial, ∂d, ∂2d, pairs = closest_point(simBorderPts, obsBorderPts, gradList, outliers=outliers)
        printstyled("        [ACCEPT] Cost decreased, accepting\n", color=:green)
        return θ_trial, d_trial, ∂d, ∂2d, cost_trial, true
    end
    
    # Try half step
    @debug "      Trying α = 0.5 (half step)..."
    θ_trial = θ_prev - 0.5 * p_damped
    val_check(θ_trial)
    model.η = [θ_trial[1]]
    scene.β = [θ_trial[2]]
    μ_list, gradList, simBorderPts, _, _, _, _, _, _, _ = simulate(model, scene, conditions)
    
    # Use cost-only version first (fast)
    d_trial_costs, _ = closest_point(simBorderPts, obsBorderPts, outliers=outliers)
    cost_trial_half = sum(d_trial_costs) / len_d_calc
    Δcost_half = cost_prev - cost_trial_half
    @debug "        cost = $(round(cost_trial_half, sigdigits=4)), Δcost = $(round(Δcost_half, sigdigits=4))"
    
    # Compute full gradients/Hessians for final result (accepted or rejected)
    d_trial, ∂d, ∂2d, pairs = closest_point(simBorderPts, obsBorderPts, gradList, outliers=outliers)
    
    if cost_trial_half < cost_prev
        printstyled("        [ACCEPT] Cost decreased, accepting\n", color=:green)
        return θ_trial, d_trial, ∂d, ∂2d, cost_trial_half, true
    else
        printstyled("        [WARNING] Neither step decreased cost. Accepting half step anyway.\n", color=:yellow)
        return θ_trial, d_trial, ∂d, ∂2d, cost_trial_half, false
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
                   outliers::Vector{Int}=Int[], line_search_method::Symbol=:armijo)
    reset_model!(model)
    model.η = [θ[1]]
    scene.β = [θ[2]]

    ηpList::Vector{Float64} = Float64[]
    βpList::Vector{Float64} = Float64[]
    cost_list::Vector{Float64} = Float64[]
    iterList::Vector{Float64} = Float64[]
    
    printstyled("Initializing simulation with η: $(round(θ[1], sigdigits=4)), β: $(round(θ[2], sigdigits=4))\n", color=:cyan)
    μ_list, gradList, simBorderPts, _, _, _, _, _, _, _ = simulate(model, scene, conditions)

    @debug begin
        path = joinpath(get_scratch_dir(), "optimization_animation")
        animate_fields(filepath=path, sim_border_nodes_2d=simBorderPts, obs_border_nodes_2d=obsBorderPts)
        "Saved optimization animation to $path"
    end

    d, ∂d, ∂2d, _ = closest_point(simBorderPts, obsBorderPts, gradList, outliers=outliers)
    totdinit::Float64 = sum(d)/length(d)

    push!(ηpList,θ[1])
    push!(βpList,θ[2])
    push!(cost_list,totdinit)
    push!(iterList,1)

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
            θ, d, ∂d, ∂2d, totd, _ = armijo_line_search(model, scene, conditions, obsBorderPts, θ_prev, p, ∂d, cost_prev, outliers=outliers)
        else
            θ, d, ∂d, ∂2d, totd, _ = backtrack_line_search(model, scene, conditions, obsBorderPts, θ_prev, p, cost_prev, outliers=outliers)
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
        @debug "Result: η = $(round(θ[1], sigdigits=4)), β = $(round(θ[2], sigdigits=4)), cost = $(round(totd, sigdigits=4))"
        @debug "Deltas: Δη/η = $(round(Δη_rel, sigdigits=3)), Δβ/β = $(round(Δβ_rel, sigdigits=3)), Δcost = $(round(c_grad, sigdigits=3)) (rel: $(round(c_grad_rel, sigdigits=3)))"

        if c_grad_rel < 1e-3 && c_grad < 1e-3 && iter > 30
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
                      outliers::Vector{Int}=Int[], λ::Float64=1e-3)
                      
    reset_model!(model)
    model.η = [θ[1]]
    scene.β = [θ[2]]

    ηpList::Vector{Float64} = Float64[]
    βpList::Vector{Float64} = Float64[]
    cost_list::Vector{Float64} = Float64[]
    iterList::Vector{Float64} = Float64[]
    
    printstyled("Initializing simulation with η: $(round(θ[1], sigdigits=4)), β: $(round(θ[2], sigdigits=4))\n", color=:cyan)
    μ_list, gradList, simBorderPts, _, _, _, _, _, _, _ = simulate(model, scene, conditions)
    d, ∂d, ∂2d, _ = closest_point(simBorderPts, obsBorderPts, gradList, outliers=outliers)
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
        d, ∂d, ∂2d, _ = closest_point(simBorderPts, obsBorderPts, gradList, outliers=outliers)
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
  - `:lm` - Levenberg-Marquardt with adaptive damping
- `line_search_method::Symbol` : Line search strategy for `:gn` method (default: `:armijo`)
  - `:armijo` - Armijo backtracking with sufficient descent condition
  - `:backtrack` - Simple backtracking (full step, then half step)

# Returns
- `stats::Dict` : Optimization results with keys "η", "β", "ηList", "βList", "cost_list", "iterList"

# Throws
- `ArgumentError` : If `method` is not `:gn` or `:lm`

# Example
```julia
# Using Gauss-Newton (default)
result_gn = fit_model(model, scene, conditions, obs_pts, θ_init)

# Using Levenberg-Marquardt
result_lm = fit_model(model, scene, conditions, obs_pts, θ_init; method=:lm)

# Using Gauss-Newton with backtracking
result_bt = fit_model(model, scene, conditions, obs_pts, θ_init; 
                      method=:gn, line_search_method=:backtrack)
```
"""
function fit_model(model::Stokes, scene::SqueezeFlow, conditions::Conditions, obsBorderPts::Vector{AbstractArray}, θ::Vector{Float64};
                   outliers::Vector{Int}=Int[], method::Symbol=:gn, line_search_method::Symbol=:armijo)
    
    # Validate method parameter
    @assert method in [:gn, :lm] "Invalid optimization method: $method. Must be :gn (Gauss-Newton) or :lm (Levenberg-Marquardt)"
    
    # Validate line_search_method
    if method == :gn
        @assert line_search_method in [:armijo, :backtrack] "Invalid line search method: $line_search_method. Must be :armijo or :backtrack for Gauss-Newton"
    end
    
    if method == :lm
        return _fit_model_LM(model, scene, conditions, obsBorderPts, θ; outliers=outliers)
    elseif method == :gn
        return _fit_model_GN(model, scene, conditions, obsBorderPts, θ; outliers=outliers, line_search_method=line_search_method)
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