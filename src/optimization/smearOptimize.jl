using LinearAlgebra
using Plots
using ArgCheck
using NearestNeighbors

function match_points(p_sim::AbstractMatrix{Float64},p_obs::AbstractMatrix{Float64})
    """
    Match simulation points to observation points using KDTree spatial indexing.
    
    Uses NearestNeighbors.jl KDTree for O(n log m) complexity instead of O(n*m) nested loop.
    Speedup: 50-1000x faster than naive O(n²) search for point clouds with 100-10000+ points.
    """
    simSize = size(p_sim,2)
    obsSize = size(p_obs,2)
    pairs = zeros(Int64, simSize, 2)
    
    p_sim_x, p_sim_y = p_sim[1,:], p_sim[2,:]
    p_obs_x, p_obs_y = p_obs[1,:], p_obs[2,:]
    
    # Build KDTree from observation points
    # KDTree expects data points as columns: (2 × obsSize) matrix
    obs_points = hcat(p_obs_x, p_obs_y)'  # Transpose to get (2 × obsSize)
    tree = KDTree(obs_points; leafsize=10)
    
    # Find nearest neighbor for each simulation point using spatial search
    for simCounter in 1:simSize
        sim_point = [p_sim_x[simCounter]; p_sim_y[simCounter]]  # Column vector for one query point
        idx, dist = knn(tree, sim_point, 1)  # Find k=1 nearest neighbor
        pairs[simCounter, :] = [simCounter, idx[1]]
    end
    
    return pairs
end

function closest_point(sim_frames::AbstractArray, obs_frames::AbstractArray; outliers::AbstractArray=[])
    # Define the cost function
    cost_list = Float64[]
    norm_cost_list = Float64[]
    pairsList = []

    # frame_counter = 1
    @argcheck length(sim_frames) == length(obs_frames) "Size of the simulation and observation scenes should be the same"
    for (frame_idx, (obs_t, sim_t)) in enumerate(zip(obs_frames, sim_frames)) # iterate over the scenes
        if frame_idx in outliers
            @info "Skipping frame $frame_idx as it is marked as an outlier."
            continue
        end
        tcost = 0
        pairs = match_points(sim_t, obs_t) # match the points using the first border
        
        pSim, qSim = sim_t[1,:], sim_t[2,:]
        pObs, qObs = obs_t[1,:], obs_t[2,:]
        
        u = [(pSim[pairs[:,1]] - pObs[pairs[:,2]]); (qSim[pairs[:,1]] - qObs[pairs[:,2]])]
        tcost = u'*u

        u_ = [(pObs[pairs[:,2]]); (qObs[pairs[:,2]])]
        denom = u_'*u_
        
        mCost = tcost/(2*length(pairs))     # 1/2m(Σ(√(xi-x_obs)^2+(yi-y_obs)^2))^2 (mean error)
        norm_cost = tcost/denom
        push!(cost_list, mCost)
        push!(pairsList,pairs)
        push!(norm_cost_list,norm_cost)
    end
    return cost_list, [pairsList, norm_cost_list]
end

function closest_point(sim_frames::AbstractArray, obs_frames::AbstractArray, dudθ::AbstractArray; outliers::AbstractArray=[])
    # Define the cost function 
    cost_list = Float64[]
    dcost_list = []
    dcost2List = []
    pairsList = []
    
    @argcheck length(sim_frames) == length(obs_frames) "Size of the simulation and observation scenes should be the same"
    @argcheck length(sim_frames) == length(dudθ) "Size of the simulation and observation scenes should be the same"

    for (frame_idx, (obs_t, sim_t, du_tdθ)) in enumerate(zip(obs_frames, sim_frames, dudθ)) # iterate over the scenes
        if frame_idx in outliers
            @info "Skipping frame $frame_idx as it is marked as an outlier."
            continue
        end
        @argcheck size(sim_t,2) == size(du_tdθ,2) "Number of the border points and the gradient points should be the same"

        mat_nan_inf_check(du_tdθ[:,:,1])
        mat_nan_inf_check(du_tdθ[:,:,2])

        nθ = size(du_tdθ,3)
        tcost = 0.0
        dtcost = zeros(Float64,1,nθ)
        dt2cost = zeros(Float64,nθ,nθ)
        
        pairs = match_points(sim_t, obs_t) # match the points using the first border
        
        pSim, qSim = sim_t[1,:], sim_t[2,:]
        dpSim, dqSim = du_tdθ[1,:,:], du_tdθ[2,:,:]
        pObs, qObs = obs_t[1,:], obs_t[2,:]

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

function init_cylinder()             
    scale = 100
    ne::Int = 4 # number of elements in the mesh for the ground truth
    FunctionClass::String = "Q2"
    camera_matrix = [[8*2048/7.07, 0.0, 2048/2] [0.0, 8*1536/5.3, 1536/2] [0.0, 0.0, 1.0]]'
    camera_pose = scale*[0 -0.25 2]'   # camera position in mm

    rList = Vector{Float64}(undef,0)
    hList = Vector{Float64}(undef,0)
    cost_list = Vector{Float64}(undef,0)
    iterList = Vector{Float64}(undef,0)

    NodeList, IEN, ID, IEN_top, IEN_bottom, IEN_side, nNodes, BorderNodes = meshgrid_cube(1, 1, 1, ne, FunctionClass=FunctionClass)

    # g_truth
    r_gt::Float64 = 0.25*scale  # radius of the cylinder in mm
    h_gt::Float64 = 0.5*scale  # height of the cylinder in mm
    NodeListCyl_gt = inflate_cylinder(NodeList, -0.5, 0.5, -0.5, 0.5, r_gt, h_gt)
    side_nodes = BorderNodes[1]
    obsBorderPts, g_SurfacePts2D = extract_borders(NodeListCyl_gt, camera_matrix, camera_pose, BorderNodesList=side_nodes, nNodes)

    # optimizer
    r = 1*scale*ones(ne)
    h = 1*scale
    NodeListCyl, ∇NodeListCyl = inflate_cylinder(NodeList, -0.5, 0.5, -0.5, 0.5, r, h, GRAD=true)
    simBorderPts, ∇BorderPts2D = extract_borders(NodeListCyl, camera_matrix, camera_pose, nNodes, BorderNodesList=side_nodes, GRAD=true, dqdθ=∇NodeListCyl, SIDES=false)

    d, ∂d, ∂2d, pairs = closest_point([simBorderPts],[obsBorderPts],[∇BorderPts2D])
    totdinit::Float64 = sum(d)/length(d)

    printstyled("Ground truth: r = $(round(r_gt, sigdigits=4)), h = $(round(h_gt, sigdigits=4))\n", color=:green)
    θ = vcat(r,h)
    iter::Int = 1
    c_grad::Float64 = 1.0
    printstyled("Initial error = $(round(totdinit, sigdigits=4))\n", color=:yellow)
    while true
        printstyled("\n================ Iteration $iter ================", color=:blue)
        println()

        t∂2d = zeros(size(∂2d[1]))
        t∂d = zeros(size(∂d[1]))
        
        println("Current parameters: r = $(round(θ[1], sigdigits=4)), h = $(round(θ[2], sigdigits=4)), cost = $(round(totdinit, sigdigits=4))")

        len_d = length(d)
        szd = 1:len_d

        for i in szd
            t∂2d = t∂2d + ∂2d[i]
            t∂d = t∂d + ∂d[i]
        end

        p = t∂2d\t∂d
        println("Newton step: [$(round(p[1], sigdigits=4)), $(round(p[2], sigdigits=4))]")
        
        α = 1
        θ = θ - α*p
        val_check(θ)
        
        r = θ[1]
        h = θ[2]

        NodeListCyl, ∇NodeListCyl = inflate_cylinder(NodeList, -0.5, 0.5, -0.5, 0.5, r, h, GRAD=true)
        simBorderPts, ∇BorderPts2D = extract_borders(NodeListCyl, camera_matrix, camera_pose, nNodes, BorderNodesList=side_nodes, GRAD=true, dqdθ=∇NodeListCyl, SIDES=false)

        d, ∂d, ∂2d, pairs = closest_point([simBorderPts],[obsBorderPts],[∇BorderPts2D])

        totd = sum(d)/len_d
        c_grad = abs(totdinit - totd)
        c_grad_rel = c_grad / (abs(totdinit) + 1e-12)
        totdinit = totd

        iter = iter + 1
        
        push!(rList,θ[1]) 
        push!(hList,θ[2])
        push!(cost_list,totd)
        push!(iterList,iter)
        println("Result: r = $(round(θ[1], sigdigits=4)), h = $(round(θ[2], sigdigits=4)), cost = $(round(totd, sigdigits=4))")
        println("Delta: Δcost = $(round(c_grad, sigdigits=3)) (rel: $(round(c_grad_rel, sigdigits=3)))")

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

function armijo_line_search(model::Stokes, scene::SqueezeFlow, conditions::Conditions, obsBorderPts::Vector{AbstractArray},
                            θ_prev::Vector{Float64}, p_damped::Vector{Float64}, ∂d_prev::Vector, cost_prev::Float64;
                            c::Float64=1e-4, max_backtracks::Int=10, outliers::Vector{Int}=Int[])
    """
    Armijo line search with sufficient descent condition.
    Accepts step if: cost(θ - α·p) ≤ cost_prev - c·α·||∇cost||²
    
    Parameters:
      c: sufficient decrease parameter (default 1e-4, range [1e-4, 0.1])
      max_backtracks: maximum backtracking iterations
    
    Returns: (θ_new, d_new, ∂d_new, ∂2d_new, cost_new, accepted)
    """
    len_d = length(∂d_prev)
    α = 1.0
    grad_prod = sum(j -> dot(∂d_prev[j], ∂d_prev[j]), 1:len_d)  # ||∇cost||²
    
    # Initialize variables for scope after loop
    θ_trial = θ_prev - α * p_damped
    d_trial = similar(∂d_prev[1])
    ∂d = similar(∂d_prev)
    ∂2d = zeros(0)  # Will be overwritten
    cost_trial = cost_prev
    
    println("      Trial α    | Cost         | Δcost        | Required Decrease | Status")
    println("      " * "-"^73)
    
    for backtrack_iter in 1:max_backtracks
        θ_trial = θ_prev - α * p_damped
        val_check(θ_trial)
        
        model.η = [θ_trial[1]]
        scene.β = [θ_trial[2]]
        μ_list, gradList, simBorderPts, splinex, spliney, pos2D = simulate(model, scene, conditions)
        d_trial, ∂d, ∂2d, pairs = closest_point(simBorderPts, obsBorderPts, gradList, outliers=outliers)
        cost_trial = sum(d_trial) / len_d
        
        # Armijo condition: cost_trial ≤ cost_prev - c·α·||∇||²
        Δcost = cost_prev - cost_trial
        sufficient_decrease = c * α * grad_prod
        
        status = Δcost ≥ sufficient_decrease ? "[ACCEPT]" : "[REJECT]"
        println("      $(round(α, sigdigits=5)) | $(round(cost_trial, sigdigits=3)) | $(round(Δcost, sigdigits=3)) | $(round(sufficient_decrease, sigdigits=3)) | $status")
        
        if cost_trial ≤ cost_prev - sufficient_decrease
            println("      " * "-"^73)
            printstyled("      [ACCEPT] Step accepted: α = $(round(α, sigdigits=5)), cost = $(round(cost_trial, sigdigits=4))\n", color=:green)
            return θ_trial, d_trial, ∂d, ∂2d, cost_trial, true
        end
        
        α *= 0.5  # Halve step size and retry
    end
    
    println("      " * "-"^73)
    printstyled("      [WARNING] Max backtracks reached. Accepting last trial (α=$(round(α, sigdigits=5)), cost=$(round(cost_trial, sigdigits=4)))\n", color=:yellow)
    return θ_trial, d_trial, ∂d, ∂2d, cost_trial, false
end

function backtrack_line_search(model::Stokes, scene::SqueezeFlow, conditions::Conditions, obsBorderPts::Vector{AbstractArray},
                               θ_prev::Vector{Float64}, p_damped::Vector{Float64}, cost_prev::Float64; outliers::Vector{Int}=Int[])
    """
    Simple backtracking line search (original method - kept for reference/comparison).
    Takes full step, then half step if needed. No principled descent guarantee.
    """
    len_d_calc = 0
    
    println("      Strategy: Try α=1.0 (full), then α=0.5 (half)")
    
    # Try full step
    println("      Trying α = 1.0 (full step)...")
    θ_trial = θ_prev - p_damped
    val_check(θ_trial)
    model.η = [θ_trial[1]]
    scene.β = [θ_trial[2]]
    μ_list, gradList, simBorderPts, splinex, spliney, pos2D = simulate(model, scene, conditions)
    d_trial, ∂d, ∂2d, pairs = closest_point(simBorderPts, obsBorderPts, gradList, outliers=outliers)
    len_d_calc = length(d_trial)
    cost_trial = sum(d_trial) / len_d_calc
    Δcost_full = cost_prev - cost_trial
    println("        cost = $(round(cost_trial, sigdigits=4)), Δcost = $(round(Δcost_full, sigdigits=4))")
    
    if cost_trial < cost_prev
        printstyled("        [ACCEPT] Cost decreased, accepting\n", color=:green)
        return θ_trial, d_trial, ∂d, ∂2d, cost_trial, true
    end
    
    # Try half step
    println("      Trying α = 0.5 (half step)...")
    θ_trial = θ_prev - 0.5 * p_damped
    val_check(θ_trial)
    model.η = [θ_trial[1]]
    scene.β = [θ_trial[2]]
    μ_list, gradList, simBorderPts, splinex, spliney, pos2D = simulate(model, scene, conditions)
    d_trial, ∂d, ∂2d, pairs = closest_point(simBorderPts, obsBorderPts, gradList, outliers=outliers)
    cost_trial_half = sum(d_trial) / len_d_calc
    Δcost_half = cost_prev - cost_trial_half
    println("        cost = $(round(cost_trial_half, sigdigits=4)), Δcost = $(round(Δcost_half, sigdigits=4))")
    
    if cost_trial_half < cost_prev
        printstyled("        [ACCEPT] Cost decreased, accepting\n", color=:green)
        return θ_trial, d_trial, ∂d, ∂2d, cost_trial_half, true
    else
        printstyled("        [WARNING] Neither step decreased cost. Accepting half step anyway.\n", color=:yellow)
        return θ_trial, d_trial, ∂d, ∂2d, cost_trial_half, false
    end
end

function fit_model(model::Stokes, scene::SqueezeFlow, conditions::Conditions, obsBorderPts::Vector{AbstractArray}, θ::Vector{Float64};
                   outliers::Vector{Int}=Int[], line_search_method::String="armijo")
    reset_model!(model)
    model.η = [θ[1]]
    scene.β = [θ[2]]

    ηpList = Vector{Float64}(undef,0)
    βpList = Vector{Float64}(undef,0)
    cost_list = Vector{Float64}(undef,0)
    iterList = Vector{Float64}(undef,0)
    
    printstyled("Initializing simulation with η: $(round(θ[1], sigdigits=4)), β: $(round(θ[2], sigdigits=4))\n", color=:cyan)
    μ_list, gradList, simBorderPts, splinex, spliney, pos2D = simulate(model, scene, conditions)
    d, ∂d, ∂2d, pairs = closest_point(simBorderPts, obsBorderPts, gradList, outliers=outliers)
    totdinit::Float64 = sum(d)/length(d)

    push!(ηpList,θ[1])
    push!(βpList,θ[2])
    push!(cost_list,totdinit)
    push!(iterList,1)

    iter::Int = 1
    c_grad::Float64 = 1.0
    printstyled("Initial error = $totdinit\n", color=:yellow)
    while true
        printstyled("\n==================== Iteration $iter ===================\n", color=:blue)

        reset_model!(model)
        t∂2d = zeros(size(∂2d[1]))
        t∂d = zeros(size(∂d[1]))

        # Current state
        if iter == 1
            printstyled("Starting optimization\n", color=:yellow)
        end
        println("Current parameters: η = $(round(θ[1], sigdigits=4)), β = $(round(θ[2], sigdigits=4)), cost = $(round(totdinit, sigdigits=4))")

        len_d = length(d)
        szd = 1:len_d

        for i in szd
            t∂2d = t∂2d + ∂2d[i]
            t∂d = t∂d + ∂d[i]
        end

        # Compute Newton step
        p = t∂2d\t∂d
        println("Newton step (undamped): [$(round(p[1], sigdigits=4)), $(round(p[2], sigdigits=4))]")
        
        # Compute Hessian diagnostics and damping
        # Parameter-specific damping: avoid global condition number damping which over-constrains
        # well-conditioned directions. Instead, only damp directions with weak Hessian curvature.
        # 
        # Each Newton step is: p_i = (∂J/∂θ_i) / (∂²J/∂θ_i²)
        # If ∂²J/∂θ_i² is small → p_i is large → needs damping
        # If ∂²J/∂θ_i² is reasonable → p_i is reasonable → no damping needed
        # 
        # This avoids the issue with global κ(H) damping: if one parameter has huge curvature,
        # it shouldn't penalize other parameters with normal curvature.
        
        H_η = abs(t∂2d[1,1])       # ∂²cost/∂η²  
        H_β = abs(t∂2d[2,2])       # ∂²cost/∂β²
        H_geom_mean = sqrt(H_η * H_β + 1e-20)  # Avoid underflow
        H_ratio = H_η / (H_β + 1e-12)
        
        # Apply clamping-based damping
        # Damping factor: α = (target_relative_step × parameter) / |step|
        # This ensures: |α × step| / parameter ≈ target_relative_step
        target_rel_step = 0.15  # Target 15% relative change
        scale_η = abs(θ[1]) + 1e-12
        scale_β = abs(θ[2]) + 1e-12

        # Relative Newton step (as % of parameter)
        rel_step_η = abs(p[1]) / scale_η
        rel_step_β = abs(p[2]) / scale_β

        # Damping: if relative step is huge, clamp it, otherwise trust Hessian
        α_η = clamp(target_rel_step / (rel_step_η + 1e-12), 0.01, 1.0)
        α_β = clamp(target_rel_step / (rel_step_β + 1e-12), 0.01, 1.0)

        println("Hessian diagnostics:")
        println("  H_η = $(round(H_η, sigdigits=3)), H_β = $(round(H_β, sigdigits=3)), H_ratio = $(round(H_ratio, sigdigits=3))")
        
        # Log damping if applied
        if α_η < 0.99 || α_η > 1.01 || α_β < 0.99 || α_β > 1.01
            printstyled("  [WARNING] Damping factors: α_η = $(round(α_η, sigdigits=3)), α_β = $(round(α_β, sigdigits=3))\n", color=:yellow)
        else
            println("  Damping factors: α_η = $(round(α_η, sigdigits=3)), α_β = $(round(α_β, sigdigits=3)) (no damping needed)")
        end
        
        # Apply damping to Newton step
        p_damped = [α_η * p[1]; α_β * p[2]]
        println("Damped step: [$(round(p_damped[1], sigdigits=4)), $(round(p_damped[2], sigdigits=4))]")
        
        # Store state before step for line search
        θ_prev = copy(θ)
        cost_prev = totdinit
        
        # ============ LINE SEARCH ============
        println("Line search ($line_search_method):")
        if lowercase(line_search_method) == "armijo"
            θ, d, ∂d, ∂2d, totd, ls_accepted = armijo_line_search(model, scene, conditions, obsBorderPts, θ_prev, p_damped, ∂d, cost_prev, outliers=outliers)
        else
            θ, d, ∂d, ∂2d, totd, ls_accepted = backtrack_line_search(model, scene, conditions, obsBorderPts, θ_prev, p_damped, cost_prev, outliers=outliers)
        end
        # ===============================
        
        c_grad = abs(cost_prev - totd)
        c_grad_rel = c_grad / (abs(cost_prev) + 1e-12)  # Relative cost change
        totdinit = totd
        
        # Compute relative parameter changes
        Δη_rel = abs(θ[1] - θ_prev[1]) / (abs(θ_prev[1]) + 1e-12)
        Δβ_rel = abs(θ[2] - θ_prev[2]) / (abs(θ_prev[2]) + 1e-12)
        Δθ_max = max(Δη_rel, Δβ_rel)

        iter = iter + 1
        
        push!(ηpList,θ[1]) 
        push!(βpList,θ[2])
        push!(cost_list,totd)
        push!(iterList,iter)
        
        # Results summary
        println("Result: η = $(round(θ[1], sigdigits=4)), β = $(round(θ[2], sigdigits=4)), cost = $(round(totd, sigdigits=4))")
        println("Deltas: Δη/η = $(round(Δη_rel, sigdigits=3)), Δβ/β = $(round(Δβ_rel, sigdigits=3)), Δcost = $(round(c_grad, sigdigits=3)) (rel: $(round(c_grad_rel, sigdigits=3)))")

        # Check convergence
        cost_converged = c_grad_rel < 1e-4  # Relative cost change < 0.01%
        
        if cost_converged
            printstyled("[CONVERGED] Relative cost change = $(round(c_grad_rel, sigdigits=3))\n", color=:green)
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
        "iterList" => iterList
    )
    return stats
end

function val_check(v::Vector{Float64})
    # Enforce physical bounds on parameters to prevent unbounded growth
    # η (viscosity): typically 1e-3 to 1e5 Pa·s
    # β (penalty parameter): can be very large in no-slip cases (1e-3 to 1e8 L/mm)
    η_min, η_max = 1e-3, 1e5
    β_min, β_max = 1e-3, 1e8  # Allow very large β for strong no-slip conditions
    
    sz = size(v,1)
    for i in 1:sz
        if i == 1  # η parameter
            if v[i] < 0 
                v[i] = abs(v[i])
                printstyled("  [WARNING] η: Negative value corrected to positive\n", color=:yellow)
            end
            if v[i] < η_min
                v[i] = η_min
                printstyled("  [WARNING] η: Below minimum ($η_min), clamped\n", color=:yellow)
            elseif v[i] > η_max
                v[i] = η_max
                printstyled("  [WARNING] η: Above maximum ($η_max), clamped\n", color=:yellow)
            end
        elseif i == 2  # β parameter (penalty for slip condition)
            if v[i] < 0 
                v[i] = abs(v[i])
                printstyled("  [WARNING] β: Negative value corrected to positive\n", color=:yellow)
            end
            if v[i] < β_min
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


function height_sample(sim_frames::Vector{AbstractArray}, obs_frames::Vector{AbstractArray})

    xSimintlst = AbstractArray[]
    ySimintlst = AbstractArray[]
    xObsintlst = AbstractArray[]
    cost_list = Float64[]
    for (obs_data, sim_data) in zip(obs_frames, sim_frames) # iterate over the scenes

        xSim, ySim = filter_points(sim_data, 2048/2)
        xObs, yObs = filter_points(obs_data, 2048/2)

        xObsint = fit_curve(borderx=xObs, bordery=yObs, samples=ySim)

        cost = 0
        iter = 1:length(ySim)
        for i in iter
            error = (xSim[i] - xObsint[i])^2
            cost += error
        end

        totPts = length(ySim)
        mCost = cost/totPts

        push!(cost_list, mCost)
        push!(xSimintlst, xSim)
        push!(ySimintlst, ySim)
        push!(xObsintlst, xObsint)
    end
    return cost_list, xObsintlst, xSimintlst, ySimintlst
end
# animate_fields(filepath = string(conditions.filepath,"/Results/images"),fields2D=borderPts2DList, pObs=splinexObs, qObs=splineyObs)

# init_cylinder()