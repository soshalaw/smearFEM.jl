using LinearAlgebra
using Plots
using ArgCheck

function match_points(pSim::AbstractMatrix{Float64},pObs::AbstractMatrix{Float64})

    simSize = size(pSim,2)
    obsSize = size(pObs,2)
    pairs = zeros(Int64,size(pSim,2),2)
    pSim, qSim = pSim[1,:], pSim[2,:]
    pObs, qObs = pObs[1,:], pObs[2,:]

    for simCounter in 1:simSize
        cost = 1e18
        closestPointIdx = 0
        for obsCounter in 1:obsSize
            error = (pSim[simCounter] - pObs[obsCounter])^2 + (qSim[simCounter] - qObs[obsCounter])^2
            if error < cost
                cost = error
                closestPointIdx = obsCounter
            end
        end
        pairs[simCounter,:] = [simCounter, closestPointIdx]
    end
    return pairs
end

function closest_point(simScene::AbstractArray, obsScene::AbstractArray)
    # Define the cost function
    costList = Float64[]
    pairsList = []
    @argcheck length(simScene) == length(obsScene) "Size of the simulation and observation scenes should be the same"
    for (obs_t, sim_t) in zip(obsScene, simScene) # iterate over the scenes
        tcost = 0
        pairs = match_points(sim_t, obs_t) # match the points using the first border
        
        pSim, qSim = sim_t[1,:], sim_t[2,:]
        pObs, qObs = obs_t[1,:], obs_t[2,:]
        
        u = [(pSim[pairs[:,1]] - pObs[pairs[:,2]]); (qSim[pairs[:,1]] - qObs[pairs[:,2]])]
        tcost = u'*u

        mCost = tcost/(2*length(pairs))     # 1/2m(Σ(√(xi-x_obs)^2+(yi-y_obs)^2))^2 (mean error)
        push!(costList, mCost)
        push!(pairsList,pairs)
    end
    return costList, pairsList
end

function closest_point(simScene::AbstractArray, obsScene::AbstractArray, dudθ::AbstractArray)

    # simScene = simScene[1:end-1] # remove the last element, which is the last time step TODO solve the issue with the last time step
    # dudθ = dudθ[1:end-1] # remove the last element, which is the last time step TODO solve the issue with the last time step
    # Define the cost function 
    costList = Float64[]
    dcostList = []
    dcost2List = []
    pairsList = []
    @argcheck length(simScene) == length(obsScene) "Size of the simulation and observation scenes should be the same"
    @argcheck length(simScene) == length(dudθ) "Size of the simulation and observation scenes should be the same"

    for (obs_t, sim_t, du_tdθ) in zip(obsScene, simScene, dudθ) # iterate over the scenes
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

        push!(costList, mCost)
        push!(dcostList,dmcost)
        push!(dcost2List,dm2cost)
        push!(pairsList,pairs)
    end
    return costList, dcostList, dcost2List, pairsList
end

function height_sample(simScene::Vector{AbstractArray}, obsScene::Vector{AbstractArray})

    xSimintlst = AbstractArray[]
    ySimintlst = AbstractArray[]
    xObsintlst = AbstractArray[]
    costList = Float64[]
    for (obsData, simData) in zip(obsScene, simScene) # iterate over the scenes

        xSim, ySim = filter_points(simData, 2048/2)
        xObs, yObs = filter_points(obsData, 2048/2)

        xObsint = fit_curve(borderx=xObs, bordery=yObs, samples=ySim)

        cost = 0
        iter = 1:length(ySim)
        for i in iter
            error = (xSim[i] - xObsint[i])^2
            cost += error
        end

        totPts = length(ySim)
        mCost = cost/totPts

        push!(costList, mCost)
        push!(xSimintlst, xSim)
        push!(ySimintlst, ySim)
        push!(xObsintlst, xObsint)
    end
    return costList, xObsintlst, xSimintlst, ySimintlst
end

function fit_model(model::Stokes, scene::SqueezeFlow, conditions::Conditions, obsBorderPts::Vector{AbstractArray}, θ::Vector{Float64})
    reset_model!(model)
    model.η = [θ[1]]
    scene.β = [θ[2]]

    ηpList = Vector{Float64}(undef,0)
    βpList = Vector{Float64}(undef,0)
    costList = Vector{Float64}(undef,0)
    iterList = Vector{Float64}(undef,0)
    
    μ_list, gradList, simBorderPts, splinex, spliney, pos2D = simulate(model, scene, conditions)
    d, ∂d, ∂2d, pairs = closest_point(simBorderPts, obsBorderPts, gradList)
    totdinit::Float64 = sum(d)
    
    push!(ηpList,θ[1])
    push!(βpList,θ[2])
    push!(costList,totdinit)
    push!(iterList,1)

    iter::Int = 1
    c_grad::Float64 = 1.0
    while true
        println("totd: ", totdinit)
        reset_model!(model)
        t∂2d = zeros(size(∂2d[1]))
        t∂d = zeros(size(∂d[1]))
        
        println("η: ", θ[1], " β: ", θ[2])

        szd = 1:length(d)

        for i in szd
            t∂2d = t∂2d + ∂2d[i]
            t∂d = t∂d + ∂d[i]
        end

        p = t∂2d\t∂d
        α = 1
        
        println("step: ", p)
        θ = θ - α*p

        val_check(θ)
        # update the model parameters
        model.η = [θ[1]]
        scene.β = [θ[2]]
        # simulate the model
        μ_list, gradList, simBorderPts, splinex, spliney, pos2D = simulate(model, scene, conditions)

        # get the cost, gradient and hessian values
        d, ∂d, ∂2d, pairs = closest_point(simBorderPts, obsBorderPts, gradList)

        totd = sum(d)
        c_grad = abs(totdinit - totd)
        totdinit = totd

        iter = iter + 1
        
        push!(ηpList,θ[1]) 
        push!(βpList,θ[2])
        push!(costList,totd)
        push!(iterList,iter)
        println("iteration: ", iter, " ratio: ", c_grad)

        if c_grad < 0.005
            break
        elseif iter ≥ 100
            break
        end
    end

    stats = Dict("η" => θ[1],
                 "β" => θ[2],
                 "ηList" => ηpList, 
                 "βList" => βpList,
                 "costList" => costList,
                 "iterList" => iterList,
                 "η" => θ[1],
                 "β" => θ[2])
    return stats
end

function val_check(v::Vector{Float64})
    sz = size(v,1)
    for i in 1:sz
        if v[i] < 0 
            println("Negative value: ", v[i])
            v[i] = abs(v[i])
        # elseif v[i] < 0
        #     println("Negative value: ", v[i])
        #     v[i] = 0.5
        end
    end
    return v
end

