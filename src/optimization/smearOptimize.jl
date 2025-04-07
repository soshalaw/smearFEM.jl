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
        cost = 1000000
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
    for (obs_t, sim_t) in zip(obsScene[1:end], simScene[1:end]) # iterate over the scenes
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

    # Define the cost function 
    costList = Float64[]
    dcostList = []
    dcost2List = []
    pairsList = []
    @argcheck length(simScene) == length(obsScene) "Size of the simulation and observation scenes should be the same"
    @argcheck length(simScene) == length(dudθ) "Size of the simulation and observation scenes should be the same"

    for (obs_t, sim_t, du_tdθ) in zip(obsScene, simScene, dudθ) # iterate over the scenes
        @argcheck size(sim_t,2) == size(du_tdθ,2) "Number of the border points and the gradient points should be the same"
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
