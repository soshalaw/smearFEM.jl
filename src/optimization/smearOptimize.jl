using LinearAlgebra
using Plots

function match_points(pSim,pObs)
    pSim, qSim = pSim[1,:], pSim[2,:]
    pObs, qObs = pObs[1,:], pObs[2,:]
    pairs = AbstractArray[]
    PointCounter = 1
    for (pi,qi) in zip(pSim, qSim)
        cost = 1000000
        closestPointIdx = 0
        ObsCounter = 1
        for (pj,qj) in zip(pObs, qObs)
            error = (pi - pj)^2 + (qi - qj)^2
            if error < cost
                cost = error
                closestPointIdx = ObsCounter
            end
            ObsCounter += 1
        end
        push!(pairs, [PointCounter, closestPointIdx])
        PointCounter += 1
    end
    return pairs
end

function closest_point(simScene, obsScene)
    # Define the cost function
    costList = Float64[]
    pairsList = []
    for (obsData, simData) in zip(obsScene[1:end], simScene[1:end]) # iterate over the scenes
        tcost = 0

        pairs = match_points(simData, obsData) # match the points using the first border

        pSim, qSim = simData[1,:], simData[2,:]
        pObs, qObs = obsData[1,:], obsData[2,:]

        for pair in pairs
            tcost += (pSim[pair[1]] - pObs[pair[2]])^2 + (qSim[pair[1]] - qObs[pair[2]])^2
        end

        mCost = tcost/length(pairs) # Σ(xi-x_obs)/n (mean error)
        push!(costList, mCost)
        push!(pairsList,pairs)
    end
    return costList, pairsList
end

function height_sample(simScene,obsScene)
    xSimintlst = AbstractArray[]
    ySimintlst = AbstractArray[]
    xObsintlst = AbstractArray[]
    costList = Float64[]
    for (obsData, simData) in zip(obsScene, simScene) # iterate over the scenes

        xSim, ySim = filter_points(simData, 2048/2)
        xObs, yObs = filter_points(obsData, 2048/2)

        # hrange = ySim[1]:1:ySim[end]

        # xSimint = fit_curve(borderx=xSim, bordery=ySim, samples=ySim)
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

function get_grad(model)

    dcMatdλ = [[ 1  1  1  0 0 0]; 
               [ 1  1  1  0 0 0]; 
               [ 1  1  1  0 0 0]; 
               [ 0  0  0  0 0 0]; 
               [ 0  0  0  0 0 0]
               [ 0  0  0  0 0 0]]  # constitutive matrix

    dcMatdμ = [[ 2  0  0  0 0 0]; 
               [ 0  2  0  0 0 0]; 
               [ 0  0  2  0 0 0]; 
               [ 0  0  0  1 0 0]; 
               [ 0  0  0  0 1 0]
               [ 0  0  0  0 0 1]]  # constitutive matrix

    dKdλ = assemble_system(model, ∂C=dcMatdλ)
    dKdμ = assemble_system(model, ∂C=dcMatdμ)

    invK = inv(K)

    dudλ = -invK*dKdλ
    dudμ = -invK*dKdμ

    return dudλ, dudμ
end


function fit_model(simborderfields, obsborderfields, model)

    opt = GradDescent(η=0.01)

    ν = 0.2
    Young = 40

    λ = Young*ν/((1+ν)*(1-2*ν))
    μ = Young/(2*(1+ν))

    dudλ, dudμ = get_grad(model)

    # TODO add simulate function to simulate the model
    # simborderfields = simulate()

    
    for i in 1:100
        x = update(opt, x, ∇f)
    end

end
