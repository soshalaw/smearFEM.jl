using LinearAlgebra
using Plots

function match_points(pSim,pObs)
    pSim, qSim = pSim[1,:], pSim[2,:]
    pObs, qObs = pObs[1], pObs[2]
    pairs = []
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

function closest_point(simScene, obsScene, pairs)
    # Define the cost function

    costList = []
    for (obsData, simData) in zip(obsScene, simScene) # iterate over the scenes
        tcost = 0
        # pSim, qSim = fit_curve(border=obsData)
        # println(size(obsData))
        pSim, qSim = simData[1,:], simData[2,:]
        pObs, qObs = obsData[1,:], obsData[2,:]

        for pair in pairs
            tcost += (pSim[pair[1]] - pObs[pair[2]])^2 + (qSim[pair[1]] - qObs[pair[2]])^2
        end
        # ObsCounter = 1
        # PointCounter = 1

        # for (pi,qi) in zip(pSim, qSim)
        #     cost = 1000000
        #     closestPointIdx = 0
        #     for (pj,qj) in zip(pObs, qObs)
        #         error = sqrt((pi - pj)^2 + (qi - qj)^2)
        #         if error < cost
        #             cost = error
        #             closestPointIdx = ObsCounter
        #         end
        #         ObsCounter += 1
        #     end
        #     push!(pairs, [PointCounter, closestPointIdx])
        #     tcost += cost
        #     PointCounter += 1
        # end
        mCost = tcost/length(pairs)
        push!(costList, mCost)
    end

    return costList
end

function height_sample(simScene,obsScene)
    xSimintlst = []
    ySimintlst = []
    xObsintlst = []

    costList = []
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

function get_grad(model,θ,u)

    λ = θ[1]
    μ = θ[2]

    model.cMat = [[ 2*μ+λ  λ    λ    0 0 0]; 
                  [  λ   2*μ+λ  λ    0 0 0]; 
                  [  λ     λ   2*μ+λ 0 0 0]; 
                  [  0     0    0    μ 0 0]; 
                  [  0     0    0    0 μ 0]
                  [  0     0    0    0 0 μ]]  # constitutive matrix

    K = assemble_system(model)    

    model.dcMatdλ = [[ 1  1  1  0 0 0]; 
               [ 1  1  1  0 0 0]; 
               [ 1  1  1  0 0 0]; 
               [ 0  0  0  0 0 0]; 
               [ 0  0  0  0 0 0]
               [ 0  0  0  0 0 0]]  # constitutive matrix

    model.dcMatdμ = [[ 2  0  0  0 0 0]; 
               [ 0  2  0  0 0 0]; 
               [ 0  0  2  0 0 0]; 
               [ 0  0  0  1 0 0]; 
               [ 0  0  0  0 1 0]
               [ 0  0  0  0 0 1]]  # constitutive matrix

    dKdθ = assemble_system(model)

    invK = inv(K)

    grad = -invK*dKdθ*u

    return grad
end


function fit_model(simborderfields, obsborderfields, model)

    opt = GradDescent(η=0.01)

    ν = 0.2
    Young = 40

    λ = Young*ν/((1+ν)*(1-2*ν))
    μ = Young/(2*(1+ν))

    θ = [λ, μ]

    # TODO add simulate function to simulate the model
    simborderfields = simulate()
    ∇f = get_grad(model, θ, simborderfields)
    
    for i in 1:100
        x = update(opt, x, ∇f)
    end

end
