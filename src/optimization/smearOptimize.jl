using LinearAlgebra
using Plots
using ArgCheck

function match_points(pSim::AbstractMatrix{Float64},pObs::AbstractMatrix{Float64})

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

function closest_point(simScene::Vector{AbstractArray}, obsScene::Vector{AbstractArray})
    # Define the cost function
    costList = Float64[]
    pairsList = []
    for (obs_t, sim_t) in zip(obsScene[1:end], simScene[1:end]) # iterate over the scenes
        tcost = 0

        pairs = match_points(simData, obsData) # match the points using the first border

        pSim, qSim = sim_t[1,:], sim_t[2,:]
        pObs, qObs = obs_t[1,:], obs_t[2,:]

        for pair in pairs
            tcost += (pSim[pair[1]] - pObs[pair[2]])^2 + (qSim[pair[1]] - qObs[pair[2]])^2
        end

        mCost = tcost/length(pairs) # Σ(xi-x_obs)/n (mean error)
        push!(costList, mCost)
        push!(pairsList,pairs)
    end
    return costList, pairsList
end

function closest_point(simScene::Vector{AbstractArray}, obsScene::Vector{AbstractArray}, dudλ::Vector{AbstractArray})
    # Define the cost function
    costList = Float64[]
    dcostList = AbstractMatrix[]
    dcost2List = AbstractMatrix[]
    pairsList = []
    @argcheck length(simScene) == length(obsScene)
    @argcheck length(simScene) == length(dudλ)

    for (obs_t, sim_t, du_tdλ) in zip(obsScene, simScene, dudλ) # iterate over the scenes
        tcost = 0.0
        dtcost = zeros(Float64,1,size(du_tdλ,2))
        dt2cost = zeros(Float64,size(du_tdλ,2),size(du_tdλ,2))
        
        pairs = match_points(sim_t, obs_t) # match the points using the first border
        
        pSim, qSim = sim_t[1,:], sim_t[2,:]
        dpSim, dqSim = du_tdλ[1,:], du_tdλ[2,:]
        pObs, qObs = obs_t[1,:], obs_t[2,:]

        for pair in pairs
            u = [(pSim[pair[1]] - pObs[pair[2]]) (qSim[pair[1]] - qObs[pair[2]])]
            du = [dpSim[pair[1]] dqSim[pair[1]]] 

            tcost =+ dot(u,u)
            dtcost =+ u*du'
            dt2cost =+ du*du'
        end

        mCost = tcost/length(pairs) # 1/m(Σ(√(xi-x_obs)^2+(yi-y_obs)^2))^2 (mean error)
        dmcost = dtcost/length(pairs) #1/m(Σ(xi-x_obs)∂x/∂θ_i +(yi-y_obs)∂y/∂θ_i)
        dm2cost = dt2cost/length(pairs) #1/m(Σ∂x2/∂2θ_i + ∂y2/∂2θ_i)

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

    K, dKdλ = assemble_system(model, GRAD=true)

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

function fit(x0, x1, y0, y1, z0, z1, ne, c1, c2, ndim::Int64, FunctionClass::String, nDof::Int64, β, CameraMatrix, endTime, tSteps, Control::String, 
    mode::String, ObsData, SIDES::Bool=false, PLOT::Bool=false, filepath=nothing)
    
    # test run with c1 = λ
    obsBorderPts = ObsData[1]

    μ_list, gradList, simBorderPts, splinex, spliney, mdl = test(x0, x1, y0, y1, z0, z1, ne, c1, c2, ndim, FunctionClass, nDof, β, CameraMatrix, 
                                                            endTime, tSteps, Control, mode=mode, SIDES=SIDES)    

    # test the closest point function
    d, ∂d, ∂d2, pairs = closest_point(simBorderPts, obsBorderPts, gradList)

    p = ∂d/∂d2
    α = 0.1

    c1 = c1 + α*p
end

# grad_check()