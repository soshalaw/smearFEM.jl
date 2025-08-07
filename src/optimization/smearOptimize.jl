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

function init_cylinder()             
    scale = 100
    ne::Int = 4 # number of elements in the mesh for the ground truth
    FunctionClass::String = "Q2"
    camera_matrix = [[8*2048/7.07, 0.0, 2048/2] [0.0, 8*1536/5.3, 1536/2] [0.0, 0.0, 1.0]]'
    camera_pose = scale*[0 -0.25 2]'   # camera position in mm

    rList = Vector{Float64}(undef,0)
    hList = Vector{Float64}(undef,0)
    costList = Vector{Float64}(undef,0)
    iterList = Vector{Float64}(undef,0)

    NodeList, IEN, ID, IEN_top, IEN_bottom, IEN_side, nNodes, BorderNodes = meshgrid_cube(1, 1, 1, ne, FunctionClass=FunctionClass)

    #g_truth
    r_gt::Float64 = 0.25*scale  # radius of the cylinder in mm
    h_gt::Float64 = 0.5*scale  # height of the cylinder in mm
    NodeListCyl_gt = inflate_cylinder(NodeList, -0.5, 0.5, -0.5, 0.5, r_gt, h_gt)
    side_nodes = BorderNodes[1]
    obsBorderPts, g_SurfacePts2D = extract_borders(NodeListCyl_gt, camera_matrix, camera_pose, side_nodes, nNodes)

    # optimizer
    r = 1*scale*ones(ne)
    h = 1*scale
    NodeListCyl, ∇NodeListCyl = inflate_cylinder(NodeList, -0.5, 0.5, -0.5, 0.5, r, h, GRAD=true)
    simBorderPts, ∇BorderPts2D = extract_borders(NodeListCyl, camera_matrix, camera_pose, side_nodes, nNodes, GRAD=true, dqdθ=∇NodeListCyl, SIDES=false)

    d, ∂d, ∂2d, pairs = closest_point([simBorderPts],[obsBorderPts],[∇BorderPts2D])
    totdinit::Float64 = sum(d)/length(d)

    println("Ground truth : r = $r_gt, h = $h_gt")
    θ = vcat(r,h)
    iter::Int = 1
    c_grad::Float64 = 1.0
    println("Initial Error: ", totdinit)
    while true

        t∂2d = zeros(size(∂2d[1]))
        t∂d = zeros(size(∂d[1]))
        
        println("r: ", θ[1], " h: ", θ[2])

        len_d = length(d)
        szd = 1:len_d

        for i in szd
            t∂2d = t∂2d + ∂2d[i]
            t∂d = t∂d + ∂d[i]
        end

        p = t∂2d\t∂d
        α = 1
        
        θ = θ - α*p

        val_check(θ)
        # update the model parameters
        r = θ[1]
        h = θ[2]

        NodeListCyl, ∇NodeListCyl = inflate_cylinder(NodeList, -0.5, 0.5, -0.5, 0.5, r, h, GRAD=true)
        simBorderPts, ∇BorderPts2D = extract_borders(NodeListCyl, camera_matrix, camera_pose, side_nodes, nNodes, GRAD=true, dqdθ=∇NodeListCyl, SIDES=false)

        d, ∂d, ∂2d, pairs = closest_point([simBorderPts],[obsBorderPts],[∇BorderPts2D])

        totd = sum(d)/len_d
        c_grad = abs(totdinit - totd)
        totdinit = totd

        iter = iter + 1
        
        push!(rList,θ[1]) 
        push!(hList,θ[2])
        push!(costList,totd)
        push!(iterList,iter)
        println("iteration $iter: steps : $p, Error = $totd, Error gradient : $c_grad")
        println(" ... ")

        if c_grad < 0.005
            break
        elseif iter ≥ 100
            break
        end

    end
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
    totdinit::Float64 = sum(d)/length(d)
    
    push!(ηpList,θ[1])
    push!(βpList,θ[2])
    push!(costList,totdinit)
    push!(iterList,1)

    iter::Int = 1
    c_grad::Float64 = 1.0
    println("Initial Error: ", totdinit)
    while true

        reset_model!(model)
        t∂2d = zeros(size(∂2d[1]))
        t∂d = zeros(size(∂d[1]))
        
        println("η: ", θ[1], " β: ", θ[2])

        len_d = length(d)
        szd = 1:len_d

        for i in szd
            t∂2d = t∂2d + ∂2d[i]
            t∂d = t∂d + ∂d[i]
        end

        p = t∂2d\t∂d
        α = 1
        
        println("iteration $iter: ", p)
        θ = θ - α*p

        val_check(θ)
        # update the model parameters
        model.η = [θ[1]]
        scene.β = [θ[2]]
        # simulate the model
        μ_list, gradList, simBorderPts, splinex, spliney, pos2D = simulate(model, scene, conditions)

        # get the cost, gradient and hessian values
        d, ∂d, ∂2d, pairs = closest_point(simBorderPts, obsBorderPts, gradList)

        totd = sum(d)/len_d
        c_grad = abs(totdinit - totd)
        totdinit = totd

        iter = iter + 1
        
        push!(ηpList,θ[1]) 
        push!(βpList,θ[2])
        push!(costList,totd)
        push!(iterList,iter)
        println("iteration $iter: steps : $p, Error = $totd, Error gradient : $c_grad")
        println(" ... ")

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
# animate_fields(filepath = string(conditions.filepath,"/Results/images"),fields2D=borderPts2DList, pObs=splinexObs, qObs=splineyObs)

init_cylinder()