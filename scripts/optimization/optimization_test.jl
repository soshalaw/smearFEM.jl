using LinearAlgebra
using ProgressMeter
using SparseArrays
using Plots

using smearFEM

"""
Match the points of the simulated and observed data

Parameters:
pSim: {Vector{Float64}} : x coordinates of the simulated data
pObs: {Vector{Float64}} : x coordinates of the observed data

Returns:
pairs: {Vector{Vector{Int}}} : pairs of the matched points
cost: {Float64} : cost of the matching
"""
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
            error = sqrt((pi - pj)^2 + (qi - qj)^2)
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

function get_gradient(simborderfields, obsborderfields, pairs, model)

    K = assemble_system(model)                 # assemble the stiffness matrix
    b = apply_boundary_conditions(model)       # apply the neumann boundary conditions
    q_d = (μ_btm*q_btm + μ_tp*q_tp)          # apply the Dirichlet boundary conditions

    K_bar = K + β*b

    invK_bar = inv(Matrix(K_bar))

    C_T = transpose(C_uc)               # transpose the constraint matrix
    K_free = C_T*K_bar*C_uc             # extract the free part of the stiffness matrix

    invK = inv(Matrix(K_free))

    q_f = invK*C_T*(-K_bar*q_d)         # solve the system of equations
    q = q_d + C_uc*q_f;                 # assemble the solution 

    # getting the ∂u/∂λ
    model.cMat = get_cMat("lame",λ = 1, μ = 0)
    ∂K = assemble_system(model)           # assemble the stiffness matrix
    ∂q = -invK_bar*∂K*q

    # getting the ∂^2u/∂λ^2
    ∂2q  = 2*invK_bar*∂K*∂q

    # post process the solution rearrange the solution into the nodal displacements
    u = [q[ID[:,1]] q[ID[:,2]] q[ID[:,3]]]'    
    ∂u = [∂q[ID[:,1]] ∂q[ID[:,2]] ∂q[ID[:,3]]]'    
    ∂2u = [∂2q[ID[:,1]] ∂2q[ID[:,2]] ∂2q[ID[:,3]]]'  

    # select border nodes
    uBorder = [u[:,BorderNodesList[1]] u[:,BorderNodesList[2]] u[:,BorderNodesList[3]]]
    ∂uBorder = [∂u[:,BorderNodesList[1]] ∂u[:,BorderNodesList[2]] ∂u[:,BorderNodesList[3]]]
    ∂2uBorder = [∂2u[:,BorderNodesList[1]] ∂2u[:,BorderNodesList[2]] ∂2u[:,BorderNodesList[3]]]

    # selection matrix
    S = zeros(size(∂2q))
    S[1,1:3:end] .= 1
    S[2,2:3:end] .= 1
    S[3,3:3:end] .= 1
end

function test(x0, x1, y0, y1, z0, z1, ne, Young, ν, ndim, FunctionClass, nDof, β, CameraMatrix, endTime, tSteps, Control; writeData=false, filepath=nothing)
    
    if writeData
        isnothing(filepath) || AssertionError("Please provide a filepath to write the data")
        set_file(filepath)
    end

    μ_tp = 0.3

    if Control == "force"
        cParam = -0.25*ones(tSteps)
    elseif Control == "displacement"
        cParam = -(μ_tp/tSteps)*ones(tSteps)
    end

    cMat = get_cMat("standard", Young = Young, ν = ν)

    μ_list, simBorderPts, simBorderNodes, splinex, spliney, mdl = simulate(x0, x1, y0, y1, z0, z1, ne, Young, ν, ndim, FunctionClass, nDof, β, CameraMatrix, endTime, tSteps, Control, cParam, cMat, writeData=writeData, filepath=filepath)

    if !writeData
        println("comparing the simulated and observed data")
        isnothing(filepath) || AssertionError("Please provide a filepath to read the data")
        # animate_fields(filepath = string(filepath,"/Results"), p=splinex, q=spliney, pObs=splinexObs, qObs=splineyObs)
        obsBorderPts, splinexObs, splineyObs = readCSV(string(filepath,"/Results/contour_data"))                         # read the observation data
        
        animate_fields(filepath = string(filepath,"/Results/cost"), BorderNodes2D=simBorderPts, p=splinex, q=spliney, pObs=splinexObs, qObs=splineyObs) # animate the fields

        pairs = match_points(simBorderPts[1], [splinexObs[1],splineyObs[1]]) # match the points using the first border

        plot_matches(simBorderPts, splinex, spliney, splinexObs, splineyObs, pairs, string(filepath,"/Results/cost"))

        # test the closest point function
        hcost, xObsintlst, xSimintlst, ySimintlst = height_sample(simBorderPts, obsBorderPts)
        plot_matches_h(xObsintlst, ySimintlst, xSimintlst, splinex, spliney, splinexObs, splineyObs, string(filepath,"/Results/cost"))

        d_cp = closest_point(simBorderPts, obsBorderPts, pairs)

        return hcost, d_cp
    else
        return 0, 0
    end
end
function write_sim_data(x0, x1, y0, y1, z0, z1, ne, Youngtst, νtst, ndim, FunctionClass, nDof, β, CameraMatrix, endTime, tSteps, Control, filename)
    writeData = true
    hcost, cpCost = test(x0, x1, y0, y1, z0, z1, ne, Youngtst, νtst, ndim, FunctionClass, nDof, β, CameraMatrix, endTime, tSteps, Control, writeData=writeData, filepath=filename)
end

function main()
    # test case 
    x0 = 0
    x1 = 1
    y0 = 0
    y1 = 1
    z0 = 0
    z1 = 1
    ne = 6
    ndim = 3
    FunctionClass = "Q1"
    nDof = ndim  # number of degree of freedom per node
    β = 100
    CameraMatrix = [[8*2048/7.07, 0.0, 2048/2] [0.0, 8*1536/5.3, 1536/2] [0.0, 0.0, 1.0]]
    endTime = 60
    tSteps = 70
    
    Control = "displacement" # "force" or "displacement"
    filepath = "/home/soshala/SMEAR-PhD/SMEAR/Data/sim_experiments/"

    Youngtst = 40
    νtst = 0.4

    filepathi = string(filepath,"experiment_single_run")
    write_sim_data(x0, x1, y0, y1, z0, z1, ne, Youngtst, νtst, ndim, FunctionClass, nDof, β, CameraMatrix, endTime, tSteps, Control,filepathi)

    Young = Youngtst
    ν = νtst-0.2

    hcost, cpCost = test(x0, x1, y0, y1, z0, z1, ne, Young, ν, ndim, FunctionClass, nDof, β, CameraMatrix, endTime, tSteps, Control, filepath=filepathi)

    plot( cpCost, label="Height Cost")
    xlabel!("Time steps")
    ylabel!("Cost")
    savefig(string(filepathi,"/Results/cost/cost_cp.png"))

    plot(hcost, label="Closest Point Cost")
    xlabel!("Time steps")
    ylabel!("Cost")
    savefig(string(filepathi,"/Results/cost/cost_height.png"))
end

main()



