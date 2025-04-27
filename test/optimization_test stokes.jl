using smearFEM
using Test

CameraMatrix = [[8*2048/7.07, 0.0, 2048/2] [0.0, 8*1536/5.3, 1536/2] [0.0, 0.0, 1.0]]'
file = "/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/sim_experiments/cost_function_test/optimization/test3"

# test case 
r = 0.5
h = 1
ne = 4
ndim = 3
FunctionClass_u = "Q2"
nDof_u = ndim  # number of degree of freedom per node
FunctionClass_p = "Q1"
nDof_p = 1  # number of degree of freedom per node

μu_tp = -0.02
μu_btm = 0
μu_side = 0

CameraMatrix = [[8*2048/7.07, 0.0, 2048/2] [0.0, 8*1536/5.3, 1536/2] [0.0, 0.0, 1.0]]'

endTime = 2
steps = 2
tSteps = endTime/steps
time = collect(range(start=tSteps,stop=endTime,step=tSteps)) # time vector

println("tStep ", tSteps)
η = 10.0
β = 100.0
Control = "force" 

Δη = 1e-8
Δβ = 1e-8

mesh_u = meshgrid_cylinder(r, h, ne, FunctionClass=FunctionClass_u)  # generate the mesh grid
mesh_p = meshgrid_cylinder(r, h, ne, FunctionClass=FunctionClass_p)  # generate the mesh grid

mdl = Stokes(ndim=ndim, mesh_u=mesh_u, nDof_u=nDof_u, mesh_p=mesh_p, nDof_p=nDof_p, η=η)
# q_tp, q_side, q_btm, C_uc = set_boundary_cond(mdl)

# A_bar = assemble_system_A(mdl)                   # assemble the stiffness matrix
# B = assemble_system_B(mdl)                   # assemble the stiffness matrix
# b = apply_boundary_conditions_stokes(mdl)           # apply the neumann boundary conditions

# # get ∂K/∂θ with finite (central) difference
# ΔKp = assemble_system(model)

# ΔcMatm = get_cMat(mode, (η-Δη), μ)
# model.cMat = ΔcMatm

# ΔKm = assemble_system(model)

# dKdη_approx = (ΔKp-ΔKm)/(2*Δη)

# iIter,jIter = size(ΔKp)
# for i in 1:iIter  
#     for j in 1:jIter
#         @test dKdη_approx[i,j] ≈ dKdη[i,j] atol=10^(-5)
#     end
# end

p, dp, model = simulate_single_tstep_stokes(r, h, ne, η, ndim, FunctionClass_u, FunctionClass_p, nDof_u, nDof_p, β, μu_tp, μu_btm, 
μu_side, GRAD=true)

Nodes = model.mesh_u.NodeList + p

# estimate dp with finite (central) difference
Δp_ηp, model = simulate_single_tstep_stokes(r, h, ne, (η+Δη), ndim, FunctionClass_u, FunctionClass_p, nDof_u, nDof_p, β, μu_tp, μu_btm, 
μu_side)
ΔNodesηp = model.mesh_u.NodeList + Δp_ηp

Δp_ηm, model = simulate_single_tstep_stokes(r, h, ne, (η-Δη), ndim, FunctionClass_u, FunctionClass_p, nDof_u, nDof_p, β, μu_tp, μu_btm, 
μu_side)
ΔNodesηm = model.mesh_u.NodeList + Δp_ηm

# estimate dp with finite (central) difference
Δp_βp, model = simulate_single_tstep_stokes(r, h, ne, η, ndim, FunctionClass_u, FunctionClass_p, nDof_u, nDof_p, (β+Δβ), μu_tp, μu_btm, 
μu_side)
ΔNodesβp = model.mesh_u.NodeList + Δp_βp

Δp_βm, model = simulate_single_tstep_stokes(r, h, ne, η, ndim, FunctionClass_u, FunctionClass_p, nDof_u, nDof_p, (β-Δβ), μu_tp, μu_btm, 
μu_side)
ΔNodesβm = model.mesh_u.NodeList + Δp_βm

dηp_approx = (ΔNodesηp-ΔNodesηm)/(2*Δη)
dβp_approx = (ΔNodesβp-ΔNodesβm)/(2*Δβ)

iIter,jIter = size(dηp_approx)        

for i in 1:iIter  
    for j in 1:jIter
        @test dp[i,j,1] ≈ dηp_approx[i,j] atol=10^(-5)
    end
end

for i in 1:iIter  
    for j in 1:jIter
        @test dp[i,j,2] ≈ dβp_approx[i,j] atol=10^(-5)
    end
end

## Testing ∇u(θ)
BorderPts2D, dudη, SurfacePts2D, ∇SurfacePts2D = extract_borders(Nodes, CameraMatrix, mdl.mesh_u.side_nodes, GRAD=true, dqdθ=dp, SIDES=false)

# estimate dudη with finite (central) difference
ΔBorderPts2Dηp, ΔSurfacePts2Dηp = extract_borders(ΔNodesηp, CameraMatrix, mdl.mesh_u.side_nodes, GRAD=false)
ΔBorderPts2Dηm, ΔSurfacePts2Dηm = extract_borders(ΔNodesηm, CameraMatrix, mdl.mesh_u.side_nodes, GRAD=false)

ΔBorderPts2Dβp, ΔSurfacePts2Dβp = extract_borders(ΔNodesβp, CameraMatrix, mdl.mesh_u.side_nodes, GRAD=false)
ΔBorderPts2Dβm, ΔSurfacePts2Dβm = extract_borders(ΔNodesβm, CameraMatrix, mdl.mesh_u.side_nodes, GRAD=false)

∇SurfacePts2dη_approx = (ΔSurfacePts2Dηp - ΔSurfacePts2Dηm)/(2*Δη)
∇SurfacePts2Dβ_approx = (ΔSurfacePts2Dβp - ΔSurfacePts2Dβm)/(2*Δβ)

iIter,jIter = size(∇SurfacePts2dη_approx)
for i in 1:iIter  
    for j in 1:jIter
        @test ∇SurfacePts2D[i,j,1] ≈ ∇SurfacePts2dη_approx[i,j] atol=10^(-4)
    end
end

for i in 1:iIter  
    for j in 1:jIter
        @test ∇SurfacePts2D[i,j,2] ≈ ∇SurfacePts2Dβ_approx[i,j] atol=10^(-4)
    end
end

dudη_approx = (ΔBorderPts2Dηp - ΔBorderPts2Dηm)/(2*Δη)
dudβ_approx = (ΔBorderPts2Dβp - ΔBorderPts2Dβm)/(2*Δβ)

iIter,jIter = size(dudη_approx)
for i in 1:iIter  
    for j in 1:jIter
        @test dudη[i,j,1] ≈ dudη_approx[i,j] atol=10^(-4)
    end
end

for i in 1:iIter  
    for j in 1:jIter
        @test dudη[i,j,2] ≈ dudβ_approx[i,j] atol=10^(-4)
    end
end

# println("dudη: ")
# display(dudη[:,:,1])
# println("dudβ: ")
# display(dudη[:,:,2])

p_gt, model = simulate_single_tstep_stokes(r, h, ne, η, ndim, FunctionClass_u, FunctionClass_p, nDof_u, nDof_p, β, μu_tp, μu_btm, 
μu_side)
Nodes_gt = model.mesh_u.NodeList + p_gt

BorderPts2D_gt = back_project(Nodes_gt, CameraMatrix)

costList, ∇d, ∇2d, pairsList = closest_point([BorderPts2D],[BorderPts2D_gt],[dudη])
ΔcostList_ηp, pairsList_ηp = closest_point([ΔBorderPts2Dηp],[BorderPts2D_gt])
ΔcostList_ηm, pairsList_ηm = closest_point([ΔBorderPts2Dηm],[BorderPts2D_gt])

ΔcostList_βp, pairsList_βp = closest_point([ΔBorderPts2Dβp],[BorderPts2D_gt])
ΔcostList_βm, pairsList_βm = closest_point([ΔBorderPts2Dβm],[BorderPts2D_gt])

∇dη_approx = (ΔcostList_ηp - ΔcostList_ηm)/(2*Δη)
∇dβ_approx = (ΔcostList_βp - ΔcostList_βm)/(2*Δβ)

println(size(∇d[1]))
iIter = size(∇d[1],1)

@test ∇d[1][1] ≈ ∇dη_approx[1] atol=10^(-4)
@test ∇d[1][2] ≈ ∇dβ_approx[1] atol=10^(-4)

conditions = Conditions(CameraMatrix=CameraMatrix)
model, scene = def_problem(r, h, ne, η, ndim, FunctionClass_u, nDof_u, FunctionClass_p, nDof_p, β, F, control, viscosity_type, sim_time, t_steps)

## testing Σ∇p(θ)
μ_list, gradList, simBorderPts, splinex, spliney, SurfacePt2D = test_stokes(x0, x1, y0, y1, z0, z1, ne, η, ndim, FunctionClass_u, nDof_u, FunctionClass_p, 
                                                                                nDof_p, β, CameraMatrix, endTime, tSteps, Control)

μ_list, simBorderPts, ΔsimBorderPts_pL, splinex, spliney, mdl, ΔSurfacePt2D_pL = test_stokes(x0, x1, y0, y1, z0, z1, ne, (η+Δη), ndim, FunctionClass_u, nDof_u, FunctionClass_p, 
                                                                                nDof_p, β, CameraMatrix, endTime, tSteps, Control)

μ_list, simBorderPts, ΔsimBorderPts_mL, splinex, spliney, mdl, ΔSurfacePt2D_mL = test_stokes(x0, x1, y0, y1, z0, z1, ne, (η-Δη), ndim, FunctionClass_u, nDof_u, FunctionClass_p, 
                                                                                nDof_p, β, CameraMatrix, endTime, tSteps, Control)

μ_list, simBorderPts, ΔsimBorderPts_pβ, splinex, spliney, mdl, ΔSurfacePt2D_mL = test_stokes(x0, x1, y0, y1, z0, z1, ne, η, ndim, FunctionClass_u, nDof_u, FunctionClass_p, 
                                                                                nDof_p, (β+Δβ), CameraMatrix, endTime, tSteps, Control)

μ_list, simBorderPts, ΔsimBorderPts_mβ, splinex, spliney, mdl, ΔSurfacePt2D_mL = test_stokes(x0, x1, y0, y1, z0, z1, ne, η, ndim, FunctionClass_u, nDof_u, FunctionClass_p, 
                                                                                nDof_p, (β-Δβ), CameraMatrix, endTime, tSteps, Control)
                                                                                    
titer = 1:length(ΔSurfacePt2D_mL)
for t in titer
    println("time: ", t)
    grad = gradList[t]

    grad_approx_η = (ΔsimBorderPts_pL[t]- ΔsimBorderPts_mL[t])/(2*Δη)

    # println("dudη: ")
    # display(grad[:,:,1])

    # println("dudη approx: ")
    # display(grad_approx_η)

    grad_approx_β = (ΔsimBorderPts_pβ[t]- ΔsimBorderPts_mβ[t])/(2*Δβ)
    # println("dudβ : ")
    # display(grad[:,:,2])

    # println("dudβ approx: ")
    # display(grad_approx_β)

    pIter,qIter = size(grad_approx_β)
    for i in 1:pIter  
        for j in 1:qIter
            @test grad[i,j,1] ≈ grad_approx_η[i,j] atol=10^(-3)
            @test grad[i,j,2] ≈ grad_approx_β[i,j] atol=10^(-3)
        end
    end
end

