using smearFEM
using Test
using LinearAlgebra

file = "/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/sim_experiments/cost_function_test/optimization/test3"

# test case 
scale = 50
r::Float64 = 0.5*scale  # radius of the cylinder in mm
h::Float64 = 1*scale  # height of the cylinder in mm
ndim::Int = 3
FunctionClass_x::String = "S2"
FunctionClass_u::String = "Q2"
nDof_u::Int = ndim  # number of degree of freedom per node
FunctionClass_p::String = "Q1"
nDof_p::Int = 1  # number of degree of freedom per node

camera_matrix = [[8*2048/7.07, 0.0, 2048/2] [0.0, 8*1536/5.3, 1536/2] [0.0, 0.0, 1.0]]'
camera_pose = scale*[0 -0.5 2.75]'   # camera position in mm

μu_tp = -0.02*scale  # top boundary condition in mm/s
μu_btm = 0
μu_side = 0

sim_time = 5.0
steps = 50.0
t_steps = sim_time/steps
time = collect(range(start=t_steps,stop=sim_time,step=t_steps)) # time vector

println("time step size", t_steps)
control = "force" 
viscosity_type = "constant" # "constant" or "bulk_viscosity"
β::Float64 = 100.0
η::Float64 = 40.0
F_ext::Float64 = 250000.0
F::Vector{Float64} = -F_ext*ones(Float64, round(Int, (sim_time/t_steps)))
# F::Float64 = 3.0
ne = 4
Δη = 1e-8
Δβ = 1e-8

filePath = "/home/soshala/SMEAR-PhD/smear-modules/smearFEM.jl/cylindergen"

mesh_x = meshgrid_cylinder(r, h, ne, FunctionClass=FunctionClass_x, filePath=filePath)  # generate the mesh grid for geometry
mesh_u = meshgrid_cylinder(r, h, ne, FunctionClass=FunctionClass_u)  # generate the mesh grid
mesh_p = meshgrid_cylinder(r, h, ne, FunctionClass=FunctionClass_p)  # generate the mesh grid

mdl = Stokes(ndim=ndim, mesh_x=mesh_x, mesh_u=mesh_u, nDof_u=nDof_u, mesh_p=mesh_p, nDof_p=nDof_p, η=[η])
# q_tp, q_side, q_btm, C_uc = set_boundary_cond(mdl)

# A_bar = assemble_system_A(mdl)                   # assemble the stiffness matrix
# B = assemble_system_B(mdl)                   # assemble the stiffness matrix
# b = apply_boundary_conditions_stokes(mdl)           # apply the neumann boundary conditions

# get ∂K/∂θ with finite (central) difference
# ΔKp = assemble_system(mdl)

# ΔcMatm = get_cMat(mode, (η-Δη), μ)
# model.cMat = ΔcMatm

# ΔKm = assemble_system(mdl)

# dKdη_approx = (ΔKp-ΔKm)/(2*Δη)

# iIter,jIter = size(ΔKp)
# for i in 1:iIter  
#     for j in 1:jIter
#         @test dKdη_approx[i,j] ≈ dKdη[i,j] atol=10^(-5)
#     end
# end

T = Matrix{Float64}(I, size(mesh_u.NodeList,2), size(mesh_u.NodeList,2))
if mesh_x.FunctionClass == "S2" && mesh_u.FunctionClass != "S2"
    T = get_nurbs_2_lagrange_proj(mesh_x.IEN, mesh_u.IEN, mesh_x.C_vol, mesh_x.NodeList, mesh_x.W)
end
T_ = T'*inv(T*T')

p_, dp_, model = simulate_single_tstep_stokes(r, h, ne, η, ndim, FunctionClass_u, FunctionClass_p, nDof_u, nDof_p, β, μu_tp, μu_btm, 
μu_side, FunctionClass_x=FunctionClass_x, GRAD=true)
p = p_*T_
dp = similar(dp_)
dp[:,:,1] = (dp_[:,:,1]*T_)*T
dp[:,:,2] = (dp_[:,:,2]*T_)*T
Nodes_ = model.mesh_x.NodeList + p
Nodes = Nodes_*T

# estimate dp with finite (central) difference
Δp_ηp_, model = simulate_single_tstep_stokes(r, h, ne, (η+Δη), ndim, FunctionClass_u, FunctionClass_p, nDof_u, nDof_p, β, μu_tp, μu_btm, 
μu_side, FunctionClass_x=FunctionClass_x)
Δp_ηp = Δp_ηp_*T_
ΔNodesηp_ = model.mesh_x.NodeList + Δp_ηp
ΔNodesηp = ΔNodesηp_*T

Δp_ηm_, model = simulate_single_tstep_stokes(r, h, ne, (η-Δη), ndim, FunctionClass_u, FunctionClass_p, nDof_u, nDof_p, β, μu_tp, μu_btm, 
μu_side, FunctionClass_x=FunctionClass_x)
Δp_ηm = Δp_ηm_*T_
ΔNodesηm_ = model.mesh_x.NodeList + Δp_ηm
ΔNodesηm = ΔNodesηm_*T

# estimate dp with finite (central) difference
Δp_βp_, model = simulate_single_tstep_stokes(r, h, ne, η, ndim, FunctionClass_u, FunctionClass_p, nDof_u, nDof_p, (β+Δβ), μu_tp, μu_btm, 
μu_side, FunctionClass_x=FunctionClass_x)
Δp_βp = Δp_βp_*T_
ΔNodesβp_ = model.mesh_x.NodeList + Δp_βp
ΔNodesβp = ΔNodesβp_*T

Δp_βm_, model = simulate_single_tstep_stokes(r, h, ne, η, ndim, FunctionClass_u, FunctionClass_p, nDof_u, nDof_p, (β-Δβ), μu_tp, μu_btm, 
μu_side, FunctionClass_x=FunctionClass_x)
Δp_βm = Δp_βm_*T_
ΔNodesβm_ = model.mesh_x.NodeList + Δp_βm
ΔNodesβm = ΔNodesβm_*T

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
BorderPts2D, dudη, SurfacePts2D, ∇SurfacePts2D = extract_borders(Nodes, camera_matrix, camera_pose, BorderNodesList = mdl.mesh_u.side_nodes, GRAD=true, dqdθ=dp, SIDES=false)

# estimate dudη with finite (central) difference
ΔBorderPts2Dηp, ΔSurfacePts2Dηp = extract_borders(ΔNodesηp, camera_matrix, camera_pose, BorderNodesList=mdl.mesh_u.side_nodes, GRAD=false)
ΔBorderPts2Dηm, ΔSurfacePts2Dηm = extract_borders(ΔNodesηm, camera_matrix, camera_pose, BorderNodesList=mdl.mesh_u.side_nodes, GRAD=false)

ΔBorderPts2Dβp, ΔSurfacePts2Dβp = extract_borders(ΔNodesβp, camera_matrix, camera_pose, BorderNodesList=mdl.mesh_u.side_nodes, GRAD=false)
ΔBorderPts2Dβm, ΔSurfacePts2Dβm = extract_borders(ΔNodesβm, camera_matrix, camera_pose, BorderNodesList=mdl.mesh_u.side_nodes, GRAD=false)

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

# println("dudη approx: ")
# display(dudη_approx)

# println("dudβ: ")
# display(dudη[:,:,2])

# println("dudβ approx: ")
# display(dudβ_approx)

p_gt_, model = simulate_single_tstep_stokes(r, h, ne, η, ndim, FunctionClass_u, FunctionClass_p, nDof_u, nDof_p, β, μu_tp, μu_btm, 
μu_side, FunctionClass_x=FunctionClass_x)
p_gt = p_gt_*T_
Nodes_gt_ = model.mesh_x.NodeList + p_gt
Nodes_gt = Nodes_gt_*T

BorderPts2D_gt = back_project(Nodes_gt, camera_matrix, camera_pose)

costList, ∇d, ∇2d, pairsList = closest_point([BorderPts2D],[BorderPts2D_gt],[dudη])
ΔcostList_ηp, pairsList_ηp = closest_point([ΔBorderPts2Dηp],[BorderPts2D_gt])
ΔcostList_ηm, pairsList_ηm = closest_point([ΔBorderPts2Dηm],[BorderPts2D_gt])

ΔcostList_βp, pairsList_βp = closest_point([ΔBorderPts2Dβp],[BorderPts2D_gt])
ΔcostList_βm, pairsList_βm = closest_point([ΔBorderPts2Dβm],[BorderPts2D_gt])

∇dη_approx = (ΔcostList_ηp - ΔcostList_ηm)/(2*Δη)
∇dβ_approx = (ΔcostList_βp - ΔcostList_βm)/(2*Δβ)

iIter = size(∇d[1],1)

@test ∇d[1][1] ≈ ∇dη_approx[1] atol=10^(-4)
@test ∇d[1][2] ≈ ∇dβ_approx[1] atol=10^(-4)

conditions = Conditions(camera_matrix=camera_matrix, camera_pose=camera_pose)
model, scene = def_problem(r, h, ne, η, ndim, FunctionClass_u, nDof_u, FunctionClass_p, nDof_p, FunctionClass_x, β, F, control, viscosity_type, 
                        sim_time, t_steps)
## testing Σ∇p(θ)
μ_list, gradList, simBorderPts, splinex, spliney, SurfacePt2D = simulate(model, scene, conditions)
reset_model!(model)
model.η = [η+Δη]
μ_list, simBorderPts, ΔsimBorderPts_pL, splinex, spliney, ΔSurfacePt2D_pL = simulate(model, scene, conditions)
reset_model!(model)
model.η = [η-Δη]
μ_list, simBorderPts, ΔsimBorderPts_mL, splinex, spliney, ΔSurfacePt2D_mL = simulate(model, scene, conditions)
reset_model!(model)
scene.β = [β+Δβ]
μ_list, simBorderPts, ΔsimBorderPts_pβ, splinex, spliney, ΔSurfacePt2D_mL = simulate(model, scene, conditions)
reset_model!(model)
scene.β = [β-Δβ]
μ_list, simBorderPts, ΔsimBorderPts_mβ, splinex, spliney, ΔSurfacePt2D_mL = simulate(model, scene, conditions)
                                                                                    
titer = 1:length(ΔSurfacePt2D_mL)
for t in titer
    println("time: ", t)
    grad = gradList[t]

    grad_approx_η = (ΔsimBorderPts_pL[t]- ΔsimBorderPts_mL[t])/(2*Δη)

    println("dudη: ")
    display(grad[:,:,1])

    println("dudη approx: ")
    display(grad_approx_η)

    grad_approx_β = (ΔsimBorderPts_pβ[t]- ΔsimBorderPts_mβ[t])/(2*Δβ)

    println("dudβ : ")
    display(grad[:,:,2])

    println("dudβ approx: ")
    display(grad_approx_β)

    pIter,qIter = size(grad_approx_η)
    for i in 1:pIter  
        for j in 1:qIter
            @test grad[i,j,1] ≈ grad_approx_η[i,j] atol=10^(-4)
            @test grad[i,j,2] ≈ grad_approx_β[i,j] atol=10^(-4)
        end
    end
end

