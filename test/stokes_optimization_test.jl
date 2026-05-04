using smearFEM
using Test
using LinearAlgebra

file = resolve_data_path("sim_experiments/cost_function_test/optimization/test3")

# test case 
r::Float64 = 25  # radius of the cylinder in mm
h::Float64 = 40.0  # height of the cylinder in mm
ndim::Int = 3
FunctionClass_x::String = "Q2"
FunctionClass_u::String = "Q2"
nDof_u::Int = ndim  # number of degree of freedom per node
FunctionClass_p::String = "Q1"
nDof_p::Int = 1  # number of degree of freedom per node

obj_pose = zeros(Float64, 4,4)
obj_pose[1,1] = -1.0
obj_pose[2,3] = -1.0
obj_pose[3,2] = -1.0
obj_pose[1:3,4] = [0.0, h/2, 150.0]
camera_matrix::Matrix{Float64} = [2.39642674e+03  0.0  1.00429248e+03; 0.0  2.40565353e+03  7.57028161e+02; 0.0  0.0 1.0;]

μu_tp = -10.0  # top boundary condition in mm/s
μu_btm = 0
μu_side = 0

sim_time = 0.5
steps = 5
t_steps = sim_time/steps
time = collect(range(start=t_steps,stop=sim_time,step=t_steps)) # time vector

println("time step size", t_steps)
control = "force" 
viscosity_type = "constant" # "constant" or "bulk_viscosity"
β::Float64 = 100.0
η::Float64 = 100.0
F_ext::Float64 = 9.813e3 * 0.85
F::Vector{Float64} = -F_ext*ones(Float64, round(Int, (sim_time/t_steps)))
# F::Float64 = 3.0
ne = 6
Δη = 1e-4
Δβ = 1e-4
error_tol = 2.5e-4
rel_error_tol = 1e-4

filePath = joinpath(@__DIR__, "..", "cylindergen")

mesh_x = meshgrid_cylinder(r, h, ne, FunctionClass=FunctionClass_x)  # generate the mesh grid for geometry
mesh_u = meshgrid_cylinder(r, h, ne, FunctionClass=FunctionClass_u)  # generate the mesh grid
mesh_p = meshgrid_cylinder(r, h, ne, FunctionClass=FunctionClass_p)  # generate the mesh grid

mdl = Stokes(ndim=ndim, mesh_x=mesh_x, mesh_u=mesh_u, nDof_u=nDof_u, mesh_p=mesh_p, nDof_p=nDof_p, η=[η])

T = Matrix{Float64}(I, size(mesh_u.NodeList,2), size(mesh_u.NodeList,2))
if mesh_x.FunctionClass == "S2" && mesh_u.FunctionClass != "S2"
    T = get_nurbs_2_lagrange_proj(mesh_x.IEN, mesh_u.IEN, mesh_x.C_vol, mesh_x.NodeList, mesh_x.W)
end
T_ = T'*inv(T*T')

p_, dp_, model = simulate_single_tstep_stokes(r, h, ne, η, ndim, FunctionClass_u, FunctionClass_p, nDof_u, nDof_p, β, μu_tp, μu_btm, 
μu_side, FunctionClass_x=FunctionClass_x, GRAD=true)
p = p_*T_
dp = similar(dp_)
dp[:,:,1] = dp_[:,:,1]*T_
dp[:,:,2] = dp_[:,:,2]*T_
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
        @test dp[i,j,1] ≈ dηp_approx[i,j] atol=error_tol rtol=rel_error_tol
    end
end

for i in 1:iIter  
    for j in 1:jIter
        @test dp[i,j,2] ≈ dβp_approx[i,j] atol=error_tol rtol=rel_error_tol
    end
end

## Testing ∇u(θ)
BorderPts2D, dudη, SurfacePts2D, ∇SurfacePts2D = extract_borders(Nodes, camera_matrix, obj_pose, BorderNodesList = mdl.mesh_u.side_nodes, GRAD=true, dqdθ=dp, SIDES=false)

# estimate dudη with finite (central) difference
ΔBorderPts2Dηp, ΔSurfacePts2Dηp = extract_borders(ΔNodesηp, camera_matrix, obj_pose, BorderNodesList=mdl.mesh_u.side_nodes, GRAD=false)
ΔBorderPts2Dηm, ΔSurfacePts2Dηm = extract_borders(ΔNodesηm, camera_matrix, obj_pose, BorderNodesList=mdl.mesh_u.side_nodes, GRAD=false)

ΔBorderPts2Dβp, ΔSurfacePts2Dβp = extract_borders(ΔNodesβp, camera_matrix, obj_pose, BorderNodesList=mdl.mesh_u.side_nodes, GRAD=false)
ΔBorderPts2Dβm, ΔSurfacePts2Dβm = extract_borders(ΔNodesβm, camera_matrix, obj_pose, BorderNodesList=mdl.mesh_u.side_nodes, GRAD=false)

∇SurfacePts2dη_approx = (ΔSurfacePts2Dηp - ΔSurfacePts2Dηm)/(2*Δη)
∇SurfacePts2Dβ_approx = (ΔSurfacePts2Dβp - ΔSurfacePts2Dβm)/(2*Δβ)

iIter,jIter = size(∇SurfacePts2dη_approx)
for i in 1:iIter  
    for j in 1:jIter
        @test ∇SurfacePts2D[i,j,1] ≈ ∇SurfacePts2dη_approx[i,j] atol=error_tol rtol=rel_error_tol
    end
end

for i in 1:iIter  
    for j in 1:jIter
        @test ∇SurfacePts2D[i,j,2] ≈ ∇SurfacePts2Dβ_approx[i,j] atol=error_tol rtol=rel_error_tol
    end
end

dudη_approx = (ΔBorderPts2Dηp - ΔBorderPts2Dηm)/(2*Δη)
dudβ_approx = (ΔBorderPts2Dβp - ΔBorderPts2Dβm)/(2*Δβ)

iIter,jIter = size(dudη_approx)
for i in 1:iIter  
    for j in 1:jIter
        @test dudη[i,j,1] ≈ dudη_approx[i,j] atol=error_tol rtol=rel_error_tol
    end
end

for i in 1:iIter  
    for j in 1:jIter
        @test dudη[i,j,2] ≈ dudβ_approx[i,j] atol=error_tol rtol=rel_error_tol
    end
end


p_gt_, model = simulate_single_tstep_stokes(r, h, ne, η, ndim, FunctionClass_u, FunctionClass_p, nDof_u, nDof_p, β, μu_tp, μu_btm, 
μu_side, FunctionClass_x=FunctionClass_x)
p_gt = p_gt_*T_
Nodes_gt_ = model.mesh_x.NodeList + p_gt
Nodes_gt = Nodes_gt_*T

BorderPts2D_gt = back_project(Nodes_gt, camera_matrix, obj_pose)

costList, ∇d, ∇2d, pairsList = closest_point([BorderPts2D],[BorderPts2D_gt],[dudη])

ΔcostList_ηp, pairsList_ηp = closest_point([ΔBorderPts2Dηp],[BorderPts2D_gt])
ΔcostList_ηm, pairsList_ηm = closest_point([ΔBorderPts2Dηm],[BorderPts2D_gt])

ΔcostList_βp, pairsList_βp = closest_point([ΔBorderPts2Dβp],[BorderPts2D_gt])
ΔcostList_βm, pairsList_βm = closest_point([ΔBorderPts2Dβm],[BorderPts2D_gt])

∇dη_approx = (ΔcostList_ηp - ΔcostList_ηm)/(2*Δη)
∇dβ_approx = (ΔcostList_βp - ΔcostList_βm)/(2*Δβ)

iIter = size(∇d[1])

@test ∇d[1][1] ≈ ∇dη_approx[1] atol=error_tol rtol=rel_error_tol
@test ∇d[1][2] ≈ ∇dβ_approx[1] atol=error_tol rtol=rel_error_tol

conditions = Conditions(camera_matrix=camera_matrix, obj_pose=obj_pose)
model, scene = def_problem(r, h, ne, η, ndim, FunctionClass_u, nDof_u, FunctionClass_p, nDof_p, FunctionClass_x, β, F, control, viscosity_type, 
                        sim_time, t_steps, GMESH_MESH=false)
## testing Σ∇p(θ)
μ_list, gradList, simBorderPts, _, _, _, Δpos3D_, _, _, _, _, gradList_3d = simulate(model, scene, conditions)
reset_model!(model)
model.η = [η+Δη]
μ_list, simBorderPts, ΔsimBorderPts_pL, _, _, _, Δpos3D_pL, _, _, _, _, _= simulate(model, scene, conditions)
reset_model!(model)
model.η = [η-Δη]
μ_list, simBorderPts, ΔsimBorderPts_mL, _, _, _, Δpos3D_mL, _, _, _, _, _ = simulate(model, scene, conditions)
reset_model!(model)
scene.β = [β+Δβ]
μ_list, simBorderPts, ΔsimBorderPts_pβ, _, _, _, Δpos3D_pβ, _, _, _, _, _ = simulate(model, scene, conditions)
reset_model!(model)
scene.β = [β-Δβ]
μ_list, simBorderPts, ΔsimBorderPts_mβ, _, _, _, Δpos3D_mβ, _, _, _, _, _ = simulate(model, scene, conditions)
                                                                                    
titer = 1:length(Δpos3D_)

for t in titer
    println("time: ", t)
    grad_3d = gradList_3d[t]

    grad_approx_η_3d = (Δpos3D_pL[t]- Δpos3D_mL[t])/(2*Δη)

    # println("dxdη: ")
    # display(grad_3d[:,:,1])

    # println("dxdη approx: ")
    # display(grad_approx_η_3d)

    grad_approx_β_3d = (Δpos3D_pβ[t]- Δpos3D_mβ[t])/(2*Δβ)

    # println("dxdβ : ")
    # display(grad_3d[:,:,2])

    # println("dxdβ approx: ")
    # display(grad_approx_β_3d)

    pIter,qIter = size(grad_approx_η_3d)
    for i in 1:pIter  
        for j in 1:qIter
            @test grad_3d[i,j,1] ≈ grad_approx_η_3d[i,j] atol=error_tol rtol=rel_error_tol
            @test grad_3d[i,j,2] ≈ grad_approx_β_3d[i,j] atol=error_tol rtol=rel_error_tol
        end
    end
end

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

    pIter,qIter = size(grad_approx_η)
    for i in 1:pIter  
        for j in 1:qIter
            @test grad[i,j,1] ≈ grad_approx_η[i,j] atol=error_tol rtol=rel_error_tol
            @test grad[i,j,2] ≈ grad_approx_β[i,j] atol=error_tol rtol=rel_error_tol
        end
    end
end

