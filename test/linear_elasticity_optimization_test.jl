using smearFEM
using Test

CameraMatrix = [[8*2048/7.07, 0.0, 2048/2] [0.0, 8*1536/5.3, 1536/2] [0.0, 0.0, 1.0]]'
file = "/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/sim_experiments/cost_function_test/optimization/test3"
# testing the gradient of the projection function
a = [1.5 -2.3 0.5]'
Δ = 1e-8

# transform point cloud wrt to camera frame 
R = [1 0 0; 0 0 1; 0 -1 0]     # rotation matrix
t = [0; -0.5; 4]               # translation vector

aTrans = R*a .+ t

∇u = ∇π(aTrans,CameraMatrix)

u = project_to(aTrans, CameraMatrix)

Δxp = [(aTrans[1]+Δ)  aTrans[2]      aTrans[3]]'
Δyp = [aTrans[1]      (aTrans[2]+Δ)  aTrans[3]]'
Δzp = [aTrans[1]      aTrans[2]      (aTrans[3]+Δ)]'

Δu_xp = project_to(Δxp, CameraMatrix)
Δu_yp = project_to(Δyp, CameraMatrix)
Δu_zp = project_to(Δzp, CameraMatrix)

Δxm = [(aTrans[1]-Δ)  aTrans[2]      aTrans[3]]'
Δym = [aTrans[1]      (aTrans[2]-Δ)  aTrans[3]]'
Δzm = [aTrans[1]      aTrans[2]      (aTrans[3]-Δ)]'

Δu_xm = project_to(Δxm, CameraMatrix)
Δu_ym = project_to(Δym, CameraMatrix)
Δu_zm = project_to(Δzm, CameraMatrix)

Δup = [Δu_xp Δu_yp Δu_zp]
Δum = [Δu_xm Δu_ym Δu_zm]

∇u_approx = (Δup-Δum)/(2*Δ)

pIter,pIter = size(∇u)
for p in 1:pIter  
    for p in 1:pIter
        @test ∇u_approx[p,p] ≈ ∇u[p,p] atol=10^(-4)
    end
end

# check the gradient of K
x0 = 0
x1 = 1
y0 = 0
y1 = 1
z0 = 0
z1 = 1
ne = 4
ndim = 3
FunctionClass = "Q2"
nDof = ndim  # number of degree of freedom per node
endTime = 10
steps = 10
tSteps = endTime/steps
Control = "force"

μ_tp = -0.03
μ_btm = 0

Youngtst = 30
νtst = 0.4
λ = Youngtst*νtst/((1+νtst)*(1-2*νtst))
μ = Youngtst/(2*(1+νtst))
β = 100
Δλ = 1e-8
Δβ = 1e-8

NodeList, IEN, ID, IEN_top, IEN_btm, BorderNodesList = meshgrid_cube(x0,x1,y0,y1,z0,z1,ne,ndim,FunctionClass=FunctionClass)  # generate the mesh grid
NodeListCylinder = inflate_cylinder(NodeList, x0, x1, y0, y1)   

mode = "lame"    
cMat = get_cMat(mode, λ, μ) # E, ν or λ, μ
dcdλ = get_cMat(mode, 1.0 , 0.0)

model = def_model("linear_elasticity", ne=ne, NodeList=NodeListCylinder, IEN=IEN, IEN_top=IEN_top, IEN_btm=IEN_btm, ndim=ndim, nDof=nDof, ID = ID,
                 FunctionClass=FunctionClass, Young=λ, ν=μ, cMat=cMat,dcMatdλ=dcdλ)

K, dKdλ = assemble_system(model, GRAD=true)  

# get ∂K/∂θ with finite (central) difference
ΔcMatp = get_cMat(mode, (λ+Δλ), μ)
model.cMat = ΔcMatp

ΔKp = assemble_system(model)

ΔcMatm = get_cMat(mode, (λ-Δλ), μ)
model.cMat = ΔcMatm

ΔKm = assemble_system(model)

dKdλ_approx = (ΔKp-ΔKm)/(2*Δλ)

iIter,jIter = size(ΔKp)
for i in 1:iIter  
    for j in 1:jIter
        @test dKdλ_approx[i,j] ≈ dKdλ[i,j] atol=10^(-5)
    end
end

p, dp, model = simulate_single_tstep(x0, x1, y0, y1, z0, z1, ne, λ, μ, ndim, FunctionClass, nDof, β, μ_tp, μ_btm, mode="lame", GRAD=true)
Nodes = model.NodeList + p

# estimate dp with finite (central) difference
Δp_Lp, model = simulate_single_tstep(x0, x1, y0, y1, z0, z1, ne, (λ+Δλ), μ, ndim, FunctionClass, nDof, β, μ_tp, μ_btm, mode="lame")
ΔNodesLp = model.NodeList + Δp_Lp

Δp_Lm, model = simulate_single_tstep(x0, x1, y0, y1, z0, z1, ne, (λ-Δλ), μ, ndim, FunctionClass, nDof, β, μ_tp, μ_btm, mode="lame")
ΔNodesLm = model.NodeList + Δp_Lm

# estimate dp with finite (central) difference
Δp_βp, model = simulate_single_tstep(x0, x1, y0, y1, z0, z1, ne, λ, μ, ndim, FunctionClass, nDof, (β+Δβ), μ_tp, μ_btm, mode="lame")
ΔNodesβp = model.NodeList + Δp_βp

Δp_βm, model = simulate_single_tstep(x0, x1, y0, y1, z0, z1, ne, λ, μ, ndim, FunctionClass, nDof, (β-Δβ), μ_tp, μ_btm, mode="lame")
ΔNodesβm = model.NodeList + Δp_βm

dLp_approx = (ΔNodesLp-ΔNodesLm)/(2*Δλ)
dβp_approx = (ΔNodesβp-ΔNodesβm)/(2*Δβ)

# display(dp[:,:,1])
# display(dLp_approx)
iIter,jIter = size(dLp_approx)        

for i in 1:iIter  
    for j in 1:jIter
        @test dp[i,j,1] ≈ dLp_approx[i,j] atol=10^(-5)
    end
end

for i in 1:iIter  
    for j in 1:jIter
        @test dp[i,j,2] ≈ dβp_approx[i,j] atol=10^(-5)
    end
end

## Testing ∇u(θ)
BorderPts2D, dudλ, SurfacePts2D, ∇SurfacePts2D = extract_borders(Nodes, CameraMatrix, BorderNodesList, GRAD=true, dqdθ=dp, SIDES=false)

# estimate dudλ with finite (central) difference
ΔBorderPts2DLp, ΔSurfacePts2DLp = extract_borders(ΔNodesLp, CameraMatrix, BorderNodesList, GRAD=false)
ΔBorderPts2DLm, ΔSurfacePts2DLm = extract_borders(ΔNodesLm, CameraMatrix, BorderNodesList, GRAD=false)

ΔBorderPts2Dβp, ΔSurfacePts2Dβp = extract_borders(ΔNodesβp, CameraMatrix, BorderNodesList, GRAD=false)
ΔBorderPts2Dβm, ΔSurfacePts2Dβm = extract_borders(ΔNodesβm, CameraMatrix, BorderNodesList, GRAD=false)

∇SurfacePts2DL_approx = (ΔSurfacePts2DLp - ΔSurfacePts2DLm)/(2*Δλ)
∇SurfacePts2Dβ_approx = (ΔSurfacePts2Dβp - ΔSurfacePts2Dβm)/(2*Δβ)

iIter,jIter = size(∇SurfacePts2DL_approx)
for i in 1:iIter  
    for j in 1:jIter
        @test ∇SurfacePts2D[i,j,1] ≈ ∇SurfacePts2DL_approx[i,j] atol=10^(-4)
    end
end

for i in 1:iIter  
    for j in 1:jIter
        @test ∇SurfacePts2D[i,j,2] ≈ ∇SurfacePts2Dβ_approx[i,j] atol=10^(-4)
    end
end

dudλ_approx = (ΔBorderPts2DLp - ΔBorderPts2DLm)/(2*Δλ)

iIter,jIter = size(dudλ_approx)
for i in 1:iIter  
    for j in 1:jIter
        @test dudλ[i,j,1] ≈ dudλ_approx[i,j] atol=10^(-4)
    end
end

p_gt, model = simulate_single_tstep(x0, x1, y0, y1, z0, z1, ne, (λ*1.3), μ, ndim, FunctionClass, nDof, β, μ_tp, μ_btm, mode="lame")
Nodes_gt = model.NodeList + p_gt

BorderPts2D_gt = back_project(Nodes_gt, CameraMatrix)

costList, ∇d, ∇2d, pairsList = closest_point([BorderPts2D],[BorderPts2D_gt],[dudλ])
ΔcostList_Lp, pairsList_Lp = closest_point([ΔBorderPts2DLp],[BorderPts2D_gt])
ΔcostList_Lm, pairsList_Lm = closest_point([ΔBorderPts2DLm],[BorderPts2D_gt])

ΔcostList_βp, pairsList_βp = closest_point([ΔBorderPts2Dβp],[BorderPts2D_gt])
ΔcostList_βm, pairsList_βm = closest_point([ΔBorderPts2Dβm],[BorderPts2D_gt])

∇dL_approx = (ΔcostList_Lp - ΔcostList_Lm)/(2*Δλ)
∇dβ_approx = (ΔcostList_βp - ΔcostList_βm)/(2*Δβ)

println(size(∇d[1]))
iIter = size(∇d[1],1)

@test ∇d[1][1] ≈ ∇dL_approx[1] atol=10^(-4)
@test ∇d[1][2] ≈ ∇dβ_approx[1] atol=10^(-4)

## testing Σ∇p(θ)
μ_list, gradList, simBorderPts, splinex, spliney, mdl, SurfacePt2D = test(x0, x1, y0, y1, z0, z1, ne, λ, μ, ndim, FunctionClass, nDof, writeData=true, filepath = file,
                                                                                    β, CameraMatrix, endTime, tSteps, Control)

μ_list, simBorderPts, ΔsimBorderPts_pL, splinex, spliney, mdl, ΔSurfacePt2D_pL = test(x0, x1, y0, y1, z0, z1, ne, (λ+Δλ), μ, ndim, FunctionClass, nDof, 
                                                                                    β, CameraMatrix, endTime, tSteps, Control)

μ_list, simBorderPts, ΔsimBorderPts_mL, splinex, spliney, mdl, ΔSurfacePt2D_mL = test(x0, x1, y0, y1, z0, z1, ne, (λ-Δλ), μ, ndim, FunctionClass, nDof, 
                                                                                    β, CameraMatrix, endTime, tSteps, Control)

μ_list, simBorderPts, ΔsimBorderPts_pβ, splinex, spliney, mdl, ΔSurfacePt2D_mL = test(x0, x1, y0, y1, z0, z1, ne, λ, μ, ndim, FunctionClass, nDof, 
                                                                                    (β+Δβ), CameraMatrix, endTime, tSteps, Control)

μ_list, simBorderPts, ΔsimBorderPts_mβ, splinex, spliney, mdl, ΔSurfacePt2D_mL = test(x0, x1, y0, y1, z0, z1, ne, λ, μ, ndim, FunctionClass, nDof, 
                                                                                (β-Δβ), CameraMatrix, endTime, tSteps, Control)
                                                                                    
global tdq = zeros(size(SurfacePt2D[1]))
global tdq_approx = zeros(size(SurfacePt2D[1]))

titer = 1:length(ΔSurfacePt2D_mL)
for t in titer
    grad = gradList[t]
    grad_approx_L = (ΔsimBorderPts_pL[t]- ΔsimBorderPts_mL[t])/(2*Δλ)
    grad_approx_β = (ΔsimBorderPts_pβ[t]- ΔsimBorderPts_mβ[t])/(2*Δβ)
    pIter,qIter = size(grad_approx_L)
    for i in 1:pIter  
        for j in 1:qIter
            @test grad[i,j,1] ≈ grad_approx_L[i,j] atol=10^(-1)
            @test grad[i,j,2] ≈ grad_approx_β[i,j] atol=10^(-1)
        end
    end
end

