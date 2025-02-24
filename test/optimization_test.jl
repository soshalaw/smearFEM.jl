using smearFEM
using Test

# testing the gradient of the projection function
a = [3.0 2.0 5.0]'
Δ = 0.00001

# transform point cloud wrt to camera frame 
R = [1 0 0; 0 0 1; 0 -1 0]     # rotation matrix
t = [0; -0.5; 4]               # translation vector

aTrans = R*a .+ t

Δx = [(aTrans[1]+Δ)  aTrans[2]     aTrans[3]]'
Δy = [aTrans[1]     (aTrans[2]+Δ)  aTrans[3]]'
Δz = [aTrans[1]      aTrans[2]    (aTrans[3]+Δ)]'

∇u = CameraMatrix*[1/aTrans[3] 0 -aTrans[1]/aTrans[3]^2; 0 1/aTrans[3] -aTrans[2]/aTrans[3]^2; 0 0 0]

u = project_to(Matrix(aTrans), CameraMatrix)

Δu_x = project_to(Matrix(Δx), CameraMatrix)
Δu_y = project_to(Matrix(Δy), CameraMatrix)
Δu_z = project_to(Matrix(Δz), CameraMatrix)

Δu = [Δu_x Δu_y Δu_z]
∇u_ = (Δu-[u u u])/Δ

pIter,qIter = size(∇u)
for p in 1:pIter  
    for q in 1:qIter
        @test ∇u_[p,q] ≈ ∇u[p,q] atol=10^(-1)
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
CameraMatrix = [[8*2048/7.07, 0.0, 2048/2] [0.0, 8*1536/5.3, 1536/2] [0.0, 0.0, 1.0]]'
endTime = 30
tSteps = 2
Control = "displacement"

μ_tp = -0.1
μ_btm = 0

Young = 30
ν = 0.3
β = 100
λ = Young*ν/((1+ν)*(1-2*ν))
μ = Young/(2*(1+ν))
Δλ = 0.00005

NodeList, IEN, ID, IEN_top, IEN_btm, BorderNodesList = meshgrid_cube(x0,x1,y0,y1,z0,z1,ne,ndim,FunctionClass=FunctionClass)  # generate the mesh grid
NodeListCylinder = inflate_cylinder(NodeList, x0, x1, y0, y1)   

mode = "lame"    
cMat = get_cMat(mode, λ, μ) # E, ν or λ, μ
dcdλ = get_cMat(mode, 1.0 , 0.0)

model = def_model("linear_elasticity", ne=ne, NodeList=NodeListCylinder, IEN=IEN, IEN_top=IEN_top, IEN_btm=IEN_btm, ndim=ndim, nDof=nDof, ID = ID,
                 FunctionClass=FunctionClass, Young=Float64(Young), ν=ν, cMat=cMat,dcMatdλ=dcdλ)

K, dKdλ = assemble_system(model, GRAD=true)  

# get ∂K/∂θ using finite differencing
ΔcMat = get_cMat(mode, (λ+Δλ), μ)
model.cMat = ΔcMat

ΔK = assemble_system(model)

dKdλ_ = (ΔK-K)/Δλ

iIter,jIter = size(ΔK)
for i in 1:iIter  
    for j in 1:jIter
        @test dKdλ_[i,j] ≈ dKdλ[i,j] atol=10^(-8)
    end
end

q, dq, model = simulate_single_tstep(x0, x1, y0, y1, z0, z1, ne, λ, μ, ndim, FunctionClass, nDof, β, μ_tp, μ_btm, mode="lame", GRAD=true)
Nodes = model.NodeList + q

Δq, model = simulate_single_tstep(x0, x1, y0, y1, z0, z1, ne, (λ+Δλ), μ, ndim, FunctionClass, nDof, β, μ_tp, μ_btm, mode="lame")
ΔNodes = model.NodeList + Δq

dq_ = (Δq - q)/Δλ

iIter,jIter = size(dq_)           
for i in 1:iIter  
    for j in 1:jIter
        @test dq[i,j] ≈ dq_[i,j] atol=10^(-8)
    end
end

BorderPts2D, dudλ = back_project(Nodes, CameraMatrix, dq)
ΔBorderPts2D = back_project(ΔNodes, CameraMatrix)

dudλ_ = (ΔBorderPts2D - BorderPts2D)/Δλ

println("estimated")
display(dudλ_)
println("analytical")
display(dudλ)

iIter,jIter = size(dudλ_)
for i in 1:iIter  
    for j in 1:jIter
        @test dudλ[i,j] ≈ dudλ_[i,j] atol=10^(-8)
    end
end
