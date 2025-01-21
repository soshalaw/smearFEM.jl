# test 3D FEM solution for 2x2x2 elements
using smearFEM
using Test

filePath = "/home/soshala/SMEAR-PhD/smear-modules/smearFEM.jl/cylindergen"

@testset "testing linear elasticity" begin
    # test case 
    x0 = 0
    x1 = 1
    y0 = 0
    y1 = 1
    z0 = 0
    z1 = 1
    H = z1 - z0
    ne = 2
    ndim = 3
    FunctionClass = "S2" #"Q1"
    nDof = ndim  # number of degree of freedom per node
    β = 1.0e-5
    Young = 40
    ν = 0.4
    μ_tp = -0.1
    μ_btm = 0
    nsub = 1

    q_, model = simulate_single_tstep(Young, ν, FunctionClass, nDof, β, μ_tp, μ_btm)

    NodeList_, IEN_list, q = eval_on_cylinder(model, q_, nsub)
    
    u_r = zeros(size(q,2))
    u_z = zeros(size(q,2))

    iter = 1:size(NodeList_, 2)

    accuzLst = []
    accurLst = []

    for i in iter
        r = sqrt(NodeList_[1,i]^2 + NodeList_[2,i]^2)
        h = NodeList_[3,i]

        if (sqrt(q[1,i]^2 + q[2,i]^2)) > 0.00001
            @test sqrt(q[1,i]^2 + q[2,i]^2) + ν*μ_tp*r/H ≈ 0 atol=10^(-5)
            @test q[3,i] - μ_tp*h/H ≈ 0 atol=10^(-5)
        end
    end
end 