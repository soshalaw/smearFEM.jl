# test 3D FEM solution for 2x2x2 elements
using smearFEM
using Test

@testset "testing linear elasticity" begin

    # test case 
    x0 = 0
    x1 = 1
    y0 = 0
    y1 = 1
    z0 = 0
    z1 = 1
    H = z1 - z0
    ne = 8
    ndim = 3
    FunctionClass = "Q1"
    nDof = ndim  # number of degree of freedom per node
    β = 1.0e-5
    Young = 40
    ν = 0.4
    μ_tp = -0.1
    μ_btm = 0

    q, model = simulate_single_tstep(x0, x1, y0, y1, z0, z1, ne, Young, ν, ndim, FunctionClass, nDof, β, μ_tp, μ_btm)

    u_r = zeros(size(q,2))
    u_z = zeros(size(q,2))

    iter = 1:size(model.NodeList, 2)
    for i in iter
        r = sqrt(model.NodeList[1,i]^2 + model.NodeList[2,i]^2)
        h = model.NodeList[3,i]

        @test sqrt(q[1,i]^2 + q[2,i]^2) + ν*μ_tp*r/H ≈ 0 atol=10^(-5)
        @test q[3,i] - μ_tp*h/H ≈ 0 atol=10^(-5)
    end

end