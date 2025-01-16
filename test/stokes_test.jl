# test 3D FEM solution for 2x2x2 elements
using smearFEM
using Test

@testset "testing the stokes model" begin

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
    FunctionClass_u = "Q2"
    FunctionClass_p = "Q1"
    nDof_u = ndim  # number of degree of freedom per node
    nDof_p = 1
    β = 1.0e-5
    ν = 1
    μu_tp = -0.1
    μu_btm = 0
    μu_side = 0

    q, model = simulate_single_tstep_stokes(x0, x1, y0, y1, z0, z1, ne, ν, ndim, FunctionClass_u, FunctionClass_p, nDof_u, nDof_p, β, μu_tp, μu_btm, 
                                            μu_side; DENSE=false)
    
    u_r = zeros(size(q,2))
    u_z = zeros(size(q,2))

    println("stokes")
    iter = 1:size(model.NodeList_u, 2)
    for i in iter
        r = sqrt(model.NodeList_u[1,i]^2 + model.NodeList_u[2,i]^2)
        h = model.NodeList_u[3,i]

        @test sqrt(q[1,i]^2 + q[2,i]^2) + 0.5*μu_tp*r/H ≈ 0 atol=10^(-5)
        @test q[3,i] - μu_tp*h/H ≈ 0 atol=10^(-5)
    end

end