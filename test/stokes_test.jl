# test 3D FEM solution for 2x2x2 elements
using smearFEM
using Test

@testset "testing the stokes model" begin

    # test case 
    r = 1   
    h = 1
    ne = 2
    ndim = 3
    FunctionClass_u = "S2"
    FunctionClass_p = "Q1"
    nDof_u = ndim  # number of degree of freedom per node
    nDof_p = 1
    β = 1.0e-5
    ν = 1.0
    μu_tp = -0.1
    μu_btm = 0
    μu_side = 0

    q, model = simulate_single_tstep_stokes(r, h, ne, ν, ndim, FunctionClass_u, FunctionClass_p, nDof_u, nDof_p, β, μu_tp, μu_btm, 
                                            μu_side)
    
    u_r = zeros(size(q,2))
    u_z = zeros(size(q,2))

    iter = 1:size(model.mesh_u.NodeList, 2)
    for i in iter
        r = sqrt(model.mesh_u.NodeList[1,i]^2 + model.mesh_u.NodeList[2,i]^2)
        h_ = model.mesh_u.NodeList[3,i]

        @test sqrt(q[1,i]^2 + q[2,i]^2) + 0.5*μu_tp*r/h ≈ 0 atol=10^(-5)
        @test q[3,i] - μu_tp*h_/h ≈ 0 atol=10^(-5)
    end

end