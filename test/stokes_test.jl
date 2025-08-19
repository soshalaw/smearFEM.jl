# test 3D FEM solution for 2x2x2 elements
using smearFEM
using Test
using Base: dirname

@testset "testing the stokes model" begin

    # test case 
    r = 1   
    h = 1
    ne = 2
    ndim = 3
    FunctionClass_x = "S2"
    FunctionClass_u = "Q2"
    FunctionClass_p = "Q1"
    nDof_u = ndim  # number of degree of freedom per node
    nDof_p = 1
    β = 1.0e-5
    ν = 1.0
    μu_tp = -0.1
    μu_btm = 0
    μu_side = 0

    q, model = simulate_single_tstep_stokes(r, h, ne, ν, ndim, FunctionClass_u, FunctionClass_p, nDof_u, nDof_p, β, μu_tp, μu_btm, 
                                            μu_side, FunctionClass_x=FunctionClass_x)
    
    u_r = zeros(size(q,2))
    u_z = zeros(size(q,2))

    NodeList_, IEN_, q_ = eval_on_cylinder(model, 1, q)

    # write_vtk("/home/soshala/SMEAR-PhD", "u", model.mesh_u.NodeList, model.mesh_u.IEN, ne, ndim, q, FunctionClass=FunctionClass_u)
    write_vtk(dirname(@__DIR__), "u", NodeList_, IEN_, ne, ndim, q_, FunctionClass=FunctionClass_u)


    # for i in iter
    #     r = sqrt(model.mesh_x.NodeList[1,i]^2 + model.mesh_x.NodeList[2,i]^2)
    #     h_ = model.mesh_x.NodeList[3,i]

    #     @test sqrt(q[1,i]^2 + q[2,i]^2) + 0.5*μu_tp*r/h ≈ 0 atol=10^(-5)
    #     @test q[3,i] - k*h_/h ≈ 0 atol=10^(-5)
    # end

end