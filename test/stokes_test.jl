# test 3D FEM solution for 2x2x2 elements
using smearFEM
using Test
using Plots
using LinearAlgebra
@testset "testing the stokes model" begin

    # test case 
    r::Float64 = 25.0  # radius of the cylinder in mm
    h::Float64 = 40.0  # height of the cylinder in mm
    ne = 4
    ndim = 3
    FunctionClass_x = "Q2"
    FunctionClass_u = "Q2"
    FunctionClass_p = "Q1"
    nDof_u = ndim  # number of degree of freedom per node
    nDof_p = 1
    β = 1e-5
    ν = 40.0
    μu_tp = -1.0
    μu_btm = 0
    μu_side = 0
    acc = 1e-3

    q, model = simulate_single_tstep_stokes(r, h, ne, ν, ndim, FunctionClass_u, FunctionClass_p, nDof_u, nDof_p, β, μu_tp, μu_btm, 
                                            μu_side, FunctionClass_x=FunctionClass_x)

    NodeList = model.mesh_x.NodeList
    iter = 1:size(NodeList, 2)

    r_list = zeros(iter)
    h_list = zeros(iter)
    for i in iter
        r_ = sqrt(NodeList[1,i]^2 + NodeList[2,i]^2)
        h_ = NodeList[3,i]

        r_list[i] = r_
        h_list[i] = h_

        @test sqrt(q[1,i]^2 + q[2,i]^2) + 0.5*μu_tp*r/h ≈ 0 atol=acc
        @test q[3,i] - μu_tp*h_/h ≈ 0 atol=acc
    end
end