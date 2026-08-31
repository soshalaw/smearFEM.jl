# test 3D FEM solution for 2x2x2 elements
using smearFEM
using Test
using Plots
using LinearAlgebra
@testset "testing the stokes model" begin

    # test case 
    r::Float64 = 25.0  # radius of the cylinder in mm
    h::Float64 = 40.0  # height of the cylinder in mm
    ne = 6
    ndim = 3
    element_shape_x::Symbol = :Hex
    basis_order_x::Int = 2  
    element_shape_u::Symbol = :Hex
    basis_order_u::Int = 2
    element_shape_p::Symbol = :Hex
    basis_order_p::Int = 1
    nDof_u = ndim  # number of degree of freedom per node
    nDof_p = 1
    β = 1e-5
    ν = 70.0
    μu_tp = -1.0
    μu_btm = 0
    μu_side = 0
    acc = 1e-3

    q, model = simulate_single_tstep_stokes(r, h, ne, ν, ndim, element_shape_u, basis_order_u, element_shape_p, basis_order_p, nDof_u, nDof_p, β, μu_tp, μu_btm,
                                            μu_side, element_shape_x=element_shape_x, basis_order_x=basis_order_x)
    
    # Get mesh data directly (no NURBS projection needed for Q2 basis)
    NodeList = model.mesh_x.NodeList
    IEN = model.mesh_u.IEN

    iter = 1:size(NodeList, 2)
    for i in iter
        r_ = sqrt(NodeList[1,i]^2 + NodeList[2,i]^2)
        h_ = NodeList[3,i]

        @test sqrt(q[1,i]^2 + q[2,i]^2) + 0.5*μu_tp*r_/h ≈ 0 atol=acc
        @test q[3,i] - μu_tp*h_/h ≈ 0 atol=acc
    end
end