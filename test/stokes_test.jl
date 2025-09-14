# test 3D FEM solution for 2x2x2 elements
using smearFEM
using Test
using Plots
using LinearAlgebra
# @testset "testing the stokes model" begin

    # test case 
    scale = 50
    r::Float64 = 0.5*scale  # radius of the cylinder in mm
    h::Float64 = 1*scale  # height of the cylinder in mm
    ne = 4
    ndim = 3
    FunctionClass_x = "S2"
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
    
    # u_r = zeros(size(q,2))
    # u_z = zeros(size(q,2))

    vis_mesh = meshgrid_cylinder(r, h, ne, FunctionClass="Q3")

    T = Matrix{Float64}(I, size(model.mesh_x.NodeList,2), size(vis_mesh.NodeList,2))

    println(vis_mesh.IEN[:,1]')

    if model.mesh_x.FunctionClass == "S2" && model.mesh_u.FunctionClass != "S2"
        T_vis = get_nurbs_2_lagrange_proj(model.mesh_x.IEN, vis_mesh.IEN, model.mesh_x.C_vol_vis, model.mesh_x.NodeList, model.mesh_x.W)
        T = get_nurbs_2_lagrange_proj(model.mesh_x.IEN, model.mesh_u.IEN, model.mesh_x.C_vol, model.mesh_x.NodeList, model.mesh_x.W)
    end

    display(T)
    T_ = T'*inv(T*T')

    display(q)
    q_ = q*T_*T
    display(q*T_)
    display(q*T_*T)
    NodeList_ = model.mesh_x.NodeList*T
    IEN_ = model.mesh_u.IEN

    NodeList_vis = model.mesh_x.NodeList*T_vis
    IEN_vis = vis_mesh.IEN

    iter = 1:size(NodeList_, 2)
    println("Control pt size: ", size(model.mesh_x.NodeList))
    println("Lagrange point mesh size: ", size(model.mesh_u.IEN))
    println("solution ",size(q))
    println("solution ",size(q_))

    r_list = zeros(iter)
    h_list = zeros(iter)
    for i in iter
        r_ = sqrt(NodeList_vis[1,i]^2 + NodeList_vis[2,i]^2)
        h_ = NodeList_vis[3,i]

        r_list[i] = r_
        h_list[i] = h_

        # @test sqrt(q_[1,i]^2 + q_[2,i]^2) + 0.5*μu_tp*r/h ≈ 0 atol=acc
        # @test q_[3,i] - μu_tp*h_/h ≈ 0 atol=acc
    end

    # @test maximum(r_list) ≈ r atol=10^(-6)
    # @test maximum(h_list) ≈ h atol=10^(-6)

    animate_fields(filepath = "/home/soshala/SMEAR-PhD/temp/", Nodes=[NodeList_] , IEN=IEN_)
    animate_fields(filepath = "/home/soshala/SMEAR-PhD/temp/cp", Nodes=[model.mesh_x.NodeList] , IEN=model.mesh_x.IEN_cp)

    set_plot(22)
    Plots.scatter3d!(NodeList_[1,:], NodeList_[2,:], NodeList_[3,:], markersize=5, label="Lagrange points")
    Plots.savefig("/home/soshala/SMEAR-PhD/temp/stokes_lagrange_pts.png")


    animate_fields(filepath = "/home/soshala/SMEAR-PhD/temp/vis", Nodes=[NodeList_vis] , IEN=IEN_vis)

    set_plot(22)
    Plots.scatter3d!(NodeList_vis[1,:], NodeList_vis[2,:], NodeList_vis[3,:], markersize=5, label="Lagrange points")
    Plots.savefig("/home/soshala/SMEAR-PhD/temp/vis_stokes_lagrange_pts.png")
    
    write_vtk("/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/sim_experiments/unit_tests/1", "u", NodeList_, IEN_, ne, ndim, q_, FunctionClass=FunctionClass_u)

    lp = model.mesh_x.NodeList * T
    n_pts = size(lp, 2)

    q_ = zeros(size(q))
    NodeList_ = zeros(size(model.mesh_x.NodeList))
    IEN_ = zeros(size(model.mesh_x.IEN))
    write_vtk("/home/soshala/SMEAR-PhD", "u", model.mesh_u.NodeList, model.mesh_u.IEN, ne, ndim, q, FunctionClass=FunctionClass_u)
    lp = get_lagrange_pts(model.mesh_x.IEN, model.mesh_u.IEN, model.mesh_x.C_vol, model.mesh_x.NodeList, model.mesh_x.W)


    
    NodeList_, IEN_, q_ = eval_on_cylinder(model, 1, q)
# end