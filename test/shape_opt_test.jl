using smearFEM
using Test
using Plots

function test_cylinder()

    scale = 50
    ne::Int = 4 # number of elements in the mesh for the ground truth
    ndim::Int = 3
    FunctionClass::String = "Q2"
    camera_matrix = [[8*2048/7.07, 0.0, 2048/2] [0.0, 8*1536/5.3, 1536/2] [0.0, 0.0, 1.0]]'
    camera_pose = scale*[0 -0.5 2.75]'   # camera position in mm
    NodeList, IEN, ID, IEN_top, IEN_bottom, IEN_side, nNodes, BorderNodes = meshgrid_cube(1, 1, 1, ne, FunctionClass=FunctionClass)
    side_node_list = BorderNodes[1]
    r = 0.5*scale
    h = 1*scale
    NodeListCyl, ∇NodeListCyl = inflate_cylinder(NodeList, -0.5, 0.5, -0.5, 0.5, r, h, GRAD=true)

    display(NodeList)
    Δr = 0.00001
    ΔNodeListCyl_rp = inflate_cylinder(NodeList, -0.5, 0.5, -0.5, 0.5, (r+Δr), h)
    ΔNodeListCyl_rm = inflate_cylinder(NodeList, -0.5, 0.5, -0.5, 0.5, (r-Δr), h)
    ∇NodeListCyl_r_approx = (ΔNodeListCyl_rp - ΔNodeListCyl_rm)/(2*Δr)

    pIter,qIter = size(∇NodeListCyl_r_approx)
    for p in 1:pIter  
        for q in 1:qIter
            @test ∇NodeListCyl_r_approx[p,q] ≈ ∇NodeListCyl[p,q,1] atol=10^(-4)
        end
    end

    Δh = 0.00001
    ΔNodeListCyl_hp = inflate_cylinder(NodeList, -0.5, 0.5, -0.5, 0.5, r, (h+Δh))
    ΔNodeListCyl_hm = inflate_cylinder(NodeList, -0.5, 0.5, -0.5, 0.5, r, (h-Δh))
    ∇NodeListCyl_h_approx = (ΔNodeListCyl_hp - ΔNodeListCyl_hm)/(2*Δh)

    pIteh,qIteh = size(∇NodeListCyl_h_approx)
    for p in 1:pIteh  
        for q in 1:qIteh
            @test ∇NodeListCyl_h_approx[p,q] ≈ ∇NodeListCyl[p,q,2] atol=10^(-4)
        end
    end

    BorderPts2D, ∇BorderPts2D = extract_borders(NodeListCyl, camera_matrix, camera_pose, side_node_list, nNodes, GRAD=true, dqdθ=∇NodeListCyl, SIDES=false)

    ΔBorderPts2D_rp, ΔSurfacePts2D_rp = extract_borders(ΔNodeListCyl_rp, camera_matrix, camera_pose, side_node_list, nNodes)
    ΔBorderPts2D_rm, ΔSurfacePts2D_rm = extract_borders(ΔNodeListCyl_rm, camera_matrix, camera_pose, side_node_list, nNodes)
    ∇BorderPts2D_r_approx = (ΔBorderPts2D_rp - ΔBorderPts2D_rm)/(2*Δr)

    pIter,qIter = size(∇BorderPts2D_r_approx)
    for p in 1:pIter  
        for q in 1:qIter
            @test ∇BorderPts2D_r_approx[p,q] ≈ ∇BorderPts2D[p,q,1] atol=10^(-5)
        end
    end

    ΔBorderPts2D_hp, ΔSurfacePts2D_hp = extract_borders(ΔNodeListCyl_hp, camera_matrix, camera_pose, side_node_list, nNodes)
    ΔBorderPts2D_hm, ΔSurfacePts2D_hm = extract_borders(ΔNodeListCyl_hm, camera_matrix, camera_pose, side_node_list, nNodes)
    ∇BorderPts2D_h_approx = (ΔBorderPts2D_hp - ΔBorderPts2D_hm)/(2*Δh)

    pIter,qIter = size(∇BorderPts2D_h_approx)
    for p in 1:pIter  
        for q in 1:qIter
            @test ∇BorderPts2D_h_approx[p,q] ≈ ∇BorderPts2D[p,q,2] atol=10^(-5)
        end
    end
end

function test_cyl_dif()

    ne::Int = 4 # number of elements in the mesh for the ground truth
    FunctionClass::String = "Q2"
    camera_matrix = [[8*2048/7.07, 0.0, 2048/2] [0.0, 8*1536/5.3, 1536/2] [0.0, 0.0, 1.0]]'
    camera_pose = scale*[0 -0.5 2.75]'   # camera position in mm
    len = 101
    cost_list = zeros(Float64,len,len)

    #g_truth
    scale = 50
    r::Float64 = 0.5*scale  # radius of the cylinder in mm
    h::Float64 = 1*scale  # height of the cylinder in mm
    NodeList_gt, IEN, ID, IEN_top, IEN_bottom, IEN_side, nNodes, BorderNodes = meshgrid_cube(1, 1, h, ne, FunctionClass=FunctionClass)
    NodeListCyl_gt = inflate_cylinder(NodeList_gt, -0.5, 0.5, -0.5, 0.5, r)
    side_nodes = BorderNodes[1]

    BorderPts2D_gt, g_SurfacePts2D = extract_borders(NodeListCyl_gt, camera_matrix, camera_pose, side_nodes, nNodes, SIDES=false)

    test_r_list = collect(range(start=(r-2),stop=(r+2),length=len))
    test_h_list = collect(range(start=(h-2),stop=(h+2),length=len))

    h_iter = 1
    for h_t in test_h_list
        r_iter = 1
        for r_t in test_r_list
            NodeList, IEN, ID, IEN_top, IEN_bottom, IEN_side, nNodes, BorderNodes = meshgrid_cube(1, 1, h_t, ne, FunctionClass=FunctionClass)
            NodeListCyl = inflate_cylinder(NodeList, -0.5, 0.5, -0.5, 0.5, r_t)
            BorderPts2D, SurfacePts2D = extract_borders(NodeListCyl, camera_matrix, camera_pose, side_nodes, nNodes, SIDES=false)

            costList, pairsList = closest_point([BorderPts2D],[BorderPts2D_gt])
            cost_list[r_iter,h_iter] = sum(costList)
            r_iter = r_iter + 1
        end
        h_iter = h_iter + 1
    end
    Plots.contourf(test_h_list,test_r_list,cost_list)
end

# test_cyl_dif()
test_cylinder()