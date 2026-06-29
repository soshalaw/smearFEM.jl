using Test
using smearFEM

# ── nodes-per-element for each (element_shape, basis_order) combo ────────────
function _npe(es::Symbol, bo::Int)
    return Dict(
        (:Hex,  1) => 8,  (:Hex,  2) => 27,
        (:Tet,  1) => 4,  (:Tet,  2) => 10,
        (:Quad, 1) => 4,  (:Quad, 2) => 9,
        (:Tri,  1) => 3,  (:Tri,  2) => 6,
        (:Line, 1) => 2,  (:Line, 2) => 3,
    )[(es, bo)]
end

# ── common checks that apply to every Gmsh mesh ──────────────────────────────
function check_mesh_basics(mesh, es::Symbol, bo::Int)
    npe = _npe(es, bo)
    @test mesh.ne > 0                        "ne must be positive"
    @test mesh.nNodes > 0                    "nNodes must be positive"
    @test size(mesh.IEN, 1) == npe           "IEN rows must equal nodes-per-element"
    @test size(mesh.IEN, 2) == mesh.ne       "IEN cols must equal ne"
    @test size(mesh.NodeList, 2) == mesh.nNodes  "NodeList cols must equal nNodes"
    # Gmsh always returns 3-row NodeList (x, y, z) even for 2-D/1-D meshes
    @test size(mesh.NodeList, 1) == 3        "NodeList must be 3×nNodes"
    @test all(1 .<= mesh.IEN .<= mesh.nNodes)  "IEN must reference valid node indices"
end

@testset "Gmsh mesh templates" begin

    mesh_dir = mktempdir()  # fresh dir for each test run; no cache contamination

    # ── 3-D cylinder ─────────────────────────────────────────────────────────
    @testset "cylinder_hex Q1" begin
        r, h = 25.0, 40.0
        m = meshgrid_cylinder(r, h;
            mesh_type=:unstructured, element_shape=:Hex, basis_order=1,
            elem_size=3.0, mesh_path=mesh_dir, ndof=1)
        check_mesh_basics(m, :Hex, 1)
        @test all(sqrt.(m.NodeList[1,:].^2 .+ m.NodeList[2,:].^2) .<= r + 1e-6)
        @test all(m.NodeList[3,:] .>= -1e-6)
        @test all(m.NodeList[3,:] .<= h + 1e-6)
        @test !isempty(m.top_nodes)
        @test !isempty(m.bottom_nodes)
        @test !isempty(m.side_nodes)

        # finer mesh must have strictly more elements
        m_fine = meshgrid_cylinder(r, h;
            mesh_type=:unstructured, element_shape=:Hex, basis_order=1,
            elem_size=5.0, mesh_path=mesh_dir, ndof=1)
        @test m_fine.ne > m.ne
    end

    @testset "cylinder_hex Q2" begin
        r, h = 25.0, 40.0
        m = meshgrid_cylinder(r, h;
            mesh_type=:unstructured, element_shape=:Hex, basis_order=2,
            elem_size=3.0, mesh_path=mesh_dir, ndof=1)
        check_mesh_basics(m, :Hex, 2)
    end

    @testset "cylinder_tet Q1" begin
        r, h = 25.0, 40.0
        m = meshgrid_cylinder(r, h;
            mesh_type=:unstructured, element_shape=:Tet, basis_order=1,
            elem_size=3.0, mesh_path=mesh_dir, ndof=1)
        check_mesh_basics(m, :Tet, 1)
        @test all(sqrt.(m.NodeList[1,:].^2 .+ m.NodeList[2,:].^2) .<= r + 1e-6)
        @test !isempty(m.top_nodes)
        @test !isempty(m.bottom_nodes)
        @test !isempty(m.side_nodes)

        m_fine = meshgrid_cylinder(r, h;
            mesh_type=:unstructured, element_shape=:Tet, basis_order=1,
            elem_size=5.0, mesh_path=mesh_dir, ndof=1)
        @test m_fine.ne > m.ne
    end

    @testset "cylinder_tet Q2" begin
        r, h = 25.0, 40.0
        m = meshgrid_cylinder(r, h;
            mesh_type=:unstructured, element_shape=:Tet, basis_order=2,
            elem_size=3.0, mesh_path=mesh_dir, ndof=1)
        check_mesh_basics(m, :Tet, 2)
    end

    # ── 3-D cuboid (hex, sharp corners) ──────────────────────────────────────
    @testset "cuboid_hex_sharp Q1" begin
        lx, ly, lz = 10.0, 10.0, 10.0
        m = meshgrid_cuboid(lx, ly, lz;
            mesh_type=:unstructured, element_shape=:Hex, basis_order=1,
            elem_size=3.0, mesh_path=mesh_dir, ndof=1)
        check_mesh_basics(m, :Hex, 1)
        @test all(m.NodeList[1,:] .>= -lx/2 - 1e-6)
        @test all(m.NodeList[1,:] .<= lx/2 + 1e-6)
        @test all(m.NodeList[3,:] .>= -1e-6)
        @test all(m.NodeList[3,:] .<= lz + 1e-6)
        @test !isempty(m.top_nodes)
        @test !isempty(m.bottom_nodes)
        @test !isempty(m.side_nodes)

        m_fine = meshgrid_cuboid(lx, ly, lz;
            mesh_type=:unstructured, element_shape=:Hex, basis_order=1,
            elem_size=5.0, mesh_path=mesh_dir, ndof=1)
        @test m_fine.ne > m.ne
    end

    # ── 3-D cuboid (hex, rounded corners) ────────────────────────────────────
    @testset "cuboid_hex_rounded Q1" begin
        lx, ly, lz = 10.0, 10.0, 10.0
        er = 1.0  # edge radius; must be < min(lx,ly)/2 = 5
        m = meshgrid_cuboid(lx, ly, lz;
            mesh_type=:unstructured, element_shape=:Hex, basis_order=1,
            elem_size=3.0, mesh_path=mesh_dir, ndof=1, edge_radius=er)
        check_mesh_basics(m, :Hex, 1)
        @test !isempty(m.top_nodes)
        @test !isempty(m.bottom_nodes)
        @test !isempty(m.side_nodes)
    end

    # ── 3-D cuboid (tet, sharp corners) ──────────────────────────────────────
    @testset "cuboid_tet_sharp Q1" begin
        lx, ly, lz = 10.0, 10.0, 10.0
        m = meshgrid_cuboid(lx, ly, lz;
            mesh_type=:unstructured, element_shape=:Tet, basis_order=1,
            elem_size=3.0, mesh_path=mesh_dir, ndof=1)
        check_mesh_basics(m, :Tet, 1)
        @test !isempty(m.top_nodes)
        @test !isempty(m.bottom_nodes)
        @test !isempty(m.side_nodes)

        m_fine = meshgrid_cuboid(lx, ly, lz;
            mesh_type=:unstructured, element_shape=:Tet, basis_order=1,
            elem_size=5.0, mesh_path=mesh_dir, ndof=1)
        @test m_fine.ne > m.ne
    end

    # ── 3-D cuboid (tet, rounded corners) ────────────────────────────────────
    @testset "cuboid_tet_rounded Q1" begin
        lx, ly, lz = 10.0, 10.0, 10.0
        er = 1.0
        m = meshgrid_cuboid(lx, ly, lz;
            mesh_type=:unstructured, element_shape=:Tet, basis_order=1,
            elem_size=3.0, mesh_path=mesh_dir, ndof=1, edge_radius=er)
        check_mesh_basics(m, :Tet, 1)
        @test !isempty(m.top_nodes)
        @test !isempty(m.bottom_nodes)
        @test !isempty(m.side_nodes)
    end

    # ── 2-D disk ──────────────────────────────────────────────────────────────
    @testset "disk_quad Q1" begin
        r = 25.0
        m = meshgrid_disk(r;
            element_shape=:Quad, basis_order=1,
            elem_size=8.0, mesh_path=mesh_dir, ndof=2)
        check_mesh_basics(m, :Quad, 1)
        @test all(sqrt.(m.NodeList[1,:].^2 .+ m.NodeList[2,:].^2) .<= r + 1e-6)
        @test !isempty(m.boundary_nodes)

        m_fine = meshgrid_disk(r;
            element_shape=:Quad, basis_order=1,
            elem_size=5.0, mesh_path=mesh_dir, ndof=2)
        @test m_fine.ne > m.ne
    end

    @testset "disk_tri Q1" begin
        r = 25.0
        m = meshgrid_disk(r;
            element_shape=:Tri, basis_order=1,
            elem_size=8.0, mesh_path=mesh_dir, ndof=2)
        check_mesh_basics(m, :Tri, 1)
        @test !isempty(m.boundary_nodes)
    end

    # ── 2-D square ────────────────────────────────────────────────────────────
    @testset "square_quad Q1" begin
        lx, ly = 50.0, 50.0
        m = meshgrid_square(lx, ly;
            mesh_type=:unstructured, element_shape=:Quad, basis_order=1,
            elem_size=10.0, mesh_path=mesh_dir, ndof=2)
        check_mesh_basics(m, :Quad, 1)
        # square_quad.geo is centered at origin: nodes in [-lx/2, lx/2] × [-ly/2, ly/2]
        @test all(m.NodeList[1,:] .>= -lx/2 - 1e-6)
        @test all(m.NodeList[1,:] .<= lx/2 + 1e-6)
        @test all(m.NodeList[2,:] .>= -ly/2 - 1e-6)
        @test all(m.NodeList[2,:] .<= ly/2 + 1e-6)
        @test !isempty(m.top_nodes)
        @test !isempty(m.bottom_nodes)

        m_fine = meshgrid_square(lx, ly;
            mesh_type=:unstructured, element_shape=:Quad, basis_order=1,
            elem_size=6.0, mesh_path=mesh_dir, ndof=2)
        @test m_fine.ne > m.ne
    end

    @testset "square_tri Q1" begin
        lx, ly = 50.0, 50.0
        m = meshgrid_square(lx, ly;
            mesh_type=:unstructured, element_shape=:Tri, basis_order=1,
            elem_size=10.0, mesh_path=mesh_dir, ndof=2)
        check_mesh_basics(m, :Tri, 1)
        @test !isempty(m.top_nodes)
        @test !isempty(m.bottom_nodes)
    end

    # ── 1-D line ──────────────────────────────────────────────────────────────
    @testset "line Q1" begin
        l = 40.0
        m = meshgrid_line(l;
            mesh_type=:unstructured, element_shape=:Line, basis_order=1,
            elem_size=8.0, mesh_path=mesh_dir, ndof=1)
        check_mesh_basics(m, :Line, 1)
        # line.geo centers on origin: nodes in [-l/2, l/2]
        @test all(abs.(m.NodeList[1,:]) .<= l/2 + 1e-6)
        @test !isempty(m.boundary_nodes)

        m_fine = meshgrid_line(l;
            mesh_type=:unstructured, element_shape=:Line, basis_order=1,
            elem_size=4.0, mesh_path=mesh_dir, ndof=1)
        @test m_fine.ne > m.ne
    end

    @testset "line Q2" begin
        l = 40.0
        m = meshgrid_line(l;
            mesh_type=:unstructured, element_shape=:Line, basis_order=2,
            elem_size=8.0, mesh_path=mesh_dir, ndof=1)
        check_mesh_basics(m, :Line, 2)
    end

end
