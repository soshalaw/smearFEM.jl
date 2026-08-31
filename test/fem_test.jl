using Test
using smearFEM

@testset "testing basis functions" begin
    lx = 2
    ly = 2
    lz = 2
    ne = 1

    # Structured mesh builders produce Line/Quad/Hex topology only.
    # :Tet is for unstructured (Gmsh) meshes and is not tested here.
    valid_shapes = Dict(1 => :Line, 2 => :Quad, 3 => :Hex)
    basis_orders = [1, 2]

    ndims = [1, 2, 3]

    for ndim in ndims
        element_shape = valid_shapes[ndim]
        for basis_order in basis_orders
            if ndim == 1
                mesh = meshgrid_line(lx; ne=ne, element_shape=element_shape, basis_order=basis_order)
                iter = 1:size(mesh.IEN, 1)
                for i in iter
                    coord = mesh.NodeList[:, mesh.IEN[i]]
                    N, dN = basis_function(coord[1], nothing, nothing, element_shape, basis_order)
                    @test findall(x->x==1, N) == [i]
                end
            elseif ndim == 2
                mesh = meshgrid_square(lx, ly; ne=ne, element_shape=element_shape, basis_order=basis_order)
                iter = 1:size(mesh.IEN, 2)
                for i in iter
                    coord = mesh.NodeList[:, mesh.IEN[i]]
                    N, dN = basis_function(coord[1], coord[2], nothing, element_shape, basis_order)
                    @test findall(x->x==1, N) == [i]
                end
            elseif ndim == 3
                mesh = meshgrid_cuboid(lx, ly, lz; ne=ne, element_shape=element_shape, basis_order=basis_order)
                iter = 1:size(mesh.IEN, 2)
                for i in iter
                    coord = mesh.NodeList[:, mesh.IEN[i]]
                    ζ = (2 * coord[3] / lz) - 1
                    N, dN = basis_function(coord[1], coord[2], ζ, element_shape, basis_order)
                    @test findall(x->x==1, N) == [i]
                end
            end
        end
    end
end

@testset "testing basis function second derivatives" begin
    N, dN, d2N = basis_function(0.0, nothing, nothing, :Line, 2; second_derivatives=true)
    @test d2N == [1.0, 1.0, -2.0]

    N, dN, d2N = basis_function(0.0, 0.0, nothing, :Quad, 1; second_derivatives=true)
    @test size(d2N) == (4, 3)
    @test d2N[:, 1] == [0.0, 0.0, 0.0, 0.0]
    @test d2N[:, 2] == [0.25, -0.25, 0.25, -0.25]
end

# test gaussian quadrature
@test gaussian_quadrature(-1,1,2) == ([-1/√3, 1/√3], [1.0, 1.0])
@test gaussian_quadrature(-1,1,3) == ([-√(3/5), 0.0, √(3/5)], [5/9, 8/9, 5/9])

# test 2D FEM solution for 2x2 elements

# test 2D FEM solution for annulus ring with 2x2 elements

