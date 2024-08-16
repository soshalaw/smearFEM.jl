using smearFEM
using Test
using Aqua

@testset "smearFEM.jl" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(smearFEM;
        ambiguities=false,
        deps_compat = false,)
    end
    
    # test mesh generation

    # test basis functions
    @test basis_function(-1)[1] == [1.0, 0.0]
    @test basis_function(1)[1] == [0.0, 1.0]

    @test basis_function(-1,-1)[1] == [1.0, 0.0, 0.0, 0.0]
    @test basis_function(1,-1)[1] == [0.0, 1.0, 0.0, 0.0]
    @test basis_function(1,1)[1] == [0.0, 0.0, 1.0, 0.0]
    @test basis_function(-1,1)[1] == [0.0, 0.0, 0.0, 1.0]

    @test basis_function(-1,-1,-1)[1] == [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    @test basis_function(1,-1,-1)[1] == [0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    @test basis_function(1,1,-1)[1] == [0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    @test basis_function(-1,1,-1)[1] == [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0]
    @test basis_function(-1,-1,1)[1] == [0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0]
    @test basis_function(1,-1,1)[1] == [0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0]
    @test basis_function(1,1,1)[1] == [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0]
    @test basis_function(-1,1,1)[1] == [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0]


    # test gaussian quadrature
    @test gaussian_quadrature(-1,1,2) == ([-1/√3, 1/√3], [1.0, 1.0])
    @test gaussian_quadrature(-1,1,3) == ([-√(3/5), 0.0, √(3/5)], [5/9, 8/9, 5/9])

    # test 2D FEM solution for 2x2 elements

    # test 2D FEM solution for annulus ring with 2x2 elements

    # test 3D FEM solution for 2x2x2 elements
    
end
