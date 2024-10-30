using smearFEM
using Test
using Aqua

@testset "smearFEM.jl" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(smearFEM;
        ambiguities=false,
        deps_compat = false,)
    end

    @testset "testing basis functions" begin
        x0 = -1
        x1 = 1
        y0 = -1
        y1 = 1
        z0 = -1
        z1 = 1
        ne = 1
    
        FunctionsClasses = ["Q1", "Q2"]
        ndims = [1, 2, 3]

        for FunctionClass in FunctionsClasses
            for ndim in ndims
                if ndim == 1
                    NodeList, IEN, BorderNodesList = meshgrid_line(x0,x1,ne;FunctionClass=FunctionClass)
                    iter = 1:size(IEN,2)
                    for i in iter
                        coord = NodeList[:,IEN[i]]
                        N, dN = basis_function(coord[1],nothing,nothing, FunctionClass)
                        @test findall(x->x==1,N)==[i]
                    end 
                elseif ndim == 2
                    NodeList, IEN, ID, IEN_top, IEN_btm, BorderNodesList = meshgrid_square(x0,x1,y0,y1,ne,ndim;FunctionClass=FunctionClass)
                    iter = 1:size(IEN,2)
                    for i in iter
                        coord = NodeList[:,IEN[i]]
                        N, dN = basis_function(coord[1],coord[2],nothing, FunctionClass)
                        @test findall(x->x==1,N)==[i]
                    end 
                elseif ndim == 3
                    NodeList, IEN, ID, IEN_top, IEN_btm, BorderNodesList = meshgrid_cube(x0,x1,y0,y1,z0,z1,ne,ndim,FunctionClass=FunctionClass)
                    iter = 1:size(IEN,2)
                    for i in iter
                        coord = NodeList[:,IEN[i]]
                        N, dN = basis_function(coord[1],coord[2],coord[3], FunctionClass)
                        @test findall(x->x==1,N)==[i]
                    end 
                end
            end
        end
    end

    # test gaussian quadrature
    @test gaussian_quadrature(-1,1,nGaussPoints=2) == ([-1/√3, 1/√3], [1.0, 1.0])
    @test gaussian_quadrature(-1,1,nGaussPoints=3) == ([-√(3/5), 0.0, √(3/5)], [5/9, 8/9, 5/9])

    # test 2D FEM solution for 2x2 elements

    # test 2D FEM solution for annulus ring with 2x2 elements

    # test 3D FEM solution for 2x2x2 elements
    
end
