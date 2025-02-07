using smearFEM
using LinearAlgebra
using Test

@testset "extraction test for Bsplines and NURBS" begin

    filePath = "/home/soshala/SMEAR-PhD/smear-modules/smearFEM.jl/cylindergen"

    CPointList, W, C_new, IEN, IEN_top, C_top_new, IEN_btm, C_btm_new, vol_BSpline, vol_NURBS, area_BSpline, area_NURBS = read_h5(string(filePath,"/cylinder.h5"),"test")

    NodeList = zeros(Float64, 3, size(IEN,2)*size(IEN,1))
    IEN_list = zeros(Int64, size(IEN))

    ndim = size(CPointList,1)
    ne = Int64(round((size(IEN,2))^(1/ndim)))

    # dummy solution
    q = ones(3,size(CPointList,2))

    ID = zeros(Int64, ndim, size(CPointList,2))
    cpiter = 1:size(CPointList,2)

    for m in cpiter
        for l in 1:ndim
            ID[l,m] = ndim*(m-1) + l
        end
    end 

    if ndim == 1
        # gaussian quadrature points for the element [-1,1] 
        ξ, w_ξ = gaussian_quadrature(-1,1)
        
        wpoints =  [w_ξ[1], w_ξ[2]]
        
        x = [ξ[1], ξ[2]]
    elseif ndim == 2
        # gaussian quadrature points for the element [-1,1]x[-1,1] 
        ξ, w_ξ = gaussian_quadrature(-1,1,nGaussPoints=3)
        η, w_η = gaussian_quadrature(-1,1,nGaussPoints=3)

        x = Float64[]
        y = Float64[]
        wpoints =  Float64[]
        
        n = 1:size(ξ,1)
        m = 1:size(η,1)
        for j in m # loop over η
            for i in n # loop over ξ
                push!(x, ξ[i])
                push!(y, η[j])
                push!(wpoints, w_ξ[i]*w_η[j])
            end
        end

    elseif ndim == 3
        # gaussian quadrature points for the element [-1,1]x[-1,1]x[-1,1] 
        ξ, w_ξ = gaussian_quadrature(-1,1,nGaussPoints=3)
        η, w_η = gaussian_quadrature(-1,1,nGaussPoints=3)
        ζ, w_ζ = gaussian_quadrature(-1,1,nGaussPoints=3)

        x = Float64[]
        y = Float64[]
        z = Float64[]
        wpoints = Float64[]
        
        l = 1:size(ζ,1)
        m = 1:size(η,1)
        n = 1:size(ξ,1)
        for k in l # loop over ζ
            for j in m # loop over η
                for i in n # loop over ξ
                    push!(x, ξ[i])
                    push!(y, η[j])
                    push!(z, ζ[k])
                    push!(wpoints, w_ξ[i]*w_η[j]*w_ζ[k])
                end
            end
        end
    end
    vol = 0

    eiter = 1:size(IEN,2)

    gpiter = 1:length(wpoints)
    for gp in gpiter
        for e in eiter 
            cIter = 1:size(C,2)
            for i in cIter
                c = C[:,:,e]
                @test 1.0 ≈ sum(c[:,i]) atol=1e-5
            end

            Be, ΔBe = basis_function(x[gp], y[gp], z[gp], C[:,:,e], W[IEN[:,e]], "S2")

            for b in Be 
                if abs(b) > 1e-10
                    @test (b > 0)
                end
            end 

            @test 1.0 ≈ sum(Be) atol=1e-5

            coords = CPointList[:,IEN[:,e]] # get the coordinates of the nodes of the element
            Jac  = coords*ΔBe # Jacobian matrix [dx/dxi dx/deta; dy/dxi dy/deta]

            w = wpoints[gp]*abs(det(Jac))

            vol += w
        end
    end
    @test vol ≈ vol_NURBS atol=1e-5
end

@testset "area test for BSplines and NURBS" begin

    filePath = "/home/soshala/SMEAR-PhD/smear-modules/smearFEM.jl/cylindergen"

    CPointList, W, C_new, IEN, IEN_top, C_top_new, IEN_btm, C_btm_new, vol_BSpline, vol_NURBS, area_BSpline, area_NURBS = read_h5(string(filePath,"/cylinder.h5"),"test")

    NodeList = zeros(Float64, 3, size(IEN,2)*size(IEN,1))
    IEN_list = zeros(Int64, size(IEN))

    ndim = size(CPointList,1)
    ne = Int64(round((size(IEN,2))^(1/ndim)))

    # dummy solution
    q = ones(3,size(CPointList,2))

    ID = zeros(Int64, ndim, size(CPointList,2))
    cpiter = 1:size(CPointList,2)

    for m in cpiter
        for l in 1:ndim
            ID[l,m] = ndim*(m-1) + l
        end
    end 

    mdl = def_model("linear_elasticity", ne=ne, NodeList=CPointList, IEN=IEN, ndim=ndim, nDof=ndim, ID = ID,
    FunctionClass="Q2", C = C, W = W)

    if mdl.ndim == 2
        # gaussian quadrature points for the element [-1,1] 
        ξ, w_ξ = gaussian_quadrature(-1,1)

        wpoints =  [w_ξ[1], w_ξ[2]]
        
        x = [ξ[1], ξ[2]]
    elseif mdl.ndim == 3            
        # gaussian quadrature points for the element [-1,1]x[-1,1] 
        ξ, w_ξ = gaussian_quadrature(-1,1,nGaussPoints=3)
        η, w_η = gaussian_quadrature(-1,1,nGaussPoints=3)

        x = Float64[]
        y = Float64[]
        wpoints =  Float64[]
        
        n = 1:size(ξ,1)
        m = 1:size(η,1)
        for j in m # loop over η
            for i in n # loop over ξ
                push!(x, ξ[i])
                push!(y, η[j])
                push!(wpoints, w_ξ[i]*w_η[j])
            end
        end
    end 
    area_top = 0
    area_btm = 0

    e_iter = 1:mdl.ne^(mdl.ndim-1)
    # integration loop
    gpiter = 1:length(wpoints)
    for gp in gpiter
        for e in e_iter
            N_top, ΔN_top = basis_function(x[gp], y[gp], C_top[:,:,e], W[IEN_top[:,e]], "S2") 
            N_btm, ΔN_btm = basis_function(x[gp], y[gp], C_btm[:,:,e], W[IEN_btm[:,e]], "S2") 

            coords_top = mdl.NodeList[:,IEN_top[:,e]] # get the coordinates of the nodes of the element
            coords_btm = mdl.NodeList[:,IEN_btm[:,e]] # get the coordinates of the nodes of the element

            dxdξ_top = coords_top*ΔN_top         # Jacobian matrix [dx/dxi dx/deta; dy/dxi dy/deta; dz/dxi dz/deta]
            dxdξ_btm = coords_btm*ΔN_btm         # Jacobian matrix [dx/dxi dx/deta; dy/dxi dy/deta; dz/dxi dz/deta]

            w_top = wpoints[gp]*norm(cross(dxdξ_top[:,1],dxdξ_top[:,2]))     # weight of the quadrature point top surface
            w_btm = wpoints[gp]*norm(cross(dxdξ_btm[:,1],dxdξ_btm[:,2]))     # weight of the quadrature point bottom surface
        
            area_top += w_top
            area_btm += w_btm
        end
    end
    @test area_btm ≈ area_top atol=1e-5
    @test area_top ≈ area_NURBS atol=1e-5
end
