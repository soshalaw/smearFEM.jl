using smearFEM
using LinearAlgebra
using Test

@testset "extraction test for Bsplines and NURBS" begin

    filePath = "/home/soshala/SMEAR-PhD/smear-modules/smearFEM.jl/cylindergen"

    CPointList, W, C, IEN, vol_BSpline, vol_NURBS = read_h5(string(filePath,"/cylinder.h5"),"test")
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

    # construct the parent element
    lPoints = zeros(27,3)
    lPoints[1,:] = [-1 , -1, -1]
    lPoints[2,:] = [1 , -1, -1]
    lPoints[3,:] = [1 , 1, -1]
    lPoints[4,:] = [-1 , 1, -1]

    lPoints[5,:] = [-1 , -1, 1]
    lPoints[6,:] = [1 , -1, 1]
    lPoints[7,:] = [1 , 1, 1]
    lPoints[8,:] = [-1 , 1, 1]
    
    lPoints[9,:] = [0 , -1, -1]
    lPoints[10,:] = [1 , 0, -1]
    lPoints[11,:] = [0 , 1, -1]
    lPoints[12,:] = [-1 , 0, -1]

    lPoints[13,:] = [0 , -1, 1]
    lPoints[14,:] = [1 , 0, 1]
    lPoints[15,:] = [0 , 1, 1]
    lPoints[16,:] = [-1 , 0, 1]

    lPoints[17,:] = [-1 , -1, 0]
    lPoints[18,:] = [1 , -1, 0]
    lPoints[19,:] = [1 , 1, 0]
    lPoints[20,:] = [-1 , 1, 0]

    lPoints[21,:] = [0 , -1, 0]
    lPoints[22,:] = [1 , 0, 0]
    lPoints[23,:] = [0 , 1, 0]
    lPoints[24,:] = [-1 , 0, 0]

    lPoints[25,:] = [0, 0, -1]
    lPoints[26,:] = [0, 0, 1]
    lPoints[27,:] = [0, 0, 0]

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
    nodeiter = 1:size(lPoints,1)

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

            vol = vol + w
        end
    end

    @test vol ≈ vol_NURBS atol=1e-5
end
