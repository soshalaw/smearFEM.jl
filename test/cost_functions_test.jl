# Contour cost functions: residual/Jacobian construction and the ContourCost types.
# Pure geometry — no FEM solve — so this is cheap enough for runtests.jl.

using smearFEM, LinearAlgebra, Random, Test

@testset "cost_functions" begin
    Random.seed!(1)
    nθ  = 2
    sim = [rand(2, 40) for _ in 1:3]
    obs = [rand(2, 55) for _ in 1:3]
    du  = [rand(2, 40, nθ) for _ in 1:3]

    # The pre-refactor inline computation, verbatim. ClosestPointCost must reproduce it
    # bit-for-bit, or every previously recorded cost and every λ calibrated against one
    # silently shifts.
    function reference(sim_frames, obs_frames, dudθ)
        cost, dcost, d2cost = Float64[], [], []
        for (obs_t, sim_t, du_tdθ) in zip(obs_frames, sim_frames, dudθ)
            pairs = match_points(sim_t, obs_t)
            pSim, qSim = sim_t[1,:], sim_t[2,:]
            dpSim, dqSim = du_tdθ[1,:,:], du_tdθ[2,:,:]
            pObs, qObs = obs_t[1,:], obs_t[2,:]
            u = [(pSim[pairs[:,1]] - pObs[pairs[:,2]]); (qSim[pairs[:,1]] - qObs[pairs[:,2]])]
            J = [dpSim[pairs[:,1],:]; dqSim[pairs[:,1],:]]
            push!(cost,   (u'*u)/(2*length(pairs)))
            push!(dcost,  (J'*u)/length(pairs))
            push!(d2cost, (J'*J)/length(pairs))
        end
        return cost, dcost, d2cost
    end

    @testset "ClosestPointCost reproduces the historical numbers" begin
        rc, rd, rh = reference(sim, obs, du)
        nc, nd, nh, _ = contour_cost(sim, obs, du)          # default cost
        @test nc == rc
        @test all(nd .== rd)
        @test all(nh .== rh)

        ec, ed, eh, _ = contour_cost(sim, obs, du; cost=ClosestPointCost())
        @test ec == nc                                        # explicit == implicit default
        @test all(ed .== nd)
        @test all(eh .== nh)
    end

    @testset "ChamferCost" begin
        cc, cd, ch, cp = contour_cost(sim, obs, du; cost=ChamferCost())
        oc, _, _, _    = contour_cost(sim, obs, du)

        @test all(isfinite, cc)
        @test all(x -> all(isfinite, x), cd)
        @test all(x -> all(isfinite, x), ch)
        @test length(cc) == length(sim)
        # Gauss-Newton Hessians must stay symmetric PSD for the Newton solve to be sane.
        @test all(h -> isapprox(h, h'; rtol=1e-12), ch)
        @test all(h -> minimum(eigvals(Symmetric(h))) > -1e-10, ch)
        # Both correspondences are returned (forward and reverse).
        @test cp[1] isa Tuple && length(cp[1]) == 2
        @test cc != oc                                        # genuinely a different cost

        # The defining property: symmetric in its arguments, where one-sided is not.
        s1, _ = contour_cost(sim, obs; cost=ChamferCost())
        s2, _ = contour_cost(obs, sim; cost=ChamferCost())
        @test isapprox(s1, s2; rtol=1e-12)
        o1, _ = contour_cost(sim, obs; cost=ClosestPointCost())
        o2, _ = contour_cost(obs, sim; cost=ClosestPointCost())
        @test !isapprox(o1, o2; rtol=1e-6)
    end

    @testset "both costs vanish on identical clouds" begin
        z = [rand(2, 30)]
        for c in (ClosestPointCost(), ChamferCost())
            zc, _ = contour_cost(z, z; cost=c)
            @test zc[1] == 0.0
        end
    end

    # The reporting counterpart of ChamferCost. Kept next to the costs deliberately: the
    # whole point of the squared metric is that it can be read against an optimization cost.
    @testset "chamfer_sq_distance_kdtree" begin
        Random.seed!(2)
        a, b = rand(40, 2), rand(55, 2)

        mp  = sum(minimum(sum((b .- permutedims(p)).^2, dims=2)) for p in eachrow(a)) / size(a, 1)
        mg  = sum(minimum(sum((a .- permutedims(q)).^2, dims=2)) for q in eachrow(b)) / size(b, 1)
        ref = 0.5 * (mp + mg)

        sq = chamfer_sq_distance_kdtree(a, b)
        @test isapprox(sq, ref; rtol=1e-12)                       # matches brute force
        @test isapprox(sq, chamfer_sq_distance_kdtree(b, a); rtol=1e-12)   # symmetric
        @test chamfer_sq_distance_kdtree(a, a) == 0.0

        # Squaring is not a post-hoc transform of the unsquared metric: the mean of squares
        # is not the square of the mean, so one cannot be recovered from the other.
        @test !isapprox(sq, chamfer_distance_kdtree(a, b)^2; rtol=1e-6)

        _, c_sq, _ = compare_pt_clouds([a], [b])
        _, c_un, _ = compare_pt_clouds([a], [b]; squared_chamfer=false)
        @test c_sq[1] == sq                                        # squared by default
        @test c_un[1] == chamfer_distance_kdtree(a, b)             # opt-out still works
    end

    @testset "outlier frames are skipped" begin
        for c in (ClosestPointCost(), ChamferCost())
            full, _    = contour_cost(sim, obs; cost=c)
            skipped, _ = contour_cost(sim, obs; outliers=[2], cost=c)
            @test length(skipped) == length(full) - 1
            @test skipped == full[[1, 3]]
        end
    end
end
