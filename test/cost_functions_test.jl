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

    @testset "ChamferCost normalizes each direction by its own point count" begin
        # Same curve, three samplings. The observed contour is the one that varies in
        # practice — a segmentation returns every boundary pixel, a projected mesh returns
        # its nodes — so the cost must not move when only that sampling does.
        circle(n, r, dx=0.0) = permutedims(hcat([dx + r*cos(t) for t in range(0, 2π, length=n+1)[1:n]],
                                                [r*sin(t) for t in range(0, 2π, length=n+1)[1:n]]))
        # Resolved finely enough that the polygon discretization is not what is being
        # measured; only the observed sampling changes, by a factor of eight.
        ring = circle(400, 100.0)
        c_sparse, _ = contour_cost([ring], [circle(800,  102.0)]; cost=ChamferCost())
        c_dense,  _ = contour_cost([ring], [circle(6400, 102.0)]; cost=ChamferCost())
        # Not exact: the reverse term measures to the simulated *vertices*, not to the edges
        # between them, so a percent of sampling dependence survives the normalization. What
        # it no longer does is drift with the reverse term's share of the residual vector,
        # which over this 8× change would have moved the pooled cost by tens of percent.
        @test isapprox(c_sparse[1], c_dense[1]; rtol=2e-2)

        # With both directions weighted equally and the sampling matched, the symmetric cost
        # collapses onto the one-directional one — which is what lets a `λ` calibrated on
        # ClosestPointCost carry over to a Chamfer fit.
        shifted = circle(400, 100.0, 3.0)
        ch, _ = contour_cost([shifted], [ring]; cost=ChamferCost())
        cp, _ = contour_cost([shifted], [ring]; cost=ClosestPointCost())
        @test isapprox(ch[1], cp[1]; rtol=1e-12)
    end

    @testset "gradients match finite differences" begin
        # The sim contour is an explicit function of θ = (radius, x-offset), so `dudθ` is
        # exact and the FEM solve is out of the picture: what is under test is only the
        # residual/Jacobian algebra. It has teeth for `ChamferCost` in particular, whose
        # per-direction weights multiply `u` and `J` — apply them to one and not the other
        # and `Jᵀu` stops being the gradient of `uᵀu`, which no cost-value test would catch.
        circle(n, r, dx) = permutedims(hcat([dx + r*cos(t) for t in range(0, 2π, length=n+1)[1:n]],
                                            [     r*sin(t) for t in range(0, 2π, length=n+1)[1:n]]))
        function sens(n)
            φ = range(0, 2π, length=n+1)[1:n]
            d = zeros(2, n, 2)
            d[1, :, 1] = cos.(φ); d[2, :, 1] = sin.(φ)   # ∂p/∂radius
            d[1, :, 2] .= 1.0                            # ∂p/∂x-offset
            return d
        end

        θ, hstep = [97.0, 3.0], 1e-5
        for c in (ClosestPointCost(), ChamferCost()),
            (nsim, nobs) in ((400, 2000), (60, 2900))    # the second: ~48 obs per sim point
            ref = circle(nobs, 100.0, 0.0)
            C(p) = contour_cost([circle(nsim, p[1], p[2])], [ref]; cost=c)[1][1]
            _, dc, _, _ = contour_cost([circle(nsim, θ[1], θ[2])], [ref], [sens(nsim)]; cost=c)

            g_fd = map(eachindex(θ)) do i
                e = zeros(length(θ)); e[i] = hstep * abs(θ[i])
                (C(θ .+ e) - C(θ .- e)) / (2e[i])
            end
            @test isapprox(vec(dc[1]), g_fd; rtol=1e-5)
        end
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
