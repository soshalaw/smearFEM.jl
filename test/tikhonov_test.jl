using smearFEM, LinearAlgebra, Test

r, h = 25.0, 40.0
ndim = 3
fem = (:Hex, 2, ndim, :Hex, 1, 1, :Hex, 2)   # u shape/order/ndof, p shape/order/ndof, x shape/order
ne = 6
sim_time, steps = 0.5, 5
t_steps = sim_time/steps
control, visc = "force", "constant"
F_ext = 9.813e3*0.85
F = -F_ext*ones(Float64, round(Int, sim_time/t_steps))

obj_pose = zeros(Float64,4,4)
obj_pose[1,1] = -1.0; obj_pose[2,3] = -1.0; obj_pose[3,2] = -1.0
obj_pose[1:3,4] = [0.0, h/2, 150.0]
cam = [2.39642674e+03 0.0 1.00429248e+03; 0.0 2.40565353e+03 7.57028161e+02; 0.0 0.0 1.0]
cond = Conditions(camera_matrix=cam, obj_pose=obj_pose)

geom = Cylinder(r,h)
θ_gt = [100.0, 100.0]

println("### generating synthetic observations at θ_gt = $θ_gt")
mdl, scn = def_problem(geom, ne, θ_gt[1], fem..., θ_gt[2], F, control, visc, sim_time, t_steps)
_, _, obsPts, _, _, _, _, _, _, _ = simulate(mdl, scn, cond)
obs = Vector{AbstractArray}(obsPts)
println("### got $(length(obs)) observation frames")

function fresh()
    m, s = def_problem(geom, ne, 70.0, fem..., 130.0, F, control, visc, sim_time, t_steps)
    return m, s
end

results = Dict{String,Any}()

for (name, kw) in (("gn",           (; method=:gn)),
                   ("tikh_λ0",      (; method=:gn_tikhonov, λ_scale=0.0)),
                   ("tikh_λ1",      (; method=:gn_tikhonov, λ_scale=1.0)),
                   ("tikh_prior",   (; method=:gn_tikhonov, λ_scale=1.0, θ_p=[95.0,105.0])),
                   ("tikh_big",     (; method=:gn_tikhonov, λ_scale=1.0e4, θ_p=[95.0,105.0])),
                   ("tikh_vec",     (; method=:gn_tikhonov, λ_scale=[1.0, 0.1])),
                   ("gn_store",     (; method=:gn, store_border_pts=true)),
                   ("tikh_store",   (; method=:gn_tikhonov, λ_scale=0.01, store_border_pts=true)))
    println("\n########## $name  $kw")
    m, s = fresh()
    st = fit_model(m, s, cond, obs, [70.0,130.0]; kw...)
    results[name] = st
    println(">>> $name : η=$(st["η"]) β=$(st["β"]) iters=$(length(st["iterList"])) " *
            "λ=$(get(st,"λ","-")) Hkeys=$(haskey(st,"H"))")
end

println("\n================ CHECKS ================")
g, z = results["gn"], results["tikh_λ0"]
@testset "tikhonov" begin
    # λ_scale = 0 must reproduce plain GN exactly (penalty machinery inert)
    @test isapprox(g["η"], z["η"]; rtol=1e-10)
    @test isapprox(g["β"], z["β"]; rtol=1e-10)
    @test all(iszero, z["λ"])
    # H and λ exposed, H symmetric PSD-ish
    for k in ("tikh_λ0","tikh_λ1","tikh_prior","tikh_big","tikh_vec")
        st = results[k]
        @test haskey(st,"H") && size(st["H"]) == (2,2)
        @test haskey(st,"λ")
        @test all(isfinite, st["H"])
        @test minimum(eigvals(Symmetric(st["H"]))) > -1e-8
    end
    # λ is one absolute weight per parameter; a scalar λ_scale broadcasts to all of them.
    @test results["tikh_λ1"]["λ"] == [1.0, 1.0]
    # ...and a vector is taken entry-wise, so η and β can be regularized differently.
    @test results["tikh_vec"]["λ"] == [1.0, 0.1]
    @test length(results["tikh_vec"]["λ"]) == 2
    # λ_scale must be a scalar or match length(θ), and must be finite and non-negative.
    m, s_ = fresh()
    @test_throws Exception fit_model(m, s_, cond, obs, [70.0,130.0]; method=:gn_tikhonov, λ_scale=[1.0,0.1,0.5])
    m, s_ = fresh()
    @test_throws Exception fit_model(m, s_, cond, obs, [70.0,130.0]; method=:gn_tikhonov, λ_scale=[-1.0,0.1])
    # a huge λ_scale must pin the answer at the prior
    @test isapprox(results["tikh_big"]["η"], 95.0;  rtol=5e-2)
    @test isapprox(results["tikh_big"]["β"], 105.0; rtol=5e-2)

    # --- stored contours ---
    for k in ("gn","tikh_λ1")                      # opt-out default stays empty
        @test isempty(results[k]["simBorderPtsList"])
    end
    for k in ("gn_store","tikh_store")
        st = results[k]
        bp = st["simBorderPtsList"]
        @test length(bp) == length(st["iterList"])          # index-aligned with the histories
        @test length(bp) == length(st["ηList"])
        @test all(f -> length(f) == length(obs), bp)        # each entry: one contour per frame
        @test !isempty(bp[1])
    end
    # the last stored contour must be the one at the returned parameters
    st = results["gn_store"]
    m2, s2 = fresh(); m2.η = [st["η"]]; s2.β = [st["β"]]
    _, _, chk, _, _, _, _, _, _, _ = simulate(m2, s2, cond)
    @test isapprox(vec(smearFEM._as_2xN(st["simBorderPtsList"][end][1])),
                   vec(smearFEM._as_2xN(chk[1])); rtol=1e-8)
end
println("θ_gt = $θ_gt")
for (k,v) in sort(collect(results), by=first)
    println(rpad(k,12), " η=", round(v["η"],sigdigits=6), " β=", round(v["β"],sigdigits=6),
            " λ=", round(get(v,"λ",NaN),sigdigits=4))
end
