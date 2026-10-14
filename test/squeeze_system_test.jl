using smearFEM
using Test
using LinearAlgebra
using SparseArrays

# Internals of the fixed-pattern solver behind `simulate` (src/solver/squeeze_system.jl).
const SQ = smearFEM

@testset "squeeze system" begin

    r, h = 25.0, 40.0
    η, β = 100.0, 100.0
    sim_time, dt = 1.0, 0.1
    nsteps = round(Int, sim_time / dt)
    F = -9518.61 * ones(nsteps)
    conditions = Conditions(camera_matrix=get_camera_matrix(), obj_pose=[150.0, 0.0, 20.0],
                            viewing_angles=[0.0, 45.0])
    problem(cparam, control) = def_problem(Cylinder(r, h), 6, η, :Hex, 2, 3, :Hex, 1, 1, :Hex, 2,
                                           β, cparam, control, "constant", sim_time, dt;
                                           GMESH_MESH=false)

    @testset "assembled system" begin
        model, scene = problem(F, "force")
        cache = BasisFunctionCache(model)
        sys = SQ._squeeze_system(model, scene)
        SQ._assemble_squeeze_system!(sys, model, cache, η, β)
        K = sys.K

        # the matrix the solver used to build with sparse products, unknowns [q_f; μ; p]
        b = SQ.apply_boundary_conditions(model, cache)
        A = η * assemble_system_A(model, cache) + β * b
        B = assemble_system_B(model, cache)
        T = [sparse(scene.C_uc) sparse(scene.q_d[:top])]
        Z = spzeros(size(B, 2), size(B, 2))
        M = [T' * A * T T' * B; B' * T Z]
        Mb = [T' * b * T spzeros(size(T, 2), size(B, 2)); spzeros(size(B, 2), size(T, 2)) Z]
        nf, np = sys.nf, sys.np
        perm = [1:nf; (nf + 2):(nf + 1 + np); nf + 1]
        @test norm(M[perm, perm] - K, 1) <= 1e-12 * norm(K, 1)
        @test norm(Mb[perm, perm] - sys.dKdβ, 1) <= 1e-12 * norm(Mb, 1)

        # residuals are evaluated as K'x
        @test norm(K - K', 1) <= 1e-14 * norm(K, 1)
    end

    @testset "homogeneity in (η, β)" begin
        # simulate derives ∂/∂η from ∂/∂β: η ∂x/∂η + β ∂x/∂β = -[u; 0] under force control,
        # [0; p; 0] under velocity control
        for control in ("force", "velocity")
            model, scene = problem(F, control)
            cache = BasisFunctionCache(model)
            sys = SQ._squeeze_system(model, scene)
            SQ._assemble_squeeze_system!(sys, model, cache, 1.0, 0.0)
            K10 = copy(sys.K)
            SQ._assemble_squeeze_system!(sys, model, cache, 0.0, 0.0)
            Ā = K10 - sys.K
            SQ._assemble_squeeze_system!(sys, model, cache, η, β)
            n, nf, np = size(sys.K, 1), sys.nf, sys.np
            rhs = zeros(n); rhs[end] = 1.0
            control == "velocity" && SQ._prescribe_plate_velocity!(sys, rhs)
            x = sys.K \ rhs
            dη, dβ = -(Ā * x), -(sys.dKdβ * x)
            if control == "velocity"
                dη[end] = 0.0; dβ[end] = 0.0
            end
            sη, sβ = sys.K \ dη, sys.K \ dβ
            expected = control == "force" ? -copy(x) : zero(x)
            expected[(nf + 1):(nf + np)] = control == "force" ? zeros(np) : x[(nf + 1):(nf + np)]
            @test norm(η * sη + β * sβ - expected) <= 1e-10 * norm(expected)
        end
    end

    @testset "force and velocity control agree" begin
        model_f, scene_f = problem(F, "force")
        out_f = simulate(model_f, scene_f, conditions)
        v = out_f[1] ./ dt                  # plate velocity of each force-controlled step
        model_v, scene_v = problem(v, "velocity")
        out_v = simulate(model_v, scene_v, conditions)
        @test out_v[1] ≈ out_f[1] rtol=1e-12
        for t in eachindex(out_f[4])
            @test out_v[4][t] ≈ out_f[4][t] rtol=1e-8          # node positions
            @test out_v[6][t] ≈ out_f[6][t] rtol=1e-6          # pressure
        end
    end

    @testset "velocity-control sensitivities vs finite differences" begin
        v = -10.0 * ones(nsteps)
        model, scene = problem(v, "velocity")
        _, _, _, pos3D, _, _, grad3D = simulate(model, scene, conditions)
        Δ = 1e-3
        fd(θset!) = begin
            θset!(+Δ); p = simulate(model, scene, conditions)[4]
            θset!(-Δ); m = simulate(model, scene, conditions)[4]
            θset!(0.0)
            [(p[t] - m[t]) / (2Δ) for t in eachindex(p)]
        end
        fd_η = fd(d -> (model.η = [η + d]))
        fd_β = fd(d -> (scene.β = [β + d]))
        # After one step the geometry does not depend on θ yet, so the sensitivities — which
        # leave out ∂K/∂x·∂x/∂θ, as they always have — are exact there; at 25 % compression
        # later steps differ from finite differences by ~1 %.
        @test grad3D[2][:, :, 1] ≈ fd_η[2] atol=1e-6 rtol=1e-4
        @test grad3D[2][:, :, 2] ≈ fd_β[2] atol=1e-6 rtol=1e-4
    end

    @testset "contour sensitivities of a rotated view vs finite differences" begin
        # the views turn the object about z; its sensitivities must turn with it
        model, scene = problem(F, "force")
        function contours()
            out = simulate(model, scene, conditions)
            return project_contours(model, out[4], conditions; grad3D=out[7])
        end
        borders, grads, _ = contours()
        Δ = 1e-3
        model.η = [η + Δ]; borders_p, _, _ = contours()
        model.η = [η - Δ]; borders_m, _, _ = contours()
        model.η = [η]
        for v in eachindex(conditions.viewing_angles)
            fd = (borders_p[2][v] - borders_m[2][v]) / (2Δ)   # border points after one step (exact there)
            @test grads[2][v][1:2, :, 1] ≈ fd atol=1e-8
        end
    end

    @testset "repeatable and prediction-only" begin
        model, scene = problem(F, "force")
        a = simulate(model, scene, conditions)
        b = simulate(model, scene, conditions)
        @test all(a[k] == b[k] for k in eachindex(a))        # the cached system keeps no state between calls
        c = simulate(model, scene, conditions; sensitivities=false)
        @test length(c[7]) == 1
        for t in eachindex(a[4])
            @test c[4][t] ≈ a[4][t] rtol=1e-8
        end
    end

    @testset "contours projected after the solve" begin
        model, scene = problem(F, "force")
        _, _, _, pos3D, _, _, grad3D = simulate(model, scene, conditions)
        borders, grads, pos2D = project_contours(model, pos3D, conditions; grad3D=grad3D)
        @test length(borders) == length(grads) == length(pos2D) == length(pos3D)
        @test all(length(b) == length(conditions.viewing_angles) for b in borders)
        @test all(all(iszero, g) for g in grads[1])      # the initial contour is θ-independent
        borders0, grads0, pos2D0 = project_contours(model, pos3D, conditions)
        @test isempty(grads0) && borders0 == borders && pos2D0 == pos2D
        c = simulate(model, scene, conditions; sensitivities=false)
        @test_throws DimensionMismatch project_contours(model, c[4], conditions; grad3D=c[7])
    end

    @testset "UMFPACK backend" begin
        # MKL Pardiso is used wherever MKL is available, so force the fallback explicitly
        model, scene = problem(F, "force")
        a = simulate(model, scene, conditions)
        sys = SQ._cached_squeeze_system(model, scene)
        sys.lu = SQ._Umfpack(nothing, nothing)
        b = simulate(model, scene, conditions)
        @test sys.lu isa SQ._Umfpack
        for t in eachindex(a[4])
            @test b[4][t] ≈ a[4][t] rtol=1e-8                       # node positions
            @test b[7][t] ≈ a[7][t] rtol=1e-6 atol=1e-12            # their sensitivities
        end
    end

    @testset "unsupported boundary drives" begin
        model, scene = problem(F, "force")
        @test_throws ArgumentError simulate(model, scene, conditions; μu_btm=1.0)
        @test_throws ArgumentError simulate(model, scene, conditions; μu_side=1.0)
    end
end
