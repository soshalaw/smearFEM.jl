using Test, LinearAlgebra
using smearFEM
using smearFEM: _close_field, _grid_centres, _field_surface

cylinder_sd(x, y, z; r=30.0, h=40.0, tz=4.685) =
    let d = (hypot(x, y) - r, max(tz - z, z - (tz + h)))
        norm(max.(d, 0)) + min(max(d...), 0)
    end

# A cylinder standing on the table at z = tz as perception leaves it: NaN more than μ outside or
# `depth` inside (the distance cap) and below the table. `extra(x, y, z)` may return a value to
# overwrite the voxel with (NaN for a hole, a negative value for a blob), or `nothing`.
function tsdf_cylinder(; r=30.0, h=40.0, tz=4.685, Δ=1.0, μ=3.0, depth=7.0, nz=53, extra=nothing)
    n = round(Int, 36 / Δ)
    grid = PerceptionGrid([0.0, 0.0, 0.0], Δ, [n, n, 0], [2n + 1, 2n + 1, nz])
    xs, ys, zs = _grid_centres(grid)
    phi = Array{Float32}(undef, length(xs), length(ys), length(zs))
    for (i, x) in enumerate(xs), (j, y) in enumerate(ys), (k, z) in enumerate(zs)
        d = cylinder_sd(x, y, z; r, h, tz)
        v = (z < tz || d > μ || d < -depth) ? NaN : d
        e = isnothing(extra) ? nothing : extra(x, y, z)
        phi[i, j, k] = isnothing(e) ? v : e
    end
    return phi, grid
end

nearest(c, v) = argmin(abs.(c .- v))

@testset "observed mesh" begin

    @testset "closing a TSDF-like cylinder" begin
        F, (px, py, pz) = _close_field(tsdf_cylinder()..., 4.685)
        @test size(F) == (75, 75, 55)
        @test !any(isnan, F)
        @test !any(iszero, F)
        @test px[1] == -37.0 && pz[1] == -1.0
        @test all(sign(F[i, j, k]) == sign(cylinder_sd(x, y, z))
                  for (i, x) in enumerate(px), (j, y) in enumerate(py), (k, z) in enumerate(pz)
                  if abs(cylinder_sd(x, y, z)) > 1.0)
    end

    @testset "a hole through the side wall is not closed" begin
        hole(x, y, z) = (21 ≤ x ≤ 35 && abs(y) ≤ 3 && 20 ≤ z ≤ 30) ? NaN : nothing
        phi, grid = tsdf_cylinder(extra=hole)
        err = try _close_field(phi, grid, 4.685); nothing catch e; e end
        @test err isa ErrorException && occursin("not closed", err.msg)
    end

    @testset "a detached lip ring is dropped" begin
        lip(x, y, z) = (32 ≤ hypot(x, y) ≤ 33 && 41 ≤ z ≤ 43) ? -0.5 :
                       (30 < hypot(x, y) < 32 && 41 ≤ z ≤ 43) ? 0.5 : nothing
        F, (px, py, pz) = @test_logs (:warn, r"spurious") match_mode=:any _close_field(tsdf_cylinder(extra=lip)..., 4.685)
        @test F[nearest(px, 32.5), nearest(py, 0.0), nearest(pz, 42.0)] > 0
        @test F[nearest(px, 0.0), nearest(py, 0.0), nearest(pz, 24.0)] < 0
    end

    @testset "a spurious blob is dropped" begin
        blob(x, y, z) = norm((x, y, z - 49)) ≤ 1.5 ? -0.5 : nothing
        F, (px, py, pz) = @test_logs (:warn, r"spurious") match_mode=:any _close_field(tsdf_cylinder(extra=blob)..., 4.685)
        @test F[nearest(px, 0.0), nearest(py, 0.0), nearest(pz, 49.0)] > 0
        @test F[nearest(px, 0.0), nearest(py, 0.0), nearest(pz, 24.0)] < 0
    end

    @testset "material at the grid side" begin
        @test_logs (:warn, r"grid boundary") match_mode=:any _close_field(tsdf_cylinder(r=36.5)..., 4.685)
    end

    @testset "bad input" begin
        phi, grid = tsdf_cylinder()
        @test_throws ArgumentError _close_field(phi, grid, -5.0)
        @test_throws ArgumentError _close_field(fill(NaN32, size(phi)), grid, 4.685)
    end

    @testset "surface of the closed field" begin
        V, T, cap = _field_surface(_close_field(tsdf_cylinder()..., 4.685)..., 4.685, 1.0)
        edges = Dict{Tuple{Int,Int},Int}()
        for t in axes(T, 2), (a, b) in ((1, 2), (2, 3), (3, 1))
            e = minmax(T[a, t], T[b, t])
            edges[e] = get(edges, e, 0) + 1
        end
        @test all(==(2), values(edges))
        @test abs(sum(det(V[:, T[:, t]]) for t in axes(T, 2))) / 6 ≈ π * 30^2 * 40 rtol=0.005
        @test count(cap) > 0
        @test all(abs.(V[3, unique(vec(T[:, cap]))] .- 4.685) .< 1e-3)
        @test all(norm(cross(V[:, T[2, t]] - V[:, T[1, t]], V[:, T[3, t]] - V[:, T[1, t]])) > 0 for t in axes(T, 2))
    end

    @testset "floating material" begin
        _, grid = tsdf_cylinder()
        xs, ys, zs = _grid_centres(grid)
        sphere = Float32[norm((x, y, z - 26)) - 10 for x in xs, y in ys, z in zs]
        @test_throws ArgumentError _field_surface(_close_field(sphere, grid, 4.685)..., 4.685, 1.0)
    end
end
