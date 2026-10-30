using Test, LinearAlgebra
using smearFEM
using smearFEM: _close_field, _grid_centres

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
end
