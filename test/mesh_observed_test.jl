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

    @testset "near-zero values make no degenerate triangles" begin
        # perception's float arithmetic leaves ±1e-12 where the surface passes through a voxel
        # centre; marching cubes then makes needles too thin for Gmsh's parametrization
        tiny(x, y, z) = abs(cylinder_sd(x, y, z)) < 1e-9 ? (isodd(round(Int, x + y + z)) ? 1e-12 : -1e-12) : nothing
        V, T, _ = _field_surface(_close_field(tsdf_cylinder(extra=tiny)..., 4.685)..., 4.685, 1.0)
        # 1 for an equilateral triangle, 0 for a degenerate one
        quality(a, b, c) = 2√3 * norm(cross(b - a, c - a)) / (sum(abs2, b - a) + sum(abs2, c - b) + sum(abs2, a - c))
        @test minimum(quality(V[:, T[1, t]], V[:, T[2, t]], V[:, T[3, t]]) for t in axes(T, 2)) > 1e-5
    end

    @testset "floating material" begin
        _, grid = tsdf_cylinder()
        xs, ys, zs = _grid_centres(grid)
        sphere = Float32[norm((x, y, z - 26)) - 10 for x in xs, y in ys, z in zs]
        @test_throws ArgumentError _field_surface(_close_field(sphere, grid, 4.685)..., 4.685, 1.0)
    end

    observed = Dict{Int,MeshgridObserved}()
    @testset "observed mesh of a cylinder (Tet$order)" for order in (1, 2)
        m = observed[order] = mesh_from_field(tsdf_cylinder()...; table_z=4.685, elem_size=6.0, basis_order=order)
        X = m.NodeList
        @test m isa MeshgridObserved && m.volume_element_shape == :Tet && size(m.IEN, 1) == (4, 10)[order]
        vols = [det([X[:, e[2]] - X[:, e[1]] X[:, e[3]] - X[:, e[1]] X[:, e[4]] - X[:, e[1]]]) / 6 for e in eachcol(m.IEN)]
        @test all(>(0), vols)
        @test sum(vols) ≈ π * 30^2 * 40 rtol=0.02
        @test all(abs.(X[3, m.bottom_nodes] .- 4.685) .< 1e-9)
        @test all(abs(cylinder_sd(X[:, a]...)) < 1.0 for a in unique(vec(m.IEN_free[1:3, :])))
        cap = sum(norm(cross(X[:, f[2]] - X[:, f[1]], X[:, f[3]] - X[:, f[1]])) / 2 for f in eachcol(m.IEN_bottom))
        @test π * 29^2 < cap < π * 30^2
    end

    @testset "table on a voxel plane" begin
        m = mesh_from_field(tsdf_cylinder(tz=0.0)...; table_z=0.0, elem_size=6.0, basis_order=1)
        @test all(abs.(m.NodeList[3, m.bottom_nodes]) .< 1e-9)
    end

    @testset "Gmsh API path" begin
        m = withenv("PATH" => "/nonexistent") do
            mesh_from_field(tsdf_cylinder()...; table_z=4.685, elem_size=6.0, basis_order=1)
        end
        @test m isa MeshgridObserved && !isempty(m.IEN_bottom) && !isempty(m.IEN_free)
    end

    @testset "volume update on the observed mesh" begin
        m = observed[2]
        X = m.NodeList
        nodes = [a for a in unique(vec(m.IEN_free[1:3, :]))
                 if abs(hypot(X[1, a], X[2, a]) - 30) < 1 && 4.685 + 5 < X[3, a] < 4.685 + 35]
        normals = reduce(hcat, [X[1, a], X[2, a], 0.0] / hypot(X[1, a], X[2, a]) for a in nodes)
        λ, μ = lame_parameters(1.0, 0.3)
        mdl = LinearElasticity(mesh=m, λ=λ, μ=μ)
        u = volume_update!(mdl, nodes, normals, fill(0.5, length(nodes)); β_table=1e-3)
        @test all(abs(dot(u[:, a], normals[:, i]) - 0.5) < 1e-8 for (i, a) in enumerate(nodes))
        @test maximum(abs.(u[3, m.bottom_nodes])) < 1e-10
    end
end
