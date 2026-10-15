using smearFEM
using LinearAlgebra
using Logging
using Test

const MESH_DIR = joinpath(pkgdir(smearFEM), "mesh_files")

# A 20 × 30 × 40 mm box on the table (z ≥ 0), Q2 hex (structured) or P1/P2 tet (Gmsh).
function test_box(shape::Symbol, order::Int)
    shape == :Hex && return meshgrid_cuboid(20.0, 30.0, 40.0; element_shape=:Hex, basis_order=order, ne=3)
    return meshgrid_cuboid(20.0, 30.0, 40.0; mesh_type=:unstructured, element_shape=shape, basis_order=order,
                           elem_size=10.0, mesh_path=MESH_DIR)
end

# Uniaxial compression between frictionless plates: u_z = ε z, u_r = -ν ε r, F = E A ε.
function uniaxial(X::Matrix{Float64}, ε::Float64, ν::Float64)
    c = [(minimum(X[1, :]) + maximum(X[1, :])) / 2, (minimum(X[2, :]) + maximum(X[2, :])) / 2, minimum(X[3, :])]
    return [-ν * ε, -ν * ε, ε] .* (X .- c)
end

@testset "linear elasticity" begin

    @testset "lame_parameters" begin
        λ, μ = lame_parameters(10.0, 0.25)
        @test λ ≈ 4.0
        @test μ ≈ 4.0
        @test_throws ArgumentError lame_parameters(1.0, 0.5)
    end

    @testset "element matrices match B C Bᵀ ($shape$order)" for (shape, order) in ((:Hex, 1), (:Hex, 2))
        mesh = meshgrid_cuboid(2.0, 3.0, 4.0; element_shape=shape, basis_order=order, ne=1)
        X = mesh.NodeList .+ 0.05 .* reshape(sin.(1:length(mesh.NodeList)), size(mesh.NodeList))   # distorted
        cache = BasisFunctionCache(mesh)
        nN, ngp = size(mesh.IEN, 1), length(cache.w)
        λ, μ = 2.0, 0.7
        KK, KL = zeros((3nN)^2), zeros((3nN)^2)
        smearFEM._element_matrices!(KK, KL, zeros(3, nN), zeros(ngp, 3, nN), zeros(ngp, 3, nN), 1, X, mesh.IEN,
                                    cache.ΔN, cache.w, λ, μ)
        C = get_c_mat(λ, μ; type="lame")
        Kref = zeros(3nN, 3nN)
        for gp in 1:ngp
            J = X[:, mesh.IEN[:, 1]] * cache.ΔN[gp]
            dN = cache.ΔN[gp] / J
            B = zeros(3nN, 6)
            B[1:3:end, 1] = dN[:, 1]; B[2:3:end, 2] = dN[:, 2]; B[3:3:end, 3] = dN[:, 3]
            B[2:3:end, 4] = dN[:, 3]; B[3:3:end, 4] = dN[:, 2]
            B[1:3:end, 5] = dN[:, 3]; B[3:3:end, 5] = dN[:, 1]
            B[1:3:end, 6] = dN[:, 2]; B[2:3:end, 6] = dN[:, 1]
            Kref .+= B * C * B' * cache.w[gp] * abs(det(J))
        end
        @test norm(reshape(KK, 3nN, 3nN) - Kref) / norm(Kref) < 1e-12
        @test norm(reshape(KL, 3nN, 3nN) - Kref) / norm(Kref) > 0.1     # L alone is not K
    end

    # The force-controlled case also guards the simplex quadrature weights: an over-integrated
    # Tet stiffens every element by the same factor, which displacement control cannot see.
    @testset "uniaxial squeeze ($shape$order)" for (shape, order) in ((:Hex, 2), (:Tet, 1), (:Tet, 2))
        mesh = test_box(shape, order)
        E, ν = 10.0, 0.3
        λ, μ = lame_parameters(E, ν)
        mdl = LinearElasticity(mesh=mesh, λ=λ, μ=μ)
        X = mesh.NodeList
        H = maximum(X[3, :]) - minimum(X[3, :])
        A = (maximum(X[1, :]) - minimum(X[1, :])) * (maximum(X[2, :]) - minimum(X[2, :]))
        δ = -0.1
        ua = uniaxial(X, δ / H, ν)
        β = 1e-9                                     # frictionless, but no free rigid-body modes

        r = simulate(mdl, [δ]; control="displacement", β=β)
        u = r.displacement[1]
        @test maximum(abs.(u[3, :] - ua[3, :])) < 1e-9
        # horizontally only up to a rigid translation, which β ≈ 0 barely pins
        uh = u[1:2, :] .- sum(u[1:2, :] - ua[1:2, :], dims=2) ./ size(u, 2)
        @test maximum(abs.(uh - ua[1:2, :])) < 1e-6
        # 1e-6: the residual tolerance (1e-9) times the conditioning of near-frictionless plates;
        # a quadrature-weight error would be a factor of ~8
        @test r.plate[1] ≈ E * A * δ / H rtol=1e-6
        @test r.pos3D[end] ≈ X + u
        @test mdl.mesh.NodeList === r.pos3D[end]

        rf = simulate(mdl, [E * A * δ / H]; control="force", β=β)
        @test rf.plate[1] ≈ δ rtol=1e-6
        @test mdl.mesh.initial_state == X            # the solve never wrote into the initial state
    end

    @testset "sensitivities match finite differences ($control)" for (control, c) in (("force", -2.0), ("displacement", -0.2))
        mesh = meshgrid_cylinder(25.0, 40.0; element_shape=:Hex, basis_order=2, ne=3)
        λ, μ = lame_parameters(10.0, 0.35)
        β = 0.5
        θ = [λ, μ, β]
        squeeze(p; sens=false) = simulate(LinearElasticity(mesh=mesh, λ=p[1], μ=p[2]), [c]; control=control, β=p[3], sensitivities=sens)
        g = squeeze(θ; sens=true).grad3D[end]
        for j in 1:3
            h = 1e-6 * θ[j]
            e = zeros(3); e[j] = h
            fd = (squeeze(θ + e).pos3D[end] - squeeze(θ - e).pos3D[end]) / 2h
            @test norm(g[:, :, j] - fd) / norm(fd) < 1e-4
        end
    end

    @testset "multi-step runs are reproducible" begin
        mesh = meshgrid_cylinder(25.0, 40.0; element_shape=:Hex, basis_order=2, ne=3)
        λ, μ = lame_parameters(10.0, 0.35)
        mdl = LinearElasticity(mesh=mesh, λ=λ, μ=μ)
        r1 = simulate(mdl, fill(-20.0, 30); β=1.0, sensitivities=true)
        simulate(mdl, fill(-5.0, 7); control="displacement", β=3.0)
        r2 = simulate(mdl, fill(-20.0, 30); β=1.0, sensitivities=true)
        @test length(r1.pos3D) == 31
        @test length(r1.grad3D) == 31
        # each call refactorizes at its first step, so earlier calls leave no trace
        @test r1.pos3D[end] == r2.pos3D[end]
        @test r1.grad3D[end] == r2.grad3D[end]
        @test all(<(0), r1.plate)
        @test issorted(cumsum(r1.plate), rev=true)
    end

    # Normal strains plus a rigid rotation about z: zero shear stress on every axis-aligned face,
    # so normal displacements on the top and sides (traction-free tangentially) and the table
    # reproduce the field exactly on any mesh.
    @testset "volume update reproduces a linear field ($shape$order)" for (shape, order) in ((:Hex, 2), (:Tet, 1), (:Tet, 2))
        mesh = test_box(shape, order)
        X = mesh.NodeList
        G = [0.01 -0.004 0.0; 0.004 -0.006 0.0; 0.0 0.0 -0.02]
        lo, hi = minimum(X, dims=2), maximum(X, dims=2)
        ua = G * (X .- [(lo[1] + hi[1]) / 2, (lo[2] + hi[2]) / 2, lo[3]])
        nodes, normals, d = Int[], Vector{Float64}[], Float64[]
        for a in axes(X, 2), i in 1:3, (bound, s) in ((lo, -1.0), (hi, 1.0))
            (i == 3 && s < 0) && continue                       # the table is added by volume_update!
            isapprox(X[i, a], bound[i]; atol=1e-9) || continue
            n = zeros(3); n[i] = s
            push!(nodes, a); push!(normals, 2n); push!(d, s * ua[i, a])   # normals need not be unit
        end
        λ, μ = lame_parameters(1.0, 0.3)
        mdl = LinearElasticity(mesh=mesh, λ=λ, μ=μ)
        u = volume_update!(mdl, nodes, reduce(hcat, normals), d)
        @test maximum(abs.(u - ua)) < 1e-12
        @test mdl.mesh.NodeList ≈ X + u
        @test mdl.mesh.initial_state == X
        # only Poisson's ratio matters
        reset_model!(mdl)
        mdl.λ, mdl.μ = lame_parameters(1e3, 0.3)
        @test volume_update!(mdl, nodes, reduce(hcat, normals), d) ≈ u atol=1e-12
    end

    @testset "volume update meets its constraints on a cylinder" begin
        mesh = meshgrid_cylinder(25.0, 40.0; element_shape=:Hex, basis_order=2, ne=3)
        X = mesh.NodeList
        λ, μ = lame_parameters(1.0, 0.3)
        mdl = LinearElasticity(mesh=mesh, λ=λ, μ=μ)
        r = sqrt.(X[1, :] .^ 2 .+ X[2, :] .^ 2)
        side = findall(r .> 25.0 - 1e-6)
        normals = vcat(X[1:2, side] ./ r[side]', zeros(1, length(side)))
        d = 0.5 .+ 0.01 .* X[3, side]
        # radial normals leave the rotation about z free: slip on the table pins it
        u = volume_update!(mdl, side, normals, d; β_table=1e-3)
        @test maximum(abs.(vec(sum(u[:, side] .* normals, dims=1)) - d)) < 1e-12
        table = unique(vec(mesh.IEN_bottom))
        @test maximum(abs.(u[3, table])) < 1e-12
        @test maximum(abs.(u[3, :])) > 0.01                    # the interior moved
    end

    # Neither sparse solver refuses a singular stiffness; both backends must say so.
    @testset "singular stiffness is reported" begin
        mesh = meshgrid_cylinder(25.0, 40.0; element_shape=:Hex, basis_order=2, ne=3)
        X = mesh.NodeList
        λ, μ = lame_parameters(1.0, 0.3)
        mdl = LinearElasticity(mesh=mesh, λ=λ, μ=μ)
        r = sqrt.(X[1, :] .^ 2 .+ X[2, :] .^ 2)
        side = findall(r .> 25.0 - 1e-6)
        normals = vcat(X[1:2, side] ./ r[side]', zeros(1, length(side)))
        d = fill(0.5, length(side))
        @test_logs (:warn, r"singular") match_mode=:any volume_update!(mdl, side, normals, d)
        reset_model!(mdl)
        @test_logs min_level=Logging.Warn volume_update!(mdl, side, normals, d; β_table=1e-3)
        @test_logs (:warn, r"singular") match_mode=:any simulate(mdl, [-0.1]; control="displacement", β=0.0)
        @test_logs min_level=Logging.Warn simulate(mdl, [-0.1]; control="displacement", β=1e-9)
    end
end
