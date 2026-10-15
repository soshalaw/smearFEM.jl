"""
    benchmark.jl

Time the linear elastic squeeze (`simulate`, per load step) and the volume update
(`volume_update!`, per call) on a Gmsh cylinder.

    julia -O3 --threads=auto --project experiments/elasticity/benchmark.jl [shape] [order] [elem_size]

Defaults: `Hex 2 6.0`, the paper-1 cylinder mesh (r = 25 mm, h = 40 mm). Run it on a quiet
machine: other load moves the numbers by tens of percent.
"""

using smearFEM
using Statistics

function main(shape::Symbol=:Hex, order::Int=2, elem_size::Float64=6.0; nsteps::Int=100)
    mesh = meshgrid_cylinder(25.0, 40.0; mesh_type=:unstructured, element_shape=shape, basis_order=order,
                             elem_size=elem_size, mesh_path=joinpath(pkgdir(smearFEM), "mesh_files"))
    λ, μ = lame_parameters(10.0, 0.35)          # kPa
    mdl = LinearElasticity(mesh=mesh, λ=λ, μ=μ)
    @info "$shape$order, element size $elem_size mm: $(mesh.nNodes) nodes, $(size(mesh.IEN, 2)) elements, $(Threads.nthreads()) threads"

    # 0.05 mm per step compresses the 40 mm cylinder by 12.5 % over 100 steps
    for (control, c) in (("displacement", -0.05), ("force", -25.0)), sensitivities in (false, true)
        simulate(mdl, fill(c, 3); control=control, β=1.0, sensitivities=sensitivities)   # compile, build the pattern
        t = minimum(@elapsed(simulate(mdl, fill(c, nsteps); control=control, β=1.0, sensitivities=sensitivities)) for _ in 1:3)
        @info "simulate, $control control, sensitivities=$sensitivities: $(round(1e3t / nsteps, digits=2)) ms/step"
    end

    # observed lateral surface pushed 0.5 mm outward; slip on the table pins the rotation about z
    reset_model!(mdl)
    X = mesh.NodeList
    side = mesh.side_nodes
    normals = vcat(X[1:2, side] ./ sqrt.(sum(abs2, X[1:2, side], dims=1)), zeros(1, length(side)))
    d = fill(0.5, length(side))
    volume_update!(mdl, side, normals, d; β_table=1e-3)
    t = [(reset_model!(mdl); @elapsed volume_update!(mdl, side, normals, d; β_table=1e-3)) for _ in 1:10]
    @info "volume_update!, $(length(side)) constrained nodes: $(round(1e3median(t), digits=2)) ms (median of 10)"
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(Symbol(get(ARGS, 1, "Hex")), parse(Int, get(ARGS, 2, "2")), parse(Float64, get(ARGS, 3, "6.0")))
end
