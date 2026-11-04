"""
    gt_comparison.jl

Paper-2 initial mesh against the ground truth of a synthetic run: the FE mesh `mesh_from_field`
builds from φ_s of one step of a reconstruction, its boundary written as an OBJ for smearPerception's
`tools/analysis/compare_mesh_gt.py`, which scores it against the run's GT mesh.

    systemd-run --user --scope -q -p MemoryMax=8G -p MemorySwapMax=0 \
        julia --project experiments/mesh_init/gt_comparison.jl <recon.h5> <out_dir> [elem_size] [table_z]

(the scope caps the memory: Gmsh has run away on bad input and taken the machine down.)
`recon.h5` comes from smearPerception's `reconstruct_bag.py <bag> <recon.h5> --gt <sim dir>
--save-field`. The table is where perception put it (the file's `table_z`: the bottom plate's top on
synthetic runs, 0 on the rig), unless `table_z` [mm] is given; files written before perception recorded
it need it. Writes `surface.msh`, `observed.msh` and `observed_<step>.obj` to `out_dir`.
"""

using smearFEM
using HDF5: h5open, attributes
using LinearAlgebra: det

"""
    export_boundary_obj(mesh::MeshgridObserved, path) -> path

Boundary of `mesh` as an OBJ [mm]: the corner nodes of its table and free faces (P2 edges are
straight, so the mid-edge nodes add nothing to the shape).
"""
function export_boundary_obj(mesh::MeshgridObserved, path::AbstractString)
    faces = hcat(mesh.IEN_bottom[1:3, :], mesh.IEN_free[1:3, :])
    nodes = unique(vec(faces))
    index = Dict(a => i for (i, a) in enumerate(nodes))
    open(path, "w") do io
        for a in nodes
            println(io, "v ", mesh.NodeList[1, a], " ", mesh.NodeList[2, a], " ", mesh.NodeList[3, a])
        end
        for f in eachcol(faces)
            println(io, "f ", index[f[1]], " ", index[f[2]], " ", index[f[3]])
        end
    end
    return path
end

function main(h5::AbstractString, out_dir::AbstractString; step::Int=0, table_z=nothing,
              elem_size::Real=3.0, basis_order::Int=2)
    data = read_perception_steps(h5)
    k = findfirst(s -> s.n == step, data.steps)
    isnothing(k) && error("$h5 has no step $step")
    if isnothing(table_z)
        table_z = h5open(h5) do f
            haskey(attributes(f), "table_z") ||
                error("$h5 has no table_z attribute (written before perception recorded it): " *
                      "pass table_z [mm]")
            1e3 * read(attributes(f)["table_z"])   # m → mm
        end
    end
    mesh = mesh_from_field(data.steps[k], data.grid; table_z, elem_size, basis_order, out_dir)
    X = mesh.NodeList
    volume = sum(det(X[:, e[2:4]] .- X[:, e[1]]) for e in eachcol(mesh.IEN)) / 6
    obj = export_boundary_obj(mesh, joinpath(out_dir, "observed_$(lpad(step, 3, '0')).obj"))
    @info "step $step, table at z = $(round(table_z, digits=3)) mm: " *
          "Tet$basis_order mesh, $(mesh.nNodes) nodes, " *
          "$(mesh.ne) tets, volume $(round(volume, digits=1)) mm³"
    @info "boundary: $obj; compare it with the GT in smearPerception (env smear_ros):\n" *
          "  python tools/analysis/compare_mesh_gt.py $obj <sim data dir> --step $step"
    return mesh
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS[1], ARGS[2]; elem_size=parse(Float64, get(ARGS, 3, "3.0")),
         table_z=length(ARGS) >= 4 ? parse(Float64, ARGS[4]) : nothing)
end
