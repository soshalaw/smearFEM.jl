using LinearAlgebra: norm
using ImageMorphology: label_components, component_lengths
using MarchingCubes: MC, march

# Paper-2 initial mesh: the FE mesh at t = 0, built from perception's fused field φ_s (a truncated
# signed distance, NaN where no camera saw the voxel). The table closes the unobserved bottom; any
# other unobserved patch is an error, not a guess.

"""
    mesh_from_field(phi, grid; table_z=0.0, elem_size, basis_order=2, out_dir=mktempdir()) -> MeshgridObserved
    mesh_from_field(step::PerceptionStep, grid; kwargs...) -> MeshgridObserved

Paper-2 initial mesh: tet mesh of the specimen perception observed, from its fused field φ_s. The
zero surface of φ_s is closed at the bottom by the table plane z = `table_z`; voxels no camera saw are
inside where the observed surface encloses them, else outside. Pieces of material apart from the
specimen are dropped, and pockets of outside enclosed by material (noise) filled, with a warning. The cap is about one voxel smaller in radius than the true
base: marching cubes bevels the bottom edge, and the bevel belongs to the free surface.

# Arguments
- `phi::AbstractArray{<:Real,3}`: φ_s (nx, ny, nz) [mm], NaN where unobserved, e.g. `step.phi`.
- `grid::PerceptionGrid`: Its grid.
- `table_z::Real`: Height of the table top [mm]: 0 on the rig, the plate top in synthetic bags.
- `elem_size::Real`: Target element size [mm].
- `basis_order::Int`: 1 or 2 (quadratic tets keep straight edges).
- `out_dir::String`: Folder for `surface.msh` (the closed surface, to inspect), `observed.geo` and
  `observed.msh`.

# Returns
- `MeshgridObserved`: Boundary groups `IEN_bottom` (on the table) and `IEN_free`.

Throws an `ArgumentError` for `table_z` outside the grid, no material, or material that does not
reach the table, and an error when the observed surface is not closed (an unobserved patch).
"""
function mesh_from_field(phi::AbstractArray{<:Real,3}, grid::PerceptionGrid; table_z::Real=0.0,
                         elem_size::Real, basis_order::Int=2, out_dir::String=mktempdir())
    elem_size > 0 || throw(ArgumentError("elem_size must be positive, got $elem_size"))
    basis_order in (1, 2) || throw(ArgumentError("basis_order must be 1 or 2, got $basis_order"))
    field, coords = _close_field(phi, grid, table_z)
    V, T, on_cap = _field_surface(field, coords, table_z, grid.delta)
    mkpath(out_dir)
    _write_surface_msh(joinpath(out_dir, "surface.msh"), V, T, on_cap)
    mesh_path = joinpath(out_dir, "observed.msh")
    rm(mesh_path; force=true)   # _get_mesh_data meshes only when the file is missing
    template = joinpath(dirname(dirname(dirname(@__DIR__))), "mesh_files", "templates", "observed_tet.geo")
    NodeList, IEN, face_IENs, nNodes, node_sets, ne = _get_mesh_data(mesh_path;
        params=Dict{String,Any}("elem_size" => Float64(elem_size)), template_path=template,
        mesh_order=basis_order, mesh_dim=3, body_group="Volume", body_dim=3,
        face_groups=["Bottom", "Free"], node_set_groups=["Bottom", "Free"])
    return MeshgridObserved(NodeList=NodeList, IEN=IEN, IEN_bottom=face_IENs[1], IEN_free=face_IENs[2],
        ID=reshape(collect(1:3nNodes), 3, nNodes), basis_order=basis_order, nNodes=nNodes, ne=ne,
        bottom_nodes=node_sets[1], free_nodes=node_sets[2], effect_elem_sze=Float64(elem_size))
end

mesh_from_field(step::PerceptionStep, grid::PerceptionGrid; kwargs...) = mesh_from_field(step.phi, grid; kwargs...)

"""
    _close_field(phi, grid, table_z) -> (field, coords)

φ_s made meshable: padded by one voxel on every side, the table plane as the bottom, every NaN voxel
decided inside or outside, one piece of material without enclosed pockets of outside, |value| ≥ 1e-3·Δ.
`field` is negative inside [mm];
`coords` are the padded grid's voxel centres along each axis [mm].
"""
function _close_field(phi::AbstractArray{<:Real,3}, grid::PerceptionGrid, table_z::Real)
    Δ = grid.delta
    xs, ys, zs = _grid_centres(grid)
    zs[1] <= table_z < zs[end] ||
        throw(ArgumentError("table_z = $table_z mm is outside the grid's z range [$(zs[1]), $(zs[end])) mm"))
    px, py, pz = (vcat(c[1] - Δ, collect(c), c[end] + Δ) for c in (xs, ys, zs))
    F = fill(NaN, length(px), length(py), length(pz))
    F[2:end-1, 2:end-1, 2:end-1] .= phi
    plane = reshape(table_z .- pz, 1, 1, :)     # signed distance to the table, positive below

    # Seal: at and below the table is outside. This also cuts the unobserved region under the
    # material off from the one below the table, through which the interior would reach the border.
    F .= ifelse.(plane .>= 0, plane, F)

    # An unobserved region reaching the border is outside; one enclosed by observed surface is the
    # material's interior, deeper than perception's distance cap.
    unobserved = isnan.(F)
    labels = label_components(unobserved)
    outside_labels = _border_labels(labels)
    enclosed = unobserved .& .!in.(labels, Ref(outside_labels))
    n_inside = count(enclosed)
    F[enclosed] .= -Δ * norm(size(F))           # below every distance inside the grid
    filled_outside = unobserved .& .!enclosed
    F[filled_outside] .= Δ
    F .= max.(F, plane)                         # the cap lies exactly on the table

    pieces = label_components(F .< 0)
    n = maximum(pieces)
    n == 0 && throw(ArgumentError("no material above the table: φ_s has no negative value"))
    if n > 1
        lengths = component_lengths(pieces)
        keep = argmax(l -> lengths[l], 1:n)
        @warn "dropped $(n - 1) spurious piece(s) of material apart from the specimen " *
              "($(sum(lengths[1:n]) - lengths[keep]) voxels)"
        F[(pieces .> 0) .& (pieces .!= keep)] .= Δ
    end

    # outside voxels the material encloses have no way out to the air: noise in an observed surface
    pockets = label_components(F .>= 0)
    inner = setdiff(unique(pockets), 0, _border_labels(pockets))
    if !isempty(inner)
        @warn "filled $(length(inner)) pocket(s) of outside enclosed by material as inside (noise)"
        F[in.(pockets, Ref(inner))] .= -Δ
    end

    ringed = _ringed_slices(F, filled_outside, pz, table_z)
    isempty(ringed) || error("observed surface is not closed: the material's interior is open to the " *
        "outside (ring-shaped cross-sections at z = $(round(first(ringed), digits=1)) … " *
        "$(round(last(ringed), digits=1)) mm); an unobserved patch of the surface, or table_z too low")

    inner = view(F, 2:size(F, 1)-1, 2:size(F, 2)-1, 2:size(F, 3)-1)
    outer_layers = (selectdim(inner, 1, 1), selectdim(inner, 1, size(inner, 1)),
                    selectdim(inner, 2, 1), selectdim(inner, 2, size(inner, 2)),
                    selectdim(inner, 3, size(inner, 3)))
    any(any(<(0), layer) for layer in outer_layers) &&
        @warn "the material reaches the grid boundary: the mesh is cut flat there"

    # a voxel centre on the surface, or within float noise of it, makes marching cubes emit needles
    # too thin for Gmsh's parametrization (quality 6e-15 ran Gmsh out of memory)
    ε = 1e-3 * Δ
    F .= ifelse.(abs.(F) .< ε, ifelse.(F .< 0, -ε, ε), F)   # -0.0 is outside, as in the labelling above
    @info "closed φ_s: $n_inside unobserved voxels filled as inside, $(count(unobserved) - n_inside) as outside"
    return F, (px, py, pz)
end

"""
    _field_surface(field, coords, table_z, delta) -> (V, T, on_cap)

Zero surface of a closed field (`_close_field`) by marching cubes: vertices `V` (3, nv) [mm],
triangles `T` (3, nt) as vertex indices, and `on_cap` (nt), the triangles on the table that close
the bottom.
"""
function _field_surface(field::Array{Float64,3}, coords::NTuple{3,Vector{Float64}}, table_z::Real, delta::Real)
    mc = MC(field, Int; x=coords[1], y=coords[2], z=coords[3])
    march(mc, 0.0)
    V = stack(mc.vertices)
    T = stack(mc.triangles)
    # exactly on the plane: Float32 rounding of φ_s moves the cap, and `_close_field`'s |φ| ≥ 1e-3·Δ
    # lifts a cap on a voxel plane by up to 1e-3·Δ
    on_table = abs.(V[3, :] .- table_z) .< 1e-2 * delta
    V[3, on_table] .= table_z
    on_cap = vec(all(on_table[T], dims=1))
    any(on_cap) || throw(ArgumentError("the material does not reach the table at z = $table_z mm; check table_z"))
    return V, T, on_cap
end

# The closed surface for Gmsh as MSH 2.2 with shared nodes, so the cap (entity 1) and the free
# surface (entity 2) meet along one curve; STL repeats the vertices of every triangle.
function _write_surface_msh(path::String, V::Matrix{Float64}, T::Matrix{Int}, on_cap::BitVector)
    open(path, "w") do io
        println(io, "\$MeshFormat\n2.2 0 8\n\$EndMeshFormat\n\$Nodes\n", size(V, 2))
        for a in axes(V, 2)
            println(io, a, " ", V[1, a], " ", V[2, a], " ", V[3, a])
        end
        println(io, "\$EndNodes\n\$Elements\n", size(T, 2))
        for t in axes(T, 2)
            entity = on_cap[t] ? 1 : 2
            println(io, t, " 2 2 ", entity, " ", entity, " ", T[1, t], " ", T[2, t], " ", T[3, t])
        end
        println(io, "\$EndElements")
    end
end

# Labels (other than 0) of the components touching the array's border.
function _border_labels(labels::AbstractArray{<:Integer})
    faces = (selectdim(labels, d, i) for d in 1:ndims(labels) for i in (1, size(labels, d)))
    return setdiff(union((unique(f) for f in faces)...), 0)
end

# z of each slice above the table where the material encloses unobserved voxels filled as outside: a
# ring around a leaked interior. Rings of observed voxels are dimples of noise under a flat top.
function _ringed_slices(F::Array{Float64,3}, filled_outside::BitArray{3}, pz::Vector{Float64}, table_z::Real)
    ringed = Float64[]
    for k in axes(F, 3)
        pz[k] > table_z || continue
        labels = label_components(view(F, :, :, k) .>= 0)
        inner = setdiff(unique(labels), 0, _border_labels(labels))
        any(in(inner), labels[view(filled_outside, :, :, k)]) && push!(ringed, pz[k])
    end
    return ringed
end
