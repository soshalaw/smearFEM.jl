using LinearAlgebra: norm
using ImageMorphology: label_components, component_lengths

# Paper-2 initial mesh: the FE mesh at t = 0, built from perception's fused field φ_s (a truncated
# signed distance, NaN where no camera saw the voxel). The table closes the unobserved bottom; any
# other unobserved patch is an error, not a guess.

"""
    _close_field(phi, grid, table_z) -> (field, coords)

φ_s made meshable: padded by one voxel on every side, the table plane as the bottom, every NaN voxel
decided inside or outside, one piece of material, no exact zeros. `field` is negative inside [mm];
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
    F[unobserved .& .!enclosed] .= Δ
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

    ringed = _ringed_slices(F, pz, table_z)
    isempty(ringed) || error("observed surface is not closed: the material's interior is open to the " *
        "outside (ring-shaped cross-sections at z = $(round(first(ringed), digits=1)) … " *
        "$(round(last(ringed), digits=1)) mm); an unobserved patch of the surface, or table_z too low")

    inner = view(F, 2:size(F, 1)-1, 2:size(F, 2)-1, 2:size(F, 3)-1)
    outer_layers = (selectdim(inner, 1, 1), selectdim(inner, 1, size(inner, 1)),
                    selectdim(inner, 2, 1), selectdim(inner, 2, size(inner, 2)),
                    selectdim(inner, 3, size(inner, 3)))
    any(any(<(0), layer) for layer in outer_layers) &&
        @warn "the material reaches the grid boundary: the mesh is cut flat there"

    # a voxel centre exactly on the surface makes marching cubes emit zero-area triangles
    F[F .== 0] .= 1e-6 * Δ
    @info "closed φ_s: $n_inside unobserved voxels filled as inside, $(count(unobserved) - n_inside) as outside"
    return F, (px, py, pz)
end

# Labels (other than 0) of the components touching the array's border.
function _border_labels(labels::AbstractArray{<:Integer})
    faces = (selectdim(labels, d, i) for d in 1:ndims(labels) for i in (1, size(labels, d)))
    return setdiff(union((unique(f) for f in faces)...), 0)
end

# z of each slice above the table where the material encloses outside voxels: a ring around a
# leaked interior, since horizontal cross-sections of a lump of material have no holes.
function _ringed_slices(F::Array{Float64,3}, pz::Vector{Float64}, table_z::Real)
    ringed = Float64[]
    for k in axes(F, 3)
        pz[k] > table_z || continue
        labels = label_components(view(F, :, :, k) .>= 0)
        isempty(setdiff(unique(labels), 0, _border_labels(labels))) || push!(ringed, pz[k])
    end
    return ringed
end
