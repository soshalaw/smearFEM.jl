using Gmsh

# Gmsh is not thread-safe: all API calls must be serialized via this lock.
const GMSH_LOCK = ReentrantLock()

# Serializes mesh file generation so concurrent threads don't race to run the
# same gmsh process on the same output file.
const GMSH_GENERATE_LOCK = ReentrantLock()

# Track whether Gmsh has been initialized. Repeated initialize/finalize cycles
# corrupt Gmsh's internal C++ heap, so we initialize once and use gmsh.clear()
# between reads instead.
const _GMSH_INITIALIZED = Ref{Bool}(false)

function _gmsh_initialize()
    if !_GMSH_INITIALIZED[]
        gmsh.initialize()
        _GMSH_INITIALIZED[] = true
    end
end

"""
    _get_elements_for_physical(name::String, dim::Int, tag_to_index::Dict{UInt64,Int}) -> Matrix{Int} or nothing

Return the element connectivity for the named physical group as a matrix of remapped
node indices, or `nothing` if the group is not found.
"""
function _get_elements_for_physical(name, dim, tag_to_index::Dict{UInt64,Int})
    phys = gmsh.model.getPhysicalGroups()

    for (d, tag) in phys
        if d == dim && gmsh.model.getPhysicalName(d, tag) == name
            entities = gmsh.model.getEntitiesForPhysicalGroup(d, tag)

            all_elems = []

            for ent in entities
                etypes, etags, enodes = gmsh.model.mesh.getElements(d, ent)

                for i in eachindex(etypes)
                    nodes = enodes[i]
                    num_elem = length(etags[i])
                    num_nodes_per_elem = length(nodes) ÷ num_elem

                    elems = reshape(nodes, num_nodes_per_elem, num_elem)'
                    elems = [tag_to_index[n] for n in elems]

                    push!(all_elems, elems)
                end
            end

            isempty(all_elems) && return nothing

            col_counts = [size(e, 2) for e in all_elems]
            if length(unique(col_counts)) > 1
                throw(ArgumentError(
                    "Physical group \"$name\" contains mixed element types " *
                    "(nodes/element: $(sort(unique(col_counts)))). " *
                    "This usually means the mesh has both hex and wedge elements due to " *
                    "incomplete quad recombination on curved surfaces. " *
                    "Delete the cached .msh file and re-run, or use element_shape=:Tet " *
                    "for geometries with rounded edges."
                ))
            end

            return vcat(all_elems...)
        end
    end

    return nothing
end

"""
    _get_nodes_from_physical(name::String, dim::Int, tag_to_index::Dict{UInt64,Int}) -> Vector{Int} or nothing

Return sorted remapped node indices belonging to the named physical group (interior
nodes only), or `nothing` if the group is not found.
"""
function _get_nodes_from_physical(name, dim, tag_to_index::Dict{UInt64,Int})
    phys = gmsh.model.getPhysicalGroups()

    for (d, tag) in phys
        if d == dim && gmsh.model.getPhysicalName(d, tag) == name
            entities = gmsh.model.getEntitiesForPhysicalGroup(d, tag)

            node_set = Set{Int}()

            for ent in entities
                nodeTags, _, _ = gmsh.model.mesh.getNodes(d, ent)

                for n in nodeTags
                    push!(node_set, tag_to_index[n])
                end
            end

            return sort(collect(node_set))
        end
    end

    return nothing
end

"""
    _get_nodes_from_physical_with_boundary(name::String, dim::Int, tag_to_index::Dict{UInt64,Int}) -> Vector{Int} or nothing

Like `_get_nodes_from_physical` but also includes nodes on the lower-dimensional
boundary entities of the group, ensuring no boundary nodes are omitted.
"""
function _get_nodes_from_physical_with_boundary(name::String, dim::Int, tag_to_index::Dict{UInt64,Int})
    phys = gmsh.model.getPhysicalGroups()
    
    for (d, tag) in phys
        if d == dim && gmsh.model.getPhysicalName(d, tag) == name
            entities = gmsh.model.getEntitiesForPhysicalGroup(d, tag)
            node_set = Set{Int}()

            for ent in entities
                # nodes on the surface
                nodeTags, _, _ = gmsh.model.mesh.getNodes(d, ent)
                for n in nodeTags
                    push!(node_set, tag_to_index[n])
                end

                # include boundary nodes
                boundaries = gmsh.model.getBoundary([(d, ent)])  # no 'oriented' keyword in Julia
                for (bdim, btag) in boundaries
                    bnodes, _, _ = gmsh.model.mesh.getNodes(bdim, btag)
                    for n in bnodes
                        push!(node_set, tag_to_index[n])
                    end
                end
            end

            return sort(collect(node_set))
        end
    end

    return nothing
end

"""
    _get_mesh_data(filePath; params, template_path, mesh_order, mesh_dim,
                  body_group, body_dim, face_groups, node_set_groups) ->
        (NodeList, IEN, face_IENs, nNodes, node_sets, ne)

Load a Gmsh mesh file and extract the connectivity arrays required by smearFEM.
Auto-generates the mesh if the `.msh` file is absent and `params` + `template_path`
are supplied.

**Thread-safety**: all Gmsh API calls are serialized via `GMSH_LOCK`.

# Keyword Arguments
- `params::Union{Dict{String,<:Any},Nothing}=nothing`: substitution dict for `_generate_mesh_geo` (required for auto-generation).
- `template_path::Union{String,Nothing}=nothing`: path to the `.geo` template (required for auto-generation).
- `mesh_order::Int=2`: element order passed to gmsh.
- `mesh_dim::Int=3`: meshing dimension passed to gmsh CLI (1, 2, or 3).
- `body_group::String="Volume"`: name of the body physical group.
- `body_dim::Int=3`: spatial dimension of the body group (3=volume, 2=surface, 1=curve).
- `face_groups::Vector{String}=["Top","Bottom","Lateral"]`: face physical group names (at `body_dim-1`).
- `node_set_groups::Vector{String}=["Top","Bottom","Lateral"]`: physical group names to extract node index sets from.

# Returns
- `NodeList::Matrix{Float64}`: 3×nNodes coordinate matrix.
- `IEN::Matrix{Int}`: nNodesPerElem×nElem body connectivity (remapped to smearFEM convention).
- `face_IENs::Vector{Union{Matrix{Int},Nothing}}`: one nNodesPerElem×nElem IEN per `face_groups` entry.
- `nNodes::Int`: total number of nodes.
- `node_sets::Vector{Union{Vector{Int},Nothing}}`: sorted node index vectors, one per `node_set_groups` entry.
- `ne::Int`: number of body elements.
"""
function _get_mesh_data(filePath::String;
                       params::Union{Dict{String,<:Any},Nothing}=nothing,
                       template_path::Union{String,Nothing}=nothing,
                       mesh_order::Int=2,
                       mesh_dim::Int=3,
                       body_group::String="Volume",
                       body_dim::Int=3,
                       face_groups::Vector{String}=["Top", "Bottom", "Lateral"],
                       node_set_groups::Vector{String}=["Top", "Bottom", "Lateral"])
    if !isfile(filePath)
        lock(GMSH_GENERATE_LOCK) do
            if !isfile(filePath)
                if isnothing(params) || isnothing(template_path)
                    error("Mesh '$filePath' not found. Provide `params` and `template_path` for auto-generation.")
                end
                geo_path = splitext(filePath)[1] * ".geo"
                @debug "Generated mesh file path: $geo_path"
                mkpath(dirname(filePath))
                _generate_mesh_geo(geo_path, template_path, params)
                _run_gmsh(geo_path, filePath, mesh_order; dim=mesh_dim) ||
                    error("gmsh failed to generate mesh at '$filePath'.")
            end
        end
    end

    lock(GMSH_LOCK) do
        _gmsh_initialize()
        gmsh.open(filePath)

        nodeTags, nodeCoords, _ = gmsh.model.mesh.getNodes()
        NodeList = reshape(nodeCoords, 3, :)
        tag_to_index = Dict(tag => i for (i, tag) in enumerate(nodeTags))
        nNodes = size(NodeList, 2)
        @debug "_get_mesh_data: NodeList size=$(size(NodeList))"

        # Body IEN
        IEN_body = vcat(_get_elements_for_physical(body_group, body_dim, tag_to_index))
        @debug "_get_mesh_data: IEN_body size=$(size(IEN_body))"
        ne = size(IEN_body, 1)

        # Face IENs (one per face_groups entry)
        face_dim = body_dim - 1
        face_IENs_raw = [_get_elements_for_physical(g, face_dim, tag_to_index) for g in face_groups]

        # Node sets
        ns_dim = body_dim - 1
        node_sets = [_get_nodes_from_physical_with_boundary(g, ns_dim, tag_to_index) for g in node_set_groups]

        # Remap body IEN node ordering to smearFEM convention
        npe = size(IEN_body, 2)
        body_map = if npe == 27        # Q2 hex
            [1, 2, 3, 4, 5, 6, 7, 8, 9, 12, 14, 10, 17, 19, 20, 18, 11, 13, 15, 16, 22, 24, 25, 23, 21, 26, 27]
        elseif npe == 10               # T2 tet
            [1, 2, 3, 4, 5, 6, 7, 8, 10, 9]
        else                           # identity (Q1 hex, T1 tet, quads, tris, lines)
            collect(1:npe)
        end
        IEN = permutedims(IEN_body[:, body_map])

        face_IENs = [isnothing(f) ? nothing : permutedims(f) for f in face_IENs_raw]

        gmsh.clear()
        return NodeList, IEN, face_IENs, nNodes, node_sets, ne
    end
end

"""
    _generate_mesh_geo(output_path, template_path, params) -> Bool

Generate a `.geo` file by substituting parameters into a template.
Each entry in `params` replaces the pattern `key = <anything>;` in the template.
Skips writing when the output file already contains identical content.

# Arguments
- `output_path::String`: destination path for the generated `.geo` file.
- `template_path::String`: path to the source `.geo` template file.
- `params::Dict{String,<:Any}`: map from Gmsh variable name to new value.
"""
function _generate_mesh_geo(output_path::String, template_path::String, params::Dict{String,<:Any})
    geo_content = read(template_path, String)
    for (key, val) in params
        geo_content = replace(geo_content, Regex("$(key)\\s*=\\s*[^;]+;") => "$key = $val;")
    end
    if isfile(output_path) && read(output_path, String) == geo_content
        @info "$(basename(output_path)) already up to date, skipping"
        return false
    end
    write(output_path, geo_content)
    @info "Generated $(basename(output_path)) from $(basename(template_path))"
    return true
end

"""
    _run_gmsh(geo_path, msh_path, mesh_order=2; dim=3) -> Bool

Invoke the `gmsh` CLI to mesh `geo_path` and write output to `msh_path`.
`dim` controls the meshing dimension (1, 2, or 3).
Returns `false` if gmsh is not on PATH or if the command fails.
"""
function _run_gmsh(geo_path::String, msh_path::String, mesh_order::Int=2; dim::Int=3)
    gmsh_path = Sys.which("gmsh")
    if gmsh_path !== nothing
        cmd = `$gmsh_path $geo_path -$dim -format msh -order $mesh_order -o $msh_path`
        @info "Running gmsh CLI: $cmd"
        try
            run(cmd)
            @info "Mesh generated: $msh_path"
            return true
        catch e
            @warn "gmsh CLI failed: $e. Falling back to Julia API."
        end
    else
        @info "gmsh CLI not found in PATH; using Julia API for mesh generation."
    end
    try
        lock(GMSH_LOCK) do
            _gmsh_initialize()
            gmsh.open(geo_path)
            gmsh.model.mesh.generate(dim)
            gmsh.model.mesh.setOrder(mesh_order)
            gmsh.write(msh_path)
            gmsh.clear()
        end
        @info "Mesh generated via Julia API: $msh_path"
        return true
    catch e
        @error "Failed to generate mesh via Julia API: $e"
        return false
    end
end
