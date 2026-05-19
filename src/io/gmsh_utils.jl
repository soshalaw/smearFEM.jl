using Gmsh

# Gmsh is not thread-safe: all API calls must be serialized via this lock.
const GMSH_LOCK = ReentrantLock()

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
    get_mesh_data(filePath; radius, height, elem_size, template_path, mesh_order) ->
        (NodeList, IEN, IEN_top, IEN_bottom, IEN_lateral, nNodes, node_sets, ne)

Load a Gmsh mesh file and extract the connectivity arrays required by smearFEM.
If the `.msh` file does not exist and `radius`, `height`, `elem_size`, and
`template_path` are all supplied, the mesh is generated automatically by calling
`generate_mesh_geo` + `run_gmsh` before loading.

**Thread-safety**: Gmsh's C API is not re-entrant. All Gmsh calls inside this
function are serialized through the module-level `GMSH_LOCK` (`ReentrantLock`),
making this function safe to call from multiple Julia threads concurrently.

# Arguments
- `filePath::String`: path to the `.msh` file.

# Keyword Arguments
- `radius::Union{Float64,Nothing}=nothing`: cylinder radius (required for auto-generation).
- `height::Union{Float64,Nothing}=nothing`: cylinder height (required for auto-generation).
- `elem_size::Union{Float64,Nothing}=nothing`: target element size (required for auto-generation).
- `template_path::Union{String,Nothing}=nothing`: path to the `.geo` template (required for auto-generation).
- `mesh_order::Int=2`: element order passed to gmsh (1 or 2).

# Returns
- `NodeList::Matrix{Float64}`: 3×nNodes coordinate matrix.
- `IEN::Matrix{Int}`: nNodesPerElem×nElem volume connectivity (Q1 or Q2 remapped).
- `IEN_top::Matrix{Int}`: connectivity for the top-face physical group.
- `IEN_bottom::Matrix{Int}`: connectivity for the bottom-face physical group.
- `IEN_lateral::Matrix{Int}`: raw (unpermuted) lateral-face connectivity.
- `nNodes::Int`: total number of nodes.
- `node_sets::Vector{Union{Vector{Int},Nothing}}`: `[nodes_lateral, nodes_bottom, nodes_top]` sorted index vectors.
- `ne::Int`: number of volume elements.
"""
function get_mesh_data(filePath::String;
                       radius::Union{Float64,Nothing}=nothing,
                       height::Union{Float64,Nothing}=nothing,
                       elem_size::Union{Float64,Nothing}=nothing,
                       template_path::Union{String,Nothing}=nothing,
                       mesh_order::Int=2)
    if !isfile(filePath)
        if any(isnothing, (radius, height, elem_size, template_path))
            error("Mesh file not found at '$filePath' and generation parameters " *
                  "(radius, height, elem_size, template_path) are not all provided.")
        end
        geo_path = splitext(filePath)[1] * ".geo"
        mkpath(dirname(filePath))
        generate_mesh_geo(something(radius), something(height), something(elem_size), geo_path, something(template_path))
        run_gmsh(geo_path, filePath, mesh_order) ||
            error("gmsh failed to generate mesh at '$filePath'.")
    end

    # All Gmsh operations must be serialized via lock due to thread-safety limitations
    lock(GMSH_LOCK) do
        gmsh.initialize()
        gmsh.open(filePath)

        # --------------------------------------------------
        # 1. NodeList
        # --------------------------------------------------
        nodeTags, nodeCoords, _ = gmsh.model.mesh.getNodes()

        NodeList = reshape(nodeCoords, 3, :)  # 3XN

        # Map: node tag → index
        tag_to_index = Dict(tag => i for (i, tag) in enumerate(nodeTags))

        nNodes = size(NodeList, 2) # number of nodes
        @debug "get_mesh_data: NodeList size=$(size(NodeList))"

        # --------------------------------------------------
        # 2. IEN arrays
        # --------------------------------------------------
        IEN_volume = vcat(_get_elements_for_physical("Cylinder", 3, tag_to_index))
        IEN_top = vcat(_get_elements_for_physical("Top", 2, tag_to_index))
        IEN_bottom = vcat(_get_elements_for_physical("Bottom", 2, tag_to_index))
        IEN_lateral = vcat(_get_elements_for_physical("Lateral", 2, tag_to_index))

        @debug "get_mesh_data: IEN_volume size=$(size(IEN_volume))"
        ne = size(IEN_volume, 1) # number of elements

        # --------------------------------------------------
        # 3. Extract node sets
        # --------------------------------------------------
        nodes_top = _get_nodes_from_physical_with_boundary("Top", 2, tag_to_index)
        nodes_bottom = _get_nodes_from_physical_with_boundary("Bottom", 2, tag_to_index)
        nodes_lateral = _get_nodes_from_physical_with_boundary("Lateral", 2, tag_to_index)

        @debug "get_mesh_data: Top nodes=$(length(nodes_top)), Bottom nodes=$(length(nodes_bottom)), Lateral nodes=$(length(nodes_lateral))"

        # --------------------------------------------------
        # 4. Remap IEN arrays
        # --------------------------------------------------
        npe = size(IEN_volume, 2)
        if npe == 27        # Q2 hex
            map = [1, 2, 3, 4, 5, 6, 7, 8, 9, 12, 14, 10, 17, 19, 20, 18, 11, 13, 15, 16, 22, 24, 25, 23, 21, 26, 27]
        elseif npe == 8     # Q1 hex
            map = collect(1:8)
        elseif npe == 10    # T2 tet (quadratic)
            map = collect(1:10)
        elseif npe == 4     # T1 tet (linear)
            map = collect(1:4)
        else
            error("Unexpected number of nodes per element: ", npe)
        end
        IEN = IEN_volume[:,map]

        # Permute IEN arrays to have nodes as columns and elements as rows
        IEN_new = permutedims(IEN)
        IEN_bottom_new = permutedims(IEN_bottom)
        IEN_top_new = permutedims(IEN_top)

        gmsh.finalize()

        return NodeList, IEN_new, IEN_top_new, IEN_bottom_new, IEN_lateral, nNodes, [nodes_lateral, nodes_bottom, nodes_top], ne
    end
end

"""
    generate_mesh_geo(radius::Float64, height::Float64, elem_size::Float64,
                      output_path::String, template_path::String) -> Bool

Generate a `.geo` file from a template by substituting `radius`, `height`, and
`elem_size_2d`. Returns `false` (and skips writing) when an existing file at
`output_path` already contains the requested parameter values.

# Arguments
- `radius::Float64`: cylinder radius to substitute into the template.
- `height::Float64`: cylinder height to substitute into the template.
- `elem_size::Float64`: characteristic element size to substitute into the template.
- `output_path::String`: destination path for the generated `.geo` file.
- `template_path::String`: path to the source `.geo` template file.
"""
function generate_mesh_geo(radius::Float64, height::Float64, elem_size::Float64, output_path::String, template_path::String)
    if isfile(output_path)
        try
            content = read(output_path, String)
            if contains(content, "radius = $radius;") &&
               contains(content, "height = $height;") &&
               contains(content, "elem_size_2d = $elem_size;")
                @info "mesh.geo already exists with correct parameters (radius=$radius, height=$height, elem_size=$elem_size), skipping generation"
                return false
            end
        catch e
            @debug "Error reading existing mesh.geo, will regenerate: $e"
        end
    end

    template_content = read(template_path, String)
    geo_content = replace(template_content,
        "radius = 25.0;"     => "radius = $radius;",
        "height = 40.0;"     => "height = $height;",
        "elem_size_2d = 10;" => "elem_size_2d = $elem_size;")
    write(output_path, geo_content)
    @info "Generated mesh.geo at $output_path from template $template_path"
    return true
end

"""
    run_gmsh(geo_path, msh_path, mesh_order=2) -> Bool

Invoke the `gmsh` CLI to mesh `geo_path` and write output to `msh_path`.
Returns `false` if gmsh is not on PATH or if the command fails.
"""
function run_gmsh(geo_path::String, msh_path::String, mesh_order::Int=2)
    gmsh_path = Sys.which("gmsh")
    if gmsh_path === nothing
        @warn "gmsh executable not found in PATH. Please install gmsh or add it to PATH."
        return false
    end

    cmd = `$gmsh_path $geo_path -3 -format msh -order $mesh_order -o $msh_path`
    @info "Running gmsh: $cmd"
    try
        run(cmd)
        @info "Mesh generated: $msh_path"
        return true
    catch e
        @error "Failed to run gmsh: $e"
        return false
    end
end
