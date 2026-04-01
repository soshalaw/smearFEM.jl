using Gmsh

# Global lock for thread-safe Gmsh operations
# Gmsh is not thread-safe, so we serialize all access via this lock
const GMSH_LOCK = ReentrantLock()

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

function get_mesh_data(filePath::String)
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
        ndim = size(NodeList, 1) # number of dimensions (should be 3 for 3D meshes)
        println("NodeList size: ", size(NodeList))

        # --------------------------------------------------
        # 2. IEN arrays
        # --------------------------------------------------
        IEN_volume = vcat(_get_elements_for_physical("Cylinder", 3, tag_to_index))
        IEN_top = vcat(_get_elements_for_physical("Top", 2, tag_to_index))
        IEN_bottom = vcat(_get_elements_for_physical("Bottom", 2, tag_to_index))
        IEN_lateral = vcat(_get_elements_for_physical("Lateral", 2, tag_to_index))

        println("IEN_volume size: ", size(IEN_volume))
        ne = size(IEN_volume, 1) # number of elements

        # --------------------------------------------------
        # 2. Extract node sets
        # --------------------------------------------------
        nodes_top = _get_nodes_from_physical_with_boundary("Top", 2, tag_to_index)
        nodes_bottom = _get_nodes_from_physical_with_boundary("Bottom", 2, tag_to_index)
        nodes_lateral = _get_nodes_from_physical_with_boundary("Lateral", 2, tag_to_index)

        println("Top nodes: ", length(nodes_top))
        println("Bottom nodes: ", length(nodes_bottom))
        println("Lateral nodes (including boundaries): ", length(nodes_lateral))


        # --------------------------------------------------
        # 2. Remap IEN Arrays
        # --------------------------------------------------
        if size(IEN_volume, 2) == 27
            map = [1, 2, 3, 4, 5, 6, 7, 8, 9, 12, 14, 10, 17, 19, 20, 18, 11, 13, 15, 16, 22, 24, 25, 23, 21, 26, 27]
        elseif size(IEN_volume, 2) == 8
             map = [1, 2, 3, 4, 5, 6, 7, 8]
        else
            error("Unexpected number of nodes per element: ", size(IEN_volume, 2))
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
