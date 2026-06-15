include("mesh_types.jl")
include("mesh_structured.jl")
include("mesh_gmsh.jl")

"""
    meshgrid_cylinder(r, h; mesh_type=:structured, element_shape=:Hex, basis_order=1,
                      ne=nothing, ndof=3, elem_size=nothing, mesh_path=nothing, template_path=nothing)

Set up a cylinder mesh, dispatching to structured (Julia) or unstructured (Gmsh) generation.

# Arguments:
- `r::Number` : radius of the cylinder
- `h::Number` : height of the cylinder
- `mesh_type::Symbol` : `:structured` (Julia) or `:unstructured` (Gmsh)
- `element_shape::Symbol` : element shape (`:Hex` or `:Tet`)
- `basis_order::Int` : polynomial order (1 or 2)
- `ne::Int` : elements per direction (`:structured` only)
- `ndof::Int` : degrees of freedom per node
- `elem_size::Float64` : target element size (`:unstructured` only)
- `mesh_path::String` : base mesh directory (`:unstructured` only)
- `template_path::String` : override `.geo` template path

# Returns:
- `MeshgridCylinder` : cylinder mesh
"""
function meshgrid_cylinder(r::T, h::U;
    mesh_type::Symbol=:structured,
    element_shape::Symbol=:Hex,
    basis_order::Int=1,
    ne::Union{Int,Nothing}=nothing,
    ndof::Int=3,
    elem_size::Union{Float64,Nothing}=nothing,
    mesh_path::Union{String,Nothing}=nothing,
    template_path::Union{String,Nothing}=nothing
) where {T<:Number, U<:Number}
    if mesh_type == :structured
        isnothing(ne) && throw(ArgumentError("`ne` is required for mesh_type=:structured"))
        mesh = _meshgrid_cylinder(r, h, ne; element_shape=element_shape, basis_order=basis_order)
        if ndof != 3
            new_ID = zeros(Int64, ndof, mesh.nNodes)
            for m in 1:mesh.nNodes, l in 1:ndof
                new_ID[l, m] = ndof*(m-1) + l
            end
            mesh.ID = new_ID
        end
        return mesh
    elseif mesh_type == :unstructured
        isnothing(mesh_path) && throw(ArgumentError("`mesh_path` is required for mesh_type=:unstructured"))
        isnothing(elem_size) && throw(ArgumentError("`elem_size` is required for mesh_type=:unstructured"))
        return _get_gmsh_cylinder(mesh_path, ndof, r, h, element_shape, basis_order, Float64(elem_size); template_path=template_path)
    else
        throw(ArgumentError("Unknown mesh_type `$mesh_type`; use :structured or :unstructured"))
    end
end

"""
    meshgrid_cube(lx, ly, lz; mesh_type=:structured, element_shape=:Hex, basis_order=1,
                 ne=nothing, ndof=3, elem_size=nothing, mesh_path=nothing, template_path=nothing,
                 edge_radius=nothing)

Set up a box mesh, dispatching to structured (Julia) or unstructured (Gmsh) generation.
Note: `IEN_front/back/left/right` are not populated for structured meshes.

# Arguments:
- `lx::Number` : length in x
- `ly::Number` : length in y
- `lz::Number` : length in z
- `mesh_type::Symbol` : `:structured` (Julia) or `:unstructured` (Gmsh)
- `element_shape::Symbol` : element shape (`:Hex` or `:Tet`)
- `basis_order::Int` : polynomial order (1 or 2)
- `ne::Int` : elements per direction (`:structured` only)
- `ndof::Int` : degrees of freedom per node
- `elem_size::Float64` : target element size (`:unstructured` only)
- `mesh_path::String` : base mesh directory (`:unstructured` only)
- `template_path::String` : override `.geo` template path
- `edge_radius::Float64` : fillet radius for the 4 vertical edges (`:unstructured` only; must satisfy `0 < edge_radius < min(lx,ly)/2`)

# Returns:
- `MeshgridCube` : box mesh
"""
function meshgrid_cube(lx::T, ly::U, lz::V;
    mesh_type::Symbol=:structured,
    element_shape::Symbol=:Hex,
    basis_order::Int=1,
    ne::Union{Int,Nothing}=nothing,
    ndof::Int=3,
    elem_size::Union{Float64,Nothing}=nothing,
    mesh_path::Union{String,Nothing}=nothing,
    template_path::Union{String,Nothing}=nothing,
    edge_radius::Union{Float64,Nothing}=nothing
) where {T<:Number, U<:Number, V<:Number}
    if mesh_type == :structured
        !isnothing(edge_radius) && throw(ArgumentError("`edge_radius` is only supported for mesh_type=:unstructured"))
        isnothing(ne) && throw(ArgumentError("`ne` is required for mesh_type=:structured"))
        NodeList, IEN, _, IEN_top, IEN_bottom, _, _, BorderNodes = _meshgrid_cube(lx, ly, lz, ne; element_shape=element_shape, basis_order=basis_order)
        total_nodes = size(NodeList, 2)
        ID = zeros(Int64, ndof, total_nodes)
        for m in 1:total_nodes, l in 1:ndofShapes
            ID[l, m] = ndof*(m-1) + l
        end
        empty_mat = Matrix{Int}(undef, 0, 0)
        return MeshgridCube(lx=lx, ly=ly, lz=lz, NodeList=NodeList, IEN=IEN,
            IEN_top=IEN_top, IEN_bottom=IEN_bottom,
            IEN_front=empty_mat, IEN_back=empty_mat, IEN_left=empty_mat, IEN_right=empty_mat,
            ID=ID, volume_element_shape=element_shape, basis_order=basis_order,
            nNodes=total_nodes, ne=ne^3,
            side_nodes=BorderNodes[1], bottom_nodes=BorderNodes[2], top_nodes=BorderNodes[3])
    elseif mesh_type == :unstructured
        isnothing(mesh_path) && throw(ArgumentError("`mesh_path` is required for mesh_type=:unstructured"))
        isnothing(elem_size) && throw(ArgumentError("`elem_size` is required for mesh_type=:unstructured"))
        return _get_gmsh_box(mesh_path, ndof, lx, ly, lz, element_shape, basis_order, Float64(elem_size);
                             template_path=template_path, edge_radius=edge_radius)
    else
        throw(ArgumentError("Unknown mesh_type `$mesh_type`; use :structured or :unstructured"))
    end
end

"""
    meshgrid_square(lx, ly; mesh_type=:structured, element_shape=:Quad, basis_order=1,
                    ne=nothing, ndof=2, elem_size=nothing, mesh_path=nothing, template_path=nothing)

Set up a square mesh, dispatching to structured (Julia) or unstructured (Gmsh) generation.
For the raw-array form use `_meshgrid_square(lx, ly, ne::Int; ...)`.
Note: `top_nodes`, `bottom_nodes`, and `IEN_sides` are not populated for structured meshes.

# Arguments:
- `lx::Number` : length in x
- `ly::Number` : length in y
- `mesh_type::Symbol` : `:structured` (Julia) or `:unstructured` (Gmsh)
- `element_shape::Symbol` : element shape (`:Quad` or `:Tri`)
- `basis_order::Int` : polynomial order (1 or 2)
- `ne::Int` : elements per direction (`:structured` only)
- `ndof::Int` : degrees of freedom per node
- `elem_size::Float64` : target element size (`:unstructured` only)
- `mesh_path::String` : base mesh directory (`:unstructured` only)
- `template_path::String` : override `.geo` template path

# Returns:
- `MeshgridSquare` : square mesh
"""
function meshgrid_square(lx::X, ly::Y;
    mesh_type::Symbol=:structured,
    element_shape::Symbol=:Quad,
    basis_order::Int=1,
    ne::Union{Int,Nothing}=nothing,
    ndof::Int=2,
    elem_size::Union{Float64,Nothing}=nothing,
    mesh_path::Union{String,Nothing}=nothing,
    template_path::Union{String,Nothing}=nothing
) where {X<:Number, Y<:Number}
    if mesh_type == :structured
        isnothing(ne) && throw(ArgumentError("`ne` is required for mesh_type=:structured"))
        NodeList, IEN, _, IEN_top, IEN_bottom, BorderNodes =
            _meshgrid_square(lx, ly, ne; element_shape=element_shape, basis_order=basis_order)
        total_nodes = size(NodeList, 2)
        ID = zeros(Int64, ndof, total_nodes)
        for m in 1:total_nodes, l in 1:ndof
            ID[l, m] = ndof*(m-1) + l
        end
        empty_mat = Matrix{Int}(undef, 0, 0)
        return MeshgridSquare(lx=lx, ly=ly, NodeList=NodeList, IEN=IEN,
            IEN_top=IEN_top, IEN_bottom=IEN_bottom,
            IEN_sides=empty_mat, ID=ID,
            volume_element_shape=element_shape, basis_order=basis_order,
            nNodes=total_nodes, ne=ne^2,
            side_nodes=BorderNodes[1], top_nodes=BorderNodes[3], bottom_nodes=BorderNodes[2])
    elseif mesh_type == :unstructured
        isnothing(mesh_path) && throw(ArgumentError("`mesh_path` is required for mesh_type=:unstructured"))
        isnothing(elem_size) && throw(ArgumentError("`elem_size` is required for mesh_type=:unstructured"))
        return _get_gmsh_square(mesh_path, ndof, lx, ly, element_shape, basis_order, Float64(elem_size); template_path=template_path)
    else
        throw(ArgumentError("Unknown mesh_type `$mesh_type`; use :structured or :unstructured"))
    end
end

"""
    meshgrid_disk(r; element_shape=:Tri, basis_order=1, ndof=2,
                  elem_size=nothing, mesh_path=nothing, template_path=nothing)

Set up a disk mesh using Gmsh (no structured equivalent).

# Arguments:
- `r::Number` : radius of the disk
- `element_shape::Symbol` : element shape (`:Tri` or `:Quad`)
- `basis_order::Int` : polynomial order (1 or 2)
- `ndof::Int` : degrees of freedom per node
- `elem_size::Float64` : target element size
- `mesh_path::String` : base mesh directory
- `template_path::String` : override `.geo` template path

# Returns:
- `MeshgridDisk` : disk mesh
"""
function meshgrid_disk(r::T;
    element_shape::Symbol=:Tri,
    basis_order::Int=1,
    ndof::Int=2,
    elem_size::Union{Float64,Nothing}=nothing,
    mesh_path::Union{String,Nothing}=nothing,
    template_path::Union{String,Nothing}=nothing
) where {T<:Number}
    isnothing(mesh_path) && throw(ArgumentError("`mesh_path` is required for disk mesh"))
    isnothing(elem_size) && throw(ArgumentError("`elem_size` is required for disk mesh"))
    return _get_gmsh_disk(mesh_path, ndof, r, element_shape, basis_order, Float64(elem_size); template_path=template_path)
end

"""
    meshgrid_line(l; mesh_type=:structured, element_shape=:Line, basis_order=1,
                  ne=nothing, ndof=1, elem_size=nothing, mesh_path=nothing, template_path=nothing)

Set up a line mesh, dispatching to structured (Julia) or unstructured (Gmsh) generation.
For the raw-array form use `_meshgrid_line(l, ne::Int; ...)`.

# Arguments:
- `l::Number` : length of the line
- `mesh_type::Symbol` : `:structured` (Julia) or `:unstructured` (Gmsh)
- `element_shape::Symbol` : element shape (`:Line`)
- `basis_order::Int` : polynomial order (1 or 2)
- `ne::Int` : number of elements (`:structured` only)
- `ndof::Int` : degrees of freedom per node
- `elem_size::Float64` : target element size (`:unstructured` only)
- `mesh_path::String` : base mesh directory (`:unstructured` only)
- `template_path::String` : override `.geo` template path

# Returns:
- `MeshgridLine` : line mesh
"""
function meshgrid_line(l::T;
    mesh_type::Symbol=:structured,
    element_shape::Symbol=:Line,
    basis_order::Int=1,
    ne::Union{Int,Nothing}=nothing,
    ndof::Int=1,
    elem_size::Union{Float64,Nothing}=nothing,
    mesh_path::Union{String,Nothing}=nothing,
    template_path::Union{String,Nothing}=nothing
) where {T<:Number}
    if mesh_type == :structured
        isnothing(ne) && throw(ArgumentError("`ne` is required for mesh_type=:structured"))
        NodeList, IEN, boundary_nodes = _meshgrid_line(l, ne; element_shape=element_shape, basis_order=basis_order)
        total_nodes = size(NodeList, 2)
        ID = zeros(Int64, ndof, total_nodes)
        for m in 1:total_nodes, dof in 1:ndof
            ID[dof, m] = ndof*(m-1) + dof
        end
        return MeshgridLine(lx=l, NodeList=NodeList, IEN=IEN, ID=ID,
            volume_element_shape=element_shape, basis_order=basis_order,
            nNodes=total_nodes, ne=ne, boundary_nodes=boundary_nodes)
    elseif mesh_type == :unstructured
        isnothing(mesh_path) && throw(ArgumentError("`mesh_path` is required for mesh_type=:unstructured"))
        isnothing(elem_size) && throw(ArgumentError("`elem_size` is required for mesh_type=:unstructured"))
        return _get_gmsh_line(mesh_path, ndof, l, element_shape, basis_order, Float64(elem_size); template_path=template_path)
    else
        throw(ArgumentError("Unknown mesh_type `$mesh_type`; use :structured or :unstructured"))
    end
end

"""
    reset_mesh!(mesh::AbstractMeshgrid)

Reset the mesh to its initial state by updating the `NodeList` to the initial state of the mesh.

# Arguments:
- `mesh::AbstractMeshgrid`: The mesh object to be reset.
"""
function reset_mesh!(mesh::T) where {T<:AbstractMeshgrid}
    mesh.NodeList = mesh.initial_state
end

"""
    update_initial_state!(mesh::AbstractMeshgrid, new_state::Matrix{Float64})

Update the initial state of the mesh and reset the current state of the mesh to the initial state.

# Arguments:
- `mesh::AbstractMeshgrid`: The mesh object whose state is to be updated.
- `new_state::Matrix{Float64}`: The new state to be assigned to the mesh.

"""
function update_initial_state!(mesh::T, new_state::Matrix{Float64}) where {T<:AbstractMeshgrid}
    mesh.initial_state = copy(new_state)  # Use `copy` to avoid unintended mutations
    reset_mesh!(mesh)  # Reset the mesh to the initial state
end
