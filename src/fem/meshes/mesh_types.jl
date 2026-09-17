# abstract type Mesh end
using ArgCheck

function _surface_element_shape(vol_shape::Symbol)
    vol_shape == :Hex  && return :Quad
    vol_shape == :Tet  && return :Tri
    vol_shape == :Quad && return :Line
    vol_shape == :Tri  && return :Line
    error("Cannot determine surface element shape for volume shape: $vol_shape")
end

"""
    MeshgridLine(; lx=0.0, NodeList, IEN, ID, volume_element_shape=:Line, basis_order=1, nNodes=0, ne=0, boundary_nodes, effect_elem_sze=nothing)

1D mesh of a segment of length `lx`. Being 1D it has no surface elements, so it carries
`boundary_nodes` instead of the per-face `IEN_*` arrays the higher-dimensional meshes have.

# Arguments
- `lx::Number`: Segment length.
- `boundary_nodes::Vector{Int}`: Nodes at the two ends.

Fields shared by every `AbstractMeshgrid`: `NodeList` (current node coordinates),
`IEN` (element connectivity), `ID` (DOF numbering), `volume_element_shape` and `basis_order`
(kept as two independent fields, never a combined "Q2"-style string), `nNodes`, `ne`,
`initial_state` (a copy of `NodeList` at construction, used to reset the mesh), and
`effect_elem_sze` (the realised element size, `nothing` for structured meshes).
"""
mutable struct MeshgridLine <: AbstractMeshgrid
    lx::Number
    NodeList::Matrix{Float64}
    IEN::Matrix{Int}
    ID::Matrix{Int}
    volume_element_shape::Symbol
    basis_order::Int
    nNodes::Int
    ne::Int
    boundary_nodes::Vector{Int}
    initial_state::Matrix{Float64}
    effect_elem_sze::Union{Float64,Nothing}

    function MeshgridLine(;
        lx::Number=0.0,
        NodeList::Matrix{Float64}=Matrix{Float64}(undef, 1, 1),
        IEN::Matrix{Int}=Matrix{Int}(undef, 2, 1),
        ID::Matrix{Int}=Matrix{Int}(undef, 1, 1),
        volume_element_shape::Symbol=:Line,
        basis_order::Int=1,
        nNodes::Int=0,
        ne::Int=0,
        boundary_nodes::Vector{Int}=Vector{Int}(),
        effect_elem_sze::Union{Float64,Nothing}=nothing
    )
        new(lx, NodeList, IEN, ID, volume_element_shape, basis_order, nNodes, ne, boundary_nodes, copy(NodeList), effect_elem_sze)
    end
end

"""
    MeshgridDisk(; r=0.0, NodeList, IEN, IEN_boundary, ID, volume_element_shape=:Quad, surface_element_shape, basis_order=1, nNodes=0, ne=0, boundary_nodes, effect_elem_sze=nothing)

2D mesh of a disk of radius `r`, the axisymmetric cross-section used in place of a full cylinder.
Its single boundary is the circumference, so it has one `IEN_boundary` rather than separate top,
bottom and side arrays.

# Arguments
- `r::Number`: Disk radius.
- `IEN_boundary::Matrix{Int}`: Connectivity of the boundary edge elements.
- `boundary_nodes::Vector{Int}`: Nodes on the circumference.
- `surface_element_shape::Symbol`: Derived from `volume_element_shape` by default.

Fields shared by every `AbstractMeshgrid`: `NodeList` (current node coordinates),
`IEN` (element connectivity), `ID` (DOF numbering), `volume_element_shape` and `basis_order`
(kept as two independent fields, never a combined "Q2"-style string), `nNodes`, `ne`,
`initial_state` (a copy of `NodeList` at construction, used to reset the mesh), and
`effect_elem_sze` (the realised element size, `nothing` for structured meshes).
"""
mutable struct MeshgridDisk <: AbstractMeshgrid
    r::Number
    NodeList::Matrix{Float64}
    IEN::Matrix{Int}
    IEN_boundary::Matrix{Int}
    ID::Matrix{Int}
    volume_element_shape::Symbol
    surface_element_shape::Symbol
    basis_order::Int
    nNodes::Int
    ne::Int
    boundary_nodes::Vector{Int}
    initial_state::Matrix{Float64}
    effect_elem_sze::Union{Float64,Nothing}

    function MeshgridDisk(;
        r::Number=0.0,
        NodeList::Matrix{Float64}=Matrix{Float64}(undef, 2, 1),
        IEN::Matrix{Int}=Matrix{Int}(undef, 3, 1),
        IEN_boundary::Matrix{Int}=Matrix{Int}(undef, 2, 1),
        ID::Matrix{Int}=Matrix{Int}(undef, 2, 1),
        volume_element_shape::Symbol=:Tri,
        surface_element_shape::Symbol=_surface_element_shape(volume_element_shape),
        basis_order::Int=1,
        nNodes::Int=0,
        ne::Int=0,
        boundary_nodes::Vector{Int}=Vector{Int}(),
        effect_elem_sze::Union{Float64,Nothing}=nothing
    )
        new(r, NodeList, IEN, IEN_boundary, ID, volume_element_shape, surface_element_shape,
            basis_order, nNodes, ne, boundary_nodes, copy(NodeList), effect_elem_sze)
    end
end

"""
    MeshgridSquare(; lx=0.0, ly=0.0, NodeList, IEN, IEN_top, IEN_bottom, IEN_sides, ID, volume_element_shape=:Quad, surface_element_shape, basis_order=1, nNodes=0, ne=0, top_nodes, bottom_nodes, side_nodes, effect_elem_sze=nothing)

2D mesh of a rectangle, the plane-strain cross-section of a squeezed block. Boundaries are split
into top, bottom and sides so the two plates can be driven independently of the free surface.

# Arguments
- `lx::Number`, `ly::Number`: Rectangle dimensions.
- `IEN_top`, `IEN_bottom`, `IEN_sides::Matrix{Int}`: Connectivity of each boundary group.
- `top_nodes`, `bottom_nodes`, `side_nodes::Vector{Int}`: Nodes in each boundary group.
- `surface_element_shape::Symbol`: Derived from `volume_element_shape` by default.

Fields shared by every `AbstractMeshgrid`: `NodeList` (current node coordinates),
`IEN` (element connectivity), `ID` (DOF numbering), `volume_element_shape` and `basis_order`
(kept as two independent fields, never a combined "Q2"-style string), `nNodes`, `ne`,
`initial_state` (a copy of `NodeList` at construction, used to reset the mesh), and
`effect_elem_sze` (the realised element size, `nothing` for structured meshes).
"""
mutable struct MeshgridSquare <: AbstractMeshgrid
    lx::Number
    ly::Number
    NodeList::Matrix{Float64}
    IEN::Matrix{Int}
    IEN_top::Matrix{Int}
    IEN_bottom::Matrix{Int}
    IEN_sides::Matrix{Int}
    ID::Matrix{Int}
    volume_element_shape::Symbol
    surface_element_shape::Symbol
    basis_order::Int
    nNodes::Int
    ne::Int
    top_nodes::Vector{Int}
    bottom_nodes::Vector{Int}
    side_nodes::Vector{Int}
    initial_state::Matrix{Float64}
    effect_elem_sze::Union{Float64,Nothing}

    function MeshgridSquare(;
        lx::Number=0.0,
        ly::Number=0.0,
        NodeList::Matrix{Float64}=Matrix{Float64}(undef, 2, 1),
        IEN::Matrix{Int}=Matrix{Int}(undef, 4, 1),
        IEN_top::Matrix{Int}=Matrix{Int}(undef, 2, 1),
        IEN_bottom::Matrix{Int}=Matrix{Int}(undef, 2, 1),
        IEN_sides::Matrix{Int}=Matrix{Int}(undef, 2, 1),
        ID::Matrix{Int}=Matrix{Int}(undef, 2, 1),
        volume_element_shape::Symbol=:Quad,
        surface_element_shape::Symbol=_surface_element_shape(volume_element_shape),
        basis_order::Int=1,
        nNodes::Int=0,
        ne::Int=0,
        top_nodes::Vector{Int}=Vector{Int}(),
        bottom_nodes::Vector{Int}=Vector{Int}(),
        side_nodes::Vector{Int}=Vector{Int}(),
        effect_elem_sze::Union{Float64,Nothing}=nothing
    )
        new(lx, ly, NodeList, IEN, IEN_top, IEN_bottom, IEN_sides, ID,
            volume_element_shape, surface_element_shape, basis_order, nNodes, ne,
            top_nodes, bottom_nodes, side_nodes, copy(NodeList), effect_elem_sze)
    end
end

"""
    MeshgridCuboid(; lx=0.0, ly=0.0, lz=0.0, NodeList, IEN, IEN_top, IEN_bottom, IEN_front, IEN_back, IEN_left, IEN_right, ID, volume_element_shape=:Hex, surface_element_shape, basis_order=1, nNodes=0, ne=0, top_nodes, bottom_nodes, side_nodes, edge_radius=nothing, effect_elem_sze=nothing)

3D mesh of a box. All six faces get their own `IEN_*` array, but the four vertical faces share a
single `side_nodes` list, since they are driven as one free surface.

# Arguments
- `lx::Number`, `ly::Number`, `lz::Number`: Box dimensions.
- `IEN_top`, `IEN_bottom`, `IEN_front`, `IEN_back`, `IEN_left`, `IEN_right::Matrix{Int}`:
  Connectivity of each face.
- `top_nodes`, `bottom_nodes`, `side_nodes::Vector{Int}`: Nodes in each boundary group.
- `edge_radius::Union{Float64,Nothing}`: Fillet radius on the vertical edges, `nothing` if sharp.
- `surface_element_shape::Symbol`: Derived from `volume_element_shape` by default.

Fields shared by every `AbstractMeshgrid`: `NodeList` (current node coordinates),
`IEN` (element connectivity), `ID` (DOF numbering), `volume_element_shape` and `basis_order`
(kept as two independent fields, never a combined "Q2"-style string), `nNodes`, `ne`,
`initial_state` (a copy of `NodeList` at construction, used to reset the mesh), and
`effect_elem_sze` (the realised element size, `nothing` for structured meshes).
"""
mutable struct MeshgridCuboid <: AbstractMeshgrid
    lx::Number
    ly::Number
    lz::Number
    NodeList::Matrix{Float64}
    IEN::Matrix{Int}
    IEN_top::Matrix{Int}
    IEN_bottom::Matrix{Int}
    IEN_front::Matrix{Int}
    IEN_back::Matrix{Int}
    IEN_left::Matrix{Int}
    IEN_right::Matrix{Int}
    ID::Matrix{Int}
    volume_element_shape::Symbol
    surface_element_shape::Symbol
    basis_order::Int
    nNodes::Int
    ne::Int
    top_nodes::Vector{Int}
    bottom_nodes::Vector{Int}
    side_nodes::Vector{Int}
    initial_state::Matrix{Float64}
    edge_radius::Union{Float64,Nothing}
    effect_elem_sze::Union{Float64,Nothing}

    function MeshgridCuboid(;
        lx::Number=0.0,
        ly::Number=0.0,
        lz::Number=0.0,
        NodeList::Matrix{Float64}=Matrix{Float64}(undef, 3, 1),
        IEN::Matrix{Int}=Matrix{Int}(undef, 8, 1),
        IEN_top::Matrix{Int}=Matrix{Int}(undef, 4, 1),
        IEN_bottom::Matrix{Int}=Matrix{Int}(undef, 4, 1),
        IEN_front::Matrix{Int}=Matrix{Int}(undef, 4, 1),
        IEN_back::Matrix{Int}=Matrix{Int}(undef, 4, 1),
        IEN_left::Matrix{Int}=Matrix{Int}(undef, 4, 1),
        IEN_right::Matrix{Int}=Matrix{Int}(undef, 4, 1),
        ID::Matrix{Int}=Matrix{Int}(undef, 3, 1),
        volume_element_shape::Symbol=:Hex,
        surface_element_shape::Symbol=_surface_element_shape(volume_element_shape),
        basis_order::Int=1,
        nNodes::Int=0,
        ne::Int=0,
        top_nodes::Vector{Int}=Vector{Int}(),
        bottom_nodes::Vector{Int}=Vector{Int}(),
        side_nodes::Vector{Int}=Vector{Int}(),
        edge_radius::Union{Float64,Nothing}=nothing,
        effect_elem_sze::Union{Float64,Nothing}=nothing
    )
        new(lx, ly, lz, NodeList, IEN, IEN_top, IEN_bottom, IEN_front, IEN_back, IEN_left, IEN_right,
            ID, volume_element_shape, surface_element_shape, basis_order, nNodes, ne, top_nodes, bottom_nodes, side_nodes, copy(NodeList), edge_radius, effect_elem_sze)
    end
end

"""
    MeshgridCylinder(; r=0.0, h=0.0, NodeList, IEN, IEN_top, IEN_bottom, IEN_sides, ID, volume_element_shape=:Hex, surface_element_shape, basis_order=1, nNodes=0, ne=0, top_nodes, bottom_nodes, side_nodes, effect_elem_sze=nothing)

3D mesh of a cylinder of radius `r` and height `h` — the default geometry for squeeze flow.
Boundaries are split into top, bottom and curved side so the two plates can be driven
independently of the free surface.

# Arguments
- `r::Number`: Cylinder radius.
- `h::Number`: Cylinder height.
- `IEN_top`, `IEN_bottom`, `IEN_sides::Matrix{Int}`: Connectivity of each boundary group.
- `top_nodes`, `bottom_nodes`, `side_nodes::Vector{Int}`: Nodes in each boundary group.
- `surface_element_shape::Symbol`: Derived from `volume_element_shape` by default.

Fields shared by every `AbstractMeshgrid`: `NodeList` (current node coordinates),
`IEN` (element connectivity), `ID` (DOF numbering), `volume_element_shape` and `basis_order`
(kept as two independent fields, never a combined "Q2"-style string), `nNodes`, `ne`,
`initial_state` (a copy of `NodeList` at construction, used to reset the mesh), and
`effect_elem_sze` (the realised element size, `nothing` for structured meshes).
"""
mutable struct MeshgridCylinder <: AbstractMeshgrid
    r::Number
    h::Number
    NodeList::Matrix{Float64}
    IEN::Matrix{Int}
    IEN_top::Matrix{Int}
    IEN_bottom::Matrix{Int}
    IEN_sides::Matrix{Int}
    ID::Matrix{Int}
    volume_element_shape::Symbol
    surface_element_shape::Symbol
    basis_order::Int
    nNodes::Int
    ne::Int
    top_nodes::Vector{Int}
    bottom_nodes::Vector{Int}
    side_nodes::Vector{Int}
    initial_state::Matrix{Float64}
    effect_elem_sze::Union{Float64,Nothing}

    function MeshgridCylinder(;
        r::Number=0.0,
        h::Number=0.0,
        NodeList::Matrix{Float64}=Matrix{Float64}(undef, 3, 1),
        IEN::Matrix{Int}=Matrix{Int}(undef, 8, 1),
        IEN_top::Matrix{Int}=Matrix{Int}(undef, 4, 1),
        IEN_bottom::Matrix{Int}=Matrix{Int}(undef, 4, 1),
        IEN_sides::Matrix{Int}=Matrix{Int}(undef, 4, 1),
        ID::Matrix{Int}=Matrix{Int}(undef, 3, 1),
        volume_element_shape::Symbol=:Hex,
        surface_element_shape::Symbol=_surface_element_shape(volume_element_shape),
        basis_order::Int=1,
        nNodes::Int=0,
        ne::Int=0,
        top_nodes::Vector{Int}=Vector{Int}(),
        bottom_nodes::Vector{Int}=Vector{Int}(),
        side_nodes::Vector{Int}=Vector{Int}(),
        effect_elem_sze::Union{Float64,Nothing}=nothing
    )
        new(r, h, NodeList, IEN, IEN_top, IEN_bottom, IEN_sides,
            ID, volume_element_shape, surface_element_shape, basis_order, nNodes, ne, top_nodes, bottom_nodes, side_nodes, copy(NodeList), effect_elem_sze)
    end
end
