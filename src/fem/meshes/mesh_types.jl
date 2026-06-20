# abstract type Mesh end
using ArgCheck

function _surface_element_shape(vol_shape::Symbol)
    vol_shape == :Hex  && return :Quad
    vol_shape == :Tet  && return :Tri
    vol_shape == :Quad && return :Line
    vol_shape == :Tri  && return :Line
    error("Cannot determine surface element shape for volume shape: $vol_shape")
end

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
