function _get_gmsh_cylinder(file_path::String, ndof::Int, r::T, h::U, element_shape::Symbol=:Hex, basis_order::Int=1, elem_size::Union{Float64,Nothing}=nothing;
                           template_path::Union{String,Nothing}=nothing) where {T<:Number,U<:Number}
    @debug "Reading cylinder mesh from $file_path"
    mesh_dir = joinpath(dirname(dirname(dirname(@__DIR__))), "mesh_files")
    if isnothing(template_path)
        templates = Dict(:Tet => "cylinder_tet.geo", :Hex => "cylinder_hex.geo")
        haskey(templates, element_shape) ||
            throw(ArgumentError("Unknown element shape $element_shape for cylinder mesh. " *
                                "Supported: $(join(keys(templates), ", "))"))
        template_path = joinpath(mesh_dir, "templates", templates[element_shape])
    end

    mesh_params = isnothing(elem_size) ? nothing :
        Dict{String,Any}("radius" => Float64(r), "height" => Float64(h), "nz" => round(Int, elem_size))

    mesh_path = joinpath(file_path, "cylinder", "$(element_shape)", "$(r)x$(h)_$(elem_size)_$(basis_order).msh")
    NodeList, IEN, face_IENs, nNodes, node_sets, ne = _get_mesh_data(mesh_path;
        params=mesh_params, template_path=template_path, mesh_order=basis_order, mesh_dim=3,
        body_group="Volume", body_dim=3,
        face_groups=["Top", "Bottom", "Lateral"],
        node_set_groups=["Lateral", "Bottom", "Top"])

    IEN_top    = face_IENs[1]
    IEN_bottom = face_IENs[2]
    nodes_lateral = something(node_sets[1], Int[])
    nodes_bottom  = something(node_sets[2], Int[])
    nodes_top     = something(node_sets[3], Int[])

    ID = zeros(Int64, ndof, nNodes)
    for m in 1:nNodes
        for l in 1:ndof
            ID[l,m] = ndof*(m-1) + l
        end
    end

    mesh = MeshgridCylinder(r=r, h=h, NodeList=NodeList, IEN=IEN, IEN_top=IEN_top, IEN_bottom=IEN_bottom, ID=ID,
            volume_element_shape=element_shape, basis_order=basis_order,
            nNodes=nNodes, ne=ne, side_nodes=nodes_lateral, top_nodes=nodes_top, bottom_nodes=nodes_bottom,
            effect_elem_sze=elem_size)
    return mesh
end


"""
    _get_gmsh_box(file_path, ndof, lx, ly, lz, element_shape, basis_order, elem_size; template_path)

Load (or auto-generate) a box mesh and return a `MeshgridCuboid`.
Physical groups expected: "Box" (volume), "Bottom", "Top", "Front", "Back", "Left", "Right" (surfaces).
"""
function _get_gmsh_box(file_path::String, ndof::Int, lx::T, ly::U, lz::V,
                      element_shape::Symbol=:Hex, basis_order::Int=1,
                      elem_size::Union{Float64,Nothing}=nothing;
                      template_path::Union{String,Nothing}=nothing,
                      edge_radius::Union{Float64,Nothing}=nothing) where {T<:Number,U<:Number,V<:Number}
    !isnothing(edge_radius) && (edge_radius >= min(lx, ly) / 2) &&
        throw(ArgumentError("`edge_radius` must be < min(lx,ly)/2 = $(min(lx,ly)/2), got $edge_radius"))

    @info "Reading box mesh from $file_path"
    effective_r = something(edge_radius, 0.0)
    mesh_dir = joinpath(dirname(dirname(dirname(@__DIR__))), "mesh_files")
    if isnothing(template_path)
        if element_shape == :Tet
            template_path = effective_r == 0.0 ?
                joinpath(mesh_dir, "templates", "box_tet_sharp.geo") :
                joinpath(mesh_dir, "templates", "box_tet.geo")
        else
            template_path = effective_r == 0.0 ?
                joinpath(mesh_dir, "templates", "box_hex_sharp.geo") :
                joinpath(mesh_dir, "templates", "box_hex.geo")
        end
    end

    mesh_params = isnothing(elem_size) ? nothing :
        Dict{String,Any}("lx" => Float64(lx), "ly" => Float64(ly), "lz" => Float64(lz),
                         "nz" => round(Int, elem_size), "r" => Float64(effective_r))

    mesh_path = joinpath(file_path, "box", "$(element_shape)", "$(lx)x$(ly)x$(lz)_$(elem_size)_$(basis_order)_r$(effective_r).msh")
    NodeList, IEN, face_IENs, nNodes, node_sets, ne = _get_mesh_data(mesh_path;
        params=mesh_params, template_path=template_path, mesh_order=basis_order, mesh_dim=3,
        body_group="Box", body_dim=3,
        face_groups=["Top", "Bottom", "Front", "Back", "Left", "Right"],
        node_set_groups=["Front", "Back", "Left", "Right", "Bottom", "Top"])

    IEN_top    = face_IENs[1]
    IEN_bottom = face_IENs[2]
    IEN_sides_ = [face_IENs[3], face_IENs[4], face_IENs[5], face_IENs[6]]

    nodes_lateral = sort(unique(vcat(
        something(node_sets[1], Int[]),
        something(node_sets[2], Int[]),
        something(node_sets[3], Int[]),
        something(node_sets[4], Int[]))))
    nodes_bottom  = something(node_sets[5], Int[])
    nodes_top     = something(node_sets[6], Int[])

    ID = zeros(Int64, ndof, nNodes)
    for m in 1:nNodes, l in 1:ndof
        ID[l,m] = ndof*(m-1) + l
    end

    empty_mat = Matrix{Int}(undef,0,0)
    return MeshgridCuboid(lx=lx, ly=ly, lz=lz, NodeList=NodeList, IEN=IEN,
        IEN_top=something(IEN_top, empty_mat),
        IEN_bottom=something(IEN_bottom, empty_mat),
        IEN_front=something(IEN_sides_[1], empty_mat),
        IEN_back=something(IEN_sides_[2], empty_mat),
        IEN_left=something(IEN_sides_[3], empty_mat),
        IEN_right=something(IEN_sides_[4], empty_mat),
        ID=ID, volume_element_shape=element_shape, basis_order=basis_order,
        nNodes=nNodes, ne=ne,
        side_nodes=nodes_lateral, top_nodes=nodes_top, bottom_nodes=nodes_bottom,
        edge_radius=edge_radius, effect_elem_sze=elem_size)
end

"""
    _get_gmsh_square(file_path, ndof, lx, ly, element_shape, basis_order, elem_size; template_path)

Load (or auto-generate) a square mesh and return a `MeshgridSquare`.
Physical groups expected: "Square" (surface), "Bottom", "Top", "Left", "Right" (curves).
"""
function _get_gmsh_square(file_path::String, ndof::Int, lx::T, ly::U,
                         element_shape::Symbol=:Quad, basis_order::Int=1,
                         elem_size::Union{Float64,Nothing}=nothing;
                         template_path::Union{String,Nothing}=nothing) where {T<:Number,U<:Number}
    @info "Reading square mesh from $file_path"
    mesh_dir = joinpath(dirname(dirname(dirname(@__DIR__))), "mesh_files")
    if isnothing(template_path)
        template_path = element_shape == :Tri ?
            joinpath(mesh_dir, "templates", "square_tri.geo") :
            joinpath(mesh_dir, "templates", "square_quad.geo")
    end

    mesh_params = isnothing(elem_size) ? nothing :
        Dict{String,Any}("lx" => Float64(lx), "ly" => Float64(ly), "elem_size" => Float64(elem_size))

    mesh_path = joinpath(file_path, "square", "$(element_shape)", "$(elem_size)_$(basis_order).msh")
    NodeList, IEN, face_IENs, nNodes, node_sets, ne = _get_mesh_data(mesh_path;
        params=mesh_params, template_path=template_path, mesh_order=basis_order, mesh_dim=2,
        body_group="Square", body_dim=2,
        face_groups=["Bottom", "Top", "Left", "Right"],
        node_set_groups=["Left", "Bottom", "Top"])

    IEN_bottom = face_IENs[1]
    IEN_top    = face_IENs[2]
    IEN_side   = let l = face_IENs[3], r = face_IENs[4]
        (!isnothing(l) && !isnothing(r)) ? hcat(l, r) :
        something(l, something(r, Matrix{Int}(undef,0,0)))
    end

    nodes_lateral = something(node_sets[1], Int[])
    nodes_bottom  = something(node_sets[2], Int[])
    nodes_top     = something(node_sets[3], Int[])

    ID = zeros(Int64, ndof, nNodes)
    for m in 1:nNodes, l in 1:ndof
        ID[l,m] = ndof*(m-1) + l
    end

    return MeshgridSquare(lx=lx, ly=ly, NodeList=NodeList, IEN=IEN,
        IEN_top=something(IEN_top, Matrix{Int}(undef,0,0)),
        IEN_bottom=something(IEN_bottom, Matrix{Int}(undef,0,0)),
        IEN_sides=IEN_side, ID=ID,
        volume_element_shape=element_shape, basis_order=basis_order,
        nNodes=nNodes, ne=ne,
        side_nodes=nodes_lateral, top_nodes=nodes_top, bottom_nodes=nodes_bottom,
        effect_elem_sze=elem_size)
end

"""
    _get_gmsh_disk(file_path, ndof, r, element_shape, basis_order, elem_size; template_path)

Load (or auto-generate) a disk mesh and return a `MeshgridDisk`.
Physical groups expected: "Disk" (surface), "Boundary" (curve).
"""
function _get_gmsh_disk(file_path::String, ndof::Int, r::T,
                       element_shape::Symbol=:Tri, basis_order::Int=1,
                       elem_size::Union{Float64,Nothing}=nothing;
                       template_path::Union{String,Nothing}=nothing) where {T<:Number}
    @info "Reading disk mesh from $file_path"
    mesh_dir = joinpath(dirname(dirname(dirname(@__DIR__))), "mesh_files")
    if isnothing(template_path)
        template_path = element_shape == :Quad ?
            joinpath(mesh_dir, "templates", "disk_quad.geo") :
            joinpath(mesh_dir, "templates", "disk_tri.geo")
    end

    mesh_params = isnothing(elem_size) ? nothing :
        Dict{String,Any}("radius" => Float64(r), "elem_size" => Float64(elem_size))

    mesh_path = joinpath(file_path, "disk", "$(element_shape)", "$(elem_size)_$(basis_order).msh")
    NodeList, IEN, face_IENs, nNodes, node_sets, ne = _get_mesh_data(mesh_path;
        params=mesh_params, template_path=template_path, mesh_order=basis_order, mesh_dim=2,
        body_group="Disk", body_dim=2,
        face_groups=["Boundary"],
        node_set_groups=["Boundary"])

    IEN_boundary   = something(face_IENs[1], Matrix{Int}(undef,0,0))
    boundary_nodes = something(node_sets[1], Int[])

    ID = zeros(Int64, ndof, nNodes)
    for m in 1:nNodes, l in 1:ndof
        ID[l,m] = ndof*(m-1) + l
    end

    return MeshgridDisk(r=r, NodeList=NodeList, IEN=IEN,
        IEN_boundary=IEN_boundary, ID=ID,
        volume_element_shape=element_shape, basis_order=basis_order,
        nNodes=nNodes, ne=ne, boundary_nodes=boundary_nodes,
        effect_elem_sze=elem_size)
end

"""
    _get_gmsh_line(file_path, ndof, len, element_shape, basis_order, elem_size; template_path)

Load (or auto-generate) a 1D line mesh and return a `MeshgridLine`.
Physical groups expected: "Line" (curve), "Left" and "Right" (points).
"""
function _get_gmsh_line(file_path::String, ndof::Int, len::T,
                       element_shape::Symbol=:Line, basis_order::Int=1,
                       elem_size::Union{Float64,Nothing}=nothing;
                       template_path::Union{String,Nothing}=nothing) where {T<:Number}
    @info "Reading line mesh from $file_path"
    mesh_dir = joinpath(dirname(dirname(dirname(@__DIR__))), "mesh_files")
    isnothing(template_path) && (template_path = joinpath(mesh_dir, "templates", "line.geo"))

    mesh_params = isnothing(elem_size) ? nothing :
        Dict{String,Any}("length" => Float64(len), "elem_size" => Float64(elem_size))

    mesh_path = joinpath(file_path, "line", "$(element_shape)", "$(elem_size)_$(basis_order).msh")
    NodeList, IEN, _, nNodes, node_sets, ne = _get_mesh_data(mesh_path;
        params=mesh_params, template_path=template_path, mesh_order=basis_order, mesh_dim=1,
        body_group="Line", body_dim=1,
        face_groups=String[],
        node_set_groups=["Left", "Right"])

    boundary_nodes = sort(vcat(something(node_sets[1], Int[]), something(node_sets[2], Int[])))

    ID = zeros(Int64, ndof, nNodes)
    for m in 1:nNodes, l in 1:ndof
        ID[l,m] = ndof*(m-1) + l
    end

    return MeshgridLine(lx=len, NodeList=NodeList, IEN=IEN, ID=ID,
        volume_element_shape=element_shape, basis_order=basis_order,
        nNodes=nNodes, ne=ne, boundary_nodes=boundary_nodes,
        effect_elem_sze=elem_size)
end
