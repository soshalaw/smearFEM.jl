"""
    get_η_power_law(t, F, R_0, H_0, η_0)

Calculate the shear viscosity using the power law viscosity model.
"""
function get_η_power_law(t::T, F::U, R_0::V, H_0::W, η_0::X) where {T<:Number,U<:Number,V<:Number,W<:Number,X<:Number}
    n::Float64 = 0.9
    K::Float64 = 100.0

    H(t) = H_0*(1+8*H_0^2*F*t/(3*π*η_0*R_0^4))^(-1/4)
    R(t) = R_0*(1+8*H_0^2*F*t/(3*π*η_0*R_0^4))^(1/8)
    H_dot(t) = 8/3*(-2*F*H(t)^3/(8*π*η_0*R(t)^4))
    γ_dot(t) = H_dot(t)/H(t)
    η(t) = K*(abs(γ_dot(t)))^(n-1)

    return η(t)
end

"""
    set_model(geom::Cylinder, ne, η, element_shape_u, basis_order_u, nDof_u, element_shape_p, basis_order_p, nDof_p, element_shape_x, basis_order_x; GMESH_MESH=true, filepath_mesh="")

Build the three meshes of the mixed Stokes discretization for a cylinder and assemble them into
a `Stokes` model. Velocity, pressure and geometry meshes come from the same geometry at
different element order; the geometry mesh always has `ndof=1`.

# Arguments
- `geom::Cylinder`: Geometry, supplying radius and height.
- `ne::Float64`: Element size for unstructured meshing; rounded to an element count for
  structured meshing.
- `η::Vector{Float64}`: Shear viscosity, either one entry or one per time step.
- `element_shape_u::Symbol`, `basis_order_u::Int`, `nDof_u::Int64`: Velocity mesh.
- `element_shape_p::Symbol`, `basis_order_p::Int`, `nDof_p::Int64`: Pressure mesh.
- `element_shape_x::Symbol`, `basis_order_x::Int`: Geometry mesh.
- `GMESH_MESH::Bool`: Mesh with Gmsh (`true`) or the structured Julia mesher (`false`).
- `filepath_mesh::String`: Root of the Gmsh mesh cache.

# Returns
- `::Stokes`: Model holding all three meshes.
"""
function set_model(geom::Cylinder, ne::Float64, η::Vector{Float64},
                   element_shape_u::Symbol, basis_order_u::Int, nDof_u::Int64,
                   element_shape_p::Symbol, basis_order_p::Int, nDof_p::Int64,
                   element_shape_x::Symbol, basis_order_x::Int;
                   GMESH_MESH::Bool=true, filepath_mesh::String="")::Stokes
    _dim = ndim(geom)
    mesh_type = GMESH_MESH ? :unstructured : :structured
    _make_mesh(es, bo, ndof) = meshgrid_cylinder(geom.r, geom.h; mesh_type=mesh_type, element_shape=es,
        basis_order=bo, ne=round(Int, ne), ndof=ndof, elem_size=ne, mesh_path=filepath_mesh)
    mesh_u = _make_mesh(element_shape_u, basis_order_u, nDof_u)
    mesh_p = _make_mesh(element_shape_p, basis_order_p, nDof_p)
    mesh_x = _make_mesh(element_shape_x, basis_order_x, 1)
    return Stokes(ndim=_dim, mesh_x=mesh_x, mesh_u=mesh_u, nDof_u=nDof_u, mesh_p=mesh_p, nDof_p=nDof_p, η=η)
end

"""
    set_model(geom::Cuboid, ne, η, ...; GMESH_MESH=true, filepath_mesh="", edge_radius=nothing)

Cuboid variant of `set_model`; see the `Cylinder` method for the shared arguments. Meshes from
`lx`, `ly`, `lz` and additionally accepts `edge_radius` to fillet the vertical edges.

# Returns
- `::Stokes`: Model holding all three meshes.
"""
function set_model(geom::Cuboid, ne::Float64, η::Vector{Float64},
                   element_shape_u::Symbol, basis_order_u::Int, nDof_u::Int64,
                   element_shape_p::Symbol, basis_order_p::Int, nDof_p::Int64,
                   element_shape_x::Symbol, basis_order_x::Int;
                   GMESH_MESH::Bool=true, filepath_mesh::String="",
                   edge_radius::Union{Float64,Nothing}=nothing)::Stokes
    _dim = ndim(geom)
    mesh_type = GMESH_MESH ? :unstructured : :structured
    _make_mesh(es, bo, ndof) = meshgrid_cuboid(geom.lx, geom.ly, geom.lz; mesh_type=mesh_type, element_shape=es,
        basis_order=bo, ne=round(Int, ne), ndof=ndof, elem_size=ne, mesh_path=filepath_mesh, edge_radius=edge_radius)
    mesh_u = _make_mesh(element_shape_u, basis_order_u, nDof_u)
    mesh_p = _make_mesh(element_shape_p, basis_order_p, nDof_p)
    mesh_x = _make_mesh(element_shape_x, basis_order_x, 1)
    return Stokes(ndim=_dim, mesh_x=mesh_x, mesh_u=mesh_u, nDof_u=nDof_u, mesh_p=mesh_p, nDof_p=nDof_p, η=η)
end

"""
    set_model(geom::Disk, ne, η, ...; GMESH_MESH=true, filepath_mesh="")

Disk variant of `set_model`; see the `Cylinder` method for the shared arguments. Errors unless
`GMESH_MESH=true` — no structured mesher exists for a disk.

# Returns
- `::Stokes`: Model holding all three meshes.
"""
function set_model(geom::Disk, ne::Float64, η::Vector{Float64},
                   element_shape_u::Symbol, basis_order_u::Int, nDof_u::Int64,
                   element_shape_p::Symbol, basis_order_p::Int, nDof_p::Int64,
                   element_shape_x::Symbol, basis_order_x::Int;
                   GMESH_MESH::Bool=true, filepath_mesh::String="")::Stokes
    _dim = ndim(geom)
    GMESH_MESH || error("Disk geometry only supports unstructured (Gmsh) mesh")
    _make_mesh(es, bo, ndof) = meshgrid_disk(geom.r; element_shape=es,
        basis_order=bo, ne=round(Int, ne), ndof=ndof, elem_size=ne, mesh_path=filepath_mesh)
    mesh_u = _make_mesh(element_shape_u, basis_order_u, nDof_u)
    mesh_p = _make_mesh(element_shape_p, basis_order_p, nDof_p)
    mesh_x = _make_mesh(element_shape_x, basis_order_x, 1)
    return Stokes(ndim=_dim, mesh_x=mesh_x, mesh_u=mesh_u, nDof_u=nDof_u, mesh_p=mesh_p, nDof_p=nDof_p, η=η)
end

"""
    set_model(geom::Square, ne, η, ...; GMESH_MESH=true, filepath_mesh="")

Square variant of `set_model`; see the `Cylinder` method for the shared arguments. Meshes from
`lx` and `ly`.

# Returns
- `::Stokes`: Model holding all three meshes.
"""
function set_model(geom::Square, ne::Float64, η::Vector{Float64},
                   element_shape_u::Symbol, basis_order_u::Int, nDof_u::Int64,
                   element_shape_p::Symbol, basis_order_p::Int, nDof_p::Int64,
                   element_shape_x::Symbol, basis_order_x::Int;
                   GMESH_MESH::Bool=true, filepath_mesh::String="")::Stokes
    _dim = ndim(geom)
    mesh_type = GMESH_MESH ? :unstructured : :structured
    _make_mesh(es, bo, ndof) = meshgrid_square(geom.lx, geom.ly; mesh_type=mesh_type, element_shape=es,
        basis_order=bo, ne=round(Int, ne), ndof=ndof, elem_size=ne, mesh_path=filepath_mesh)
    mesh_u = _make_mesh(element_shape_u, basis_order_u, nDof_u)
    mesh_p = _make_mesh(element_shape_p, basis_order_p, nDof_p)
    mesh_x = _make_mesh(element_shape_x, basis_order_x, 1)
    return Stokes(ndim=_dim, mesh_x=mesh_x, mesh_u=mesh_u, nDof_u=nDof_u, mesh_p=mesh_p, nDof_p=nDof_p, η=η)
end

"""
    set_model(geom::Segment, ne, η, ...; GMESH_MESH=true, filepath_mesh="")

Segment variant of `set_model`; see the `Cylinder` method for the shared arguments. Always
meshes structurally, so `GMESH_MESH` and `filepath_mesh` are ignored.

# Returns
- `::Stokes`: Model holding all three meshes.
"""
function set_model(geom::Segment, ne::Float64, η::Vector{Float64},
                   element_shape_u::Symbol, basis_order_u::Int, nDof_u::Int64,
                   element_shape_p::Symbol, basis_order_p::Int, nDof_p::Int64,
                   element_shape_x::Symbol, basis_order_x::Int;
                   GMESH_MESH::Bool=true, filepath_mesh::String="")::Stokes
    _dim = ndim(geom)
    _make_mesh(es, bo, ndof) = meshgrid_line(geom.l; element_shape=es,
        basis_order=bo, ne=round(Int, ne), ndof=ndof)
    mesh_u = _make_mesh(element_shape_u, basis_order_u, nDof_u)
    mesh_p = _make_mesh(element_shape_p, basis_order_p, nDof_p)
    mesh_x = _make_mesh(element_shape_x, basis_order_x, 1)
    return Stokes(ndim=_dim, mesh_x=mesh_x, mesh_u=mesh_u, nDof_u=nDof_u, mesh_p=mesh_p, nDof_p=nDof_p, η=η)
end

"""
    def_problem(geom::Cylinder, ne, η_0, element_shape_u, basis_order_u, nDof_u, element_shape_p, basis_order_p, nDof_p, element_shape_x, basis_order_x, β, cParam, control, viscosity_type, sim_time, t_steps; viscosity_model="power_law", GMESH_MESH=true, mesh_path=...)

Wire up a complete squeeze-flow problem for a cylinder: build the viscosity history, the model
and its boundary conditions, and pack them into a `SqueezeFlow` scenario.

With `viscosity_type == "bulk_viscosity"` the viscosity is resolved per time step, from the
power-law model when `viscosity_model == "power_law"` and as a constant fill otherwise. Any
other `viscosity_type` yields a single constant entry.

# Arguments
- `geom::Cylinder`: Geometry, supplying radius and height.
- `ne::Z`: Element size, passed through to `set_model`.
- `η_0::V`: Reference shear viscosity.
- `element_shape_u`, `basis_order_u`, `nDof_u`: Velocity mesh.
- `element_shape_p`, `basis_order_p`, `nDof_p`: Pressure mesh.
- `element_shape_x`, `basis_order_x`: Geometry mesh.
- `β::Y`: Boundary slip/friction parameter.
- `cParam::Vector{Float64}`: Control history, one entry per time step. A force here is in
  kg*mm/s^2, not newtons. Logs an error if shorter than the time array.
- `control::String`: Driving mode, e.g. constant velocity or constant force.
- `viscosity_type::String`: `"bulk_viscosity"` for a per-step history, anything else for a
  constant.
- `sim_time::W`: Total simulated time in seconds.
- `t_steps::X`: Time step size in seconds.
- `viscosity_model::String`: Law used when `viscosity_type == "bulk_viscosity"`.
- `GMESH_MESH::Bool`: Mesh with Gmsh rather than the structured mesher.
- `mesh_path::String`: Root of the Gmsh mesh cache.

# Returns
- `stokes::Stokes`: The assembled model.
- `squeeze::SqueezeFlow`: The scenario driving it.
"""
function def_problem(geom::Cylinder, ne::Z, η_0::V,
                    element_shape_u::Symbol, basis_order_u::Int, nDof_u::Int64,
                    element_shape_p::Symbol, basis_order_p::Int, nDof_p::Int64,
                    element_shape_x::Symbol, basis_order_x::Int,
                    β::Y, cParam::Vector{Float64}, control::String, viscosity_type::String,
                    sim_time::W, t_steps::X;
                    viscosity_model::String="power_law", GMESH_MESH::Bool=true,
                    mesh_path::String = joinpath(dirname(dirname(@__DIR__)), "mesh_files")) where {V<:Number,W<:Number,X<:Number,Y<:Number,Z<:Number}

    time = collect(Float64, range(start=t_steps, stop=sim_time, step=t_steps))
    len_t::Int = length(time)
    @info "Simulation time: $sim_time, Time step: $t_steps, Number of time steps: $(round(Int, sim_time/t_steps))"
    @info "Length of time array: $(len_t)"

    if length(cParam) < len_t
        @error "Length of the Force vector ($(length(cParam))) is less than length of time array ($(length(time)))"
    end

    η = if viscosity_type == "bulk_viscosity"
            if viscosity_model == "power_law"
                @info "Using power law viscosity model"
                get_η_power_law.(time, -cParam[1:len_t], geom.r, geom.h, η_0)
            else
                fill(Float64(η_0), len_t)
            end
        else
            [Float64(η_0)]
        end

    stokes = set_model(geom, float(ne), η, element_shape_u, basis_order_u, nDof_u,
                       element_shape_p, basis_order_p, nDof_p, element_shape_x, basis_order_x;
                       filepath_mesh=mesh_path, GMESH_MESH=GMESH_MESH)
    q_tp, q_side, q_btm, C_uc = set_boundary_cond(stokes)
    squeeze = SqueezeFlow([β], [q_tp, q_btm, q_side], C_uc, control, sim_time, t_steps, viscosity_type, cParam)
    return stokes, squeeze
end

"""
    def_problem(geom::Cuboid, ne, η_0, ...; viscosity_model="power_law", GMESH_MESH=true, mesh_path=...)

Cuboid variant of `def_problem`; see the `Cylinder` method for the shared arguments. The
power-law viscosity uses the cross-sectional diagonal `sqrt(lx^2 + ly^2)` as the equivalent
radius and `lz` as the height, and `geom.edge_radius` is forwarded to `set_model`.

# Returns
- `stokes::Stokes`: The assembled model.
- `squeeze::SqueezeFlow`: The scenario driving it.
"""
function def_problem(geom::Cuboid, ne::Z, η_0::V,
                    element_shape_u::Symbol, basis_order_u::Int, nDof_u::Int64,
                    element_shape_p::Symbol, basis_order_p::Int, nDof_p::Int64,
                    element_shape_x::Symbol, basis_order_x::Int,
                    β::Y, cParam::Vector{Float64}, control::String, viscosity_type::String,
                    sim_time::W, t_steps::X;
                    viscosity_model::String="power_law", GMESH_MESH::Bool=true,
                    mesh_path::String = joinpath(dirname(dirname(@__DIR__)), "mesh_files")) where {V<:Number,W<:Number,X<:Number,Y<:Number,Z<:Number}

    edge_radius = geom.edge_radius
    time = collect(Float64, range(start=t_steps, stop=sim_time, step=t_steps))
    len_t::Int = length(time)
    @info "Simulation time: $sim_time, Time step: $t_steps, Number of time steps: $(round(Int, sim_time/t_steps))"
    @info "Length of time array: $(len_t)"

    if length(cParam) < len_t
        @error "Length of the Force vector ($(length(cParam))) is less than length of time array ($(length(time)))"
    end

    η = if viscosity_type == "bulk_viscosity"
        if viscosity_model == "power_law"
            @info "Using power law viscosity model"
            get_η_power_law.(time, -cParam[1:len_t], sqrt(geom.lx^2 + geom.ly^2), geom.lz, η_0)
        else
            fill(Float64(η_0), len_t)
        end
    else
        [Float64(η_0)]
    end

    stokes = set_model(geom, float(ne), η, element_shape_u, basis_order_u, nDof_u,
                       element_shape_p, basis_order_p, nDof_p, element_shape_x, basis_order_x;
                       filepath_mesh=mesh_path, GMESH_MESH=GMESH_MESH, edge_radius=edge_radius)
    q_tp, q_side, q_btm, C_uc = set_boundary_cond(stokes)
    squeeze = SqueezeFlow([β], [q_tp, q_btm, q_side], C_uc, control, sim_time, t_steps, viscosity_type, cParam)
    return stokes, squeeze
end
