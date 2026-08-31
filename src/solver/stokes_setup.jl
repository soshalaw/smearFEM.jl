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
    squeeze = SqueezeFlow(stokes, [β], [q_tp, q_side, q_btm], C_uc, control, sim_time, t_steps, viscosity_type, cParam)
    return stokes, squeeze
end

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
    squeeze = SqueezeFlow(stokes, [β], [q_tp, q_side, q_btm], C_uc, control, sim_time, t_steps, viscosity_type, cParam)
    return stokes, squeeze
end
