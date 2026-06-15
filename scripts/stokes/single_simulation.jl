using smearFEM
using Dates
using Plots, Plots.PlotMeasures
using LaTeXStrings
using StatsPlots


function const_vel(r::Float64, h::Float64, ndim::Int, element_shape_u::Symbol, basis_order_u::Int, nDof_u::Int, element_shape_p::Symbol, basis_order_p::Int,
                    nDof_p::Int, element_shape_x::Symbol, basis_order_x::Int, ne::Int,
                    camera_matrix::AbstractMatrix{Float64}, obj_pose::AbstractMatrix{Float64};
                    geometry::Symbol=:cylinder, edge_radius::Union{Float64,Nothing}=nothing)

    control::String = "force"            # "force" or "velocity"
    viscosity_type::String = "constant"  # "constant" or "bulk_viscosity"

    sim_time::Float64 = 20.0
    t_steps::Float64 = 1.0

    β_gt = 100.0
    η_gt = 70.0  # viscosity in (kg/(mm⋅s))

    F_ext = 9.813e3  # force applied in kg.mm/s^2 (N)
    F::Vector{Float64} = -F_ext * ones(Float64, round(Int, sim_time / t_steps))

    model, scene = def_problem(r, h, ne, η_gt, ndim, element_shape_u, basis_order_u, nDof_u,
                               element_shape_p, basis_order_p, nDof_p, element_shape_x, basis_order_x,
                               β_gt, F, control, viscosity_type, sim_time, t_steps;
                               geometry=geometry, edge_radius=edge_radius)

    filepath = resolve_data_path("sim_experiments/single_simulation/FEM/Stokes/$control/$viscosity_type/test_sim")

    write_sim_data(model, scene, camera_matrix, obj_pose, filepath)
end

function main()
    r::Float64 = 25.0   # side half-length for cube / radius for cylinder (mm)
    h::Float64 = 40.0   # height (mm)
    ndim::Int = 3
    element_shape_x::Symbol = :Hex
    basis_order_x::Int = 2
    element_shape_u::Symbol = :Hex
    basis_order_u::Int = 2
    element_shape_p::Symbol = :Hex
    basis_order_p::Int = 1
    nDof_p::Int = 1
    nDof_u::Int = ndim

    ne::Int = 6

    geometry::Symbol = :cube # :cylinder or :cube
    edge_radius::Float64 = 1.0   # fillet radius on vertical edges (mm)

    camera_matrix = [[8*2048/7.07, 0.0, 2048/2] [0.0, 8*1536/5.3, 1536/2] [0.0, 0.0, 1.0]]'
    obj_pose = Float64.([-1.0 0.0 0.0 0.0; 0.0 0.0 -1.0 20.0; 0.0 -1.0 0.0 150; 0.0 0.0 0.0 1.0])

    const_vel(r, h, ndim, element_shape_u, basis_order_u, nDof_u, element_shape_p, basis_order_p,
              nDof_p, element_shape_x, basis_order_x, ne, camera_matrix, obj_pose;
              geometry=geometry, edge_radius=edge_radius)
end

main()

