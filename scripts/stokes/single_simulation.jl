using smearFEM
using Dates
using Plots, Plots.PlotMeasures
using LaTeXStrings
using StatsPlots


function const_vel(r::Float64, h::Float64,
                    element_shape_u::Symbol, basis_order_u::Int, nDof_u::Int,
                    element_shape_p::Symbol, basis_order_p::Int, nDof_p::Int,
                    element_shape_x::Symbol, basis_order_x::Int, ne::Int,
                    camera_matrix::AbstractMatrix{Float64}, obj_pose::Vector{Float64}, z_angle_list::Vector{Float64})

    control::String = "force"
    viscosity_type::String = "constant"
    sim_time::Float64 = 20.0
    t_steps::Float64 = 1.0
    β_gt = 100.0
    η_gt = 70.0
    F_ext = 9.813e3
    F::Vector{Float64} = -F_ext * ones(Float64, round(Int, sim_time / t_steps))

    model, scene = def_problem(Cylinder(r, h), ne, η_gt,
                               element_shape_u, basis_order_u, nDof_u,
                               element_shape_p, basis_order_p, nDof_p,
                               element_shape_x, basis_order_x,
                               β_gt, F, control, viscosity_type, sim_time, t_steps)

    filepath = resolve_data_path("sim_experiments/single_simulation/FEM/Stokes/$control/$viscosity_type/test_sim")
    write_sim_data(model, scene, camera_matrix, obj_pose, z_angle_list, filepath)
end

function const_vel(lx::Float64, ly::Float64, lz::Float64,
                    element_shape_u::Symbol, basis_order_u::Int, nDof_u::Int,
                    element_shape_p::Symbol, basis_order_p::Int, nDof_p::Int,
                    element_shape_x::Symbol, basis_order_x::Int, ne::Int,
                    camera_matrix::AbstractMatrix{Float64}, obj_pose::Vector{Float64}, z_angle_list::Vector{Float64};
                    edge_radius::Union{Float64,Nothing}=nothing)

    control::String = "force"
    viscosity_type::String = "constant"
    sim_time::Float64 = 20.0
    t_steps::Float64 = 1.0
    β_gt = 100.0
    η_gt = 70.0
    F_ext = 9.813e3
    F::Vector{Float64} = -F_ext * ones(Float64, round(Int, sim_time / t_steps))

    model, scene = def_problem(Cuboid(lx, ly, lz), ne, η_gt,
                               element_shape_u, basis_order_u, nDof_u,
                               element_shape_p, basis_order_p, nDof_p,
                               element_shape_x, basis_order_x,
                               β_gt, F, control, viscosity_type, sim_time, t_steps;
                               edge_radius=edge_radius)

    filepath = resolve_data_path("sim_experiments/single_simulation/FEM/Stokes/$control/$viscosity_type/test_sim")
    write_sim_data(model, scene, camera_matrix, obj_pose, z_angle_list, filepath)
end

function main()
    element_shape_x::Symbol = :Hex
    basis_order_x::Int = 2
    element_shape_u::Symbol = :Hex
    basis_order_u::Int = 2
    element_shape_p::Symbol = :Hex
    basis_order_p::Int = 1
    nDof_p::Int = 1
    nDof_u::Int = 3
    ne::Int = 6

    camera_matrix = [[8*2048/7.07, 0.0, 2048/2] [0.0, 8*1536/5.3, 1536/2] [0.0, 0.0, 1.0]]'
    obj_pose = [150.0, 0.0, 20.0]
    z_angle_list = [0.0, 30.0, 60.0]

    lx::Float64 = 50.0
    ly::Float64 = 50.0
    lz::Float64 = 40.0
    edge_radius::Float64 = 1.0

    const_vel(lx, ly, lz, element_shape_u, basis_order_u, nDof_u, element_shape_p, basis_order_p,
              nDof_p, element_shape_x, basis_order_x, ne, camera_matrix, obj_pose, z_angle_list;
              edge_radius=edge_radius)
end

main()
