using Test
using smearFEM

@testset "testing the stokes model" begin

    # test case 
    r::Float64 = 25   # radius of the cylinder in mm
    h::Float64 = 40.0     # height of the cylinder in mm
    ndim::Int = 3
    element_shape_x::Symbol = :Tet
    basis_order_x::Int = 2  
    element_shape_u::Symbol = :Tet
    basis_order_u::Int = 2
    element_shape_p::Symbol = :Tet
    basis_order_p::Int = 1
    nDof_p::Int = 1                 # number of degree of freedom per node
    nDof_u::Int = ndim              # number of degree of freedom per node

    β = 1e-5
    η = 70.0

    ne::Int = 11 # number of elements in the mesh for the ground truth

    camera_matrix = get_camera_matrix()
    obj_pose = Float64.([-1.0 0.0 0.0 0.0; 0.0 0.0 -1.0 20.0; 0.0 -1.0 0.0 150; 0.0 0.0 0.0 1.0])
    acc = 1e-3

    control::String = "force"            # "force" or "velocity"
    viscosity_type::String = "constant"  # "constant" or "bulk_viscosity"

    sim_time::Float64 = 20.0           # simulation time in seconds
    t_steps::Float64 = 1              # number of time steps
    steps::Float64 = sim_time/t_steps

    F_ext = 9.813e3 # force applied to the cylinder in kg.mm/s^2 (N)
    F::Vector{Float64} = -F_ext*ones(Float64, round(Int, (sim_time/t_steps))) # force applied to the cylinder in N
    

    model_tet, scene_tetq = def_problem(r, h, ne, η, ndim, element_shape_u, basis_order_u, nDof_u, element_shape_p, basis_order_p, nDof_p, 
                    element_shape_x, basis_order_x, β, F, control, viscosity_type, sim_time, t_steps, 
                    mesh_path=joinpath(dirname(@__DIR__), "mesh_files"))

    conditions = Conditions(camera_matrix=camera_matrix, obj_pose=obj_pose)
    
    output, _ = simulate(model_tet, scene_tetq, conditions)
    
    element_shape_x = :Hex
    basis_order_x = 2  
    element_shape_u = :Hex
    basis_order_u = 2
    nDof_u = ndim  # number of degree of freedom per node
    element_shape_p = :Hex
    basis_order_p = 1

    model_hex, scene_hexq = def_problem(r, h, ne, η, ndim, element_shape_u, basis_order_u, nDof_u, element_shape_p, basis_order_p, nDof_p, 
                element_shape_x, basis_order_x, β, F, control, viscosity_type, sim_time, t_steps, 
                mesh_path=joinpath(dirname(@__DIR__), "mesh_files"))

    output_hex, _ = simulate(model_hex, scene_hexq, conditions)

    l_iter = 1:length(output)
    
    for i in l_iter
        @test abs(output[i]) - abs(output_hex[i]) ≈ 0 atol=acc
    end

end