using smearFEM
using Dates
using Plots, Plots.PlotMeasures
using LaTeXStrings
using StatsPlots


function const_vel(r::Float64, h::Float64, ndim::Int, FunctionClass_u::String, nDof_u::Int, FunctionClass_p::String, nDof_p::Int, FunctionClass_x::String, ne::Int, 
                    camera_matrix::AbstractMatrix{Float64}, obj_pose::AbstractMatrix{Float64})

  control::String = "force"            # "force" or "velocity"
  viscosity_type::String = "constant"  # "constant" or "bulk_viscosity"

  sim_time::Float64 = 20.0           # simulation time in seconds
  t_steps::Float64 = 0.1              # number of time steps
  steps::Float64 = sim_time/t_steps
    
  β_gt = 500.0 # penalty parameter for the ground truth
  η_gt =  70.0 # viscosity in (kg/(mm⋅s))

if β_gt <= 1.0
    F_ext = 9.813e3*0.85
elseif β_gt == 10.0
    F_ext = 9.813e3*0.93
elseif β_gt == 50.0
    F_ext = 9.813e3*0.955
elseif β_gt == 100.0
    F_ext = 9.813e3*0.97
elseif β_gt == 500.0
    F_ext = 9.813e3*0.99
elseif β_gt == 1e3
    F_ext = 9.813e3*1.0
elseif β_gt == 5e3
    F_ext = 9.813e3*1.01 # force applied to the cylinder in kg.mm/s^2 (N)

elseif β_gt == 2e3
    F_ext = 9.813e3*0.995 # force applied to the cylinder in kg.mm/s^2 (N)
elseif β_gt == 1e4
    F_ext = 9.813e3*0.85*700 # force applied to the cylinder in kg.mm/s^2 (N)
elseif β_gt == 1e5
    F_ext = 9.813e3*0.85*7e3 # force applied to the cylinder in kg.mm/s^2 (N)
elseif β_gt == 1e10
    F_ext = 9.813e3*0.85*7e8 # force applied to the cylinder in kg.mm/s^2 (N)
end

  F_ext = 9.813e3 # force applied to the cylinder in kg.mm/s^2 (N)
  F::Vector{Float64} = -F_ext*ones(Float64, round(Int, (sim_time/t_steps))) # force applied to the cylinder in N
#   F::Vector{Float64} = -0.6*ones(Float64, round(Int, (sim_time/t_steps))) # force applied to the cylinder in N

  model, scene = def_problem(r, h, ne, η_gt, ndim, FunctionClass_u, nDof_u, FunctionClass_p, nDof_p, FunctionClass_x, β_gt, F, control, 
                            viscosity_type, sim_time, t_steps)
  filepath = string("/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/sim_experiments/single_simulation/FEM/Stokes/$control/$viscosity_type/test_sim")
#   filepath = string("/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/ground_truth/sim_data/Stokes/force/constant/Q2_6/5")

  write_sim_data(model, scene, camera_matrix, obj_pose, filepath)
end

function main()

  r::Float64 = 25   # radius of the cylinder in mm
  h::Float64 = 40.0     # height of the cylinder in mm
  ndim::Int = 3
  FunctionClass_x::String = "Q2"
  FunctionClass_u::String = "Q2"
  nDof_u::Int = ndim              # number of degree of freedom per node
  FunctionClass_p::String = "Q1"
  nDof_p::Int = 1                 # number of degree of freedom per node

  ne::Int = 6 # number of elements in the mesh for the ground truth

  camera_matrix = [[8*2048/7.07, 0.0, 2048/2] [0.0, 8*1536/5.3, 1536/2] [0.0, 0.0, 1.0]]'
  obj_pose = Float64.([-1.0 0.0 0.0 0.0; 0.0 0.0 -1.0 20.0; 0.0 -1.0 0.0 150; 0.0 0.0 0.0 1.0])


  const_vel(r, h, ndim, FunctionClass_u, nDof_u, FunctionClass_p, nDof_p, FunctionClass_x, ne, camera_matrix, obj_pose)
  # bulk_vel(r, h, ndim, FunctionClass_u, nDof_u, FunctionClass_p, nDof_p, ne, camera_matrix, obj_pose)


end

main()

