using smearFEM
using Dates
using Plots, Plots.PlotMeasures
using LaTeXStrings
using StatsPlots

function const_vel(r::Float64, h::Float64, ndim::Int, FunctionClass_u::String, nDof_u::Int, FunctionClass_p::String, nDof_p::Int, ne::Int, 
                    camera_matrix::AbstractMatrix{Float64}, camera_pose::AbstractMatrix{Float64})

  control::String = "force"            # "force" or "velocity"
  viscosity_type::String = "constant"  # "constant" or "bulk_viscosity"

  # simulation parameters for the ground truth
  sim_time::Float64 = 20.0            # simulation time in seconds
  steps::Float64 = 200.0              # number of time steps
  t_steps::Float64 = sim_time/steps

  F_ext::Float64 = 1500.0
  F::Vector{Float64} = -F_ext*ones(Float64, round(Int, (sim_time/t_steps))) # force applied to the cylinder in N

  gt_β = 15.0
  gt_η = 40.0

  model, scene = def_problem(r, h, ne, gt_η, ndim, FunctionClass_u, nDof_u, FunctionClass_p, nDof_p, gt_β, F, control, viscosity_type, sim_time, t_steps)
  filepath = string("/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/sim_experiments/single_simulation/fem_runs/Stokes/",control,"_gt/",viscosity_type,"/test_sim")
  write_sim_data(model, scene, camera_matrix, camera_pose, filepath)

end

function main()

  scale = 100
  r::Float64 = 0.25*scale    # radius of the cylinder in mm
  h::Float64 = 0.5*scale     # height of the cylinder in mm
  ndim::Int = 3
  FunctionClass_x::String = "S2"
  FunctionClass_u::String = "Q2"
  nDof_u::Int = ndim              # number of degree of freedom per node
  FunctionClass_p::String = "Q1"
  nDof_p::Int = 1                 # number of degree of freedom per node

  ne::Int = 4 # number of elements in the mesh for the ground truth

  camera_matrix = [[8*2048/7.07, 0.0, 2048/2] [0.0, 8*1536/5.3, 1536/2] [0.0, 0.0, 1.0]]'
  camera_pose = scale*[0 -0.25 2]'   # camera position in mm

  const_vel(r, h, ndim, FunctionClass_u, nDof_u, FunctionClass_p, nDof_p, ne, camera_matrix, camera_pose)
  # bulk_vel(r, h, ndim, FunctionClass_u, nDof_u, FunctionClass_p, nDof_p, ne, camera_matrix, camera_pose)

end

main()

