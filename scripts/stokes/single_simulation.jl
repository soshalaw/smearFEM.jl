using smearFEM
using Dates
using Plots, Plots.PlotMeasures
using LaTeXStrings
using StatsPlots

function const_vel(r::Float64, h::Float64, ndim::Int, FunctionClass_u::String, nDof_u::Int, FunctionClass_p::String, nDof_p::Int, FunctionClass_x::String, ne::Int, 
                    camera_matrix::AbstractMatrix{Float64}, camera_pose::AbstractMatrix{Float64})

  control::String = "force"            # "force" or "velocity"
  viscosity_type::String = "constant"  # "constant" or "bulk_viscosity"

  sim_time::Float64 = 20.0           # simulation time in seconds
  steps::Float64 = 20.0              # number of time steps
  t_steps::Float64 = sim_time/steps

  gt_β = 50.0
  gt_η = 40.0

  if gt_β <= 1.0
    F_ext::Float64 = 1500.0
  else  
    F_ext::Float64 = 1500.0*gt_β
  end
  F::Vector{Float64} = -F_ext*ones(Float64, round(Int, (sim_time/t_steps))) # force applied to the cylinder in N

  model, scene = def_problem(r, h, ne, gt_η, ndim, FunctionClass_u, nDof_u, FunctionClass_p, nDof_p, FunctionClass_x, gt_β, F, control, 
                            viscosity_type, sim_time, t_steps)
  filepath = string("/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/sim_experiments/single_simulation/fem_runs/Stokes/",control,"_gt/",viscosity_type,"/test_sim/",FunctionClass_x)
  write_sim_data(model, scene, camera_matrix, camera_pose, filepath)
  
end

function main()

  scale = 50
  r::Float64 = 0.5*scale    # radius of the cylinder in mm
  h::Float64 = 1*scale     # height of the cylinder in mm
  ndim::Int = 3
  FunctionClass_x::String = "Q2"
  FunctionClass_u::String = "Q2"
  nDof_u::Int = ndim              # number of degree of freedom per node
  FunctionClass_p::String = "Q1"
  nDof_p::Int = 1                 # number of degree of freedom per node

  ne::Int = 8 # number of elements in the mesh for the ground truth

  camera_matrix = [[8*2048/7.07, 0.0, 2048/2] [0.0, 8*1536/5.3, 1536/2] [0.0, 0.0, 1.0]]'
  camera_pose = scale*[0 -0.5 2.75]'   # camera position in mm

  const_vel(r, h, ndim, FunctionClass_u, nDof_u, FunctionClass_p, nDof_p, FunctionClass_x, ne, camera_matrix, camera_pose)
  # bulk_vel(r, h, ndim, FunctionClass_u, nDof_u, FunctionClass_p, nDof_p, ne, camera_matrix, camera_pose)

end

main()

