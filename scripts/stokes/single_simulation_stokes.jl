using smearFEM
using Dates

function main()
  scale = 100
  r::Float64 = 0.25*scale  # radius of the cylinder in mm
  h::Float64 = 0.5*scale  # height of the cylinder in mm
  ndim::Int = 3
  FunctionClass_u::String = "Q2"
  nDof_u::Int = ndim  # number of degree of freedom per node
  FunctionClass_p::String = "Q1"
  nDof_p::Int = 1  # number of degree of freedom per node

  camera_matrix = [[8*2048/7.07, 0.0, 2048/2] [0.0, 8*1536/5.3, 1536/2] [0.0, 0.0, 1.0]]'
  camera_pose = scale*[0 -0.25 2]'   # camera position in mm

  control::String = "force" # "force" or "velocity"
  viscosity_type::String = "constant" # "constant" or "bulk_viscosity"

  # simulation parameters for the ground truth
  sim_time::Float64 = 30.0# simulation time in seconds
  steps::Float64 = 30.0 # number of time steps
  t_steps::Float64 = sim_time/steps

  gt_β::Float64 = 100.0
  gt_η::Float64 = 40.0
  F::Float64 = 250000.0
  ne::Int = 4

  β_list = [100 1000 1e5]
  η_list = [30 40 50]
  i::Int = 1
  for gt_β in β_list
    for gt_η in η_list

      model, scene = def_problem(r, h, ne, gt_η, ndim, FunctionClass_u, nDof_u, FunctionClass_p, nDof_p, gt_β, F, control, viscosity_type, sim_time, t_steps)
      filepath = string("/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/sim_experiments/single_simulation/fem_runs/Stokes/",control,"_gt/exp_$i")
      write_sim_data(model, scene, camera_matrix, camera_pose, filepath)

      
      i += 1
    end
  end
end

fields = main()
