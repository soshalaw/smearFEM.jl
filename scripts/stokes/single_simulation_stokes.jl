using smearFEM
using Dates

function main()
  # test case 
  r = 0.5
  h = 1
  ne = 4
  ndim = 3
  FunctionClass_u = "Q2"
  nDof_u = ndim  # number of degree of freedom per node
  FunctionClass_p = "Q1"
  nDof_p = 1  # number of degree of freedom per node

  CameraMatrix = [[8*2048/7.07, 0.0, 2048/2] [0.0, 8*1536/5.3, 1536/2] [0.0, 0.0, 1.0]]'

  sim_time = 10.0
  steps = 10.0
  t_steps = sim_time/steps
  viscosity_type = "constant" # "constant" or "linear"

  β = 100.0
  η = 40.0
  F = 3.0
  control = "force" # "force" or "displacement"
  dateTime = Dates.now()

  model, scene = def_problem(r, h, ne, η, ndim, FunctionClass_u, nDof_u, FunctionClass_p, nDof_p, β, F, control, viscosity_type, sim_time, t_steps)
  
  # filepath = string("/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/sim_experiments/single_simulation/fem_runs/Linear_Elasticity/",mode,"/",Control,"/",Date(dateTime),"/",Time(dateTime),"/")
  # filepath = string("/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/sim_experiments/single_simulation/fem_runs/Linear_Elasticity/",mode,"/",Control,"/",Date(dateTime),"/09:41:12.027/")
  filepath = string("/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/sim_experiments/single_simulation/fem_runs/Stokes/",control,"_gt/2025-03-07/09:41:12.027")
  write_sim_data(model, scene, CameraMatrix, filepath)
end

fields = main()