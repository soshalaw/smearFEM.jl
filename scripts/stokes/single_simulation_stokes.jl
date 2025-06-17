using smearFEM
using Dates

function main()
  # test case 
  scale = 100
  r = 0.2*scale  # radius of the cylinder in mm
  h = 0.5*scale  # height of the cylinder in mm
  ne = 10  # number of elements in each direction
  ndim = 3
  FunctionClass_u = "Q2"
  nDof_u = ndim  # number of degree of freedom per node
  FunctionClass_p = "Q1"
  nDof_p = 1  # number of degree of freedom per node

  camera_matrix = [[8*2048/7.07, 0.0, 2048/2] [0.0, 8*1536/5.3, 1536/2] [0.0, 0.0, 1.0]]'
  camera_pose = scale*[0 -0.25 2]'   # camera position in mm

  sim_time = 60.0 # simulation time in seconds
  steps = 120.0 # number of time steps 
  t_steps = sim_time/steps
  viscosity_type = "constant" # "constant" or "linear"

  β = 100.0
  η = 40.0
  F = 30000.0
  control = "force" # "force" or "displacement"
  dateTime = Dates.now()

  model, scene = def_problem(r, h, ne, η, ndim, FunctionClass_u, nDof_u, FunctionClass_p, nDof_p, β, F, control, viscosity_type, sim_time, t_steps)
  
  # filepath = string("/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/sim_experiments/single_simulation/fem_runs/Linear_Elasticity/",mode,"/",Control,"/",Date(dateTime),"/",Time(dateTime),"/")
  # filepath = string("/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/sim_experiments/single_simulation/fem_runs/Linear_Elasticity/",mode,"/",Control,"/",Date(dateTime),"/09:41:12.027/")
  filepath = string("/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/sim_experiments/single_simulation/fem_runs/Stokes/",control,"_gt")
  write_sim_data(model, scene, camera_matrix, camera_pose, filepath)
end

fields = main()
