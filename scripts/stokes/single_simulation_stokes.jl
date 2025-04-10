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

  endTime = 15
  steps = 15
  tSteps = endTime/steps

  β = 100
  η = 20
  Control = "force" # "force" or "displacement"
  dateTime = Dates.now()

  # filepath = string("/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/sim_experiments/single_simulation/fem_runs/Linear_Elasticity/",mode,"/",Control,"/",Date(dateTime),"/",Time(dateTime),"/")
  # filepath = string("/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/sim_experiments/single_simulation/fem_runs/Linear_Elasticity/",mode,"/",Control,"/",Date(dateTime),"/09:41:12.027/")
  filepath = string("/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/sim_experiments/single_simulation/fem_runs/Stokes/",Control,"_gt/2025-03-07/09:41:12.027")
  write_sim_data_stokes(r, h, ne, η, ndim, FunctionClass_u, nDof_u, FunctionClass_p, nDof_p, β, CameraMatrix, endTime, tSteps, Control,
                        filepath)
end

fields = main()