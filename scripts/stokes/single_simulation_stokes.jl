using smearFEM
using Dates
using Plots, Plots.PlotMeasures
using LaTeXStrings
using StatsPlots

function const_vel(r::Float64, h::Float64, ndim::Int, FunctionClass_u::String, nDof_u::Int, FunctionClass_p::String, nDof_p::Int, ne::Int, 
                    camera_matrix::AbstractMatrix{Float64}, camera_pose::AbstractMatrix{Float64})

  control::String = "force" # "force" or "velocity"
  viscosity_type::String = "constant" # "constant" or "bulk_viscosity"

  # simulation parameters for the ground truth
  sim_time::Float64 = 30.0# simulation time in seconds
  steps::Float64 = 30.0 # number of time steps
  t_steps::Float64 = sim_time/steps

  gt_β::Float64 = 100.0
  gt_η::Float64 = 40.0

  # F::Float64 = 150000000.0 force applied to the cylinder for β = 1e5
  # F::Float64 = 1500000.0 force applied to the cylinder for β = 1e3
  # F::Float64 = 150000.0 force applied to the cylinder for β = 1e2

  F_ext::Float64 = 1500.0
  F_::Vector{Float64} = -F_ext*ones(Float64, round(Int, (sim_time/t_steps))) # force applied to the cylinder in N

  β_list = [100.0 1000.0 1e5]
  η_list = [30.0 40.0 50.0]
  i::Int = 1

  for gt_β::Float64 in β_list
    F = F_*gt_β
    for gt_η::Float64 in η_list

      model, scene = def_problem(r, h, ne, gt_η, ndim, FunctionClass_u, nDof_u, FunctionClass_p, nDof_p, gt_β, F, control, viscosity_type, sim_time, t_steps)
      filepath = string("/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/sim_experiments/single_simulation/fem_runs/Stokes/",control,"_gt/",viscosity_type,"/exp_$i")
      write_sim_data(model, scene, camera_matrix, camera_pose, filepath)

      i += 1
    end
  end
end

function bulk_vel(r::Float64, h::Float64, ndim::Int, FunctionClass_u::String, nDof_u::Int, FunctionClass_p::String, nDof_p::Int, ne::Int, 
                    camera_matrix::AbstractMatrix{Float64}, camera_pose::AbstractMatrix{Float64})

  control::String = "force" # "force" or "velocity"
  viscosity_type::String = "bulk_viscosity" # "constant" or "bulk_viscosity"

  # simulation parameters for the ground truth
  sim_time::Float64 = 90.0
  steps::Float64 = 90.0
  t_steps::Float64 = sim_time/steps

  F_ext::Float64 = 50000.0 # force applied to the cylinder in N
  F::Vector{Float64} = -F_ext*ones(Float64, round(Int, (sim_time/t_steps))) # force applied to the cylinder in N

  β_list = [100.0]
  η_list = [60.0]
  i::Int = 1

  for gt_β::Float64 in β_list
    # F = F_*gt_β
    for gt_η::Float64 in η_list

      model, scene = def_problem(r, h, ne, gt_η, ndim, FunctionClass_u, nDof_u, FunctionClass_p, nDof_p, gt_β, F, control, viscosity_type, sim_time, t_steps)
      filepath = string("/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/sim_experiments/single_simulation/fem_runs/Stokes/",control,"_gt/",viscosity_type,"/exp_$i")
      write_sim_data(model, scene, camera_matrix, camera_pose, filepath)

      set_file(string(filepath,"/Results/plots"))
      plt_η = set_plot(22)
      t_windows = collect(range(start=t_steps, stop=sim_time, step=t_steps))
      Plots.plot!(model.η, t_windows, label=L"\eta(t)", dpi=400)
      Plots.xlabel!("Time (s)")
      Plots.ylabel!(L"η")
      Plots.savefig(string(filepath,"/Results/plots/η.pdf"))

      i += 1
    end
  end
end

function main()

  scale = 100
  r::Float64 = 0.25*scale  # radius of the cylinder in mm
  h::Float64 = 0.5*scale  # height of the cylinder in mm
  ndim::Int = 3
  FunctionClass_u::String = "Q2"
  nDof_u::Int = ndim  # number of degree of freedom per node
  FunctionClass_p::String = "Q1"
  nDof_p::Int = 1  # number of degree of freedom per node

  ne::Int = 2 # number of elements in the mesh for the ground truth

  camera_matrix = [[8*2048/7.07, 0.0, 2048/2] [0.0, 8*1536/5.3, 1536/2] [0.0, 0.0, 1.0]]'
  camera_pose = scale*[0 -0.25 2]'   # camera position in mm

  # const_vel(r, h, ndim, FunctionClass_u, nDof_u, FunctionClass_p, nDof_p, ne, camera_matrix, camera_pose)
  bulk_vel(r, h, ndim, FunctionClass_u, nDof_u, FunctionClass_p, nDof_p, ne, camera_matrix, camera_pose)

end

main()

