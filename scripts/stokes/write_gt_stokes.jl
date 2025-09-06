using smearFEM
using Dates
using Plots, Plots.PlotMeasures
using LaTeXStrings
using StatsPlots

function const_vel(r::Float64, h::Float64, ndim::Int, FunctionClass_u::String, nDof_u::Int, FunctionClass_p::String, nDof_p::Int, ne::Int, 
                    camera_matrix::AbstractMatrix{Float64}, camera_pose::AbstractMatrix{Float64}, β_list::Vector{Float64}, η_list::Vector{Float64})

  control::String = "force" # "force" or "velocity"
  viscosity_type::String = "constant" # "constant" or "bulk_viscosity"

  # simulation parameters for the ground truth
  sim_time::Float64 = 10.0# simulation time in seconds
  steps::Float64 = 100.0 # number of time steps
  t_steps::Float64 = sim_time/steps

  # F::Float64 = 150000000.0 force applied to the cylinder for β = 1e5
  # F::Float64 = 1500000.0 force applied to the cylinder for β = 1e3
  # F::Float64 = 150000.0 force applied to the cylinder for β = 1e2

  F_ext::Float64 = 1500.0
  F_::Vector{Float64} = -F_ext*ones(Float64, round(Int, (sim_time/t_steps))) # force applied to the cylinder in N

  # β_list = [100.0]
  # η_list = [50.0]
  i::Int = 1

  for gt_β::Float64 in β_list
    if gt_β < 100
      F = F_*100
    else
      F = F_*gt_β
    end
    for gt_η::Float64 in η_list

      model, scene = def_problem(r, h, ne, gt_η, ndim, FunctionClass_u, nDof_u, FunctionClass_p, nDof_p, "Q2", gt_β, F, control, viscosity_type, 
                    sim_time, t_steps)
      filepath = string("/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/sim_experiments/ground_truth/fem_runs/Stokes/$control/$viscosity_type/$i")
      write_sim_data(model, scene, camera_matrix, camera_pose, filepath)

      i += 1
    end
  end
end

function bulk_vel(r::Float64, h::Float64, ndim::Int, FunctionClass_u::String, nDof_u::Int, FunctionClass_p::String, nDof_p::Int, ne::Int, 
                    camera_matrix::AbstractMatrix{Float64}, camera_pose::AbstractMatrix{Float64}, β_list::Vector{Float64}, η_list::Vector{Float64})

  control::String = "force" # "force" or "velocity"
  viscosity_type::String = "bulk_viscosity" # "constant" or "bulk_viscosity"

  # simulation parameters for the ground truth
  sim_time::Float64 = 50.0
  steps::Float64 = 500.0
  t_steps::Float64 = sim_time/steps

  F_ext::Float64 = 200.0 # force applied to the cylinder in N
  F_::Vector{Float64} = -F_ext*ones(Float64, round(Int, (sim_time/t_steps))) # force applied to the cylinder in N

  i::Int = 1

  for gt_β::Float64 in β_list
    if gt_β <= 100
      F = F_*100
    else
      F = F_*gt_β
    end

    for gt_η::Float64 in η_list

      model, scene = def_problem(r, h, ne, gt_η, ndim, FunctionClass_u, nDof_u, FunctionClass_p, nDof_p, "Q2", gt_β, F, control, viscosity_type, 
              sim_time, t_steps)
      filepath = string("/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/sim_experiments/ground_truth/FEM/Stokes/$control/$viscosity_type/$i")
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

  scale = 50
  r::Float64 = 0.5*scale  # radius of the cylinder in mm
  h::Float64 = 1*scale  # height of the cylinder in mm
  ndim::Int = 3
  FunctionClass_u::String = "Q2"
  nDof_u::Int = ndim  # number of degree of freedom per node
  FunctionClass_p::String = "Q1"
  nDof_p::Int = 1  # number of degree of freedom per node

  ne::Int = 2 # number of elements in the mesh for the ground truth

  camera_matrix = [[8*2048/7.07, 0.0, 2048/2] [0.0, 8*1536/5.3, 1536/2] [0.0, 0.0, 1.0]]'
  camera_pose = scale*[0 -0.5 2.75]'   # camera position in mm

  # β_list = [10.0, 50.0, 100.0, 400.0, 700.0, 1e3, 4e3, 7e3, 1e4, 9e4, 1e5, 1e6]
  # η_list = [60.0]

  β_list = [10.0, 50.0, 100.0, 1e3]
  η_list = [60.0]
  
  const_vel(r, h, ndim, FunctionClass_u, nDof_u, FunctionClass_p, nDof_p, ne, camera_matrix, camera_pose, β_list, η_list)
  bulk_vel(r, h, ndim, FunctionClass_u, nDof_u, FunctionClass_p, nDof_p, ne, camera_matrix, camera_pose, β_list, η_list)

end

main()

