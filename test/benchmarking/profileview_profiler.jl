"""
    profileview_profiler.jl
    
Uses ProfileView to generate a flamegraph of the Stokes solver execution.
Identifies bottlenecks by visualizing call stack and time spent in each function.
"""

using smearFEM
using Profile

# Test case parameters
r::Float64 = 25.0
h::Float64 = 40.0
ne = 6
ndim = 3
FunctionClass_u = "Q2"
nDof_u = ndim
FunctionClass_p = "Q1"
nDof_p = 1
FunctionClass_x = "Q2"

camera_matrix = [[8*2048/7.07, 0.0, 2048/2] [0.0, 8*1536/5.3, 1536/2] [0.0, 0.0, 1.0]]'
camera_pose = Float64.([-1.0 0.0 0.0 0.0; 0.0 0.0 -1.0 20.0; 0.0 -1.0 0.0 150; 0.0 0.0 0.0 1.0])

sim_time = 10.0
steps = 100.0
t_steps = sim_time/steps

β = 0.007
η = 40.0
F_ext::Float64 = 2
F::Vector{Float64} = -F_ext*ones(Float64, round(Int, (sim_time/t_steps)))

control = "force"
viscosity_type = "constant"
filepath = string("/home/soshala/tmp/profile_view/")

println("\n" * "="^80)
println("PROFILING WITH PROFILESVG (FLAMEGRAPH)")
println("="^80)
println("Cylinder: r=$(r) mm, h=$(h) mm")
println("Mesh elements: ne=$(ne)")
println("Will profile write_sim_data (includes full simulation)")
println("")

# Initialize model
println("Initializing model...")
model, scene = def_problem(r, h, ne, η, ndim, FunctionClass_u, nDof_u, FunctionClass_p, nDof_p, FunctionClass_x, β, F, control, viscosity_type, sim_time, t_steps)

conditions = Conditions(ANIMATE=false, WRITECONTOUR=false, RENDER=false, WRITEVTK=false, 
                       camera_matrix=camera_matrix, obj_pose=camera_pose, filepath=filepath)

println("Model initialized")
println("  - Velocity nodes: $(model.mesh_u.nNodes)")
println("  - Pressure nodes: $(model.mesh_p.nNodes)")
println("  - Elements: $(model.mesh_x.ne)")
println("")

println("Starting profiling...")
println("This will profile the full simulation with 100 timesteps.")
println("")

# Warm-up: compile/JIT code paths
println("Warming up (first call - may be slow due to JIT)...")
@time write_sim_data(model, scene, camera_matrix, camera_pose, filepath, ANIMATE=false, WRITECONTOUR=false, RENDER=false, WRITEVTK=false)

println("")
println("Warm-up complete. Now profiling with Profile...")
println("")

# Profile the second run
Profile.clear()
@profile begin
    write_sim_data(model, scene, camera_matrix, camera_pose, filepath, ANIMATE=false, WRITECONTOUR=false, RENDER=false, WRITEVTK=false)
end

println("Profiling complete!")
println("")

# Save profile to file for analysis
profile_file = "/home/soshala/tmp/profile_analysis.txt"
println("Saving profile results to: $profile_file")
println("")

open(profile_file, "w") do io
    # Print flat format with good depth
    println(io, "="^80)
    println(io, "PROFILING RESULTS (FLAT FORMAT)")
    println(io, "="^80)
    println(io, "")
    Profile.print(io, format=:flat, maxdepth=30, mincount=10)
end

println("Profile saved. Now displaying top bottlenecks...")
println("")
println("="^80)
println("KEY BOTTLENECKS (from flat profile)")
println("="^80)
println("")

# Print a brief summary
Profile.print(format=:flat, maxdepth=20, mincount=10)

rm(filepath, recursive=true, force=true)
println("\nCleaned up temporary files.")
