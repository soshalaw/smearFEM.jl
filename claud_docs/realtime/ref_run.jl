# usage: julia --project=<repo> ref_run.jl <out.jls> [constant|bulk_viscosity] [force|velocity]
# force: F = -9518.61 kg·mm/s² on the top plate; velocity: plate velocity -1 mm/s.
using smearFEM, Serialization, LinearAlgebra

out = ARGS[1]
vtype = length(ARGS) > 1 ? ARGS[2] : "constant"
control = length(ARGS) > 2 ? ARGS[3] : "force"
mesh_path = "/home/soshala/SMEAR-PhD/smear-modules/smearFEM.jl/mesh_files"

sim_time, dt = 4.0, 0.1
F = (control == "force" ? -9518.61 : -1.0) * ones(round(Int, sim_time / dt))
model, scene = def_problem(Cylinder(25.0, 40.0), 6, 100.0,
                           :Hex, 2, 3, :Hex, 1, 1, :Hex, 2,
                           100.0, F, control, vtype, sim_time, dt; mesh_path=mesh_path)
cond = Conditions(camera_matrix=get_camera_matrix(), obj_pose=[150.0, 0.0, 20.0],
                  viewing_angles=[0.0, 30.0, 60.0])

simulate(model, scene, cond)   # warm-up
t = @elapsed res = simulate(model, scene, cond)
println("$control/$vtype: $(round(1e3 * t / length(F), digits=1)) ms/step")
serialize(out, res)
