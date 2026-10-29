# Per-call overhead vs steady per-step cost of simulate: times runs of 2, 10 and 40 steps and
# fits t = a + b·n. Also reports the slowest single step of the 40-step run.
# usage: julia -O3 --threads=auto --project=<repo> step_cost.jl [constant|bulk_viscosity] [force|velocity]
using smearFEM, Printf, LinearAlgebra
haskey(ENV, "BLAS_THREADS") && BLAS.set_num_threads(parse(Int, ENV["BLAS_THREADS"]))
vtype = length(ARGS) > 0 ? ARGS[1] : "constant"
control = length(ARGS) > 1 ? ARGS[2] : "force"
mesh_path = joinpath(dirname(dirname(@__DIR__)), "mesh_files")
cparam = control == "force" ? -9518.61 : -1.0       # plate force [kg·mm/s²] or velocity [mm/s]
cond = Conditions(camera_matrix=get_camera_matrix(), obj_pose=[150.0, 0.0, 20.0],
                  viewing_angles=[0.0, 30.0, 60.0])

function run(nsteps)
    model, scene = def_problem(Cylinder(25.0, 40.0), 6, 100.0, :Hex, 2, 3, :Hex, 1, 1, :Hex, 2,
                               100.0, cparam * ones(nsteps), control, vtype, nsteps * 0.1, 0.1;
                               mesh_path=mesh_path)
    simulate(model, scene, cond)                    # warm-up, fills the per-model cache
    minimum(@elapsed(simulate(model, scene, cond)) for _ in 1:2)
end
ns = [2, 10, 40]
ts = run.(ns)
b = (ts[3] - ts[2]) / (ns[3] - ns[2])
a = ts[3] - b * ns[3]
@printf("%s/%s: per step %.1f ms, per call %.0f ms  (runs: %s ms)\n", control, vtype, 1e3b, 1e3a,
        join(round.(Int, 1e3 .* ts), " / "))
