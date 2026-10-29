# Sampling profile of simulate (force control) on the benchmark case, flat by function.
# usage: julia -O3 --threads=auto --project=<repo> profile_flat.jl
using smearFEM, Profile
mesh_path = joinpath(dirname(dirname(@__DIR__)), "mesh_files")
F = -9518.61 * ones(40)
model, scene = def_problem(Cylinder(25.0, 40.0), 6, 100.0, :Hex, 2, 3, :Hex, 1, 1, :Hex, 2,
                           100.0, F, "force", "constant", 4.0, 0.1; mesh_path=mesh_path)
cond = Conditions(camera_matrix=get_camera_matrix(), obj_pose=[150.0, 0.0, 20.0],
                  viewing_angles=[0.0, 30.0, 60.0])
simulate(model, scene, cond)
Profile.clear()
Profile.init(n=10^7, delay=0.0005)
@profile simulate(model, scene, cond)
Profile.print(format=:flat, sortedby=:count, mincount=40, noisefloor=1.0, threads=1, C=true, maxdepth=200)
