#!/usr/bin/env python3
"""Copy the package to <out_dir> and put lap timers into simulate's force loop.

usage: python3 instrument.py <repo> <out_dir>
then:  julia -O3 --threads=auto --project=<out_dir> instrument_run.jl [constant|bulk_viscosity]

Each marker is matched by the code text it follows, so this survives line-number churn; it
fails loudly if a marker is not found exactly once.
"""
import os, re, shutil, sys

repo, out = sys.argv[1], sys.argv[2]
shutil.rmtree(out, ignore_errors=True)
os.makedirs(out)
for f in ("Project.toml", "Manifest.toml"):
    shutil.copy(os.path.join(repo, f), out)
shutil.copytree(os.path.join(repo, "src"), os.path.join(out, "src"))
os.symlink(os.path.join(repo, "mesh_files"), os.path.join(out, "mesh_files"))

path = os.path.join(out, "src", "solver", "stokes_solver.jl")
src = open(path).read()
force = src.index('    sys = _cached_squeeze_system(mdl, scene)')
end = src.index('    if conditions.WRITEVTK')
head, loop, tail = src[:force], src[force:end], src[end:]

# (text the lap goes after, lap index); laps time the code since the previous lap
marks = [
    ("_assemble_squeeze_system!(sys, mdl, cache, η_t, β_t)\n", 1),
    ("iter == 1 && _refactor!(sys)\n", 2),
    ("_scale_unknowns!(_push_history!(z_hist, z)[1], sys, η_t^(1 - a), η_t^(-a))\n", 3),
    ("nodeList_cached = nodeList_cached + motion\n", 4),
    ("dmdθ_out = @views cat(dNodeList_dη,dNodeList_dβ,dims=3)\n", 5),
    ("GRAD=true, dqdθ=dmdθ_out, angles=rot_angle_cached, FIT=conditions.WRITECONTOUR)\n", 6),
    ('next!(pr, showvalues = [(:iterations,iter),(:time,string(t," seconds"))])\n', 7),
]
for text, k in marks:
    n = loop.count(text)
    if n != 1:
        sys.exit(f"marker {k} found {n} times: {text!r}")
    loop = loop.replace(text, text + f"        _lap!({k})\n")
loop = loop.replace("    for t in time\n", "    for t in time\n        LAST[] = time_ns(); STEP0[] = LAST[]; LASTB[] = Base.gc_bytes(); T_START .= TIMES\n", 1)
loop = loop.replace("_lap!(7)\n", "_lap!(7); push!(STEPS, (time_ns() - STEP0[]) / 1e6); push!(STEP_SECT, TIMES .- T_START)\n", 1)
open(path, "w").write(head + loop + tail)

with open(os.path.join(out, "src", "solver", "squeeze_system.jl"), "a") as fh:
    fh.write("\nconst TIMES = zeros(8)\nconst ALLOCS = zeros(8)\nconst LAST = Ref(UInt64(0))\nconst LASTB = Ref(Int64(0))\n"
             "const STEP0 = Ref(UInt64(0))\nconst STEPS = Float64[]\nconst T_START = zeros(8)\nconst STEP_SECT = Vector{Vector{Float64}}()\n"
             "_lap!(k) = (now = time_ns(); b = Base.gc_bytes(); TIMES[k] += (now - LAST[]) / 1e6; ALLOCS[k] += (b - LASTB[]) / 1e6;"
             " LAST[] = now; LASTB[] = b; nothing)\n")
print("instrumented copy at", out)
