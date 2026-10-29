# Per-section timing of simulate's force loop; run against a copy made by instrument.py.
using smearFEM, LinearAlgebra, Logging
const S = smearFEM
haskey(ENV, "BLAS_THREADS") && BLAS.set_num_threads(parse(Int, ENV["BLAS_THREADS"]))
struct L <: AbstractLogger; m::Vector{String}; end
Logging.min_enabled_level(::L) = Logging.Debug
Logging.shouldlog(::L, level, _module, args...) = _module === smearFEM
Logging.handle_message(l::L, level, msg, args...; kw...) = push!(l.m, string(msg))
Logging.catch_exceptions(::L) = true
mesh_path = "/home/soshala/SMEAR-PhD/smear-modules/smearFEM.jl/mesh_files"
vtype = length(ARGS) > 0 ? ARGS[1] : "constant"
control = length(ARGS) > 1 ? ARGS[2] : "force"
nsteps = length(ARGS) > 2 ? parse(Int, ARGS[3]) : 40
F = (control == "force" ? -9518.61 : -1.0) * ones(nsteps)
model, scene = def_problem(Cylinder(25.0, 40.0), 6, 100.0, :Hex, 2, 3, :Hex, 1, 1, :Hex, 2,
                           100.0, F, control, vtype, nsteps * 0.1, 0.1; mesh_path=mesh_path)
cond = Conditions(camera_matrix=get_camera_matrix(), obj_pose=[150.0, 0.0, 20.0], viewing_angles=[0.0, 30.0, 60.0])
simulate(model, scene, cond)
S.TIMES .= 0; S.ALLOCS .= 0; empty!(S.STEPS); empty!(S.STEP_SECT)
lg = L(String[])
tm = @timed with_logger(() -> simulate(model, scene, cond), lg)
t = tm.time
names = ["assembly", "refactor(step1)", "z+w solve", "state update", "sens update", "2D data", "push+progress"]
for (k, nm) in enumerate(names)
    println(rpad(nm, 16), round(S.TIMES[k] / nsteps, digits=2), " ms/step, ", round(S.ALLOCS[k] / nsteps, digits=2), " MB/step")
end
println("loop sum ", round(sum(S.TIMES) / nsteps, digits=2), " ms/step; total simulate ", round(1e3t / nsteps, digits=2), " ms/step")
st = S.STEPS
using Statistics
println("per step: step 1 ", round(st[1], digits=1), ", steps 2-5 ", round.(st[2:5], digits=1),
        ", steps ≥ 6: median ", round(median(st[6:end]), digits=1), " max ", round(maximum(st[6:end]), digits=1),
        " at step ", 5 + argmax(st[6:end]))
slow = findall(>(20), st)
println("steps over 20 ms: ", slow)
corr = [parse(Int, match(r"(\d+) corrections", m)[1]) for m in lg.m if occursin("corrections", m)]
waits = [parse(Float64, match(r"waited ([\d.]+) ms", m)[1]) for m in lg.m if occursin("swapped in", m)]
println("corrections/solve mean ", round(mean(corr), digits=2), ", refreshes started ", count(m -> occursin("started in the background", m), lg.m),
        ", swaps ", length(waits), " (waits ms: ", waits, "), synchronous refactors ", count(m -> occursin("synchronous", m), lg.m),
        ", GC ", round(1e3 * tm.gctime, digits=0), " ms in total")
rounds = [parse(Int, match(r"(\d+) corrections", m)[1]) for m in lg.m if occursin("corrections", m)]
println("step ms / correction rounds:")
println(join(["$k:$(round(st[k], digits=1))/$(k <= length(rounds) ? rounds[k] : -1)" for k in eachindex(st)], "  "))
# per step: correction rounds and whether it refactorized (the log has one "initial residual"
# line per step, emitted before any refactorization of that step's solve)
steplog = Vector{Vector{String}}()
for m in lg.m
    occursin("initial residual", m) && push!(steplog, String[])
    isempty(steplog) || push!(steplog[end], m)
end
refac = [any(occursin("refactorization", x) for x in st_) for st_ in steplog]
println("refactor steps (within a solve): ", findall(refac))
println("refactorization times (ms): ", [parse(Float64, match(r", ([\d.]+) ms", m)[1]) for m in lg.m if occursin("refactorization after", m)])
slowest = sortperm(st, rev=true)[1:min(8, end)]
println("slowest steps, section ms [asm, refactor1, solve, state, sens, 2D, push]:")
for k in sort(slowest)
    println("  step $k: $(round(st[k], digits=1)) ms = ", round.(S.STEP_SECT[k][1:7], digits=1))
end
