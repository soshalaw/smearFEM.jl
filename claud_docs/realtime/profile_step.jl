# Per-component timing of one force-control step on the benchmark case.
# usage: julia -O3 --threads=auto --project=<repo> profile_step.jl [constant|bulk_viscosity]
using smearFEM, LinearAlgebra, SparseArrays, Logging, Statistics

const S = smearFEM
vtype = length(ARGS) > 0 ? ARGS[1] : "constant"
mesh_path = joinpath(dirname(dirname(@__DIR__)), "mesh_files")

sim_time, dt = 4.0, 0.1
F = -9518.61 * ones(round(Int, sim_time / dt))
model, scene = def_problem(Cylinder(25.0, 40.0), 6, 100.0,
                           :Hex, 2, 3, :Hex, 1, 1, :Hex, 2,
                           100.0, F, "force", vtype, sim_time, dt; mesh_path=mesh_path)
cond = Conditions(camera_matrix=get_camera_matrix(), obj_pose=[150.0, 0.0, 20.0],
                  viewing_angles=[0.0, 30.0, 60.0])

# per-step wall times and refactor count from the solver's own @debug lines
struct StepLog <: AbstractLogger
    msgs::Vector{String}
end
Logging.min_enabled_level(::StepLog) = Logging.Debug
Logging.shouldlog(::StepLog, level, _module, args...) = _module === smearFEM
Logging.handle_message(l::StepLog, level, msg, args...; kw...) = push!(l.msgs, string(msg))
Logging.catch_exceptions(::StepLog) = true

simulate(model, scene, cond)                              # warm-up
logger = StepLog(String[])
t_total = @elapsed with_logger(logger) do
    simulate(model, scene, cond)
end
n = length(F)
nfact = count(m -> occursin(r", 1 steps since factorization", m), logger.msgs)
corr = [parse(Int, match(r"(\d+) corrections", m)[1]) for m in logger.msgs if occursin("corrections", m)]
println("simulate: $(round(1e3 * t_total / n, digits=1)) ms/step over $n steps, ",
        "$nfact factorizations, corrections per solve: mean $(round(mean(corr), digits=2)), max $(maximum(corr))")

# components, at the initial geometry
S.reset_model!(model)
cache = S.BasisFunctionCache(model)
t_sys = @elapsed sys = S._squeeze_system(model, scene)   # uncached build
t_asm = minimum(@elapsed(S._assemble_squeeze_system!(sys, model, cache, 100.0, 100.0)) for _ in 1:5)
t_lu = @elapsed sys.F = lu(sys.K)
t_relu = minimum(@elapsed(lu!(sys.F, sys.K)) for _ in 1:3)
x = rand(size(sys.K, 1)); X = rand(size(sys.K, 1), 2)
t_solve1 = minimum(@elapsed(sys.F \ x) for _ in 1:10)
t_solve2 = minimum(@elapsed(sys.F \ X) for _ in 1:10)
t_matvec = minimum(@elapsed(sys.K * x) for _ in 1:10)
L, U = sys.F.L, sys.F.U
println("n = $(size(sys.K, 1)), nnz(K) = $(nnz(sys.K)), nnz(L+U) = $(nnz(L) + nnz(U))")

X0 = model.mesh_x.NodeList
side = model.mesh_u.side_nodes
cm, op, h = get_camera_matrix(), [150.0, 0.0, 20.0], 40.0
dq = rand(size(X0)..., 2)
t_2d_grad = minimum(@elapsed(S._get_2D_data(X0, cm, op, h; BorderNodesList=side, GRAD=true, dqdθ=dq, angles=[0.0, 30.0, 60.0])) for _ in 1:3)
t_2d = minimum(@elapsed(S._get_2D_data(X0, cm, op, h; BorderNodesList=side, angles=[0.0, 30.0, 60.0])) for _ in 1:3)

ms(t) = round(1e3 * t, digits=2)
println("""
  _squeeze_system (uncached build)  : $(ms(t_sys)) ms
  assembly                          : $(ms(t_asm)) ms
  first lu (symbolic + numeric)     : $(ms(t_lu)) ms
  lu! refactor (numeric)            : $(ms(t_relu)) ms
  F \\ x, 1 rhs / 2 rhs              : $(ms(t_solve1)) / $(ms(t_solve2)) ms
  K * x                             : $(ms(t_matvec)) ms
  _get_2D_data GRAD=true / false    : $(ms(t_2d_grad)) / $(ms(t_2d)) ms   (3 views, $(length(side)) side nodes)""")
