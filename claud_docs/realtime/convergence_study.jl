# How good are the initial guess and a stale factorization along the benchmark trajectory?
# Records K_t and the exact solution of every step, then reports
#   - relative initial residual for extrapolation of order 0, 1, 2;
#   - contraction factor of one iterative-refinement step with the factorization of step t-s;
#   - GMRES vs iterative refinement iteration counts with a stale factorization.
# usage: julia -O3 --project=<repo> convergence_study.jl
using smearFEM, LinearAlgebra, SparseArrays, Printf, Statistics
const S = smearFEM
mesh_path = joinpath(dirname(dirname(@__DIR__)), "mesh_files")

nsteps, dt = 40, 0.1
F = -9518.61 * ones(nsteps)
model, scene = def_problem(Cylinder(25.0, 40.0), 6, 100.0, :Hex, 2, 3, :Hex, 1, 1, :Hex, 2,
                           100.0, F, "force", "constant", nsteps * dt, dt; mesh_path=mesh_path)
S.reset_model!(model)
cache = S.BasisFunctionCache(model)
sys = S._squeeze_system(model, scene)
n, nf = size(sys.K, 1), sys.nf
C_uc = sparse(scene.C_uc); q_top = vec(Matrix(scene.q_d[:top])); ID = model.mesh_u.ID
X = copy(model.mesh_x.NodeList)

Ks, sols = SparseMatrixCSC{Float64,Int}[], Vector{Float64}[]
for t in 1:nsteps
    model.mesh_x.NodeList = X
    S._assemble_squeeze_system!(sys, model, cache, 100.0, 100.0)
    r = zeros(n); r[end] = F[t]
    sol = lu(sys.K) \ r
    push!(Ks, copy(sys.K)); push!(sols, sol)
    q = C_uc * sol[1:nf] + sol[end] * q_top
    global X = X + dt * vcat(q[ID[1, :]]', q[ID[2, :]]', q[ID[3, :]]')
end
rhs(t) = (r = zeros(n); r[end] = F[t]; r)
relres(t, x) = norm(rhs(t) - Ks[t] * x) / norm(rhs(t))

println("initial residual, mean over steps 4..40:")
for (name, g) in (("order 0", t -> sols[t-1]),
                  ("order 1", t -> 2sols[t-1] - sols[t-2]),
                  ("order 2", t -> 3sols[t-1] - 3sols[t-2] + sols[t-3]))
    @printf("  %-8s %.2e\n", name, mean(relres(t, g(t)) for t in 4:nsteps))
end

# one IR step from the order-1 guess, factor from step t-s
println("contraction of one refinement step vs factorization age s (t = 20, 40):")
for s in (1, 2, 5, 10, 19)
    ρ = map((20, 40)) do t
        Fk = lu(Ks[t-s]); x = 2sols[t-1] - sols[t-2]
        r0 = rhs(t) - Ks[t] * x
        x1 = x + Fk \ r0
        norm(rhs(t) - Ks[t] * x1) / norm(r0)
    end
    @printf("  s = %2d   ρ = %.2e, %.2e\n", s, ρ...)
end

# right-preconditioned GMRES (no restart) vs iterative refinement to a relative residual tol
function ir_count(K, Fk, b, x; tol)
    r = b - K * x; k = 0
    while norm(r) > tol * norm(b) && k < 50
        x = x + Fk \ r; r = b - K * x; k += 1
    end
    k
end
function gmres_count(K, Fk, b, x; tol)
    r = b - K * x; β = norm(r); k = 0
    β <= tol * norm(b) && return 0
    V = [r / β]; H = zeros(51, 50)
    for j in 1:50
        w = K * (Fk \ V[j])
        for i in 1:j
            H[i, j] = dot(V[i], w); w -= H[i, j] * V[i]
        end
        H[j+1, j] = norm(w); push!(V, w / H[j+1, j])
        e1 = zeros(j + 1); e1[1] = β
        y = H[1:j+1, 1:j] \ e1
        norm(H[1:j+1, 1:j] * y - e1) <= tol * norm(b) && return j
    end
    50
end
println("iterations to tol from the order-1 guess (IR / GMRES), t = 40:")
for tol in (1e-8, 1e-10), s in (1, 5, 19, 39)
    t = 40; Fk = lu(Ks[t-s]); x = 2sols[t-1] - sols[t-2]
    @printf("  tol %.0e  s = %2d   IR %d   GMRES %d\n", tol, s,
            ir_count(Ks[t], Fk, rhs(t), x; tol=tol), gmres_count(Ks[t], Fk, rhs(t), x; tol=tol))
end
