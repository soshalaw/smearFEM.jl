using LinearAlgebra
using ProgressMeter
using SparseArrays
using Plots

using smearFEM

function main()
    x0 = 0
    x1 = 1
    y0 = 0
    y1 = 1
    z0 = 0
    z1 = 1
    H = z1 - z0
    ne = 8
    ndim = 3
    FunctionClass = "S2" #"Q1"
    nDof = ndim  # number of degree of freedom per node
    β = 1.0e-5
    Young = 40
    ν = 0.4
    μ_tp = -0.1
    μ_btm = 0

    q, model = simulate_single_tstep(x0, x1, y0, y1, z0, z1, ne, Young, ν, ndim, FunctionClass, nDof, β, μ_tp, μ_btm)
end

main()


