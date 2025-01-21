using LinearAlgebra
using ProgressMeter
using SparseArrays
using Plots

using smearFEM

β = 1.0e5
Young = 40
ν = 0.4
μ_tp = -0.1
μ_btm = 0
nsub = 1
FunctionClass = "S2" #"Q1"
nDof = 3
filePath = "/home/soshala/SMEAR-PhD/smear-modules/smearFEM.jl/cylindergen"
ne = 1
ndim = 3

q_, model = simulate_single_tstep(Young, ν, FunctionClass, nDof, β, μ_tp, μ_btm)

NodeList_, IEN_list, q = eval_on_cylinder(model, q_, nsub)

write_vtk(filePath, "q", NodeList_, IEN_list, ne, ndim, q, ID=model.ID, FunctionClass="Q2")