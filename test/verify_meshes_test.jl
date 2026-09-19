using Test
using smearFEM

# The same squeeze-flow problem solved on a Tet and a Hex discretization must give the same
# answer. Neither `stokes_test.jl` (Hex only) nor `mesh_templates_test.jl` (mesh generation
# only) covers this, and it is the standing guard on the simplex quadrature: `_basis_tet`/
# `_basis_tri` define an equilateral/regular master simplex, so the weights must sum to *its*
# measure rather than the right-simplex values most published tables give. When that drifts,
# every Tet volume integral is off by 6√2 and every Tri surface integral by 2/√3.
#
# The β sweep is what makes both halves observable. At β→0 the boundary term β·b vanishes and
# only the Tet volume rule is exercised; at large β that surface term dominates, exercising the
# Tri rule (Tet meshes carry Tri faces, Hex meshes carry Quad).

const VM_R, VM_H      = 25.0, 40.0   # cylinder radius and height, mm
const VM_NDIM         = 3
const VM_NE           = 11
const VM_ETA          = 70.0
const VM_F_EXT        = 9.813e3      # kg·mm/s², not newtons
const VM_T_STEPS      = 1.0
# Relative, not absolute: step magnitudes vary with β, and Tet/Hex cannot agree better than
# their differing boundary faceting allows. A convergence sweep (2026-09-19) measured the
# genuine discretization gap at 0.79% (ne=4) falling to 0.23-0.31% (ne=10) — it converges, and
# vanishes entirely at β→0 where the surface term drops out. 2% therefore clears the real
# spread by ~2.5x while staying ~40x below the 87% signal of a simplex-quadrature mismatch.
# Hex is the denominator: for equal DOFs it is the more accurate discretization.
const VM_RTOL         = 0.02

# (β, sim_time). The β→0 case keeps the full 20 s because it is the only one that would catch
# Tet and Hex *drifting apart* over time. The two surface-dominated cases run short: a weight
# mismatch shows up as a ~88% discrepancy on the very first step, so extra steps cost minutes
# and add nothing.
const VM_CASES = ((1e-5, 20.0), (100.0, 5.0), (1e4, 5.0))

function _vm_solve(shape::Symbol, β::Float64, sim_time::Float64)
    F = -VM_F_EXT * ones(Float64, round(Int, sim_time / VM_T_STEPS))
    model, scene = def_problem(Cylinder(VM_R, VM_H), VM_NE, VM_ETA,
                               shape, 2, VM_NDIM,   # velocity mesh
                               shape, 1, 1,         # pressure mesh
                               shape, 2,            # geometry mesh
                               β, F, "force", "constant", sim_time, VM_T_STEPS,
                               mesh_path=joinpath(dirname(@__DIR__), "mesh_files"))
    conditions = Conditions(camera_matrix=get_camera_matrix(),
                            obj_pose=Float64.([-1.0 0.0 0.0 0.0;
                                                0.0 0.0 -1.0 20.0;
                                                0.0 -1.0 0.0 150.0;
                                                0.0 0.0 0.0 1.0]))
    return simulate(model, scene, conditions)[1]
end

@testset "Tet and Hex discretizations agree" begin
    # β→0 isolates the volume rule; the larger two make the surface rule dominant.
    for (β, sim_time) in VM_CASES
        @testset "β = $β" begin
            out_tet = _vm_solve(:Tet, β, sim_time)
            out_hex = _vm_solve(:Hex, β, sim_time)
            @test length(out_tet) == length(out_hex)
            for i in eachindex(out_tet)
                @test abs(out_tet[i] - out_hex[i]) <= VM_RTOL * abs(out_hex[i])
            end
        end
    end
end
