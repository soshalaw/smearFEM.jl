using LinearAlgebra
using SparseArrays

# The object squeezed between two rigid horizontal plates: the bottom plate holds u_z = 0, the
# top plate moves every top node by one plate displacement — the last unknown, μ — and both
# plates carry Navier slip β on the tangential components. Under force control μ is solved for
# with the plate force on its row; under displacement control its row and column become the
# identity and the column moves to the right-hand side.

# Relative-residual tolerances of the refined solves. The state accumulates into the geometry
# over the steps; the sensitivities are first-order information only.
const STATE_TOL = 1e-9
const SENS_TOL = 1e-8

mutable struct SqueezePlates
    sys::ElasticSystem
    plate_col::Vector{Int}       # nzval range of the plate unknown's column
    plate_row::Vector{Int}       # nzval indices of its row outside that column
end

# Unknowns: every component of every node, except z on the bottom plate (eliminated) and z on
# the top plate (all mapped to the plate unknown, the last one). Plate nodes come from the face
# connectivity: the meshes' `top_nodes`/`bottom_nodes` leave out the nodes shared with the side.
function _squeeze_plates(mesh::AbstractMeshgrid)
    top = unique(vec(mesh.IEN_top))
    on_plate = falses(mesh.nNodes)
    on_plate[top] .= true
    on_plate[vec(mesh.IEN_bottom)] .= true
    dof = zeros(Int, 3, mesh.nNodes)
    k = 0
    for a in 1:mesh.nNodes, i in 1:3
        (i == 3 && on_plate[a]) || (dof[i, a] = (k += 1))
    end
    dof[3, top] .= k + 1
    sys = _elastic_system(mesh, dof; slip_faces=[mesh.IEN_top, mesh.IEN_bottom])
    K, n = sys.K, k + 1
    row = [p for c in 1:(n - 1) for p in K.colptr[c]:(K.colptr[c+1] - 1) if K.rowval[p] == n]
    return SqueezePlates(sys, collect(K.colptr[n]:(K.colptr[n+1] - 1)), row)
end

# The pattern, gather lists and symbolic factorization depend only on the mesh topology, so
# each mesh keeps one system across `simulate` calls (a fit calls `simulate` many times).
const _SQUEEZE_PLATES = WeakKeyDict{AbstractMeshgrid,SqueezePlates}()

_cached_squeeze_plates(mesh::AbstractMeshgrid) = get!(() -> _squeeze_plates(mesh), _SQUEEZE_PLATES, mesh)

# Displacement control: the plate displacement is given, so its row and column of K become the
# identity and the column moves to the right-hand side, written to `rhs` for a unit plate
# displacement. The pattern is unchanged, so both control modes share the factorization.
function _prescribe_plate!(sq::SqueezePlates, rhs::Vector{Float64})
    K = sq.sys.K
    n = size(K, 1)
    fill!(rhs, 0.0)
    for p in sq.plate_col
        r = K.rowval[p]
        r != n && (rhs[r] = -K.nzval[p])
        K.nzval[p] = r == n ? 1.0 : 0.0
    end
    K.nzval[sq.plate_row] .= 0.0
    rhs[n] = 1.0
    return rhs
end

"""
    simulate(mdl::LinearElasticity, cParam; control="force", β=0.0, sensitivities=false)

Squeeze the model's mesh between two rigid horizontal plates, one load increment per entry of
`cParam`: each increment is solved on the current geometry (small strain), and the mesh is then
moved by it (updated Lagrangian). The bottom plate is fixed in z, the top plate moves rigidly
in z, and both carry Navier slip `β` on the tangential components. The mesh starts from its
initial state (`reset_model!`) and is left at the final one.

One step costs one assembly and, typically, a refinement round or two against a factorization
reused from earlier steps (refactorized when it stops converging fast).

# Arguments
- `mdl::LinearElasticity`: Model; needs a 3D mesh with `IEN_top` and `IEN_bottom` (cylinder
  or cuboid).
- `cParam::Vector{Float64}`: Per-step increment of the top-plate force (kg·mm/s², the z
  component; negative compresses) under force control, or of the top-plate displacement (mm)
  under displacement control.
- `control::String`: `"force"` or `"displacement"`.
- `β::Real`: Navier slip coefficient of both plates, in kPa/mm. `β = 0` leaves the horizontal
  rigid-body modes free, so the stiffness is singular; use a small positive value for
  frictionless plates.
- `sensitivities::Bool`: Also return the nodal sensitivities to (λ, μ, β), accumulated over the
  steps with the geometry held fixed within each step (the approximation the Stokes solver
  uses too). Exact at the first step; against finite differences of the full run, about 0.5 %
  off after 0.6 % compression and 3 % after 12.5 % (25 steps, test cylinder).

# Returns
A named tuple:
- `plate::Vector{Float64}`: Per step, the top-plate displacement increment (force control) or
  the plate force increment it took (displacement control).
- `displacement::Vector{Matrix{Float64}}`: Per step, the nodal displacement increments (3 × nNodes).
- `pos3D::Vector{Matrix{Float64}}`: Nodal positions, the initial state first.
- `grad3D::Vector{Array{Float64,3}}`: Per state, ∂positions/∂(λ, μ, β) as (3, nNodes, 3); empty
  without `sensitivities`.
"""
function simulate(mdl::LinearElasticity, cParam::Vector{Float64}; control::String="force",
                  β::Real=0.0, sensitivities::Bool=false)
    control in ("force", "displacement") || throw(ArgumentError("control must be \"force\" or \"displacement\", got \"$control\""))
    reset_model!(mdl)
    mesh = mdl.mesh
    λ, μ = mdl.λ, mdl.μ
    IEN::Matrix{Int} = mesh.IEN
    X::Matrix{Float64} = mesh.NodeList
    displacement_control = control == "displacement"

    sq = _cached_squeeze_plates(mesh)
    sys = sq.sys
    n = size(sys.K, 1)
    col = zeros(length(sq.plate_col))     # the plate column before it is prescribed: its reaction

    # Solve for the response to a unit plate force (or displacement) and scale by the control
    # value: the unit response depends on the geometry alone, so it extrapolates well even when
    # a measured control signal is noisy.
    rhs = zeros(n)
    rhs[n] = 1.0
    z = zeros(n)
    z_hist = [zeros(n) for _ in 1:4]
    # ∂z/∂λ and ∂z/∂β; ∂z/∂μ follows from homogeneity
    w = [zeros(n) for _ in 1:2]
    w_hist = [[zeros(n) for _ in 1:4] for _ in 1:2]
    rhs_w = zeros(n)

    u = zeros(3, mesh.nNodes)
    du = [zeros(3, mesh.nNodes) for _ in 1:3]
    dX = zeros(3, mesh.nNodes, 3)
    plate = Float64[]
    displacement = Matrix{Float64}[]
    pos3D = Matrix{Float64}[X]
    grad3D = sensitivities ? Array{Float64,3}[copy(dX)] : Array{Float64,3}[]

    for (step, c) in enumerate(cParam)
        _assemble!(sys, X, IEN, λ, μ, β; sensitivities=sensitivities)
        if displacement_control
            copyto!(col, view(sys.K.nzval, sq.plate_col))
            _prescribe_plate!(sq, rhs)
        end
        # a factorization left over from an earlier call would make results depend on it
        step == 1 && _refactor!(sys)

        k = min(step - 1, 4)
        _refined_solve!(sys, _extrapolate!(z, z_hist, k), rhs; tol=STATE_TOL)
        _push_history!(z_hist, z)

        if sensitivities
            # K ∂z/∂θ = -(∂K/∂θ) z; under displacement control the plate row stays fixed
            for (j, (wi, hist)) in enumerate(zip(w, w_hist))
                if j == 1
                    _symv!(rhs_w, sys.K, sys.nz_L, z, -1.0)     # ∂K/∂λ = L
                else
                    mul!(rhs_w, sys.dKdβ, z, -1.0, 0.0)          # ∂K/∂β = b
                end
                displacement_control && (rhs_w[n] = 0.0)
                _refined_solve!(sys, _extrapolate!(wi, hist, k), rhs_w; tol=SENS_TOL)
                _push_history!(hist, wi)
            end
        end

        _nodal!(u, z, sys.dof)
        u .*= c
        push!(plate, displacement_control ? c * dot(col, view(z, view(sys.K.rowval, sq.plate_col))) : c * z[n])
        push!(displacement, copy(u))
        X = X + u                 # a new array: the mesh's initial state may alias NodeList
        push!(pos3D, X)

        if sensitivities
            # K(sθ) = s K(θ), so scaling every parameter by s scales the unit response by 1/s
            # under force control and leaves it unchanged under displacement control:
            # λ ∂z/∂λ + μ ∂z/∂μ + β ∂z/∂β = -z or 0. That gives ∂z/∂μ without a third solve.
            _nodal!(du[1], w[1], sys.dof)
            _nodal!(du[3], w[2], sys.dof)
            _nodal!(du[2], z, sys.dof)
            s = displacement_control ? 0.0 : 1.0
            @. du[2] = -(s * du[2] + λ * du[1] + β * du[3]) / μ
            for j in 1:3
                @views dX[:, :, j] .+= c .* du[j]
            end
            push!(grad3D, copy(dX))
        end
    end
    mesh.NodeList = X
    return (plate=plate, displacement=displacement, pos3D=pos3D, grad3D=grad3D)
end
