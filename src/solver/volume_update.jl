using LinearAlgebra
using SparseArrays

# Volume update (paper 2, at each window start): move the FE mesh so its surface meets the
# observed surface, the interior following a compressible pseudo-elastic solve on the current
# geometry. Boundary conditions: a prescribed normal displacement u·n = d at observed and
# rigid-contact surface nodes (tangential motion free), u_z = 0 on the table, traction-free
# elsewhere.
#
# A constrained node is rotated into a frame whose first axes span its constraint normals, so
# each constraint prescribes one unknown, and a prescribed unknown becomes an identity row. Both
# keep the sparsity pattern of the unconstrained stiffness — every node-pair block of it is a
# full 3×3 — so one symbolic analysis serves every window, whichever nodes are observed.

mutable struct VolumeSystem
    sys::ElasticSystem
    blocks::Matrix{Int}          # (9, nblocks) nzval index of each entry of a node-pair block, column-major
    block_nodes::Matrix{Int}     # (2, nblocks) row node and column node of each block
    mirror::Vector{Int}          # nzval index of each entry's transpose
    diag::Vector{Int}            # nzval index of each diagonal entry
end

function _volume_system(mesh::AbstractMeshgrid)
    nN = mesh.nNodes
    sys = _elastic_system(mesh, reshape(collect(1:3nN), 3, nN); slip_faces=[mesh.IEN_bottom])
    K = sys.K
    blocks = Int[]; block_nodes = Int[]
    for b in 1:nN, p in K.colptr[3b-2]:(K.colptr[3b-1] - 1)
        r = K.rowval[p]
        r % 3 == 1 || continue
        a = div(r, 3) + 1
        append!(blocks, (_nz_index(K, 3(a-1) + i, 3(b-1) + j) for j in 1:3 for i in 1:3))
        append!(block_nodes, (a, b))
    end
    mirror = [_nz_index(K, c, K.rowval[p]) for c in 1:size(K, 2) for p in K.colptr[c]:(K.colptr[c+1] - 1)]
    diag = [_nz_index(K, k, k) for k in 1:size(K, 1)]
    return VolumeSystem(sys, reshape(blocks, 9, :), reshape(block_nodes, 2, :), mirror, diag)
end

const _VOLUME_SYSTEMS = WeakKeyDict{AbstractMeshgrid,VolumeSystem}()

_cached_volume_system(mesh::AbstractMeshgrid) = get!(() -> _volume_system(mesh), _VOLUME_SYSTEMS, mesh)

# Orthonormal frame Q (columns) of a node whose first m axes span its constraint normals, and the
# prescribed displacement along each of them: u·n_k = d_k for every constraint k, u = Q ũ, gives
# ũ_1..ũ_m. Returns m. A constraint (numerically) parallel to earlier ones is dropped.
function _constraint_frame!(Q::AbstractMatrix{Float64}, vals::AbstractVector{Float64},
                            constraints::Vector{Tuple{Vector{Float64},Float64}})
    m = 0
    for (n, d) in constraints
        m == 3 && break
        v, dv = copy(n), d
        for l in 1:m
            c = dot(n, view(Q, :, l))
            v .-= c .* view(Q, :, l)
            dv -= c * vals[l]
        end
        nv = norm(v)
        if nv < 1e-6
            @warn "dropping a constraint parallel to another on the same node (residual $dv mm)"
            continue
        end
        m += 1
        Q[:, m] .= v ./ nv
        vals[m] = dv / nv
    end
    if m == 1
        e = zeros(3)
        e[argmin(abs.(view(Q, :, 1)))] = 1.0
        Q[:, 2] .= normalize(cross(Q[:, 1], e))
    end
    m <= 2 && (Q[:, 3] .= cross(Q[:, 1], Q[:, 2]))
    return m
end

# K ← R'KR for the block-diagonal R of the node frames (identity where `rotated` is false).
function _rotate_blocks!(nz::Vector{Float64}, vs::VolumeSystem, Q::Array{Float64,3}, rotated::BitVector)
    B, T = zeros(3, 3), zeros(3, 3)
    @inbounds for k in axes(vs.blocks, 2)
        a, b = vs.block_nodes[1, k], vs.block_nodes[2, k]
        (rotated[a] || rotated[b]) || continue
        for i in 1:9
            B[i] = nz[vs.blocks[i, k]]
        end
        if rotated[b]                       # B ← B Q_b
            for j in 1:3, i in 1:3
                T[i, j] = B[i, 1]*Q[1, j, b] + B[i, 2]*Q[2, j, b] + B[i, 3]*Q[3, j, b]
            end
            B .= T
        end
        if rotated[a]                       # B ← Q_a' B
            for j in 1:3, i in 1:3
                T[i, j] = Q[1, i, a]*B[1, j] + Q[2, i, a]*B[2, j] + Q[3, i, a]*B[3, j]
            end
            B .= T
        end
        for i in 1:9
            nz[vs.blocks[i, k]] = B[i]
        end
    end
    return nz
end

"""
    volume_update!(mdl::LinearElasticity, nodes, normals, d; β_table=0.0)

Volume update of a window start: move the model's mesh so that each node in `nodes` is
displaced by `d` along its unit normal (`u·n = d`, tangential motion free), with `u_z = 0` on
the table (every node of the mesh's `IEN_bottom` faces) and the rest of the surface
traction-free. The interior
follows a compressible linear elastic solve on the current geometry; with displacement
constraints only, the result depends on Poisson's ratio alone, not on the stiffness scale.

Constraints on observed and on rigid-contact nodes take the same form; a node may appear more
than once (e.g. an observed node on the table). The sparsity pattern and symbolic factorization
are cached per mesh, so a call costs one assembly, one numeric factorization and a solve.

# Arguments
- `mdl::LinearElasticity`: Model; its mesh is moved (`NodeList` replaced) and needs
  `IEN_bottom` (cylinder or cuboid).
- `nodes::Vector{Int}`: Constrained nodes.
- `normals::Matrix{Float64}`: (3, length(nodes)) normal of each constraint; normalized here.
- `d::Vector{Float64}`: Displacement of each node along its unit normal, in mm.
- `β_table::Real`: `[not decided]` Navier slip on the table, in kPa/mm. With the default 0 the
  table is frictionless, and the constraints must remove the horizontal rigid-body modes
  themselves. They do not on an axisymmetric surface with radial normals (the rotation about z
  stays free): the factorization then warns that the stiffness is singular, and the result
  carries an arbitrary rotation.

# Returns
- `u::Matrix{Float64}`: Nodal displacements (3 × nNodes).
"""
function volume_update!(mdl::LinearElasticity, nodes::Vector{Int}, normals::Matrix{Float64},
                        d::Vector{Float64}; β_table::Real=0.0)
    size(normals) == (3, length(nodes)) || throw(DimensionMismatch("normals must be 3 × $(length(nodes)), got $(size(normals))"))
    length(d) == length(nodes) || throw(DimensionMismatch("d has $(length(d)) entries for $(length(nodes)) nodes"))
    mesh = mdl.mesh
    nN::Int = mesh.nNodes          # mesh is abstractly typed: keep the loops below concrete
    X::Matrix{Float64} = mesh.NodeList
    vs = _cached_volume_system(mesh)
    sys, K = vs.sys, vs.sys.K
    _assemble!(sys, X, mesh.IEN, mdl.λ, mdl.μ, β_table)

    constraints = [Tuple{Vector{Float64},Float64}[] for _ in 1:nN]
    # the table nodes from the face connectivity: `bottom_nodes` leaves out those on the side
    for a in unique(vec(mesh.IEN_bottom))
        push!(constraints[a], ([0.0, 0.0, 1.0], 0.0))
    end
    for (k, a) in enumerate(nodes)
        n = normals[:, k]
        push!(constraints[a], (n ./ norm(n), d[k]))
    end

    Q = zeros(3, 3, nN)
    vals = zeros(3, nN)
    nfix = zeros(Int, nN)
    rotated = falses(nN)
    for a in 1:nN
        isempty(constraints[a]) && continue
        nfix[a] = _constraint_frame!(view(Q, :, :, a), view(vals, :, a), constraints[a])
        rotated[a] = true
    end
    _rotate_blocks!(K.nzval, vs, Q, rotated)

    # prescribed unknowns: move their columns to the right-hand side, then make them identity rows
    rhs = zeros(size(K, 1))
    for a in 1:nN, m in 1:nfix[a]
        c, g = 3(a-1) + m, vals[m, a]
        for p in K.colptr[c]:(K.colptr[c+1] - 1)
            rhs[K.rowval[p]] -= K.nzval[p] * g
        end
    end
    for a in 1:nN, m in 1:nfix[a]
        c = 3(a-1) + m
        for p in K.colptr[c]:(K.colptr[c+1] - 1)
            K.nzval[p] = 0.0
            K.nzval[vs.mirror[p]] = 0.0
        end
    end
    for a in 1:nN, m in 1:nfix[a]
        c = 3(a-1) + m
        K.nzval[vs.diag[c]] = 1.0
        rhs[c] = vals[m, a]
    end

    _refactor!(sys)
    x = zeros(size(K, 1))
    _refined_solve!(sys, x, rhs; tol=STATE_TOL)

    u = reshape(x, 3, nN)
    for a in 1:nN
        rotated[a] && (u[:, a] = view(Q, :, :, a) * u[:, a])
    end
    mesh.NodeList = X + u           # a new array: the mesh's initial state may alias NodeList
    return u
end
