using Sockets
using JSON
using HDF5
using Interpolations: BSpline, Linear, extrapolate, interpolate!, scale   # explicit: Distributions also exports `scale`

# Link to the perception process (smearPerception, Python): perception listens on a Unix-domain
# socket and sends one message per FE step; this side connects and reads. The messages are
# documented in smearPerception's `src/smear_perception/fe_link.py`; keep the two in step. The
# offline form is the HDF5 file of its `scripts/reconstruct_bag.py --save-field`.
#
# Perception works in SI units with C-order arrays. This side converts once, on read, to the FE's
# units (mm; forces in kg·mm/s²). The fields are permuted back to perception's axis order, so
# `phi[i, j, k]` is the same voxel on both sides; points (P, 3) arrive as (3, P), and the 4×4 probe
# pose is transposed back.

const PERCEPTION_LINK_PATH = "/tmp/smear_perception.sock"
const LINK_MAGIC = b"SMEARLNK"
const LINK_VERSION = UInt32(2)

"""
    PerceptionGrid

The voxel grid of the perception fields, in mm: voxel (i, j, k) (0-based, perception's array
index) is centred at `x0 + delta * ((i, j, k) - offset)` and stored at `phi[i+1, j+1, k+1]`.
"""
struct PerceptionGrid
    x0::Vector{Float64}       # centre of the bottom face [mm]
    delta::Float64            # voxel size [mm]
    offset::Vector{Int}       # array index of x0
    shape::Vector{Int}        # (nx, ny, nz), perception's order
end

"""
    PerceptionStep

What perception reports for one FE step, in FE units (mm, kg·mm/s²).

- `n`: global step index; `t_ns` = t0 + n·Δt [ns].
- `n_views`: cameras that delivered a frame.
- `phi`: fused field φ_s, (nx, ny, nz), NaN where unobserved.
- `count`: contributing cameras N(x), same dims.
- `points`, `normals`: (3, P) zero crossings of φ_s and their unit normals (out of the material).
- `contact`: 1 if |F| exceeded perception's contact threshold, 0 if not, -1 unknown (no threshold
  set, or no F/T sample).
- `wrench`: [Fx, Fy, Fz, Mx, My, Mz] at t_n in the F/T sensor frame [kg·mm/s², kg·mm²/s²], at every
  step; NaN without a sample.
- `probe_pose`: 4×4 pose of the probe (probe → world) [mm]; NaN while it is not tracked.
"""
mutable struct PerceptionStep{M<:AbstractMatrix{Float32}}
    n::Int64
    t_ns::Int64
    n_views::Int
    phi::Array{Float32,3}
    count::Array{UInt8,3}
    points::M
    normals::M
    contact::Int8
    wrench::Vector{Float64}
    probe_pose::Matrix{Float64}
end

function _perception_grid(x0, delta, offset, shape)
    return PerceptionGrid(1e3 .* Float64.(x0), 1e3 * Float64(delta), Int.(offset), Int.(shape))
end

function _to_fe_units!(st::PerceptionStep)
    st.phi .*= 1f3
    st.points .*= 1f3
    @views st.wrench[1:3] .*= 1e3     # N → kg·mm/s²
    @views st.wrench[4:6] .*= 1e6     # N·m → kg·mm²/s²
    @views st.probe_pose[1:3, 4] .*= 1e3
    return st
end

"""
    read_perception_steps(path)

Read a camera reconstruction written by smearPerception's `scripts/reconstruct_bag.py` (run with
`--save-field`): the offline form of the perception link.

# Arguments
- `path::AbstractString`: The HDF5 file.

# Returns
- `(; grid::PerceptionGrid, steps::Vector{PerceptionStep}, t0_ns, dt_ns)`, in FE units.
"""
function read_perception_steps(path::AbstractString)
    h5open(path, "r") do f
        haskey(f, "phi") || error("$path has no phi: run reconstruct_bag.py with --save-field")
        haskey(f, "n") || error("$path has no step records (written before 2026-10-23): run reconstruct_bag.py again")
        a = attributes(f)
        m = read(a["grid_m"])
        grid = _perception_grid(read(a["grid_x0"]), read(a["grid_delta"]), [m[1], m[2], 0],
                                read(a["grid_shape"]))
        n, t_ns, n_views, index = read(f["n"]), read(f["t_ns"]), read(f["n_views"]), read(f["index"])
        contact, wrench, pose = read(f["contact"]), read(f["wrench"]), read(f["probe_pose"])
        points, normals = read(f["points"]), read(f["normals"])
        phi, count = read(f["phi"]), read(f["count"])
        steps = map(eachindex(t_ns)) do s
            start, len = index[1, s], index[2, s]
            cols = start+1:start+len
            _to_fe_units!(PerceptionStep(n[s], t_ns[s], Int(n_views[s]),
                                         permutedims(phi[:, :, :, s], (3, 2, 1)), permutedims(count[:, :, :, s], (3, 2, 1)),
                                         points[:, cols], normals[:, cols], contact[s], wrench[:, s],
                                         permutedims(pose[:, :, s])))
        end
        return (; grid, steps, t0_ns=Int64(read(a["t0_ns"])), dt_ns=Int64(read(a["dt_ns"])))
    end
end

"""
    PerceptionLink

A perception stream being read: the run's constants and the buffers every step is read into. See
`open_perception_link`.
"""
mutable struct PerceptionLink
    io::IO
    grid::PerceptionGrid
    t0_ns::Int64
    dt_ns::Int64
    meta::Dict{String,Any}      # the hello as sent (SI units)
    step::PerceptionStep        # reused for every step
    points::Matrix{Float32}     # capacity behind step.points and step.normals
    normals::Matrix{Float32}
    pose::Matrix{Float64}       # the pose as sent (transposed)
    phi::Array{Float32,3}       # the fields as sent (axes reversed), permuted into step
    count::Array{UInt8,3}
end

# Short reads mean the stream ended: perception stopped without END.
function _read_bytes(io::IO, n::Integer)
    b = read(io, n)
    length(b) == n || throw(EOFError())
    return b
end

"""
    PerceptionLink(io)

Start reading a perception stream from `io` (a connected socket, or a file holding the same
bytes): read its hello message.
"""
function PerceptionLink(io::IO)
    _read_bytes(io, 8) == LINK_MAGIC || error("not a perception link")
    version = read(io, UInt32)
    version == LINK_VERSION || error("perception link version $version, expected $LINK_VERSION")
    meta = JSON.parse(String(_read_bytes(io, read(io, UInt32))))
    meta["units"] == "m" || error("perception link in $(meta["units"]), expected m")
    grid = _perception_grid(meta["grid_x0"], meta["grid_delta"], meta["grid_offset"], meta["grid_shape"])
    nx, ny, nz = grid.shape
    points, normals = Matrix{Float32}(undef, 3, 0), Matrix{Float32}(undef, 3, 0)
    step = PerceptionStep(Int64(0), Int64(0), 0, Array{Float32}(undef, nx, ny, nz), Array{UInt8}(undef, nx, ny, nz),
                          view(points, :, 1:0), view(normals, :, 1:0), Int8(-1), fill(NaN, 6), fill(NaN, 4, 4))
    return PerceptionLink(io, grid, meta["t0_ns"], meta["dt_ns"], meta, step, points, normals, zeros(4, 4),
                          Array{Float32}(undef, nz, ny, nx), Array{UInt8}(undef, nz, ny, nx))
end

"""
    open_perception_link(path=PERCEPTION_LINK_PATH; timeout=Inf)

Connect to the perception process (smearPerception's `fe_link.StepPublisher`, e.g.
`reconstruct_bag.py --link`). Start order does not matter: this retries until perception listens.

# Arguments
- `path::AbstractString`: The socket file perception listens on.
- `timeout::Real`: Seconds to keep retrying.

# Returns
- `PerceptionLink`, once perception has sent the hello (with its first step); read the steps with
  `next_step!`.
"""
function open_perception_link(path::AbstractString=PERCEPTION_LINK_PATH; timeout::Real=Inf)
    deadline = time() + timeout
    while true
        sock = try
            connect(path)
        catch e
            e isa Base.IOError || rethrow()
            time() > deadline && error("no perception link at $path")
            sleep(0.05)
            continue
        end
        return PerceptionLink(sock)
    end
end

"""
    next_step!(link)

The next step from perception, in order; blocks until it arrives. No step is skipped: perception
waits when the FE falls behind.

# Arguments
- `link::PerceptionLink`: From `open_perception_link`.

# Returns
- `PerceptionStep` (FE units), or `nothing` once perception has finished, so
  `while (st = next_step!(link)) !== nothing` ends. It is the link's own buffer, overwritten by the
  next call (`deepcopy` it to keep it): reading allocates almost nothing, since the GC pauses
  (~40–250 ms) that per-step allocation caused stalled the link. Throws `EOFError` if perception
  stopped without finishing (a crash).
"""
function next_step!(link::PerceptionLink)
    io, st = link.io, link.step
    tag = _read_bytes(io, 4)
    tag == b"END." && return nothing
    tag == b"STEP" || error("perception link out of step: $(repr(String(tag))) where a message should start")
    st.n = read(io, Int64)
    st.t_ns = read(io, Int64)
    st.n_views = read(io, Int32)
    st.contact = read(io, Int8)
    n_points = Int(read(io, Int64))
    read!(io, st.wrench)
    read!(io, link.pose)
    permutedims!(st.probe_pose, link.pose, (2, 1))
    read!(io, link.phi)
    read!(io, link.count)
    permutedims!(st.phi, link.phi, (3, 2, 1))
    permutedims!(st.count, link.count, (3, 2, 1))
    if n_points > size(link.points, 2)    # rare: room for twice as many
        link.points = Matrix{Float32}(undef, 3, 2n_points)
        link.normals = Matrix{Float32}(undef, 3, 2n_points)
    end
    st.points, st.normals = view(link.points, :, 1:n_points), view(link.normals, :, 1:n_points)
    read!(io, st.points)
    read!(io, st.normals)
    return _to_fe_units!(st)
end

"""
    close_link(link)

Close the connection to perception; perception stops at its next step (its send fails).
"""
close_link(link::PerceptionLink) = close(link.io)

"""
    sample_field(phi, grid, points)

Trilinear interpolation of a perception field at points (e.g. FE surface nodes), by
Interpolations.jl; the same rule as smearPerception's `voxels.sample_field` (SciPy).

# Arguments
- `phi::Array{<:Real,3}`: Field (nx, ny, nz), e.g. `step.phi` [mm]; wrapped, not copied.
- `grid::PerceptionGrid`: Its grid.
- `points::AbstractMatrix{<:Real}`: (3, n) points in the world frame [mm].

# Returns
- `Vector{Float64}` of n values; NaN outside the voxel centres or next to an undefined voxel (a NaN
  corner makes the result NaN).
"""
function sample_field(phi::Array{<:Real,3}, grid::PerceptionGrid, points::AbstractMatrix{<:Real})
    centres = [grid.x0[a] .+ grid.delta .* ((0:grid.shape[a]-1) .- grid.offset[a]) for a in 1:3]
    # interpolate! because linear B-splines need no prefiltering: no copy of the field per step.
    itp = extrapolate(scale(interpolate!(phi, BSpline(Linear())), centres...), NaN)
    return [Float64(itp(points[1, p], points[2, p], points[3, p])) for p in axes(points, 2)]
end
