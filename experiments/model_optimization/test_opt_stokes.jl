using LinearAlgebra
using ProgressMeter
using SparseArrays
using NearestNeighbors

using smearFEM
using StatsPlots
using Distributions
using Dates
using Plots
using LaTeXStrings
using DelimitedFiles
using Printf
using Statistics

using ArgCheck
using Base.Threads
using Random
using Colors

using Dates
using Plots.PlotMeasures

include("../ParallelExecution.jl")
using .ParallelExecution

global def_orange = RGB(245/255,118/255,0)
global def_blue = RGB(5/255,79/255,185/255)
global def_red = RGB(196/255,70/255,1/255)
global def_green = RGB(2/255,147/255,86/255)
global end_obs_win = 20.1

# PLOT_CONFIG: default geometry seeding the globals below at load time.
# PLOT_PRESETS: per-(data_type[, viscosity_type]) overrides applied by
# set_plot_config() — the place to edit sizes per data type.
# for 1/2 linewidth
const PLOT_CONFIG = Dict(
    :font_size => 12,
    :plot_height => 320,
    :plot_width => 480,
    :left_margin => 1pt,
    :right_margin => 5pt,
    :top_margin => 1pt,
)

# for 1/3 linewidth
# const PLOT_CONFIG = Dict(:font_size => 12, :plot_height => 360, :plot_width => 330,
#                           :left_margin => -6pt, :right_margin => 10pt, :top_margin => 0pt)

# for 1/4 linewidth
# const PLOT_CONFIG = Dict(:font_size => 10, :plot_height => 300, :plot_width => 239,
#                           :left_margin => -6pt, :right_margin => 0pt, :top_margin => -1pt)

global fs::Int = PLOT_CONFIG[:font_size]
global plt_height::Int = PLOT_CONFIG[:plot_height]
global plt_width::Int = PLOT_CONFIG[:plot_width]
global plt_lft_margin = PLOT_CONFIG[:left_margin]
global plt_right_margin = PLOT_CONFIG[:right_margin]
global plt_top_margin = PLOT_CONFIG[:top_margin]

# Fixed 1/2-linewidth geometry for the multi-window η/β plots (plt_η, plt_β),
# independent of the data_type/viscosity_type preset applied by set_plot_config().
global eta_beta_plt_width::Int = 330
global eta_beta_plt_height::Int = 360

# Per-(data_type[, viscosity_type]) geometry + y-limit presets, applied by
# set_plot_config(). Keyed by data_type alone except for "synthetic", which
# is further split by viscosity_type since its geometry differs.
const PLOT_PRESETS = Dict(
    "physical" => Dict(
        :font_size => 12, :plot_height => 360, :plot_width => 330,
        :left_margin => 1pt, :right_margin => 5pt, :top_margin => 1pt,
        :y_lims_h_norm => (0.8, 1.05), :y_lims_rel_error => (-0.05, 20),
    ),
    ("synthetic", "constant") => Dict(
        :font_size => 12, :plot_height => 350, :plot_width => 477,
        :left_margin => 1pt, :right_margin => 5pt, :top_margin => 1pt,
        :y_lims_h_norm => (0.97, 1.03), :y_lims_rel_error => (-0.1, 3.0),
    ),
    ("synthetic", "bulk_viscosity") => Dict(
        :font_size => 12, :plot_height => 360, :plot_width => 330,
        :left_margin => -2pt, :right_margin => 1pt, :top_margin => 0pt,
        :y_lims_h_norm => (0.97, 1.031), :y_lims_rel_error => (-0.1, 3.0),
    ),
    "simulated" => Dict(
        :font_size => 12, :plot_height => 350, :plot_width => 477,
        :left_margin => 1pt, :right_margin => 5pt, :top_margin => 1pt,
        :y_lims_h_norm => (0.995, 1.005), :y_lims_rel_error => (-0.05, 0.1),
    ),
)

# for synthetic data (default preset until set_plot_config() is called)
global y_lims_h_norm = PLOT_PRESETS[("synthetic", "constant")][:y_lims_h_norm]
global y_lims_rel_error = PLOT_PRESETS[("synthetic", "constant")][:y_lims_rel_error]

# Colorblind-friendly 20-color palette, spaced across hue/brightness.
global color_palette = [
    RGB(230/255, 159/255, 0/255),       # 1. orange (warm, bright)
    RGB(0/255, 114/255, 178/255),       # 2. medium blue (cool)
    RGB(204/255, 0/255, 0/255),         # 3. red (warm)
    RGB(0/255, 158/255, 115/255),       # 4. green (cool, mid-tone)
    RGB(204/255, 41/255, 130/255),      # 5. magenta/pink (warm)
    RGB(0/255, 0/255, 153/255),         # 6. dark blue (cool, dark)
    RGB(213/255, 94/255, 0/255),        # 7. vermillion (warm, dark)
    RGB(75/255, 0/255, 130/255),        # 8. indigo (cool, dark)
    RGB(255/255, 127/255, 80/255),      # 9. coral (warm, light)
    RGB(86/255, 180/255, 233/255),      # 10. sky blue (cool, light)
    RGB(178/255, 140/255, 51/255),      # 11. tan/brown (neutral, warm)
    RGB(0/255, 150/255, 200/255),       # 12. cyan (cool, bright)
    RGB(153/255, 0/255, 76/255),        # 13. deep magenta (warm, dark)
    RGB(51/255, 153/255, 51/255),       # 14. forest green (cool, dark)
    RGB(200/255, 50/255, 100/255),      # 15. deep rose (warm)
    RGB(0/255, 120/255, 100/255),       # 16. dark teal (cool, dark)
    RGB(102/255, 0/255, 204/255),       # 17. purple (cool, bright)
    RGB(150/255, 75/255, 0/255),        # 18. dark brown (neutral, dark)
    RGB(220/255, 80/255, 150/255),      # 19. light pink (warm, light)
    RGB(89/255, 89/255, 89/255)         # 20. dark gray (neutral)
]

# JSON stores matrices column-major flattened; reshape back to 3x3.
_read_camera_matrix(sim_params) = reshape(Array(float.(sim_params["camera_matrix"])), 3, 3)

# Synthetic pose stores only a camera position ([x,y,z] from column 4, permuted
# to [3,1,2]); the orientation is synthesised downstream.
function _read_obj_pose(sim_params)
    pose = Array(float.(sim_params["obj_pose"]))
    return length(pose) == 16 ? reshape(pose, 4, 4)[[3,1,2], 4] : pose
end

# Physical pose is marker-measured: keep the full 4×4 (carries the camera tilt).
_read_obj_pose_matrix(sim_params) = reshape(Array(float.(sim_params["obj_pose"])), 4, 4)

# Physical ground truth records its data_type; synthetic does not.
_obj_pose_for(sim_params) = get(sim_params, "data_type", "") == "physical" ?
    _read_obj_pose_matrix(sim_params) : _read_obj_pose(sim_params)

# Path helpers for the per-experiment `data/` and `plots/` subdirectories.
datapath(p, parts...) = joinpath(p, "data", parts...)
plotpath(p, parts...) = joinpath(p, "plots", parts...)

# Levels of the `collect_experiment_groups` tree that post_analysis_* aggregates over.
# A group must hold exactly one dt_*/view_* leaf; more is a dimension these do not plot.
_elem_size_folders(groups) = unique(g.elem_size_folder for g in groups)
_sim_time_folders(groups, elem_size_folder) =
    unique(g.sim_time_folder for g in groups if g.elem_size_folder == elem_size_folder)
_noise_folders(groups, elem_size_folder, sim_time_folder) =
    unique(g.noise_folder for g in groups if g.elem_size_folder == elem_size_folder &&
                                             g.sim_time_folder == sim_time_folder)

function _find_group(groups, elem_size_folder, sim_time_folder, noise_folder)
    matches = filter(g -> g.elem_size_folder == elem_size_folder &&
                          g.sim_time_folder == sim_time_folder &&
                          g.noise_folder == noise_folder, groups)
    isempty(matches) && return nothing
    return only(matches)
end

function _single_leaf(group)
    if isempty(group.leaves)
        @warn "No dt_*/view_* leaf found, skipping" group.step_path
        return nothing
    elseif length(group.leaves) > 1
        error("Expected a single dt_*/view_* leaf under $(group.step_path), found $(length(group.leaves)): " *
              join((joinpath(l.step_folder, l.view_folder) for l in group.leaves), ", "))
    end
    return group.leaves[1]
end

# Default figure geometry (the globals seeded from PLOT_CONFIG at load time).
default_plot() = _fig(margins=:all)

"""
    _ensure_scalar_float(x) -> Float64

Coerce `x` to a `Float64`, taking the first element if `x` is a collection
rather than a plain `Number` (defensive against JSON-backed values that may
come back as length-1 arrays instead of scalars).

# Arguments
- `x`: a `Number`, or a collection whose first element is a `Number`.

# Returns
- `Float64`: `x` (or its first element) as a `Float64`.
"""
function _ensure_scalar_float(x)
    if isa(x, Number)
        return Float64(x)
    end
    try
        c = collect(x)
        if length(c) >= 1
            return Float64(c[1])
        else
            error("_ensure_scalar_float: container is empty")
        end
    catch err
        error("_ensure_scalar_float: cannot convert input to scalar Float64: $err")
    end
end

"""
    interp_z_at(x, y, xg, yg, Z) -> Float64

Bilinearly interpolate `Z` (defined on the grid `xg` × `yg`, with rows
indexing `yg` and columns indexing `xg`) at the point `(x, y)`, clamping
`(x, y)` to the grid extents. Tolerates `Z` being 1D, transposed, or of a
non-`Float64` element type (e.g. JSON-backed arrays) by coercing/reshaping it
to the expected `(length(yg), length(xg))` orientation first.

# Arguments
- `x::Real`, `y::Real`: point to interpolate at.
- `xg::AbstractVector`, `yg::AbstractVector`: grid coordinates `Z` is defined on.
- `Z::AbstractMatrix`: surface values on the `xg` × `yg` grid.

# Returns
- `Float64`: the interpolated value at `(x, y)`.
"""
function interp_z_at(x::Real, y::Real, xg::AbstractVector, yg::AbstractVector, Z::AbstractMatrix)
    # Coerce grid vectors and Z to concrete Float64 arrays/matrix and
    # normalize shapes. This avoids errors when Z is an Adjoint, 1D Vec,
    # or contains Union element types (JSON-backed arrays).
    xg_arr = Float64.(collect(xg))
    yg_arr = Float64.(collect(yg))
    nx = length(xg_arr); ny = length(yg_arr)

    Zmat = try
        Array(Z) |> x -> Float64.(x)
    catch err
        error("Failed to coerce Z for interpolation: $err")
    end

    # If Zmat ended up 1D (or different shape), try to reshape/transpose to
    # match (ny, nx) where rows -> y and cols -> x.
    if ndims(Zmat) == 1
        if length(Zmat) == nx*ny
            Zmat = reshape(Zmat, ny, nx)
        else
            error("interp_z_at: Z is 1D with length=$(length(Zmat)) which is not nx*ny=$(nx*ny)")
        end
    elseif ndims(Zmat) == 2
        if size(Zmat,1) == ny && size(Zmat,2) == nx
        elseif size(Zmat,1) == nx && size(Zmat,2) == ny
            Zmat = Zmat'
        elseif prod(size(Zmat)) == nx*ny
            Zmat = reshape(vec(Zmat), ny, nx)
        else
            @warn "interp_z_at: Z matrix shape $(size(Zmat)) does not match (ny,nx)=($(ny),$(nx)). Attempting to continue but interpolation may be wrong."
        end
    else
        error("interp_z_at: Z has unexpected number of dimensions: $(ndims(Zmat))")
    end

    # clamp x,y to grid extents
    xcl = clamp(x, first(xg_arr), last(xg_arr))
    ycl = clamp(y, first(yg_arr), last(yg_arr))

    ix = searchsortedfirst(xg_arr, xcl)
    if ix == 1
        i0 = 1; i1 = 1; tx = 0.0
    elseif ix > nx
        i0 = nx; i1 = nx; tx = 0.0
    else
        i1 = ix; i0 = max(1, ix-1)
        x0 = xg_arr[i0]; x1 = xg_arr[i1]
        tx = x1==x0 ? 0.0 : (xcl - x0)/(x1 - x0)
    end

    iy = searchsortedfirst(yg, ycl)
    if iy == 1
        j0 = 1; j1 = 1; ty = 0.0
    elseif iy > ny
        j0 = ny; j1 = ny; ty = 0.0
    else
        j1 = iy; j0 = max(1, iy-1)
        y0 = yg_arr[j0]; y1 = yg_arr[j1]
        ty = y1==y0 ? 0.0 : (ycl - y0)/(y1 - y0)
    end

    # values at corners: note Z rows->y (j index), cols->x (i index)
    z00 = float(Zmat[j0, i0]); z10 = float(Zmat[j0, i1])
    z01 = float(Zmat[j1, i0]); z11 = float(Zmat[j1, i1])

    # bilinear interpolation
    z0 = (1-tx)*z00 + tx*z10
    z1 = (1-tx)*z01 + tx*z11
    z = (1-ty)*z0 + ty*z1
    return z
end

"""
    interp_z_at(xs, ys, xg, yg, Z) -> Vector{Float64}

Vectorized form of `interp_z_at(x, y, xg, yg, Z)`: bilinearly
interpolate `Z` at each `(xs[i], ys[i])` pair.

# Arguments
- `xs::AbstractVector{<:Real}`, `ys::AbstractVector{<:Real}`: points to interpolate at.
- `xg::AbstractVector`, `yg::AbstractVector`: grid coordinates `Z` is defined on.
- `Z::AbstractMatrix`: surface values on the `xg` × `yg` grid.

# Returns
- `Vector{Float64}`: the interpolated value at each `(xs[i], ys[i])`.
"""
function interp_z_at(xs::AbstractVector{<:Real}, ys::AbstractVector{<:Real}, xg::AbstractVector, yg::AbstractVector, Z::AbstractMatrix)
    if length(xs) == 1 && length(ys) == 1
        return interp_z_at(xs[1], ys[1], xg, yg, Z)
    end
    # If paired vectors, compute element-wise
    if length(xs) == length(ys)
        return [interp_z_at(xi, yi, xg, yg, Z) for (xi, yi) in zip(xs, ys)]
    end
    error("interp_z_at: expecting scalar x,y or vectors of equal length")
end

"""
    set_time_window(time_step_len, data; method="linear", window_size=10.0, global_window_sz=30.0) -> (time_windows, windows, data_ranges, t_windows)

Split `data` (indexed along its first dimension by time) into successive time
windows whose cumulative end time grows according to `method`
(`"linear"`: `window_size*iter`; `"quadratic"`: `window_size*iter^2`;
`"exponential"`: `window_size*exp(3*(iter-1))`), stopping once the window end
reaches the end of `data` or `global_window_sz`.

# Arguments
- `time_step_len::Float64`: time step length of `data`.
- `data::AbstractArray`: array indexed along its first dimension by time.

# Keyword Arguments
- `method::String`: window growth schedule, `"linear"`, `"quadratic"`, or
  `"exponential"` (default: `"linear"`).
- `window_size::Float64`: base window size (default: `10.0`).
- `global_window_sz::Float64`: overall time budget the windows are truncated
  to (default: `30.0`).

# Returns
- `time_windows::Vector{Float64}`: duration of each window.
- `windows::Vector{AbstractArray}`: the slice of `data` in each window.
- `data_ranges::Vector{AbstractArray}`: the index range (relative to `data`)
  of each window.
- `t_windows::Vector{Float64}`: cumulative end time of each window.
"""
function set_time_window(time_step_len::Float64, data::AbstractArray; method::String="linear", window_size::Float64=10.0, global_window_sz::Float64=50.0)
    windows::Vector{AbstractArray} = Vector{AbstractArray}()
    time_windows::Vector{Float64} = Vector{Float64}()
    data_ranges::Vector{AbstractArray} = Vector{AbstractArray}()
    t_windows::Vector{Float64} = Vector{Float64}()

    function get_t_window(window_size::Float64, step_len::Float64, iter::Int, method::String)::Float64
        t_window = 0.0
        if method == "linear"
            t_window = round(window_size*iter, digits=1)
        elseif method == "quadratic"
            t_window = round(window_size*iter^2, digits=1)
        elseif method == "exponential"
            if iter == 1
                t_window = round(window_size*exp(0.5*(iter-1)), digits=1)
            else 
                δt = window_size*round(exp(3*(iter-1)), digits=1)
                if δt < window_size
                    @warn "Computed time window increment $δt is less than minimum Window $window_size; using minimum."
                    δt = window_size
                end
                t_window = round(window_size*exp(3*(iter-1)), digits=1)
            end
        end
        return t_window
    end

    iter::Int = 1
    start_point::Int = 1
    t_window_prev::Float64 = 0.0
    t_window_end::Float64 = get_t_window(window_size, time_step_len, iter, method)
    end_point::Int = round(Int,t_window_end*time_step_len)+1

    END_FLAG::Bool = false
    while true
            # If the computed end_point reaches or exceeds the data length, adjust.
            # Keep the original requested value for clearer messages.
            if end_point >= size(data, 1) || t_window_end >= global_window_sz
                if end_point == size(data, 1)
                    @info "Reached end of data at point $end_point."
                elseif end_point > size(data, 1)
                    requested_end = end_point
                    # adjust the time window end to the last available sample
                    t_window_end = round((size(data, 1)-1)/time_step_len, digits=1)
                    end_point = size(data, 1)
                    @warn "Requested end point $requested_end exceeds data size $(size(data,1)); adjusting to end of data (end_point=$end_point). Adjusted end time to $t_window_end seconds."
                elseif t_window_end >= global_window_sz
                    requested_end = t_window_end
                    # adjust the end point to match the global window size
                    t_window_end = global_window_sz
                    end_point = round(Int,t_window_end*time_step_len)+1
                    @warn "Requested time window end $requested_end seconds exceeds global window size $global_window_sz seconds; adjusting to global window size (end_point=$end_point)."
                end
                END_FLAG = true
            end

        data_range = start_point:end_point
        data_range_ = start_point:(end_point-1)

        @debug "Data frame : $data_range"
        @debug "time windows from : $t_window_prev to $t_window_end"
        @debug "time Window : $(t_window_end - t_window_prev) seconds"
        @debug "data length : $(size(data[data_range], 1))"
        @debug "data range size : $(length(data_range_))"
        @debug "----------"
        t_window_size = round(t_window_end - t_window_prev, digits=1)
        push!(time_windows, t_window_size)
        push!(windows, data[data_range])
        push!(data_ranges, data_range_)
        push!(t_windows, t_window_end)

        if END_FLAG == true
            break
        end
        iter = iter + 1
        start_point = end_point
        t_window_prev = t_window_end
        t_window_end = get_t_window(window_size, time_step_len, iter, method)
        end_point = round(Int,t_window_end*time_step_len)+1
    end
    return time_windows, windows, data_ranges, t_windows
end

"""
    _get_borders(data_type, filepath_gt, exp_path, num_exp_points; view_folder="view_1") -> (obs_border_pt_lst, sim_border_pt_lst, gt_Splinex, gt_Spliney, splinex, spliney)

Load and align the observed (ground-truth) and simulated 2D border/contour
data for one experiment, truncated to the first `num_exp_points` frames.

# Arguments
- `data_type::String`: `"synthetic"`, `"physical"`, or `"simulated"` — selects
  where the ground-truth contour data is read from under `filepath_gt`.
- `filepath_gt::String`: ground-truth data root.
- `exp_path::String`: experiment result directory containing the simulated
  `2D_border_points`.
- `num_exp_points::Int`: number of leading time frames to keep.

# Keyword Arguments
- `view_folder::String`: camera-view subdirectory to read ground-truth data
  from (default: `"view_1"`).

# Returns
- `obs_border_pt_lst`, `sim_border_pt_lst`: observed and simulated border
  points, truncated to `num_exp_points` frames.
- `gt_Splinex`, `gt_Spliney`, `splinex`, `spliney`: matching spline-sampled
  x/y coordinates for the observed and simulated contours.
"""
function _get_borders(data_type::String, filepath_gt::String, exp_path::String, num_exp_points::Int; view_folder::String="view_1")

    if data_type  == "synthetic"
        ObsDataList, splinexObs, splineyObs = read_csv(datapath(filepath_gt,"img_data",view_folder,"contour_data"))
        @info "Data type $data_type Reading synthetic ground truth contour data of $(length(ObsDataList)) time steps"
    elseif data_type == "physical"
        # Physical ground truth is single-camera and stores contours flat (no view_* level),
        # matching what `optimize` reads for data_type == "physical".
        ObsDataList, splinexObs, splineyObs = read_csv(datapath(filepath_gt,"img_data","contour_data"))
        # The video starts before the press engages, so index by `compressed_frames` exactly
        # as `optimize` does. Without this, frame k is scored against simulation step k and
        # lands ~100 frames early — every contour metric downstream is measured against
        # footage of an unloaded specimen.
        meta = read_json(datapath(filepath_gt,"video_metadata"))
        ObsDataList = ObsDataList[meta["compressed_frames"]]
        splinexObs = splinexObs[meta["compressed_frames"]]
        splineyObs = splineyObs[meta["compressed_frames"]]
        @info "Data type $data_type Reading physical contour data of $(length(ObsDataList)) time steps"
    elseif data_type == "simulated"
        ObsDataList, splinexObs, splineyObs = read_csv(datapath(filepath_gt,"sim_data",view_folder,"contour_data"))
        @info "Data type : $data_type Reading simulated ground truth contour data from $(datapath(filepath_gt,"sim_data",view_folder,"contour_data")) of $(length(ObsDataList)) time steps"
    else
        error("Unknown data type: $data_type")
    end

    obs_border_pt_lst, gt_Splinex, gt_Spliney, pd = add_noise(ObsDataList, nFactor=0.0)
    # write_2d_data numbers angles within one simulation; experiments sim a single angle.
    sim_border_pt_lst, splinex, spliney = read_csv(datapath(exp_path,"sim_data","view_1","2D_border_points_est"))
    @info "Reading simulated contour data from $(datapath(exp_path,"sim_data","view_1","2D_border_points_est")) of $(length(sim_border_pt_lst)) time steps"

    @assert length(ObsDataList) >= num_exp_points "Not enough observation border points: have $(length(ObsDataList)), need at least $num_exp_points"
    @assert length(sim_border_pt_lst) >= num_exp_points "Not enough simulation border points: have $(length(sim_border_pt_lst)), need at least $num_exp_points"

    obs_border_pt_lst = obs_border_pt_lst[1:num_exp_points, :]
    sim_border_pt_lst = sim_border_pt_lst[1:num_exp_points, :]

    gt_Splinex = gt_Splinex[1:num_exp_points, :]
    gt_Spliney = gt_Spliney[1:num_exp_points, :]

    splinex = splinex[1:num_exp_points, :]
    spliney = spliney[1:num_exp_points, :]

    return obs_border_pt_lst, sim_border_pt_lst, gt_Splinex, gt_Spliney, splinex, spliney
end

"""
    iteration_contour_metrics(data_type, filepath_gt, exp_path; view_folder, sim_view, write_results)
        -> Dict

Score every optimizer iterate against the observed contours.

Reads the per-iteration simulated contours written by `optimize` under
`<exp_path>/data/opt_data/iter_<n>/<sim_view>/2D_border_points` (produced when
`fit_model` is called with `store_border_pts=true`) and compares each against the
ground-truth contours with `compare_pt_clouds`, giving Hausdorff and squared-Chamfer
distances per frame, per iteration. The Chamfer figure is squared so it is directly
comparable to the optimizer's `ChamferCost`; Hausdorff and closest-point stay in pixel
units, so the three series only share an axis because each is normalized by its own first
iterate before plotting.

Iterate `n` is index-aligned with `stats["iterList"]`/`ηList`/`βList`, so `iter_1`
is the initial guess and `iter_2` the first accepted step. Frames are truncated to
the shorter of the observed and simulated sequences.

# Arguments
- `data_type::String`: `"synthetic"`, `"physical"` or `"simulated"` — selects where
  the ground-truth contours are read from under `filepath_gt` (same sources as
  [`_get_borders`](@ref)).
- `filepath_gt::String`: ground-truth data root.
- `exp_path::String`: experiment result directory containing `data/opt_data`.

# Keyword Arguments
- `view_folder::String`: camera-view subdirectory for the ground-truth contours
  (default: `"view_1"`).
- `sim_view::String`: camera-view subdirectory for the simulated contours
  (default: `"view_1"`).
- `write_results::Bool`: also write the per-iteration series to
  `<exp_path>/data/opt_data/metrics/` (default: `true`).
- `plot_results::Bool`: also plot the time-averaged (frame-averaged) Chamfer,
  Hausdorff and closest-point distances against iteration, together on
  `<exp_path>/plots/contour_metrics_iter.pdf` (default: `true`). Each curve is
  normalized by its own value at the first iterate, so all three start at 1 and the
  axis reads as the fraction of the initial cost remaining; this affects the plot
  only, and the returned `Dict` and the written CSVs keep the raw px distances.
  Uses a log y-axis when every value is positive, falling back to linear otherwise.

# Returns
- `Dict` with `"iters"` (iteration indices) and, per iteration,
  `"hausdorff_mean"`, `"hausdorff_max"`, `"chamfer_mean"`, `"chamfer_max"` (squared),
  `"closest_pt_mean"`, plus the per-frame series `"hausdorff_frames"` and
  `"chamfer_frames"` (one vector per iteration).
"""
function iteration_contour_metrics(data_type::String, filepath_gt::String, exp_path::String;
                                   view_folder::String="view_1", sim_view::String="view_1",
                                   write_results::Bool=true, plot_results::Bool=true)

    if data_type == "synthetic"
        ObsDataList, _, _ = read_csv(datapath(filepath_gt,"img_data",view_folder,"contour_data"))
    elseif data_type == "physical"
        # Physical ground truth is single-camera and stores contours flat (no view_* level).
        ObsDataList, _, _ = read_csv(datapath(filepath_gt,"img_data","contour_data"))
    elseif data_type == "simulated"
        ObsDataList, _, _ = read_csv(datapath(filepath_gt,"sim_data",view_folder,"contour_data"))
    else
        error("Unknown data type: $data_type")
    end
    @info "Read $(length(ObsDataList)) ground-truth contour frames for data type $data_type"

    opt_root = datapath(exp_path, "opt_data")
    isdir(opt_root) || throw(SystemError("No opt_data directory at $opt_root"))

    # "iter_10" sorts before "iter_2" lexicographically, so order by the parsed index.
    iter_dirs = filter(d -> occursin(r"^iter_[0-9]+$", d), readdir(opt_root))
    isempty(iter_dirs) && throw(ArgumentError(
        "No iter_* directories under $opt_root — was the fit run with store_border_pts=true?"))
    iters = sort(parse.(Int, replace.(iter_dirs, "iter_" => "")))

    hausdorff_mean = Float64[]; hausdorff_max = Float64[]
    chamfer_mean   = Float64[]; chamfer_max   = Float64[]
    closest_pt_mean = Float64[]
    hausdorff_frames = Vector{Vector{Float64}}()
    chamfer_frames   = Vector{Vector{Float64}}()

    for n in iters
        sim_pts, _, _ = read_csv(datapath(exp_path,"opt_data","iter_$n",sim_view,"2D_border_points"))
        nf = min(length(sim_pts), length(ObsDataList))
        nf > 0 || throw(ArgumentError("No overlapping frames for iteration $n"))
        h, c, cp = compare_pt_clouds(sim_pts[1:nf], ObsDataList[1:nf])

        push!(hausdorff_frames, h);      push!(chamfer_frames, c)
        push!(hausdorff_mean, mean(h));  push!(hausdorff_max, maximum(h))
        push!(chamfer_mean, mean(c));    push!(chamfer_max, maximum(c))
        push!(closest_pt_mean, mean(cp))
        @info "iter $n ($nf frames): chamfer(mean)=$(round(mean(c), sigdigits=4)), " *
              "hausdorff(mean)=$(round(mean(h), sigdigits=4)), hausdorff(max)=$(round(maximum(h), sigdigits=4))"
    end

    metrics = Dict{String,Any}(
        "iters"            => iters,
        "hausdorff_mean"   => hausdorff_mean,
        "hausdorff_max"    => hausdorff_max,
        "chamfer_mean"     => chamfer_mean,
        "chamfer_max"      => chamfer_max,
        "closest_pt_mean"  => closest_pt_mean,
        "hausdorff_frames" => hausdorff_frames,
        "chamfer_frames"   => chamfer_frames,
    )

    if write_results
        for k in ("iters","hausdorff_mean","hausdorff_max","chamfer_mean","chamfer_max","closest_pt_mean")
            write_csv(datapath(exp_path,"opt_data","metrics",k), metrics[k])
        end
        # per-frame series: one file per iteration
        for (i, n) in enumerate(iters)
            write_csv(datapath(exp_path,"opt_data","metrics","hausdorff_frames","iter_$n"), hausdorff_frames[i])
            write_csv(datapath(exp_path,"opt_data","metrics","chamfer_frames","iter_$n"), chamfer_frames[i])
        end
        @info "Wrote iteration contour metrics to $(datapath(exp_path,"opt_data","metrics"))"
    end

    if plot_results
        set_file(plotpath(exp_path))

        # Normalize each metric by its own value at the first iterate (the initial guess),
        # so every curve starts at 1 and the plot reads as "fraction of the initial cost
        # remaining". This is presentation-only: `metrics` and the written CSVs keep the
        # raw distances in px. A non-positive or non-finite reference cannot be divided
        # out meaningfully, so that series is left in raw units instead.
        function normalize_to_initial(v::Vector{Float64}, name::String)
            ref = isempty(v) ? 0.0 : v[1]
            if !(isfinite(ref) && ref > 0)
                @warn "Initial $name is $ref; plotting $name unnormalized"
                return v
            end
            return v ./ ref
        end
        chamfer_rel     = normalize_to_initial(chamfer_mean, "Chamfer")
        hausdorff_rel   = normalize_to_initial(hausdorff_mean, "Hausdorff")
        closest_pt_rel  = normalize_to_initial(closest_pt_mean, "closest-point")

        # Metrics decay over several decades as the fit converges, so prefer a log axis;
        # an exact zero (or a single iterate) would make that invalid.
        use_log = all(>(0), chamfer_rel) && all(>(0), hausdorff_rel) && all(>(0), closest_pt_rel)
        yscale = use_log ? :log10 : :identity

        # The three metrics can sit almost on top of each other (a near-uniform contour
        # offset makes mean and max nearest-neighbour distance nearly equal), so vary the
        # line style as well as the colour to keep them separable.
        metric_plt = default_plot()
        Plots.plot!(metric_plt, iters, chamfer_rel, label=L"\mathrm{Chamfer}^2", marker=1,
                    linestyle=:solid, yscale=yscale, xminorgrid=:false,
                    legend=:outerbottom, legend_column=3)
        Plots.plot!(metric_plt, iters, hausdorff_rel, label=L"\mathrm{Hausdorff}", marker=1,
                    linestyle=:dash, yscale=yscale, legend=:outerbottom, legend_column=3)
        Plots.plot!(metric_plt, iters, closest_pt_rel, label=L"\mathrm{Closest\;point}", marker=1,
                    linestyle=:dot, yscale=yscale, legend=:outerbottom, legend_column=3)
        _label!(metric_plt, L"\mathrm{Iterations}", L"d^{[\imath]}/d^{[1]}")
        Plots.savefig(metric_plt, plotpath(exp_path,"contour_metrics_iter.pdf"))
        @info "Wrote $(plotpath(exp_path,"contour_metrics_iter.pdf"))"
    end

    return metrics
end

"""
    get_surface_mosd(surface_list, cam_pose, height; tol=1e-3) -> Vector{Float64}

Max/Mean Of Surface Depth (MOSD): the depth (z, along the camera's
viewing axis) of the furthest surface points, per frame.

Each surface is transformed into the camera frame via `project_to_camera_frame`,
then reduced with `_max_band_mean` so that flat faces/edges facing the
camera (e.g. a cubic mesh) are averaged over rather than reduced to a single
mesh vertex.

# Arguments
- `surface_list::AbstractArray`: list of 3×N surface point matrices, one per frame.
- `cam_pose::AbstractArray`: camera pose used to build the camera frame.
- `height::Float64`: current specimen height, passed through to `project_to_camera_frame`.

# Keyword Arguments
- `tol::Float64`: band width (in the same units as the mesh) around the max
  depth over which points are averaged (default: `1e-3`).

# Returns
- `mosd::Vector{Float64}`: one MOSD value per frame.
"""
function get_surface_mosd(surface_list::AbstractArray, cam_pose::AbstractArray, height::Float64; tol::Float64=1e-3)
    mosd = Vector{Float64}(undef, length(surface_list))

    for (i, sim_surface) in enumerate(surface_list)
        transformed_sim_surface = project_to_camera_frame(sim_surface, cam_pose, height)
        mosd[i] = _max_band_mean(transformed_sim_surface[3, :], tol)
    end

    return mosd
end

"""
    _max_band_mean(z, tol) -> Float64

Mean of all entries of `z` within `tol` of `maximum(z)`.

Used instead of a hard `maximum` so that a flat edge/face (many points tying
for the max within meshing/floating-point tolerance) doesn't collapse to a
single, arbitrarily-chosen point, and doesn't produce a MOSD that jumps
discretely between mesh nodes as the geometry deforms.

# Arguments
- `z::AbstractVector`: values to reduce.
- `tol::Float64`: band width around `maximum(z)` over which entries are averaged.

# Returns
- `Float64`: mean of all entries of `z` within `tol` of `maximum(z)`.
"""
function _max_band_mean(z::AbstractVector, tol::Float64)
    z_max = maximum(z)
    band = z[z .≥ z_max - tol]
    return sum(band) / length(band)
end

"""
    _handle_worker_error(err, i, params)

Log a failed parallel worker task (index `i`, its `params`, and the caught
exception `err` with backtrace) via `@error`, for use in a `catch` block.

# Arguments
- `err`: the caught exception.
- `i::Int`: index of the failed task.
- `params`: parameters passed to the failed task, logged for diagnosis.

# Returns
None.
"""
function _handle_worker_error(err, i::Int, params)
    bt = catch_backtrace()
    @error "Task $i failed" params exception=(err, bt)
end

"""
    run_window_predictions(model, scene, conditions, est_ηpList, est_βpList,
                           data_ranges_, time_windows, F, h)
        -> (h_pred_vec, pos3D_pred_vec, surface_pts_3D_pred_vec,
            fields_pred_vec, pos2D_pred_vec, borderPts_pred_vec)

Per-window forward prediction for a bulk-viscosity experiment, shared by
`optimize` (called inline, with the window estimates still in memory) and
`predict` (called from disk after reading the estimates back).

For each time window `ti` the within-window estimate `(est_ηpList, est_βpList)`
is simulated forward. From the second window on, the *previous* window's
estimate is additionally simulated forward into the current window — the
`k → k+1` prediction: parameters identified in window `k` are used to predict
the dynamics of window `k+1`. Window 1 has no predecessor, so its estimation
run stands in for its prediction.

`model`/`scene` are the (constant-viscosity) simulation model and are mutated in
place: `reset_model!` before each run returns to the committed reference, and
`update_model!` after the estimation run commits the deformed configuration so
the next window starts from the current specimen height (`_h` carries the end
height forward).

# Arguments
- `model`, `scene`, `conditions`: constant-viscosity simulation problem and observation conditions.
- `est_ηpList`, `est_βpList`: per-frame estimated η/β (constant within each window).
- `data_ranges_`: per-window frame index ranges (`data_ranges_[ti]`).
- `time_windows`: per-window simulation time.
- `F`: control parameter (force/displacement) time series.
- `h`: initial specimen height.

# Returns
`vcat`-reduced trajectories over all windows: predicted height, 3D points, 3D
surface points, motion fields, 2D points, and 2D border points.
"""
function run_window_predictions(model, scene, conditions, est_ηpList, est_βpList,
                                data_ranges_, time_windows, F, h)
    h_pred_list = AbstractArray[]
    borderPts_pred_list = AbstractArray[]
    fields_list_pred_list = AbstractArray[]
    pos2D_pred_list = AbstractArray[]
    pos3D_pred_list = AbstractArray[]
    surface_pts_3D_pred_list = AbstractArray[]

    _h = h
    for ti::Int in 1:(size(data_ranges_, 1))

        data_range_ = data_ranges_[ti]
        scene.sim_time = time_windows[ti]

        η = est_ηpList[data_range_][1] # first value in the window; the rest are identical
        β = est_βpList[data_range_][1]
        scene.cParam = F[data_range_]

        printstyled("Time window: $(ti), time frames: $(scene.sim_time)\n"; color = :blue)
        @debug "Data frame : $(data_range_)"

        if ti > 1
            @info "Predicting dynamics in time window $(ti) from window $(ti-1)'s estimate..."
            reset_model!(model)
            data_range_prev = data_ranges_[ti-1]
            η_pred = est_ηpList[data_range_prev][1]
            β_pred = est_βpList[data_range_prev][1]

            model.η = [η_pred]
            scene.β = [β_pred]
            pred_μ_list, gradList, pred_simBorderPts, fields_pred, surface_pts_3D_pred, pos2D_pred, pos3D_pred, _, _, _ = simulate(model, scene, conditions)
            pred_h = get_height(pred_μ_list, _h)

            push!(h_pred_list, pred_h)
            push!(borderPts_pred_list, pred_simBorderPts)
            push!(fields_list_pred_list, fields_pred)
            push!(pos2D_pred_list, pos2D_pred)
            push!(pos3D_pred_list, pos3D_pred)
            push!(surface_pts_3D_pred_list, surface_pts_3D_pred)
        end

        reset_model!(model)
        model.η = [η]
        scene.β = [β]
        est_μ_list, gradList, est_simBorderPts, fields_est, surface_pts_3D_est, pos2D_est, pos3D_est, _, _, _ = simulate(model, scene, conditions)

        h_est = get_height(est_μ_list, _h)
        if ti == 1
            push!(h_pred_list, h_est)
            push!(borderPts_pred_list, est_simBorderPts)
            push!(fields_list_pred_list, fields_est)
            push!(pos2D_pred_list, pos2D_est)
            push!(pos3D_pred_list, pos3D_est)
            push!(surface_pts_3D_pred_list, surface_pts_3D_est)
        end
        update_model!(model)
        _h = h_est[end] # carry the end height into the next window
    end

    h_pred_vec              = reduce(vcat, h_pred_list)
    pos3D_pred_vec          = reduce(vcat, pos3D_pred_list)
    surface_pts_3D_pred_vec = reduce(vcat, surface_pts_3D_pred_list)
    fields_pred_vec         = reduce(vcat, fields_list_pred_list)
    pos2D_pred_vec          = reduce(vcat, pos2D_pred_list)
    borderPts_pred_vec      = reduce(vcat, borderPts_pred_list)

    return h_pred_vec, pos3D_pred_vec, surface_pts_3D_pred_vec, fields_pred_vec, pos2D_pred_vec, borderPts_pred_vec
end

"""
    write_window_predictions(dest_data_dir, h_pred_vec, pos3D_pred_vec,
                             surface_pts_3D_pred_vec, fields_pred_vec,
                             pos2D_pred_vec, borderPts_pred_vec)

Write the per-window prediction trajectories returned by
`run_window_predictions` into `dest_data_dir` (an experiment's `data`
directory), matching the layout consumed by `replot`.
"""
function write_window_predictions(dest_data_dir, h_pred_vec, pos3D_pred_vec,
                                  surface_pts_3D_pred_vec, fields_pred_vec,
                                  pos2D_pred_vec, borderPts_pred_vec)
    write_csv(joinpath(dest_data_dir, "pred_h"), h_pred_vec)
    write_data(joinpath(dest_data_dir, "sim_data", "3D_points_pred"), pos3D_pred_vec)
    write_data(joinpath(dest_data_dir, "sim_data", "3D_surface_points_pred"), surface_pts_3D_pred_vec)
    write_data(joinpath(dest_data_dir, "sim_data", "motion_fields_pred "), fields_pred_vec)
    write_2d_data(joinpath(dest_data_dir, "sim_data", "2D_points_pred"), pos2D_pred_vec)
    write_2d_data(joinpath(dest_data_dir, "sim_data", "2D_border_points_pred"), borderPts_pred_vec)
end

"""
    optimize(exp_params)

Run the full parameter-estimation pipeline for one experiment configuration:
build the ground-truth (or physical) model and observation data described by
`exp_params`, fit viscosity `η` (and slip `β`) to the observed contours via
Gauss-Newton (`fit_model`), then re-simulate with the estimated parameters and
write all results (fitted parameters, cost/iteration history, 2D/3D fields,
heights) to `exp_params["filepath_res"]`.

Supports `exp_params["data_type"] ∈ ("simulated"|"synthetic", "physical")`,
`viscosity_type ∈ ("constant", "bulk_viscosity")` (the latter optionally
windowed in time via `exp_params["mode"]`), and multiple camera view angles
(`exp_params["z_angle_list"]` in the ground-truth data).

# Arguments
- `exp_params::Dict`: experiment configuration (data/viscosity type, mesh,
  camera, and output-path parameters).

# Returns
None. All outputs are side effects written to disk.
"""
function optimize(exp_params::Dict)

    η_gt::Vector{Float64} = [0.0]
    β_gt::Vector{Float64} = [0.0]
    F_ext::Float64 = 0.0
    sim_time_gt::Float64 = 0.0
    t_steps_gt::Float64 = 0.0
    steps_exp::Float64 = 0.0
    outlier_frames::Vector{Int} = Int[]
    time_window_type::String = ""

    ndim::Int = 3
    nDof_p::Int = 1
    nDof_u::Int = ndim

    element_shape_u = Symbol(exp_params["element_shape_u"])
    basis_order_u   = Int(exp_params["basis_order_u"])
    element_shape_p = Symbol(exp_params["element_shape_p"])
    basis_order_p   = Int(exp_params["basis_order_p"])
    element_shape_x = Symbol(exp_params["element_shape_x"])
    basis_order_x   = Int(exp_params["basis_order_x"])

    _mesh_path_kw(d) = haskey(d, "mesh_path") ? (mesh_path=d["mesh_path"],) : (;)

    # Optimization method and the `fit_model` keyword arguments that go with it, supplied
    # by the driver (`optimize_sim`/`optimize_syn`/`optimize_real`) through `exp_params`.
    # Absent keys mean plain Gauss-Newton with no extra kwargs, which is what every
    # experiment written before 2026-08-26 used.
    #
    # `exp_params` reaches here either as native Julia (in-process or via
    # `run_parallel_tasks`) or round-tripped through `experiment_parameters.json`, so
    # coerce the types a JSON round-trip changes: Symbols come back as `String`, and a
    # `λ_scale` written as `1` comes back `Int`, which `fit_model`'s `::Float64` would
    # reject with a MethodError.
    opt_method::Symbol = Symbol(get(exp_params, "opt_method", :gn))

    # Which contour cost drives the fit. `:closest_point` is the historical default, so an
    # experiment that does not name one behaves exactly as it did before the cost became
    # configurable.
    cost_function::Symbol = Symbol(get(exp_params, "cost_function", :closest_point))
    contour_cost_fn = _contour_cost(cost_function)

    function _opt_kwargs(d)
        kw = Dict{Symbol,Any}(Symbol(k) => v for (k, v) in get(d, "opt_kwargs", Dict{String,Any}()))
        haskey(kw, :line_search_method) && (kw[:line_search_method] = Symbol(kw[:line_search_method]))
        haskey(kw, :λ_scale)            && (kw[:λ_scale] = float(kw[:λ_scale]))
        haskey(kw, :θ_p)                && (kw[:θ_p] = Vector{Float64}(kw[:θ_p]))
        # Γ takes a spec rather than a matrix: `:relative` (drop it, so the optimizer's own
        # default diag(1 ./ |θ|) applies) or `:identity`. A raw nested array still works as
        # an escape hatch — JSON hands matrices back row-wise, hence the permutedims.
        if haskey(kw, :Γ)
            spec = kw[:Γ]
            if spec isa Symbol || spec isa String
                if Symbol(spec) === :relative
                    delete!(kw, :Γ)
                elseif Symbol(spec) === :identity
                    kw[:Γ] = Matrix{Float64}(I, 2, 2)
                else
                    error("Unknown Γ spec :$spec — use :relative, :identity, or a matrix")
                end
            else
                kw[:Γ] = permutedims(reduce(hcat, [Vector{Float64}(r) for r in spec]))
            end
        end
        return (; kw...)
    end
    opt_kw = _method_kwargs(opt_method, (; _opt_kwargs(exp_params)..., cost=contour_cost_fn))
    @info "Optimization method: $opt_method" cost_function opt_kw

    _fem = (element_shape_u, basis_order_u, nDof_u, element_shape_p, basis_order_p, nDof_p, element_shape_x, basis_order_x)

    iterList::Vector{Float64} = Vector{Float64}()
    costList::Vector{Float64}  = Vector{Float64}()
    ηpList::Vector{Float64} = Vector{Float64}()
    βpList::Vector{Float64} = Vector{Float64}()

    # simulation parameters for the ground truth
    start_time = Dates.now()
    filepath_gt::String = exp_params["filepath_gt"]
    filepath_res::String = exp_params["filepath_res"]

    obj_pose::AbstractArray{Float64} = zeros(Float64, 3)
    camera_matrix::AbstractArray = exp_params["camera_matrix"]

    sim_time_exp::Float64 = exp_params["sim_time_exp"]
    t_steps_exp::Float64 = if haskey(exp_params, "dt")
        steps_exp = sim_time_exp / Float64(exp_params["dt"])
        Float64(exp_params["dt"])
    elseif haskey(exp_params, "steps_exp")
        steps_exp = Float64(exp_params["steps_exp"])
        sim_time_exp / steps_exp
    else
        steps_exp = sim_time_exp / 0.1
    end

    ne_exp::Int = exp_params["ne_exp"]

    r::Float64 = 0.0
    h::Float64 = 0.0
    data_type ::String = exp_params["data_type"] # "simulated" or "physical" or "real"
    geometry::Union{Symbol, Nothing} = nothing
    z_angles::Vector{Float64} = [0.0]
    edge_radius::Union{Float64,Nothing} = nothing
    noiseLevel::Float64 = 0.0
    viscosity_model::String = "" # "constant" or "bulk_viscosity"

    obs_data_angles = AbstractArray[]
    splinex_angles = AbstractArray[]
    spliney_angles = AbstractArray[]

    if data_type == "simulated" || data_type == "synthetic"

        WRITE_GT = exp_params["WRITE_GT"] 
        noiseLevel = exp_params["noise_level"]
        outlier_frames = Int[1]

        sim_params = read_json(datapath(filepath_gt,"sim_params"))

        r = sim_params["r"]
        h = sim_params["h"]
        geometry = get!(sim_params, "geometry", :cylinder)
        z_angles = get!(sim_params, "z_angle_list", [0.0])
        edge_radius = get!(sim_params, "edge_radius", nothing)
        gt_viscosity_type = sim_params["viscosity_type"]
        F = Array(float.(sim_params["cParam"]))

        sim_time_gt = sim_params["simulation_time"]
        t_steps_gt = sim_params["time_steps"]

        camera_matrix = _read_camera_matrix(sim_params)
        obj_pose = _read_obj_pose(sim_params)

        control = sim_params["control_type"]

        if WRITE_GT == true
            @info "Writing ground truth gt data to with $ne_exp elements to $exp_path"
            write_gt_data(exp_params)
        end

        if gt_viscosity_type == "bulk_viscosity"
            time_window_type = get(exp_params, "time_windows", "linear")
            if haskey(sim_params, "model_type") && sim_params["model_type"] == "carreau"
                viscosity_model = sim_params["model_type"]
            else
                viscosity_model = "power_law"
                η_gt = sim_params["η"]
                β_gt = sim_params["β"]
            end
        else
            η_gt = sim_params["η"]
            β_gt = sim_params["β"]
        end

        _shared = (η_gt[1], element_shape_u, basis_order_u, nDof_u, element_shape_p, basis_order_p, nDof_p, element_shape_x, basis_order_x, β_gt[1], F, control, gt_viscosity_type, sim_time_gt, t_steps_gt)
        model_gt, scene_exp = if geometry === :cylinder
                                def_problem(Cylinder(r, h), ne_exp, _shared...; _mesh_path_kw(exp_params)...)
                            elseif geometry === :cube
                                lx = Float64(get(exp_params, "lx", r))
                                ly = Float64(get(exp_params, "ly", r))
                                lz = Float64(get(exp_params, "lz", h))
                                def_problem(Cuboid(lx, ly, lz, edge_radius), ne_exp, _shared...; _mesh_path_kw(exp_params)...)
                            else
                                error("Unsupported geometry type: $geometry")
                            end

        for i::Int in 1:length(z_angles)
            if data_type == "synthetic"
                ObsDataList, splinexObs, splineyObs = read_csv(datapath(filepath_gt,"img_data","view_$i","contour_data"))
            elseif data_type == "simulated"
                ObsDataList, splinexObs, splineyObs = read_csv(datapath(filepath_gt,"sim_data","view_$i","contour_data"))  
            end

            push!(obs_data_angles, ObsDataList)
            push!(splinex_angles, splinexObs)
            push!(spliney_angles, splineyObs)
        end

    elseif data_type == "physical"
        r = exp_params["r"]  # radius of the cylinder in mm
        h = exp_params["h"]  # height of the cylinder in mm

        control::String = exp_params["control"]
        gt_viscosity_type::String = exp_params["viscosity_type"] # "constant" or "bulk_viscosity"

        t_obs, t_top_plt, t_btm_plt = read_perception_data(datapath(filepath_gt, "sequence.hdf5"))

        _ObsDataList, _splinexObs, _splineyObs = read_csv(datapath(filepath_gt, "img_data", "contour_data"))
        meta_data = read_json(datapath(filepath_gt, "video_metadata"))

        frame_rate::Float64 = round(meta_data["frame_rate"], digits=1)
        frame_width::Int = meta_data["frame_width"]
        frame_height::Int = meta_data["frame_height"]
        compression_frames = meta_data["compressed_frames"]

        sim_time_gt = length(compression_frames)/frame_rate  # seconds
        steps_exp = sim_time_exp*frame_rate
        t_steps_exp = 1/frame_rate
        t_steps_gt = t_steps_exp
        time_window_type = get(exp_params, "time_windows", "quadratic")

        F_ext = exp_params["F_ext"]
        @debug "Applied force: $F_ext N"

        F = -F_ext*ones(Float64, round(Int, sim_time_gt*frame_rate)) # force applied to the cylinder in N

        # Keep the measured 4×4 as-is: its rotation carries a ~1-2° camera tilt, worth
        # tens of pixels at this standoff, which a synthesised camera frame cannot express.
        obj_pose = get_pose(t_obs)

        push!(obs_data_angles, _ObsDataList[compression_frames])
        push!(splinex_angles, _splinexObs[compression_frames])
        push!(spliney_angles, _splineyObs[compression_frames])

        function get_gt_height(t_top_plt::AbstractArray, t_btm_plt::AbstractArray, h_gt_init::Float64)
            n_frames = size(t_top_plt,2)
            Δh_gt = zeros(Float64, n_frames)
            h_init = abs(t_top_plt[2,1] - t_btm_plt[2,1])
            for i::Int in 1:n_frames
                Δh_gt[i] = h_init - abs(t_top_plt[2,i] - t_btm_plt[2,i])
            end

            h_gt = h_gt_init .- Δh_gt

            return h_gt, Δh_gt
        end

        h_gt , _ = get_gt_height(t_top_plt[1:3,4,compression_frames], t_btm_plt[1:3,4,compression_frames], 38.5)

        sim_data = Dict("camera_matrix" => camera_matrix,
                    "obj_pose" => obj_pose,
                    "time_steps" => t_steps_exp,
                    "simulation_time" => sim_time_gt,
                    "cParam" => F,
                    "r" => r,
                    "h" => h,
                    "control_type" => control,
                    "viscosity_type" => gt_viscosity_type,
                    "data_type" => data_type)

        write_json(datapath(filepath_gt,"sim_params"), sim_data)
        write_csv(datapath(filepath_gt,"h"), h_gt)

        valid_frames, _outlier_frames = detect_outlier_observations(_ObsDataList[compression_frames])
        outlier_frames = vcat(outlier_frames, _outlier_frames)
        geometry = :cylinder
    else
        error("data_type should be either simulated or physical")
    end

    geom_exp = if geometry === :cylinder
        Cylinder(r, h)
    elseif geometry === :cube
        lx = Float64(get(exp_params, "lx", r))
        ly = Float64(get(exp_params, "ly", r))
        lz = Float64(get(exp_params, "lz", h))
        Cuboid(lx, ly, lz, edge_radius)
    else
        error("Unsupported geometry type: $geometry")
    end

    if sim_time_gt < sim_time_exp
        @warn "Ground truth simulation time $sim_time_gt is less than experimental simulation time $sim_time_exp , switching to ground truth simulation time"
        sim_time_exp = sim_time_gt
    end

    if t_steps_gt > t_steps_exp
        @warn "time resolution of the ground truth $t_steps_gt is larger than the experimental $t_steps_exp, switching to ground truth resolution"
        t_steps_exp = t_steps_gt
    end

    for (i,z_angle) in enumerate(z_angles)

        printstyled("Processing view $i with z_angle = $z_angle degrees\n"; color = :blue)
        if data_type == "physical"  || viscosity_model == "carreau"
            η_start = exp_params["η_start"]
            β_start = exp_params["β_start"]
        else
            dev::Float64 = 0.3

            dev_η::Float64 = dev*η_gt[1]
            η_start::Float64 = abs(η_gt[1] - dev_η)

            if β_gt[1] >= 200.0 # setting the slip to partial slip for no slip cases
                β_start = 50.0
            elseif β_gt[1] <= 1 # setting the slip to partial slip for free slip cases
                β_start = 10.0
            else
                dev_β::Float64 = dev*β_gt[1]
                β_start::Float64 = abs(β_gt[1] - dev_β)
            end
        end

        θ::Vector{Float64} = [η_start, β_start]

        ObsDataList = obs_data_angles[i]
        splinexObs = splinex_angles[i]
        splineyObs = spliney_angles[i]

        exp_params["z_angle"] = z_angle
        if data_type != "physical"
            exp_params["num_ne"] = model_gt.mesh_u.ne
        end
        exp_path = joinpath(dirname(filepath_res), "view_$i", basename(filepath_res))
        write_json(datapath(exp_path,"experiment_parameters"), exp_params)

        if gt_viscosity_type == "constant"

            printstyled("Ground truth η: $(η_gt), ground truth β: $(β_gt)\n"; color = :green)
            _range = collect(range(start = 1, stop = (round(Int,sim_time_exp/t_steps_gt)+1), step = round(Int,t_steps_exp/t_steps_gt)))
            @debug "Step progress: $(round(Int,sim_time_exp/t_steps_gt)+1) / $(round(Int,t_steps_exp/t_steps_gt)) / $(length(_range))"
            @info "Considering from frame $(first(_range)) to frame $(last(_range)) in the observations"
            ObsDataList = ObsDataList[_range] # align the observation points with the simulation time
            printstyled("Observation data length: $(length(ObsDataList))\n"; color = :blue)

            model, scene = def_problem(geom_exp, ne_exp, η_gt[1], _fem..., β_gt[1], F, control, gt_viscosity_type, sim_time_exp, t_steps_exp; _mesh_path_kw(exp_params)...)
            est_model, est_scene = def_problem(geom_exp, ne_exp, η_start, _fem..., β_start, F, control, gt_viscosity_type, sim_time_gt, t_steps_exp; _mesh_path_kw(exp_params)...)
            conditions = Conditions(camera_matrix=camera_matrix, obj_pose=obj_pose, filepath=exp_path, ANIMATE=false, viewing_angles=[z_angle])

            gt_h_ = readdlm(datapath(filepath_gt,"h.csv"), ',', Float64)
            gt_h = gt_h_[1:(round(Int,sim_time_gt/t_steps_exp)+1)]

            if noiseLevel == 0.0

                obs_border_pt_lst, nSplinex, nSpliney, pd = add_noise(ObsDataList, nFactor=0.0)
                stats = fit_model(model, scene, conditions, obs_border_pt_lst, θ; outliers=outlier_frames,
                                  method=opt_method, store_border_pts=true, opt_kw...)

                iterList = stats["iterList"]
                costList = stats["cost_list"]
                ηpList = stats["ηList"]
                βpList = stats["βList"]

                η = stats["η"]
                β = stats["β"]

                printstyled("Estimated η : $(η), estimated β: $(β)\n"; color = :green)

                η_accuracy = (1-abs((η_gt[1]-η)/η_gt[1]))*100
                β_accuracy = (1-abs((β_gt[1]-β)/β_gt[1]))*100
                printstyled("η accuracy: $(η_accuracy) %\n"; color = :green)
                printstyled("β accuracy: $(β_accuracy) %\n"; color = :green)

                set_file(plotpath(exp_path))

                write_json(datapath(exp_path,"stats"), stats)
                write_csv(datapath(exp_path,"η"), ηpList)
                write_csv(datapath(exp_path,"β"), βpList)
                write_csv(datapath(exp_path,"gt_h"), gt_h)
                write_csv(datapath(exp_path,"cost_iter"), costList)

                for (n, iter_border_pts) in enumerate(stats["simBorderPtsList"])
                    write_2d_data(datapath(exp_path,"opt_data","iter_$n","2D_border_points"), iter_border_pts)
                end

                reset_model!(est_model)
                est_model.η = [η]
                est_scene.β = [β]

                conditions = Conditions(camera_matrix=camera_matrix, obj_pose=obj_pose, filepath=exp_path, ANIMATE=true, viewing_angles=[z_angle])
                est_μ_list, gradList, borderPts2DList, fields, surface_pts_3D, pos2D, pos3D, _, _, _ = simulate(est_model, est_scene, conditions)
                est_h = get_height(est_μ_list, h)

                write_csv(datapath(exp_path,"est_h"), est_h)
                write_data(datapath(exp_path,"sim_data","3D_point_est"), pos3D)
                write_data(datapath(exp_path,"sim_data","motion_fields_est "), fields)
                write_data(datapath(exp_path,"sim_data","3D_surface_points_est"), surface_pts_3D)

                write_2d_data(datapath(exp_path,"sim_data","2D_surface_points_est"), pos2D)
                write_2d_data(datapath(exp_path,"sim_data","2D_border_points_est"), borderPts2DList)

                if maximum(ηpList) > η+dev_η
                    ηStop = maximum(ηpList)*1.1
                else
                    ηStop = η+dev_η
                end

                if minimum(ηpList) < η-dev_η
                    if minimum(ηpList) <= 0.0
                        η_start = 0.1
                    else
                        η_start = minimum(ηpList)*0.9
                    end
                else
                    η_start = η-dev_η
                end

            else
                n_samples = 10
                η_pred = zeros(Float64, n_samples)
                β_pred = zeros(Float64, n_samples)
                costnList = Vector{Vector{Float64}}(undef, n_samples)
                iternList = Vector{AbstractVector}(undef, n_samples)
                est_h_list = Matrix{Float64}(undef, n_samples, round(Int,sim_time_gt/t_steps_exp)+1)
                ANIMATED = false

                # Pre-allocate containers for batch file I/O (collect phase)
                sim_params_list = Vector{Dict}(undef, n_samples)
                pos2D_list = Vector{Any}(undef, n_samples)
                pos3D_list = Vector{Any}(undef, n_samples)
                fields_list = Vector{Any}(undef, n_samples)
                borderPts2D_list = Vector{Any}(undef, n_samples)
                η_steps_list = Vector{Vector}(undef, n_samples)
                β_steps_list = Vector{Vector}(undef, n_samples)
                cost_steps_list = Vector{Vector}(undef, n_samples)
                iter_steps_list = Vector{Vector}(undef, n_samples)

                set_file(plotpath(exp_path))
                for n::Int in 1:n_samples
                    obs_border_pt_lst, nSplinex, nSpliney, pd = add_noise(ObsDataList, nFactor=noiseLevel)

                    stats = fit_model(model, scene, conditions, obs_border_pt_lst, θ; outliers=outlier_frames,
                                      method=opt_method, opt_kw...)

                    iterList = stats["iterList"]
                    costList = stats["cost_list"]
                    ηpList = stats["ηList"]
                    βpList = stats["βList"]

                    η = stats["η"]
                    β = stats["β"]

                    reset_model!(est_model)
                    est_model.η = [η]
                    est_scene.β = [β]

                    est_μ_list, gradList, borderPts2DList, fields, surface_pts_3D, pos2D, pos3D, _, _, _ = simulate(est_model, est_scene, conditions)

                    est_h = get_height(est_μ_list, h)

                    η_accuracy = (1-abs((η_gt[1]-η)/η_gt[1]))*100
                    β_accuracy = (1-abs((β_gt[1]-β)/β_gt[1]))*100
                    printstyled("η accuracy: $(η_accuracy) %\n"; color = :green)
                    printstyled("β accuracy: $(β_accuracy) %\n"; color = :green)

                    η_pred[n] = η
                    β_pred[n] = β
                    est_h_list[n,:] = est_h
                    costnList[n] = costList
                    iternList[n] = iterList

                    # Collect phase: Store all data in memory instead of writing immediately
                    sim_params = Dict("gt_η" => η_gt,
                    "gt_β" => β_gt,
                    "η" => η,
                    "β" => β,
                    "η_accuracy" => η_accuracy,
                    "β_accuracy" => β_accuracy)

                    sim_params_list[n] = sim_params
                    cost_steps_list[n] = costList
                    η_steps_list[n] = ηpList
                    β_steps_list[n] = βpList
                    iter_steps_list[n] = iterList
                    pos2D_list[n] = pos2D
                    pos3D_list[n] = pos3D
                    fields_list[n] = fields
                    borderPts2D_list[n] = borderPts2DList

                end

                # Write phase: Batch write all collected data at once (after loop completes)
                @info "Beginning batch file write phase ($n_samples runs)..."
                for n::Int in 1:n_samples
                    write_json(datapath(exp_path,"stats","run_$n"), sim_params_list[n])
                    write_csv(datapath(exp_path,"opt_data","cost_steps","run_$n"), cost_steps_list[n])
                    write_csv(datapath(exp_path,"opt_data","eta_steps","run_$n"), η_steps_list[n])
                    write_csv(datapath(exp_path,"opt_data","beta_steps","run_$n"), β_steps_list[n])
                    write_csv(datapath(exp_path,"opt_data","iter","run_$n"), iter_steps_list[n])
                    write_csv(datapath(exp_path,"opt_data","est_height","run_$n"), est_h_list[n,:])
                    write_2d_data(datapath(exp_path,"sim_data","2D_points","run_$n"), pos2D_list[n])
                    write_data(datapath(exp_path,"sim_data","3D_points","run_$n"), pos3D_list[n])
                    write_data(datapath(exp_path,"sim_data","motion_fields ","run_$n"), fields_list[n])
                    write_2d_data(datapath(exp_path,"sim_data","2D_border_points","run_$n"), borderPts2D_list[n])
                end

                write_csv(datapath(exp_path,"eta_est"), η_pred)
                write_csv(datapath(exp_path,"beta_est"), β_pred)
                write_csv(datapath(exp_path,"h_est"), est_h_list)

                @info "Data writing completed. Results saved to $exp_path"
            end

        elseif gt_viscosity_type == "bulk_viscosity"

            av_η::Float64 = 0.0
            obs_border_pt_lst, nSplinex, nSpliney, pd = add_noise(ObsDataList, nFactor=0.0)
            mode::String = exp_params["mode"] # "single_window" or "multiple_window" or "full_time"
            sim_time_window::Float64 = 30.0 # time window size for optimization in seconds
            if data_type != "physical"
                model_gt.η = η_gt
            end

            viscosity_type = "constant"

            conditions = Conditions(camera_matrix=camera_matrix, obj_pose=obj_pose, viewing_angles=[z_angle])
            model, scene = def_problem(geom_exp, ne_exp, η_gt[1], _fem..., β_gt[1], F, control, viscosity_type, sim_time_exp, t_steps_exp; _mesh_path_kw(exp_params)...)

            set_file(plotpath(exp_path))

            time_windows, windows, data_ranges_, t_windows = set_time_window(1/t_steps_exp, obs_border_pt_lst, method=time_window_type, window_size=sim_time_exp)
            _, splinexObs_win, _, _ = set_time_window(1/t_steps_exp, splinexObs, method=time_window_type, window_size=sim_time_exp)
            _, splineyObs_win, _, _ = set_time_window(1/t_steps_exp, splineyObs, method=time_window_type, window_size=sim_time_exp)
            @debug "Time windows: $(time_windows)"
            obs_time = sum(time_windows)

            if obs_time < sim_time_gt
                @warn "Observation time frame $obs_time is less than preset ground truth time frame $sim_time_gt, switching to observation time frame"
                sim_time_gt = obs_time
            end

            if obs_time < sim_time_exp
                @warn "Observation time frame $obs_time is less than experimental simulation time frame $sim_time_exp, switching to observation time frame"
                sim_time_exp = obs_time
            end

            data_pt_len = round(Int,obs_time/t_steps_exp)
            est_ηpList = Vector{Float64}(undef,data_pt_len)
            avg_ηList = Vector{Float64}(undef,data_pt_len)
            est_βpList = Vector{Float64}(undef,data_pt_len)
            cost_list = AbstractArray[]
            pred_h_list = AbstractArray[]

            if mode == "single_window"
                @info "Optimizing over a single time window"
                ti = 1
                data_range_ = data_ranges_[ti]
                scene.sim_time = time_windows[ti]
                if data_type == "physical"
                    _F = -F_ext*ones(Float64, round(Int, scene.sim_time*frame_rate)) # force applied to the cylinder in N
                    scene.cParam = _F
                else
                    scene.cParam = F[data_range_]
                end
                @info "Time window $(time_windows[ti])"

                @debug "Data frame : $(data_range_)"
                @debug "Time frame : $(scene.sim_time)"

                printstyled("Time window: $(ti), time frames: $(scene.sim_time)\n"; color = :blue)

                obs_border_pt_lst_t = windows[ti] # align the observation points with the simulation time
                splinexObs_t = splinexObs_win[ti]
                splineyObs_t = splineyObs_win[ti]

                @debug "observation Window size: $(size(obs_border_pt_lst_t,1)) seconds"

                @debug "Time frame : $data_range_"
                printstyled("Time window: $(ti)\n"; color = :green)

                if data_type != "physical"
                    η_gt_ = model_gt.η[data_range_]
                    av_η = mean(η_gt_)
                    avg_ηList[data_range_] .= av_η
                    printstyled("Average ground truth η in the window: $(av_η), ground truth β: $(β_gt)\n"; color = :green)
                end

                stats = fit_model(model, scene, conditions, obs_border_pt_lst_t, θ; outliers=outlier_frames,
                                  method=opt_method, opt_kw...)

                est_ηpList[data_range_] .= stats["η"]
                est_βpList[data_range_] .= stats["β"]

                θ[1] = stats["η"]
                θ[2] = stats["β"]

                update_model!(model)

                iterList = stats["iterList"]
                costList = stats["cost_list"]
                ηpList = stats["ηList"]
                βpList = stats["βList"]

                viscosity_type = "bulk_viscosity"
                est_model, est_scene = def_problem(geom_exp, ne_exp, η_gt[1], _fem..., β_gt[1], F[data_range_], control, viscosity_type, sim_time_exp, t_steps_exp; _mesh_path_kw(exp_params)...)
                est_model.η = est_ηpList[data_range_] 
                est_scene.β = est_βpList[data_range_]
                est_μ_list, gradList, borderPts2DList, fields_est, _, pos2D_est, pos3D_est, _, _, _ = simulate(est_model, est_scene, conditions)
                est_h_list = get_height(est_μ_list, h)

            else
                @debug "Number of time windows: $(length(windows))"

                # Separate constant-viscosity model for the k -> k+1 prediction,
                # advanced in lockstep so `_h` carries height across windows without
                # disturbing the fitting `model`.
                pred_model, pred_scene = def_problem(geom_exp, ne_exp, η_gt[1], _fem..., β_gt[1], F, control, "constant", sim_time_exp, t_steps_exp; _mesh_path_kw(exp_params)...)
                _h_pred = h
                h_pred_list = AbstractArray[]
                borderPts_pred_list = AbstractArray[]
                fields_list_pred_list = AbstractArray[]
                pos2D_pred_list = AbstractArray[]
                pos3D_pred_list = AbstractArray[]
                surface_pts_3D_pred_list = AbstractArray[]

                for ti::Int in 1:length(windows)
                    data_range_ = data_ranges_[ti]
                    data_range_prev = ti > 1 ? data_ranges_[ti-1] : 1:data_range_[1]

                    scene.sim_time = time_windows[ti]
                    if data_type == "physical"
                        _F = -F_ext*ones(Float64, round(Int, scene.sim_time*frame_rate)) # force applied to the cylinder in N
                        scene.cParam = _F
                    else
                        scene.cParam = F[data_range_]
                    end
                    obs_border_pt_lst_t = windows[ti] # align the observation points with the simulation time

                    printstyled("Time window: $(ti), time frames: $(scene.sim_time)\n"; color = :blue)
                    @debug "Data frame : $(data_range_)"
                    @debug "Time frame : $(scene.sim_time)"
                    @info "Time window $(t_windows[ti])"
                    @debug "observation Window size: $(size(obs_border_pt_lst_t,1)) seconds"

                    if data_type != "physical" && viscosity_model != "carreau"
                        η_gt_ = model_gt.η[data_range_]
                        av_η = mean(η_gt_)
                        avg_ηList[data_range_] .= av_η
                        printstyled("Average ground truth η in the window: $(av_η), ground truth β: $(β_gt)\n"; color = :green)
                    end

                    if data_type == "physical" && ti == 2
                        @debug "Second time window, applying correction factor to the parameters for better convergence..."
                        θ[1] = stats["η"]*4
                        θ[2] = stats["β"]
                    end

                    @info "Fitting model in time window $(ti)..."
                    stats = fit_model(model, scene, conditions, obs_border_pt_lst_t, θ; outliers=outlier_frames,
                                  method=opt_method, opt_kw...)

                    est_ηpList[data_range_] .= stats["η"]
                    est_βpList[data_range_] .= stats["β"]
                    push!(cost_list, stats["cost_list"])

                    θ[1] = stats["η"]
                    θ[2] = stats["β"]

                    update_model!(model)

                    # k -> k+1 prediction: this window's dynamics from the previous
                    # window's estimate (window 1 uses its own). Same control as the fit.
                    pred_scene.sim_time = time_windows[ti]
                    pred_scene.cParam = scene.cParam

                    if ti > 1
                        @info "Predicting dynamics in time window $(ti) from window $(ti-1)'s estimate..."
                        reset_model!(pred_model)
                        pred_model.η = [est_ηpList[data_range_prev][1]]
                        pred_scene.β = [est_βpList[data_range_prev][1]]
                        pred_μ_list, _, pred_simBorderPts, fields_pred, surface_pts_3D_pred, pos2D_pred, pos3D_pred, _, _, _ = simulate(pred_model, pred_scene, conditions)
                        pred_h = get_height(pred_μ_list, _h_pred)

                        push!(h_pred_list, pred_h)
                        push!(borderPts_pred_list, pred_simBorderPts)
                        push!(fields_list_pred_list, fields_pred)
                        push!(pos2D_pred_list, pos2D_pred)
                        push!(pos3D_pred_list, pos3D_pred)
                        push!(surface_pts_3D_pred_list, surface_pts_3D_pred)
                    end

                    reset_model!(pred_model)
                    pred_model.η = [est_ηpList[data_range_][1]]
                    pred_scene.β = [est_βpList[data_range_][1]]
                    pred_est_μ_list, _, pred_est_simBorderPts, fields_pred_est, surface_pts_3D_pred_est, pos2D_pred_est, pos3D_pred_est, _, _, _ = simulate(pred_model, pred_scene, conditions)
                    h_pred_est = get_height(pred_est_μ_list, _h_pred)
                    if ti == 1
                        push!(h_pred_list, h_pred_est)
                        push!(borderPts_pred_list, pred_est_simBorderPts)
                        push!(fields_list_pred_list, fields_pred_est)
                        push!(pos2D_pred_list, pos2D_pred_est)
                        push!(pos3D_pred_list, pos3D_pred_est)
                        push!(surface_pts_3D_pred_list, surface_pts_3D_pred_est)
                    end
                    update_model!(pred_model)
                    _h_pred = h_pred_est[end] # carry the end height into the next window
                end

                @info "Completed all time windows."

                viscosity_type = "bulk_viscosity"
                est_model, est_scene = def_problem(geom_exp, ne_exp, η_start, _fem..., β_start, F, control, viscosity_type, sim_time_gt, t_steps_gt; viscosity_model=viscosity_model, _mesh_path_kw(exp_params)...)
                est_model.η = est_ηpList
                est_scene.β = est_βpList

                write_csv(datapath(exp_path, "est_η"), est_ηpList)
                write_csv(datapath(exp_path, "est_β"), est_βpList)

                write_csv(datapath(exp_path, "avg_η"), avg_ηList)
                write_csv(datapath(exp_path, "window_data","time_windows"), time_windows)
                write_csv(datapath(exp_path, "window_data","t_windows"), t_windows)
                write_csv(datapath(exp_path, "window_data","data_ranges"), data_ranges_)
                write_csv(datapath(exp_path, "window_data","windows_sizes"), windows)
                write_csv(datapath(exp_path, "window_data","cost_windows"), cost_list)

                est_μ_list, gradList, simBorderPts, fields_est, surface_pts_3D_est, pos2D_est, pos3D_est, _, _, _ = simulate(est_model, est_scene, conditions)
                est_h_list = get_height(est_μ_list, h)

                # Write the per-window predictions accumulated in the loop above.
                write_window_predictions(datapath(exp_path),
                    reduce(vcat, h_pred_list),
                    reduce(vcat, pos3D_pred_list),
                    reduce(vcat, surface_pts_3D_pred_list),
                    reduce(vcat, fields_list_pred_list),
                    reduce(vcat, pos2D_pred_list),
                    reduce(vcat, borderPts_pred_list))

                if data_type != "physical" && viscosity_model != "carreau"
                    gt_μ_list, gradList, borderPts2DList_gt, fields_gt, surface_pts_3D_gt, pos2D_gt, pos3D_gt, _, _, _ = simulate(model_gt, scene_exp, conditions)
                    gt_h_list = get_height(gt_μ_list, h)
                    write_csv(datapath(exp_path, "η_gt"), model_gt.η)
                    write_csv(datapath(exp_path, "β_gt"), β_gt)
                    write_csv(datapath(exp_path, "gt_h"), gt_h_list)

                    write_2d_data(datapath(exp_path, "sim_data","2D_surface_points_gt"), pos2D_gt)
                    write_2d_data(datapath(exp_path, "sim_data","2D_border_points_gt"), borderPts2DList_gt)

                    write_data(datapath(exp_path, "sim_data","3D_points_gt"), pos3D_gt)
                    write_data(datapath(exp_path, "sim_data","motion_fields_gt "), fields_gt)
                    write_data(datapath(exp_path, "sim_data","3D_surface_points_gt"), surface_pts_3D_gt)
                end
            end

            write_csv(datapath(exp_path, "est_h"), est_h_list)

            write_data(datapath(exp_path, "sim_data","3D_points_est"), pos3D_est)
            write_data(datapath(exp_path, "sim_data","motion_fields_est"), fields_est)
            write_data(datapath(exp_path, "sim_data","3D_surface_points_est"), surface_pts_3D_est)
            write_2d_data(datapath(exp_path, "sim_data","2D_surface_points_est"), pos2D_est)
            write_2d_data(datapath(exp_path, "sim_data","2D_border_points_est"), simBorderPts)
        end
        end_time = Dates.now()
        write_time_log(start_time, end_time, exp_params; dest_dir=joinpath(exp_path, "logs"))
        @info "Data writing completed. Results saved to $exp_path"
    end
end

"""
    predict(filepath, filepath_gt)

Out-of-sample forward prediction for bulk-viscosity experiments already fit by
`optimize`/`replot`-style time windowing. For every leaf
directory under `filepath` (excluding `view_1`, walking
`elem_size/simtime/noise/dt/view/window`), re-simulate each time window using
the η/β estimated for that window (and, from the second window onward, also
simulate the *previous* window's estimate forward for comparison), then write
the resulting predicted height, 3D points, motion fields, and 2D
border/surface points back into that experiment's `data` directory.

Only applies when `sim_params["viscosity_type"] == "bulk_viscosity"`; a no-op
otherwise.

# Arguments
- `filepath::String`: root of the result tree to walk and update in place.
- `filepath_gt::String`: ground-truth data root (supplies `sim_params`).

# Returns
None. All outputs are side effects written to disk.
"""
function predict(filepath, filepath_gt)

    sim_params = read_json(datapath(filepath_gt,"sim_params")) 

    r::Float64 = sim_params["r"]
    h::Float64 = sim_params["h"]

    viscosity_type::String = sim_params["viscosity_type"]
    F = Array(float.(sim_params["cParam"]))

    sim_time::Float64 = sim_params["simulation_time"]
    t_steps::Float64 = sim_params["time_steps"]

    camera_matrix::AbstractArray = _read_camera_matrix(sim_params)
    obj_pose::AbstractArray{Float64} = _obj_pose_for(sim_params)
    control::String = sim_params["control_type"]
    edge_radius::Union{Float64, Nothing} = get(sim_params, "edge_radius", nothing)

    viscosity_model::String = ""
    if viscosity_type == "bulk_viscosity" && haskey(sim_params, "model_type") && sim_params["model_type"] == "carreau"
        viscosity_model = sim_params["model_type"]
    end

    ndim::Int = 3
    nDof_p::Int = 1
    nDof_u::Int = ndim

    geometry = get(sim_params, "geometry", :cylinder)
    _mesh_path_kw(d) = haskey(d, "mesh_path") ? (mesh_path=d["mesh_path"],) : (;)
    geom = if geometry === :cube
        @info "Using cuboid geometry with dimensions $r x $r x $h"
        lx = Float64(get(sim_params, "lx", r))
        ly = Float64(get(sim_params, "ly", r))
        lz = Float64(get(sim_params, "lz", h))
        Cuboid(lx, ly, lz, edge_radius)
    elseif geometry === :cylinder
        @info "Using cylinder geometry with radius $r and height $h"
        Cylinder(r, h)
    else
        error("Unsupported geometry type: $geometry")
    end

    gt_h_::Matrix{Float64} = Matrix{Float64}(undef,0,0)

    if viscosity_type == "bulk_viscosity"
        viscosity_type = "constant"

        # Each leaf is a dt_/view_ directory; predictions live one level deeper,
        # under per-window subdirectories of the view.
        for group in collect_experiment_groups(filepath; method=method)
            for leaf in group.leaves
                view_path = leaf.exp_path
                for window_dir in readdir(view_path)
                    win_exp_path = joinpath(view_path, window_dir)

                    @debug "Processing window: $win_exp_path"
                    exp_params = read_json(datapath(win_exp_path,"experiment_parameters"))

                    data_type::String = exp_params["data_type"]
                    ne_exp::Int = exp_params["ne_exp"]
                    z_angle::Float64 = exp_params["z_angle"]
                    element_shape_u = Symbol(exp_params["element_shape_u"])
                    basis_order_u   = Int(exp_params["basis_order_u"])
                    element_shape_p = Symbol(exp_params["element_shape_p"])
                    basis_order_p   = Int(exp_params["basis_order_p"])
                    element_shape_x = Symbol(exp_params["element_shape_x"])
                    basis_order_x   = Int(exp_params["basis_order_x"])

                    est_ηpList = readdlm(datapath(win_exp_path,"est_η.csv"), ',', Float64)
                    est_βpList = readdlm(datapath(win_exp_path,"est_β.csv"), ',', Float64)
                    est_h_list = readdlm(datapath(win_exp_path,"est_h.csv"), ',', Float64)

                    data_ranges_ = get_time_windows(datapath(win_exp_path,"window_data","data_ranges.csv"))
                    time_windows = readdlm(datapath(win_exp_path,"window_data","time_windows.csv"),',',Float64)

                    conditions = Conditions(camera_matrix=camera_matrix, obj_pose=obj_pose, viewing_angles=[z_angle])
                    _fem = (element_shape_u, basis_order_u, nDof_u, element_shape_p, basis_order_p, nDof_p, element_shape_x, basis_order_x)
                    model, scene = def_problem(geom, ne_exp, 1.0, _fem..., 1.0, F, control, viscosity_type, sim_time, t_steps; _mesh_path_kw(exp_params)...)

                    h_pred_vec, pos3D_pred_vec, surface_pts_3D_pred_vec, fields_pred_vec, pos2D_pred_vec, borderPts_pred_vec =
                        run_window_predictions(model, scene, conditions, est_ηpList, est_βpList, data_ranges_, time_windows, F, h)

                    write_window_predictions(datapath(win_exp_path),
                        h_pred_vec, pos3D_pred_vec, surface_pts_3D_pred_vec,
                        fields_pred_vec, pos2D_pred_vec, borderPts_pred_vec)
                end
            end
        end
    end
    return
end

"""
    set_plot_config(data_type, viscosity_type)

Set the module-global plot styling variables (`fs`, `plt_height`, `plt_width`,
margins, and the `y_lims_h_norm`/`y_lims_rel_error` axis limits used
throughout `replot`/`post_analysis_*`) from `PLOT_PRESETS` (defined near the
top of this file), keyed by `data_type ∈ ("physical", "synthetic",
"simulated")` — `"synthetic"` is further split by `viscosity_type` since it
uses a narrower plot for `"bulk_viscosity"`. Falls back to `PLOT_CONFIG`
(with a warning) for an unrecognized `data_type`/`viscosity_type` combo.

# Arguments
- `data_type::String`: `"physical"`, `"synthetic"`, or `"simulated"`.
- `viscosity_type::String`: `"constant"` or `"bulk_viscosity"`; only
  consulted when `data_type == "synthetic"`.

# Returns
None. The effect is entirely via the mutated globals.
"""
function set_plot_config(data_type::String, viscosity_type::String)
    global fs, plt_height, plt_width, plt_lft_margin, plt_right_margin, plt_top_margin
    global y_lims_h_norm, y_lims_rel_error

    preset_key = data_type == "synthetic" ? (data_type, viscosity_type) : data_type
    preset = get(PLOT_PRESETS, preset_key, nothing)
    if preset === nothing
        @warn "Unknown data_type/viscosity_type combo '$data_type'/'$viscosity_type', using PLOT_CONFIG default"
        return
    end

    fs = preset[:font_size]
    plt_height = preset[:plot_height]
    plt_width = preset[:plot_width]
    plt_lft_margin = preset[:left_margin]
    plt_right_margin = preset[:right_margin]
    plt_top_margin = preset[:top_margin]

    y_lims_h_norm = preset[:y_lims_h_norm]
    y_lims_rel_error = preset[:y_lims_rel_error]
end

"""
    _palette(i) -> colour

The `i`-th plot colour, wrapping around `color_palette`. `_exp_color` takes the experiment
directory name instead, so a given experiment keeps one colour across every figure.
"""
_palette(i::Integer) = color_palette[(i - 1) % length(color_palette) + 1]
_exp_color(dir::AbstractString) = _palette(parse(Int, dir))

"""
    _fig(; legend_column=1, margins=:sides, sz=(plt_width, plt_height)) -> plot

An empty figure at this file's standard size.

`margins` names the three geometries actually in use — `:all` (left, right and top),
`:sides`, or `:none` — rather than repeating the margin constants at every call. They are
not interchangeable: which one a figure uses changes its layout.
"""
function _fig(; legend_column::Int=1, margins::Symbol=:sides, sz=(plt_width, plt_height))
    margins === :all && return set_plot(fs; sz=sz, legend_column=legend_column,
                                        left_margin=plt_lft_margin, right_margin=plt_right_margin,
                                        top_margin=plt_top_margin)
    margins === :sides && return set_plot(fs; sz=sz, legend_column=legend_column,
                                          left_margin=plt_lft_margin, right_margin=plt_right_margin)
    return set_plot(fs; sz=sz, legend_column=legend_column)
end

"""
    _label!(p, xlabel, ylabel; xlims=nothing, ylims=nothing) -> p

Label a plot's axes and optionally set their limits. Kept separate from `_fig` because the
margin geometry is a per-figure choice while labelling is uniform.
"""
function _label!(p, xlabel, ylabel; xlims=nothing, ylims=nothing)
    Plots.xlabel!(p, xlabel)
    Plots.ylabel!(p, ylabel)
    isnothing(xlims) || Plots.xlims!(p, xlims...)
    isnothing(ylims) || Plots.ylims!(p, ylims...)
    return p
end

"""
    _align_windowed(series, data_ranges, n) -> Vector{Float64}

Map a per-window prediction series onto the `n`-step time axis.

`pred_h` is not a time series: it is the per-window predictions concatenated, and each window
contributes one sample *more* than its own length — the step that carries into the next
window. Indexing it against time directly therefore drifts by one step per window boundary
(0.7 s by the last window of a standard eight-window fit), which silently corrupts anything
computed from it.

Windows overlap by that one step: window `k` ends on the index window `k+1` starts on, and
the later window wins. Unwritten steps stay `NaN` rather than 0, so a gap cannot be mistaken
for a measurement. Mirrors the `range_`/`range_pred` bookkeeping `replot` uses to plot the
same quantity.
"""
function _align_windowed(series::AbstractVector, data_ranges, n::Int)
    aligned = fill(NaN, n)
    pred_start = 1
    for ti in 1:size(data_ranges, 1)
        w = data_ranges[ti]
        rng      = w[1]:min(w[end] + 1, n)
        rng_pred = ti == 1 ? (w[1]:(w[end] + 1)) : (pred_start:(pred_start + length(w)))
        pred_start = rng_pred[end] + 1
        k = min(length(rng), length(rng_pred), max(0, length(series) - rng_pred[1] + 1))
        k > 0 || continue
        aligned[rng[1:k]] = series[rng_pred[1:k]]
    end
    return aligned
end

"""
    _windowed_series(series, data_ranges, n) -> Vector{Tuple{UnitRange{Int},Vector{Float64}}}

Split a per-window prediction into `(time indices, values)`, one entry per window.

Each window runs to `w[end] + 1` — its own final value, the step that carries into the next
window — so consecutive windows both hold a value at that index. A single time-indexed
vector cannot represent that: [`_align_windowed`](@ref) has to let one of them win, which is
fine for statistics but wrong for drawing, because it makes window `k` end on window `k+1`'s
value and the two appear joined. Trimming the range instead just loses window `k`'s last
point. Keeping the windows apart is the only representation that is right on both counts.
"""
function _windowed_series(series::AbstractVector, data_ranges, n::Int)
    out = Tuple{UnitRange{Int},Vector{Float64}}[]
    pred_start = 1
    for ti in 1:size(data_ranges, 1)
        w = data_ranges[ti]
        rng      = w[1]:min(w[end] + 1, n)
        rng_pred = ti == 1 ? (w[1]:(w[end] + 1)) : (pred_start:(pred_start + length(w)))
        pred_start = rng_pred[end] + 1
        k = min(length(rng), length(rng_pred), max(0, length(series) - rng_pred[1] + 1))
        k > 0 || continue
        push!(out, (rng[1]:rng[k], Float64.(series[rng_pred[1:k]])))
    end
    return out
end

"""
    _stat_fig(groups, tgrid, outpath; xlabel, ylabel, kwargs...) -> plt

A statistical figure: for each group, every run drawn thin and its replicate mean ± 95% CI
on top.

`groups` is an iterable of `(runs, stats, colour, linestyle, label)`. Estimation and
prediction go on one figure — solid red and dashed blue, as the per-experiment figures draw
them — because the two are read against each other, and separate axes make that comparison
by eye impossible.

Runs are drawn at their final width rather than thinned afterwards the way [`_band!`](@ref)
does for shared plots: this figure exists only for the band, so there is no plain version
whose widths need preserving.
"""
function _stat_fig(groups, tgrid, outpath; xlabel, ylabel, xlims=nothing, ylims=nothing,
                   hline=nothing, vlines=nothing, legend_column::Int=3, run_lw=0.6)
    plt = _fig(legend_column=legend_column)
    isnothing(hline) || Plots.hline!(plt, [hline], linestyle=:dash, label=false, color=:black)
    # Window boundaries: a windowed fit restarts at each one, so a step in the curves there
    # is the method, not the material.
    isnothing(vlines) || for t in vlines
        Plots.vline!(plt, [t], color=:gray, linestyle=:dash, label=false)
    end
    _label!(plt, xlabel, ylabel; xlims=xlims, ylims=ylims)

    groups = collect(groups)
    for (gi, (runs, st, col, ls, lbl)) in enumerate(groups)
        (isnothing(st) || isempty(runs)) && continue
        last_group = gi == length(groups)
        band_label = string(lbl, L"\;\mathrm{mean}\pm95\%\;\mathrm{CI}")

        if eltype(runs) <: AbstractVector{<:Tuple}
            # Per-window segments: each window is its own series, and its band is computed
            # from that window's values across runs — not sliced out of a shared vector,
            # which would take the neighbouring window's value at the boundary.
            for (wi, _) in enumerate(runs[1])
                cols = [r[wi][2] for r in runs if wi <= length(r)]
                isempty(cols) && continue
                m = minimum(length, cols)
                M = hcat((c[1:m] for c in cols)...)
                rng = runs[1][wi][1]
                rr = rng[1]:rng[m]
                for c in cols
                    Plots.plot!(plt, tgrid[rr], c[1:m]; lw=run_lw, color=col, linestyle=ls,
                                linealpha=0.45, label=false)
                end
                Plots.plot!(plt, tgrid[rr], vec(mean(M; dims=2));
                            ribbon=replicate_ci(vec(std(M; dims=2)), size(M, 2)),
                            lw=1, color=col, linestyle=ls,
                            label=(wi == 1 ? band_label : false))
            end
            last_group && Plots.savefig(plt, outpath)
        else
            for v in runs
                k = min(length(v), length(tgrid))
                Plots.plot!(plt, tgrid[1:k], v[1:k]; lw=run_lw, color=col, linestyle=ls,
                            linealpha=0.45, label=false)
            end
            _band!(plt, st, tgrid, outpath; color=col, linestyle=ls, thin_runs=false,
                   save=last_group, label=band_label)
        end
    end
    return plt
end

"""
    _method_kwargs(method, kw) -> NamedTuple

Drop the entries of `kw` that `method`'s optimizer does not accept, warning for each.

`opt_kwargs` is one configuration bag shared by every method, so it naturally accumulates
keys that only some of them take — `λ_scale`/`θ_p`/`Γ` belong to `:gn_tikhonov`, `λ` to
`:lm`. Forwarding those to `:gn` is a `MethodError` that kills the run partway through a
long sweep. Warning and dropping keeps a stray key visible without being fatal; the accepted
set is read from the optimizer itself, so it cannot drift out of date.
"""
function _method_kwargs(method::Symbol, kw)
    target = method === :lm           ? smearFEM._fit_model_LM :
             method === :gn_tikhonov  ? smearFEM._fit_model_GN_tikhonov :
                                        smearFEM._fit_model_GN
    accepted = Set(Base.kwarg_decl(first(methods(target))))
    dropped = [k for k in keys(kw) if !(k in accepted)]
    isempty(dropped) || @warn "Ignoring keyword(s) that :$method does not accept" dropped
    return (; (k => v for (k, v) in pairs(kw) if k in accepted)...)
end

"""
    _contour_cost(name) -> ContourCost

Map a configuration symbol to the cost object `fit_model` takes.

Experiments are configured with symbols so they survive the JSON round-trip through
`experiment_parameters.json`, which cannot carry a Julia struct.
"""
function _contour_cost(name::Symbol)
    name === :closest_point && return ClosestPointCost()
    name === :chamfer       && return ChamferCost()
    error("Unknown cost function :$name — use :closest_point or :chamfer")
end

"""
    _run_dir(method, cost) -> String

The result-path segment identifying how a fit was produced.

`:closest_point` keeps the bare method name because every result written before the cost
became configurable used it — so existing trees stay valid and only a new cost gets a new
directory. Any other cost is suffixed, e.g. `gn_chamfer`.
"""
_run_dir(method, cost) =
    Symbol(cost) === :closest_point ? string(method) : string(method, "_", cost)

"""
    _method_dir(base, method) -> String

Post-analysis output directory for one optimizer. `gn` and `gn_tikhonov` results live under
the same experiment, so writing their plots to one directory would silently overwrite one
with the other. `nothing` keeps the flat layout for a tree that has only ever seen one.
"""
_method_dir(base::String, method) = isnothing(method) ? base : joinpath(base, method)

"""
    _skip_incomplete(exp_path) -> Bool

Whether `exp_path` holds no results — an interrupted fit — warning if so, for use as
`_skip_incomplete(p) && continue` in a walk over `ExpLeaf`s.

`optimize` writes `experiment_parameters` before it starts fitting, so that file proves only
that a run *started*; the test is that something else was written after it. Without this a
single dead run aborts the whole pass, and with several optimizers under one experiment it
would take every other method's output down with it.

Both layouts count: a single-window fit writes `<leaf>/data`, a windowed one
`<leaf>/<window>/data`.
"""
function _skip_incomplete(exp_path::String)
    has(d) = isdir(d) && any(f -> !startswith(f, "experiment_parameters"), readdir(d))
    ok = has(datapath(exp_path)) ||
         (isdir(exp_path) && any(sub -> has(datapath(joinpath(exp_path, sub))), readdir(exp_path)))
    ok || @warn "Incomplete run, skipping" exp_path
    return !ok
end

"""
    _band!(plt, st, tgrid, outpath; color, label, lw, save) -> plt

Draw a replicate mean ± 95% CI ribbon over the individual runs already on `plt`, and save it.

Defaults to `def_orange`, the one scheme colour with no series meaning attached — red,
blue and green already stand for estimation, prediction and ground truth, and every
`color_palette` entry is in use by a per-run curve underneath. Reserving orange for the
aggregate keeps the band from reading as one more run, or as a series semantic it does not
carry. Pass `color=` to override, as the MOSD figure does for its two curves.

`thin_runs` drops the curves already on `plt` to `run_lw` first, so the mean reads against
them. This happens here rather than where the runs are drawn because the same plot object is
also saved without a band: thinning at the source would change that figure too. Pass
`thin_runs=false` when layering a second band, or the first band's mean gets thinned as well.

Pass `save=false` to layer several bands on one figure and write it once.
"""
function _band!(plt, st, tgrid, outpath; color=def_orange, lw=1, save::Bool=true,
                thin_runs::Bool=true, run_lw=0.6, linestyle=:solid, segments=nothing,
                label=L"\mathrm{mean}\;\pm\;95\%\;\mathrm{CI}")
    # Only the data curves. `hline!`/`vline!` reference lines and the annotation arrows are
    # 3-point series, and thinning those would fade the very guides the band is read against.
    thin_runs && for ser in plt.series_list
        length(ser[:y]) > 3 && (ser[:linewidth] = run_lw)
    end
    k = min(st.n_steps, length(tgrid))
    ci = replicate_ci(st.pointwise_sd[1:k], st.n_runs)
    if isnothing(segments)
        Plots.plot!(plt, tgrid[1:k], st.pointwise_mean[1:k]; ribbon=ci, lw=lw, color=color,
                    linestyle=linestyle, label=label)
    else
        # One series per window so the mean breaks at each boundary like the runs beneath it.
        for (si, seg) in enumerate(segments)
            r = seg[1]:min(seg[end], k)
            isempty(r) && continue
            Plots.plot!(plt, tgrid[r], st.pointwise_mean[r]; ribbon=ci[r], lw=lw, color=color,
                        linestyle=linestyle, label=(si == 1 ? label : false))
        end
    end
    save && Plots.savefig(plt, outpath)
    return plt
end

"""
    _result_methods(filepath_res, dirs, avoid_dirs) -> Vector{Union{Nothing,String}}

Which optimizers have results under a result tree, read off the `dt_*/<method>/view_*` path
level. Returns `[nothing]` for a pre-2026-08-26 tree with no method level, so it is
post-processed once into the flat layout it already uses.
"""
function _result_methods(filepath_res::String, dirs, avoid_dirs)
    found = Set{String}()
    flat = false
    for dir in dirs
        dir in avoid_dirs && continue
        root = joinpath(filepath_res, dir)
        isdir(root) || continue
        for (parent, subdirs, _) in walkdir(root)
            startswith(basename(parent), "dt_") || continue
            for sd in subdirs
                startswith(sd, "view_") ? (flat = true) : push!(found, sd)
            end
        end
    end
    methods = Union{Nothing,String}[m for m in sort(collect(found))]
    flat && push!(methods, nothing)
    return isempty(methods) ? Union{Nothing,String}[nothing] : methods
end

"""
    replot(filepath, filepath_gt)

Regenerate all comparison plots (and some derived post-processed data) for
every fitted experiment found under `filepath`, using the ground-truth
`sim_params`/height/contour data in `filepath_gt`. Does not re-run any
optimization; it only reads previously written `optimize`/`predict` outputs
and produces plots/derived data from them. Iterates leaves grouped by
`collect_experiment_groups`.

For `viscosity_type == "constant"`: plots the observed-vs-simulated contour,
closest-point distance error, estimated-vs-ground-truth height (and its
error), η/β vs. iteration, and cost vs. iteration. When `noise_level == 0`,
additionally post-processes the stored ``(η, β)`` cost-surface grid — locating
its minimum, computing the local Hessian and its principal (steepest/flattest)
directions, and plotting cost-surface slices along them — and writes the
result to `data/direction_analysis`/`data/slice_data`. When `noise_level != 0`
(repeated noisy-sample runs), instead plots the covariance ellipse and
error-banded height/relative-error/normalized-error series across samples.

For `viscosity_type == "bulk_viscosity"`: for each time-windowed
sub-experiment, plots per-window estimated vs. ground-truth η and β,
height and height error/normalized error, closest-point
distance error, and (via `compare_pt_clouds`) Hausdorff/Chamfer/
closest-point 3D surface distances — each comparing the within-window
estimate against the forward-predicted trajectory from the previous window,
with vertical lines marking window boundaries.

# Arguments
- `filepath::String`: root of the result tree to walk and regenerate plots for.
- `filepath_gt::String`: ground-truth data root (supplies `sim_params`, height,
  and contour data).

# Returns
None. All outputs are side effects (PDFs under each experiment's `plots/`
directory, JSON/CSV under `data/`).
"""
function replot(filepath, filepath_gt; method::Union{Nothing,String}=nothing)
    sim_params = read_json(datapath(filepath_gt,"sim_params"))

    r::Float64 = sim_params["r"]
    h::Float64 = sim_params["h"]

    viscosity_type::String = sim_params["viscosity_type"]
    F = Array(float.(sim_params["cParam"]))

    sim_time::Float64 = sim_params["simulation_time"]
    t_steps::Float64 = sim_params["time_steps"]

    camera_matrix::Matrix{Float64} = _read_camera_matrix(sim_params)
    obj_pose::AbstractArray{Float64} = _obj_pose_for(sim_params)

    control::String = sim_params["control_type"]

    viscosity_model::String = ""
    if viscosity_type == "bulk_viscosity" && haskey(sim_params, "model_type") && sim_params["model_type"] == "carreau"
            viscosity_model = sim_params["model_type"]
    end

    ndim::Int = 3
    nDof_p::Int = 1
    nDof_u::Int = ndim

    gt_h_::Matrix{Float64} = Matrix{Float64}(undef,0,0)

    for group in collect_experiment_groups(filepath; method=method)
        for leaf in group.leaves
            view_folder = leaf.view_folder

            exp_path = leaf.exp_path

            _skip_incomplete(exp_path) && continue
            if viscosity_type == "constant"
                @debug "Comparing experiments in: $exp_path"

                η_gt = float(sim_params["η"][1])
                β_gt = float(sim_params["β"][1])

                exp_params = read_json(datapath(exp_path,"experiment_parameters"))

                noise_level = exp_params["noise_level"]
                sim_time_exp = exp_params["sim_time_exp"]   
                ne_exp = exp_params["ne_exp"]
                data_type = exp_params["data_type"]
                control = exp_params["control"]
                num_exp_points::Int = round(Int, 20/t_steps)

                printstyled(num_exp_points, " experimental points will be used for optimization\n"; color = :blue)

                if data_type != "physical"
                    gt_h_ = readdlm(datapath(filepath_gt,"h.csv"), ',', Float64)
                end
                if data_type  == "synthetic"
                    # Ground-truth contours live under a view_* level, as `_get_borders` reads them.
                    ObsDataList, splinexObs, splineyObs = read_csv(datapath(filepath_gt,"img_data",view_folder,"contour_data"))
                    @info "Data type $data_type Reading synthetic contour data of $(length(ObsDataList)) time steps"
                elseif data_type == "simulated"
                    ObsDataList, splinexObs, splineyObs = read_csv(datapath(filepath_gt,"sim_data",view_folder,"contour_data"))
                    @info "Data type $data_type Reading simulated contour data from $(datapath(filepath_gt,"sim_data",view_folder,"contour_data")) of $(length(ObsDataList)) time steps"
                else
                    error("Unknown data type: $data_type")
                end

                time = collect(Float64, range(start=0, stop=sim_time, step=t_steps))
                gt_h = gt_h_[1:(round(Int, sim_time/ t_steps))+1]

                conditions = Conditions(camera_matrix=camera_matrix, obj_pose=obj_pose)
                if noise_level == 0.0
                    obs_border_pt_lst, sim_border_pt_lst, nSplinex, nSpliney, splinex, spliney = _get_borders(data_type, filepath_gt, exp_path, num_exp_points; view_folder=view_folder)

                    est_η = readdlm(datapath(exp_path,"η.csv"), ',', Float64)
                    est_β = readdlm(datapath(exp_path,"β.csv"), ',', Float64)
                    est_h = readdlm(datapath(exp_path,"est_h.csv"), ',', Float64)
                    cost_iter = readdlm(datapath(exp_path,"cost_iter.csv"), ',', Float64)

                    stats = read_json(datapath(exp_path,"stats"))
                    η = stats["η"]
                    β = stats["β"]

                    d_est, pairs = contour_cost(sim_border_pt_lst, obs_border_pt_lst)

                    if sim_time_exp == 5.0
                        cont_y_min = 350
                    elseif sim_time_exp == 10.0
                        cont_y_min = 400
                    elseif sim_time_exp == 20.0
                        cont_y_min = 420
                    else
                        cont_y_min = 480
                    end

                    contour_plt = default_plot()
                    Plots.plot!(contour_plt, [], label=false, legend=:outerbottom, legend_column=2, aspect_ratio = :equal)
                    Plots.plot!(contour_plt, nSplinex[end], nSpliney[end], label="Ground truth", color=:red)
                    Plots.yflip!(contour_plt, true)
                    Plots.xlims!(contour_plt, 480, 1520)
                    Plots.ylims!(contour_plt, cont_y_min, 1200)
                    _label!(contour_plt, L"x\;\mathrm{[px]}", L"y\;\mathrm{[px]}")
                    Plots.savefig(contour_plt, plotpath(exp_path,"contour_comparison.pdf"))

                    contour_plt_zoom = set_plot(fs, sz=(1200, plt_height), left_margin=plt_lft_margin, right_margin=plt_right_margin, top_margin=plt_top_margin)
                    Plots.plot!(contour_plt_zoom, [], label=false, legend=:outerbottom, legend_column=2, aspect_ratio = :equal)
                    Plots.plot!(contour_plt_zoom, nSplinex[end], nSpliney[end], label=false, color=:red)
                    Plots.xlims!(contour_plt_zoom, 1000, 1520)
                    Plots.ylims!(contour_plt_zoom, cont_y_min, 1200)
                    Plots.xticks!(contour_plt_zoom, 1100:200:1520)
                    Plots.yflip!(contour_plt_zoom, true)
                    _label!(contour_plt_zoom, L"x\;\mathrm{[px]}", L"y\;\mathrm{[px]}")
                    Plots.savefig(contour_plt_zoom, plotpath(exp_path,"contour_comparison_zoom.pdf"))

                    costList = stats["cost_list"]
                    iterList = stats["iterList"]

                    plt_cnt_error = default_plot()
                    Plots.plot!(plt_cnt_error, time[1:(length(d_est)-1)], d_est[2:end], label="Closest point distance error", legend=:outerbottom, legend_column=2)
                    _label!(plt_cnt_error, L"\mathrm{Time\;[s]}", L"\mathrm{Closest\;Point\;Distance\;[px]}"; xlims=(0, end_obs_win), ylims=(0, max(maximum(d_est[2:end])*1.1, 0.5)))
                    Plots.savefig(plt_cnt_error, plotpath(exp_path,"closest_point_distance_error.pdf"))

                    est_h = vec(Float64.(collect(est_h)))
                    gt_h = vec(Float64.(collect(gt_h)))

                    n_time = min(length(time), length(est_h), length(gt_h)) #, (end_obs_win/t_steps + 1))
                    if n_time < length(time) || n_time < length(est_h) || n_time < length(gt_h)
                        @warn "Time and height vectors have mismatched lengths: time=$(length(time)), est_h=$(length(est_h)), gt_h=$(length(gt_h)). Truncating to $n_time samples for plotting."
                    end

                    h_plt = default_plot()
                    Plots.plot!(h_plt, time[1:n_time], est_h[1:n_time], label="Estimation", legend=:outerbottom, legend_column=2)
                    Plots.plot!(h_plt, time[1:n_time], gt_h[1:n_time], label="Ground truth", legend=:outerbottom, legend_column=2)
                    _label!(h_plt, L"\mathrm{Time\;[s]}", L"\mathrm{Height\;[mm]}"; xlims=(0, end_obs_win))
                    Plots.savefig(h_plt, plotpath(exp_path,"h_est.pdf"))

                    error_plt = default_plot()
                    Plots.plot!(error_plt, time[1:n_time], abs.(est_h[1:n_time] .- gt_h[1:n_time]), label="Estimation error", legend=:outerbottom, legend_column=1)
                    _label!(error_plt, L"\mathrm{Time\;[s]}", "Height Error [mm]"; xlims=(0, end_obs_win))
                    Plots.savefig(error_plt, plotpath(exp_path,"h_est_error.pdf"))

                    η_plt = default_plot()
                    Plots.plot!(η_plt, est_η, label=L"\eta_{\mathrm{est}}", marker=1, legend=:outerbottom, legend_column=2)
                    Plots.hline!([η_gt], label=L"\eta_{\mathrm{gt}}", legend=:outerbottom, legend_column=2)
                    _label!(η_plt, L"\mathrm{Iterations}", latexstring("\$\\eta\$ [kPa s]"))
                    Plots.savefig(η_plt, plotpath(exp_path,"η.pdf"))

                    β_plt = default_plot()
                    Plots.plot!(β_plt, est_β, label=L"\beta_{\mathrm{est}}", marker=1, legend=:outerbottom, legend_column=2)
                    Plots.hline!(β_plt, [β_gt], label=L"\beta_{\mathrm{gt}}",legend=:outerbottom, legend_column=2)
                    _label!(β_plt, L"\mathrm{Iterations}", latexstring("\$\\beta\$ [MPa s m\$^{-1}\$]"))
                    Plots.savefig(β_plt, plotpath(exp_path,"β.pdf"))

                    cost_plt = default_plot()
                    Plots.plot!(cost_plt, iterList, costList, label="Cost", marker=1, yscale=:log10, xminorgrid = :false,legend=:outerbottom, legend_column=1)
                    _label!(cost_plt, L"\mathrm{Iterations}", L"\mathrm{Cost\;[px]}")
                    Plots.savefig(cost_plt, plotpath(exp_path,"cost_steps.pdf"))

                    cost_log_plt = default_plot()
                    Plots.plot!(cost_log_plt, iterList, costList, label="Cost", marker=1, yscale=:log10, xscale=:log10, legend=:outerbottom, legend_column=1)
                    _label!(cost_log_plt, L"\mathrm{Iterations}", L"\mathrm{Cost\;[px]}")
                    Plots.savefig(cost_log_plt, plotpath(exp_path,"cost_steps_log.pdf"))

                    # Contour metrics per optimizer iterate, available only when the fit stored
                    # its iterates (`store_border_pts=true`). `opt_data/` also holds the
                    # noise-sweep series, so test for `iter_*` specifically; result directories
                    # written before this existed are skipped rather than failing the replot.
                    opt_root = datapath(exp_path, "opt_data")
                    has_iterates = isdir(opt_root) &&
                                   any(d -> occursin(r"^iter_[0-9]+$", d), readdir(opt_root))
                    if has_iterates
                        try
                            iteration_contour_metrics(data_type, filepath_gt, exp_path;
                                                      view_folder=view_folder)
                        catch err
                            @warn "Iteration contour metrics failed for $exp_path" exception=err
                        end
                    else
                        @debug "No stored optimizer iterates under $opt_root — skipping iteration contour metrics"
                    end

                    # Coerce contour parameters to concrete arrays (be defensive against scalar values)
                    function _to_vector(x)
                        if isa(x, Number)
                            return Float64[x]
                        elseif isa(x, AbstractVector)
                            return Float64.(collect(x))
                        else
                            try
                                return Float64.(collect(x))
                            catch err
                                error("Cannot coerce input to vector: $err")
                            end
                        end
                    end

                    function _to_matrix(x)
                        if isa(x, AbstractMatrix)
                            return Float64.(collect(x))
                        else
                            try
                                rows = collect(x)

                                # Handle a collection/set of row vectors by stacking to a matrix
                                if !isempty(rows)
                                    row_arrays = map(rows) do r
                                        if isa(r, AbstractVector)
                                            Float64.(collect(r))
                                        elseif isa(r, AbstractMatrix) && (size(r,1) == 1 || size(r,2) == 1)
                                            vec(Float64.(collect(r)))
                                        else
                                            nothing
                                        end
                                    end

                                    if all(!isnothing, row_arrays)
                                        row_arrays = Vector{Vector{Float64}}(row_arrays)
                                        ncols = length(row_arrays[1])
                                        if any(length(row) != ncols for row in row_arrays)
                                            error("Cannot coerce input to matrix: inconsistent row lengths")
                                        end
                                        return reduce(vcat, permutedims.(row_arrays))
                                    end
                                end

                                return Float64.(rows)
                            catch err
                                error("Cannot coerce input to matrix: $err")
                            end
                        end
                    end

                    try 
                        contour_plot_params = read_json(datapath(exp_path,"contour_plot_params")) 
                        ηList = _to_vector(contour_plot_params["η_list"])
                        βList = _to_vector(contour_plot_params["β_list"])
                        CostMat = _to_matrix(contour_plot_params["cost_mat"])

                        # --- Post-process: compute Hessian-based principal directions from CostMat ---
                            # Ensure CostMat shape matches (length(ηList), length(βList))
                            nx = length(ηList); ny = length(βList)
                            sz = size(CostMat)
                            if ndims(CostMat) == 1 && length(CostMat) == nx*ny
                                CostMat = reshape(CostMat, nx, ny)
                            elseif ndims(CostMat) == 2
                                if sz[1] == ny && sz[2] == nx
                                    # CostMat is (ny, nx) => transpose to (nx, ny)
                                    CostMat = CostMat'
                                elseif sz[1] == nx && sz[2] == ny
                                    # already in expected orientation
                                else
                                    @warn "CostMat shape $(sz) does not match (nx,ny)=($(nx),$(ny)). Attempting to reshape if possible."
                                    if prod(sz) == nx*ny
                                        CostMat = reshape(vec(CostMat), nx, ny)
                                    end
                                end
                            end
                            # Prepare interpolation grid (rows->β, cols->η)
                            Z_for_interp = CostMat'

                            # locate minimum in the cost matrix (η-index, β-index)
                            minval, minidx = findmin(CostMat)
                            i0 = Int(minidx[1]); j0 = Int(minidx[2])
                            # guard against out-of-range indices
                            ni = length(ηList); nj = length(βList)
                            if i0 < 1 || i0 > ni || j0 < 1 || j0 > nj
                                error("Minimum index from CostMat (i0=$i0, j0=$j0) is out of bounds for ηList/βList lengths (", ni, ",", nj, ")")
                            end
                            η0 = ηList[i0]; β0 = βList[j0]

                            # choose interior index for finite-difference Hessian computation
                            ii = clamp(i0, 2, max(2, ni-1))
                            jj = clamp(j0, 2, max(2, nj-1))
                            if ii != i0 || jj != j0
                                @warn "Minimum on boundary; using nearby interior index (ii,jj)=($(ii),$(jj)) for Hessian"
                            end

                            # grid spacings (assume near-uniform)
                            dη = (ni > 1) ? (ηList[min(2,ni)] - ηList[1]) : 1.0
                            dβ = (nj > 1) ? (βList[min(2,nj)] - βList[1]) : 1.0

                            C = CostMat
                            # centered finite differences for second derivatives
                            d2C_dη2 = (C[ii+1,jj] - 2C[ii,jj] + C[ii-1,jj]) / (dη^2)
                            d2C_dβ2 = (C[ii,jj+1] - 2C[ii,jj] + C[ii,jj-1]) / (dβ^2)
                            d2C_dηdβ = (C[ii+1,jj+1] - C[ii+1,jj-1] - C[ii-1,jj+1] + C[ii-1,jj-1]) / (4.0*dη*dβ)

                            H = [d2C_dη2 d2C_dηdβ; d2C_dηdβ d2C_dβ2]
                            ev = eigen(Symmetric(H))
                            evals = ev.values
                            evecs = ev.vectors

                            # steepest = largest eigenvalue, flattest = smallest
                            idx_max = argmax(evals); idx_min = argmin(evals)
                            v_steep = evecs[:, idx_max] / norm(evecs[:, idx_max])
                            v_flat  = evecs[:, idx_min] / norm(evecs[:, idx_min])

                            # param-space sampling along directions (normalized) centered at minimum
                            span_η = maximum(ηList) - minimum(ηList)
                            span_β = maximum(βList) - minimum(βList)
                            diag_span = sqrt(span_η^2 + span_β^2)
                            tvec = collect(range(-0.6*diag_span, stop=0.6*diag_span, length=401))

                            etas_steep = η0 .+ tvec .* v_steep[1]
                            betas_steep = β0 .+ tvec .* v_steep[2]
                            etas_flat = η0 .+ tvec .* v_flat[1]
                            betas_flat = β0 .+ tvec .* v_flat[2]

                            # keep samples inside the original grid extents
                            minη, maxη = minimum(ηList), maximum(ηList)
                            minβ, maxβ = minimum(βList), maximum(βList)
                            mask_steep = (etas_steep .>= minη) .& (etas_steep .<= maxη) .& (betas_steep .>= minβ) .& (betas_steep .<= maxβ)
                            mask_flat  = (etas_flat  .>= minη) .& (etas_flat  .<= maxη) .& (betas_flat  .>= minβ) .& (betas_flat  .<= maxβ)

                            t_steep = tvec[mask_steep]; t_flat = tvec[mask_flat]
                            etas_steep = etas_steep[mask_steep]; betas_steep = betas_steep[mask_steep]
                            etas_flat  = etas_flat[mask_flat];  betas_flat  = betas_flat[mask_flat]

                            # interpolate cost along these samples
                            zs_steep = interp_z_at(etas_steep, betas_steep, ηList, βList, Z_for_interp)
                            zs_flat  = interp_z_at(etas_flat,  betas_flat,  ηList, βList, Z_for_interp)

                            # Contour with overlaid directions and minimum (use est_η/est_β as iteration path)
                            # ensure estimation paths are 1D Float64 vectors
                            est_η = vec(Float64.(collect(est_η)))
                            est_β = vec(Float64.(collect(est_β)))

                            @info "Post-process grid sizes: len(ηList)=$(length(ηList)), len(βList)=$(length(βList)), CostMat_size=$(size(CostMat))"
                            # Compute adaptive clims from CostMat to represent all values
                            cost_min, cost_max = extrema(CostMat)
                            cost_max = min(cost_max, 1e3) # avoid zero max
                            cost_clims = (0, cost_max)
                            plt_dirs = default_plot()
                            Plots.contour!(plt_dirs, ηList, βList, CostMat', color=:turbo, fill=false, levels=100, legend=:outerbottom, legend_column=3, clims=cost_clims)
                            Plots.plot!(plt_dirs, est_η, est_β, label="Estimation", ms=:4, m=:x, color=:red, legend=:outerbottom, legend_column=3, markerstrokewidth=0.1)
                            Plots.plot!(plt_dirs, etas_steep, betas_steep, label = "Steepest dir", color=:black, legend=:outerbottom, legend_column=3)
                            Plots.plot!(plt_dirs, etas_flat,  betas_flat,  label = "Flattest dir",  lw=1, color=:magenta, legend=:outerbottom, legend_column=3)
                            Plots.scatter!(plt_dirs, [η0], [β0], label="Minimum Cost", ms=:5, m=:star5, color=:black, legend=:outerbottom, legend_column=3, markerstrokewidth=0.1)
                            Plots.scatter!(plt_dirs, [η_gt], [β_gt], label="Ground truth", ms=:5, m=:star5, color=def_red, legend=:outerbottom, legend_column=3, markerstrokewidth=0.1)
                            _label!(plt_dirs, L"\eta\;\mathrm{[kPa\, s]}", L"\beta\;\mathrm{[Pa\, s \, m]}")
                            Plots.savefig(plt_dirs, plotpath(exp_path, "cost_surface_with_directions.pdf"))

                            plt_cont = default_plot()
                            Plots.contour!(plt_cont, ηList, βList, CostMat', color=:turbo, fill=false, levels=100, legend=:outerbottom, legend_column=3, clims=cost_clims)
                            Plots.plot!(plt_cont, est_η, est_β, label="Estimation", ms=:4, m=:x, color=:red, legend=:outerbottom, legend_column=3, markerstrokewidth=0.1)
                            Plots.scatter!(plt_cont, [η0], [β0], label="Minimum Cost", ms=:5, m=:star5, color=:black, legend=:outerbottom, legend_column=3, markerstrokewidth=0.1)
                            Plots.scatter!(plt_cont, [η_gt], [β_gt], label="Ground truth", ms=:5, m=:star5, color=def_red, legend=:outerbottom, legend_column=3, markerstrokewidth=0.1)
                            _label!(plt_cont, L"\eta\;\mathrm{[kPa\, s]}", L"\beta\;\mathrm{[Pa\, s \, m]}")
                            Plots.savefig(plt_cont, plotpath(exp_path, "cost_surface.pdf"))

                            plt_cont = default_plot()
                            Plots.contourf!(plt_cont, ηList, βList, CostMat', color=:turbo, fill=false, levels=100, legend=:outerbottom, legend_column=3, clims=cost_clims)
                            Plots.plot!(plt_cont, est_η, est_β, label="Estimation", ms=:4, m=:x, color=:red, legend=:outerbottom, legend_column=3, markerstrokewidth=0.1)
                            Plots.scatter!(plt_cont, [η0], [β0], label="Minimum Cost", ms=:5, m=:star5, color=:black, legend=:outerbottom, legend_column=3, markerstrokewidth=0.1)
                            Plots.scatter!(plt_cont, [η_gt], [β_gt], label="Ground truth", ms=:5, m=:star5, color=def_red, legend=:outerbottom, legend_column=3, markerstrokewidth=0.1)
                            _label!(plt_cont, L"\eta\;\mathrm{[kPa\, s]}", L"\beta\;\mathrm{[Pa\, s \, m]}")
                            Plots.savefig(plt_cont, plotpath(exp_path, "cost_surface_iter.pdf"))

                            # 2D slices: cost vs distance along the two directions
                            plt_slices = default_plot()
                            if length(t_steep) > 0 && length(zs_steep) == length(t_steep)
                                Plots.plot!(plt_slices, t_steep, zs_steep, label = "Steepest direction", color=:black, legend=:outerbottom, legend_column=2)
                            else
                                @warn "Skipping steep slice plot: empty or mismatched lengths: $(length(t_steep)) vs $(length(zs_steep))"
                            end
                            if length(t_flat) > 0 && length(zs_flat) == length(t_flat)
                                Plots.plot!(plt_slices, t_flat,  zs_flat,  label = "Flattest direction",  lw=1, color=:gray, legend=:outerbottom, legend_column=2)
                            else
                                @warn "Skipping flat slice plot: empty or mismatched lengths: $(length(t_flat)) vs $(length(zs_flat))"
                            end
                            Plots.vline!(plt_slices, [0.0], color=:blue, linestyle=:dash, label="Minimum")
                            _label!(plt_slices, L"\mathrm{Distance\;from\;minimum\;[px]}", L"\mathrm{Cost}")
                            Plots.savefig(plt_slices, plotpath(exp_path, "cost_slices_along_directions.pdf"))

                            Plots.ylims!(plt_slices, 0, 10)
                            Plots.savefig(plt_slices, plotpath(exp_path, "cost_slices_along_directions_zoomed.pdf"))

                            # save analysis metadata
                            dir_info = Dict("η_min"=>η0, "β_min"=>β0, "Hessian"=>H, "evals"=>evals, "v_steep"=>v_steep, "v_flat"=>v_flat)
                            slice_data = Dict(
                                "steep"=>Dict("t"=>t_steep, "etas"=>etas_steep, "betas"=>betas_steep, "zs"=>zs_steep),
                                "flat"=>Dict("t"=>t_flat, "etas"=>etas_flat, "betas"=>betas_flat, "zs"=>zs_flat)
                            )
                            write_json(datapath(exp_path, "direction_analysis"), dir_info)
                            write_json(datapath(exp_path, "slice_data"), slice_data)
                    catch err
                        @warn "Failed to post-process contour parameters and cost surface: $err"
                    end
                else
                    # try
                        η_pred = readdlm(datapath(exp_path,"eta_est.csv"), ',', Float64) # estimated η values per sample [n x n_iter]
                        β_pred = readdlm(datapath(exp_path,"beta_est.csv"), ',', Float64) # estimated β values per sample [n x n_iter]
                        h_pred = readdlm(datapath(exp_path,"h_est.csv"), ',', Float64) # estimated height values per sample [n x sim_time]

                        exp_params = read_json(datapath(exp_path,"experiment_parameters"))

                        sim_time_exp = exp_params["sim_time_exp"]
                        obs_border_pt_lst, nSplinex, nSpliney, pd = add_noise(ObsDataList, nFactor=noise_level)

                        covarience_plt = default_plot()
                        plot_covariance!(covarience_plt, η_pred[:,1], β_pred[:,1])
                        _label!(covarience_plt, latexstring("\$\\eta\$ [kPa s]"), latexstring("\$\\beta\$ [MPa s m\$^{-1}\$]"))
                        Plots.savefig(covarience_plt, plotpath(exp_path,"covariance.pdf"))

                        n_time = min(length(time), size(h_pred, 2), length(gt_h)) #, (end_obs_win/t_steps + 1))
                        if n_time < length(time) || n_time < size(h_pred, 2) || n_time < length(gt_h)
                            @warn "Time and height vectors have mismatched lengths: time=$(length(time)), est_h=$(size(h_pred, 2)), gt_h=$(length(gt_h)). Truncating to $n_time samples for plotting."
                        end

                        h_plot = _fig(margins=:all, legend_column=2)
                        Plots.plot!(h_plot, time[1:n_time], gt_h[1:n_time], label="Ground truth")
                        StatsPlots.errorline!(h_plot, time[1:n_time], h_pred[:,1:n_time], label="Estimation")
                        _label!(h_plot, L"\mathrm{Time\;[s]}", L"\mathrm{Height\;[mm]}"; xlims=(0, end_obs_win))
                        Plots.savefig(h_plot, plotpath(exp_path,"h_est_noisy.pdf"))

                        # Add inset subplot for zoomed region (10-20 seconds)
                        t_min, t_max = 5.0, 5.1
                        mask_zoom = (time .>= t_min) .& (time .<= t_max)
                        if sum(mask_zoom) > 0
                            idx_zoom = findall(mask_zoom)
                            # Build a single combined figure with two subplots and one shared legend
                            # Use a 2-row layout: top row 2 subplots, bottom row small height for a centered legend
                            plt_layout = @layout [a b; c{0.12h}]

                            # Build a single combined figure and draw into subplots 1,2 and 3
                            plt_combined = set_subplot(fs, sz=(3500, 350), layout=plt_layout)

                            # Full time series on left subplot (subplot=1)
                            Plots.plot!(plt_combined, time, gt_h, subplot=1, label="", color=:blue)
                            try
                                StatsPlots.errorline!(plt_combined, time[1:n_time], h_pred[:,1:n_time], subplot=1, label="", color=def_red)
                            catch e
                                @debug "errorline! failed with h_pred as-is, retrying transposed: $e"
                                StatsPlots.errorline!(plt_combined, time[1:n_time], h_pred[:,1:n_time]', subplot=1, label="", color=def_red)
                            end
                            Plots.xlabel!(plt_combined, L"\mathrm{Time\;[s]}", subplot=1)
                            Plots.ylabel!(plt_combined, L"\mathrm{Height\;[mm]}", subplot=1)

                            # Zoomed region on right subplot (subplot=2)
                            Plots.plot!(plt_combined, time[idx_zoom], gt_h[idx_zoom], subplot=2, label="", color=:blue)
                            try
                                StatsPlots.errorline!(plt_combined, time[idx_zoom], h_pred[:, idx_zoom], subplot=2, label="", color=def_red)
                            catch e
                                @debug "errorline! failed with h_pred as-is, retrying without leading colon: $e"
                                StatsPlots.errorline!(plt_combined, time[idx_zoom], h_pred[idx_zoom], subplot=2, label="", color=def_red)
                            end
                            Plots.xlabel!(plt_combined, L"\mathrm{Time\;[s]}", subplot=2)
                            Plots.xlims!(plt_combined, t_min, t_max, subplot=2)

                            # Bottom subplot (subplot=3) — draw a single centered annotation as deterministic legend
                            Plots.plot!(plt_combined, [], [], subplot=3, label="Ground truth", framestyle=:none, legend=:outerbottom, color=:blue, legend_column=2, background_color=:transparent)
                            Plots.plot!(plt_combined, [], [], subplot=3, label="Estimation", framestyle=:none, legend=:outerbottom, color=def_red, background_color=:transparent)
                            Plots.xlims!(plt_combined, 0.0, 1.0, subplot=3)
                            Plots.ylims!(plt_combined, 0.0, 1.0, subplot=3)

                            # Save combined figure
                            Plots.savefig(plt_combined, plotpath(exp_path, "h_est_noisy_zoomed.pdf"))
                        end

                        error = abs.(h_pred[:,1:n_time]' .- gt_h[1:n_time]) ./ gt_h[1:n_time]
                        h_error_plot = default_plot()
                        StatsPlots.errorline!(h_error_plot, time[1:n_time], error', label="Estimation error", legend=:outerbottom, legend_column=1)
                        _label!(h_error_plot, L"\mathrm{Time\;[s]}", L"\mathrm{Relative\;Height\;Error}"; xlims=(0, end_obs_win))
                        Plots.savefig(h_error_plot, plotpath(exp_path,"h_rel_error_noisy.pdf"))

                        h_norm = h_pred[:,1:n_time]' ./ gt_h[1:n_time]
                        h_normalized_plot = default_plot()
                        StatsPlots.errorline!(h_normalized_plot, time[1:n_time], h_norm', label="Estimation error", legend=:outerbottom, legend_column=1)
                        _label!(h_normalized_plot, L"\mathrm{Time\;[s]}", L"h_{\mathrm{est}}/h_{\mathrm{gt}}"; xlims=(0, end_obs_win))
                        Plots.savefig(h_normalized_plot, plotpath(exp_path,"h_normalized_rel_error_noisy.pdf"))
                end
            elseif viscosity_type == "bulk_viscosity"
                local η_gt = Vector{Float64}(undef, 0)
                local β_gt = Vector{Float64}(undef, 0)

                window_dirs = readdir(exp_path)
                for window_dir in window_dirs
                    if window_dir == "Results" || window_dir == "post_analysis_window" || window_dir == "single_window" || window_dir == "post_analysis_noise"
                        @debug "Skipping directory: $window_dir"
                        continue
                    end
                    win_exp_path = joinpath(exp_path, window_dir)

                    @debug "Processing window: $win_exp_path"
                    exp_params = read_json(datapath(win_exp_path,"experiment_parameters"))
                    sim_time_exp = exp_params["sim_time_exp"]
                    data_type = exp_params["data_type"]
                    ne = exp_params["ne_exp"]

                    gt_h_ = readdlm(datapath(filepath_gt,"h.csv"), ',', Float64)
                    if data_type != "physical"
                        noise_level = exp_params["noise_level"]
                    else
                        noise_level = 0.0
                    end

                    est_ηpList = readdlm(datapath(win_exp_path,"est_η.csv"), ',', Float64)
                    est_βpList = readdlm(datapath(win_exp_path,"est_β.csv"), ',', Float64)
                    cost_init = readdlm(datapath(win_exp_path,"window_data","cost_windows.csv"), ',', '\n')[1,1]
                    # cost_init = 1.0

                    est_h_list = readdlm(datapath(win_exp_path,"est_h.csv"), ',', Float64)
                    pred_h_list = readdlm(datapath(win_exp_path,"pred_h.csv"))
                    avg_ηList = similar(est_ηpList)

                    data_ranges_ = get_time_windows(datapath(win_exp_path,"window_data","data_ranges.csv"))
                    t_windows = readdlm(datapath(win_exp_path,"window_data","t_windows.csv"),',',Float64)
                    time_windows = readdlm(datapath(win_exp_path,"window_data","time_windows.csv"),',',Float64)

                    @debug "Time windows: $(time_windows)"
                    obs_time = sum(time_windows)
                    effective_sim_time = sim_time

                    if obs_time != effective_sim_time
                        @warn "Observation time frame $obs_time is less than preset ground truth time frame $effective_sim_time, updating time frame"
                        effective_sim_time = min(obs_time, effective_sim_time)
                        obs_time = effective_sim_time
                    end

                    if obs_time < sim_time_exp
                        @warn "Observation time frame $obs_time is less than experimental simulation time frame $sim_time_exp, switching to observation time frame"
                        sim_time_exp = obs_time
                    end

                    data_point_len = round(Int, obs_time/t_steps)
                    obs_border_pt_lst, sim_border_pt_lst, nSplinex, nSpliney, splinex, spliney = _get_borders(data_type, filepath_gt, win_exp_path, data_point_len+1; view_folder=view_folder)
                    pred_border_pt_lst, _, _ = read_csv(datapath(win_exp_path,"sim_data","view_1","2D_border_points_pred"))

                    if data_type != "physical"
                        η_gt = float.(sim_params["η"])
                        β_gt = float.(sim_params["β"])
                        η_gt = η_gt[1:data_point_len]
                        avg_ηList = readdlm(datapath(win_exp_path,"avg_η.csv"), ',', Float64)
                    end

                    gt_h = gt_h_[1:(data_point_len+1)]
                    est_h_list = est_h_list[1:data_point_len+1, :]

                    t_full_h = collect(range(start=0, stop=effective_sim_time, step=t_steps))

                    # Physical data has no ground-truth mesh, so there is no sim_data/surface_nodes to
                    # compare against; the surface/MOSD comparison only applies to simulated ground truth.
                    if data_type != "physical"
                        pred_surface_pt_lst, _, _ = read_csv(datapath(win_exp_path,"sim_data","3D_surface_points_pred"))
                        gt_surface_pt_lst, _, _ = read_csv(datapath(filepath_gt,"sim_data","surface_nodes"))
                        est_surface_pt_lst, _, _ = read_csv(datapath(win_exp_path,"sim_data","3D_surface_points_est"))

                        # The ground-truth run can span more frames than the observation window,
                        # so trim it to the observed range like gt_h above.
                        gt_surface_pt_lst = gt_surface_pt_lst[1:(data_point_len+1)]

                        mosd_est = get_surface_mosd(est_surface_pt_lst, obj_pose, h)
                        mosd_pred = get_surface_mosd(pred_surface_pt_lst, obj_pose, h)
                        mosd_gt = get_surface_mosd(gt_surface_pt_lst, obj_pose, h)
                        plt_surface_mosd = _fig(margins=:all, legend_column=2)
                        for ti::Int in 1:(size(data_ranges_, 1))
                            range_ = collect(range(start=data_ranges_[ti][1], stop=(data_ranges_[ti][end]+1), step=1))
                            range_pred = if ti === 1
                                collect(range(start=data_ranges_[ti][1], stop=(data_ranges_[ti][end]+1), step=1))
                            else
                                collect(range(start=pred_range_start, stop=(pred_range_start+size(data_ranges_[ti], 1)), step=1))
                            end
                            t = t_windows[ti]
                            Plots.vline!(plt_surface_mosd, [t], color=:gray, linestyle=:dash, label=false)
                            if ti === 1
                                Plots.plot!(plt_surface_mosd, t_full_h[range_], mosd_pred[range_pred]./mosd_gt[range_], label="Prediction", color=def_blue, linestyle=:dash)
                            else
                                Plots.plot!(plt_surface_mosd, t_full_h[range_], mosd_pred[range_pred]./mosd_gt[range_], label=false, color=def_blue, linestyle=:dash)
                            end
                            pred_range_start = range_pred[end] + 1
                        end
                        Plots.plot!(plt_surface_mosd, t_full_h, mosd_est./mosd_gt, label="Estimation", color=def_red)
                        _label!(plt_surface_mosd, L"\mathrm{Time\;[s]}", L"\mathrm{Relative\;MOSD}"; xlims=(0, end_obs_win), ylims=(y_lims_h_norm))
                        Plots.savefig(plt_surface_mosd, plotpath(win_exp_path,"surface_area_qoi.pdf"))

                        dc_surface_est, dh_surface_est, dcp_surface_est = compare_pt_clouds(est_surface_pt_lst, gt_surface_pt_lst)
                        plt_surface_error_dc = _fig(margins=:all, legend_column=2)
                        plt_surface_error_dh = _fig(margins=:all, legend_column=2)
                        plt_surface_error_dcp = _fig(margins=:all, legend_column=2)
                        pred_range_start = 0
                        for ti::Int in 1:(size(data_ranges_, 1))
                            range_ = collect(range(start=data_ranges_[ti][1], stop=(data_ranges_[ti][end]+1), step=1))
                            range_pred = if ti === 1
                                            collect(range(start=data_ranges_[ti][1], stop=(data_ranges_[ti][end]+1), step=1))
                                        else
                                            collect(range(start=pred_range_start, stop=(pred_range_start+size(data_ranges_[ti], 1)), step=1))
                                        end
                            dc_surface_pred, dh_surface_pred, dcp_surface_pred = compare_pt_clouds(pred_surface_pt_lst[range_pred], gt_surface_pt_lst[range_])
                            t = t_windows[ti]
                            Plots.vline!(plt_surface_error_dc, [t], color=:gray, linestyle=:dash, label=false)
                            Plots.vline!(plt_surface_error_dh, [t], color=:gray, linestyle=:dash, label=false)
                            Plots.vline!(plt_surface_error_dcp, [t], color=:gray, linestyle=:dash, label=false)
                            if ti === 1
                                Plots.plot!(plt_surface_error_dc, t_full_h[range_], dc_surface_pred, label="Prediction error", color=def_blue, linestyle=:dash)
                                Plots.plot!(plt_surface_error_dh, t_full_h[range_], dh_surface_pred, label="Prediction error", color=def_blue, linestyle=:dash)
                                Plots.plot!(plt_surface_error_dcp, t_full_h[range_], dcp_surface_pred, label="Prediction error", color=def_blue, linestyle=:dash)
                            else
                                Plots.plot!(plt_surface_error_dc, t_full_h[range_], dc_surface_pred, label=false, color=def_blue, linestyle=:dash)
                                Plots.plot!(plt_surface_error_dh, t_full_h[range_], dh_surface_pred, label=false, color=def_blue, linestyle=:dash)
                                Plots.plot!(plt_surface_error_dcp, t_full_h[range_], dcp_surface_pred, label=false, color=def_blue, linestyle=:dash)
                            end
                            pred_range_start = range_pred[end] + 1
                        end
                        Plots.plot!(plt_surface_error_dc, t_full_h, dc_surface_est, label="Estimation error", color=def_red)
                        _label!(plt_surface_error_dc, L"\mathrm{Time\;[s]}", L"\mathrm{Hausdorff\;Distance\;[mm]}"; xlims=(0, end_obs_win))
                        Plots.savefig(plt_surface_error_dc, plotpath(win_exp_path,"surface_point_distance_error_dh.pdf"))

                        Plots.plot!(plt_surface_error_dh, t_full_h, dh_surface_est, label="Estimation error", color=def_red)
                        _label!(plt_surface_error_dh, L"\mathrm{Time\;[s]}", L"\mathrm{Chamfer\;Distance\;[mm]}"; xlims=(0, end_obs_win))
                        Plots.savefig(plt_surface_error_dh, plotpath(win_exp_path,"surface_point_distance_error_dc.pdf"))

                        Plots.plot!(plt_surface_error_dcp, t_full_h, dcp_surface_est, label="Estimation error", color=def_red)
                        _label!(plt_surface_error_dcp, L"\mathrm{Time\;[s]}", L"\mathrm{Closest\;Point\;Distance\;[mm]}"; xlims=(0, end_obs_win))
                        Plots.savefig(plt_surface_error_dcp, plotpath(win_exp_path,"surface_point_distance_error_dcp.pdf"))
                    end

                    d_est, _ = contour_cost(sim_border_pt_lst, obs_border_pt_lst)

                    plt_cnt_error = _fig(margins=:all, legend_column=2)
                    pred_range_start = 0
                    for ti::Int in 1:(size(data_ranges_, 1))
                        range_ = collect(range(start=data_ranges_[ti][1], stop=(data_ranges_[ti][end]+1), step=1))
                        range_pred = if ti === 1
                                        collect(range(start=data_ranges_[ti][1], stop=(data_ranges_[ti][end]+1), step=1))
                                    else
                                        collect(range(start=pred_range_start, stop=(pred_range_start+size(data_ranges_[ti], 1)), step=1))
                                    end
                        d_pred, _ = contour_cost(pred_border_pt_lst[range_pred], obs_border_pt_lst[range_])
                        t = t_windows[ti]
                        Plots.vline!(plt_cnt_error, [t], color=:gray, linestyle=:dash, label=false)
                        if ti === 1
                            Plots.plot!(plt_cnt_error, t_full_h[range_], d_pred./cost_init, label="Prediction", color=def_blue, linestyle=:dash)
                        else
                            Plots.plot!(plt_cnt_error, t_full_h[range_], d_pred./cost_init, label=false, color=def_blue, linestyle=:dash)
                        end
                        pred_range_start = range_pred[end] + 1
                    end
                    Plots.plot!(plt_cnt_error, t_full_h, d_est./cost_init, label="Estimation", color=def_red)
                    _label!(plt_cnt_error, L"\mathrm{Time\;[s]}", "Relative Cost"; xlims=(0, end_obs_win), ylims=(0, 3))
                    Plots.savefig(plt_cnt_error, plotpath(win_exp_path,"closest_point_distance_error.pdf"))

                    t_full = collect(range(start=t_steps, stop=effective_sim_time, step=t_steps))

                    plt_η = set_plot(fs, sz=(eta_beta_plt_width, eta_beta_plt_height), legend_column=4, left_margin = plt_lft_margin, right_margin = plt_right_margin, top_margin = plt_top_margin)
                    Plots.plot!(plt_η, [], label=false, legend=:outerbottom, left_margin = plt_lft_margin, right_margin = plt_right_margin, top_margin = plt_top_margin)
                    for ti::Int in 1:(size(data_ranges_, 1))
                        t = t_windows[ti]
                        Plots.vline!(plt_η, [t], color=:gray, linestyle=:dash, label=false)
                    end
                    t_prev = 0.1
                    prev_η = 0.0
                    for ti::Int in 1:(size(data_ranges_, 1))
                        t = t_windows[ti]
                        data_range_ = data_ranges_[ti]
                        t_win = collect(range(start=t_prev, stop=t, step=t_steps))
                        if ti == 1
                            Plots.plot!(plt_η, t_win, est_ηpList[data_range_], label=L"\eta_{\mathrm{est}}(t)", color=def_red)
                            Plots.plot!(plt_η, [], label=L"\eta_{\mathrm{pred}}(t)", color=def_blue, linestyle=:dash)
                        else
                            Plots.plot!(plt_η, t_win, est_ηpList[data_range_], color=def_red, label=false)
                            Plots.plot!(plt_η, t_win, prev_η*ones(length(t_win)), label=false, color=def_blue, linestyle=:dash)
                        end
                        prev_η = est_ηpList[data_range_[end]]
                        t_prev = t+t_steps
                    end
                    if data_type != "physical"
                        Plots.plot!(plt_η, t_full, η_gt, label=L"\eta_{\mathrm{gt}}(t)", color=def_green)
                    end
                    _label!(plt_η, L"\mathrm{Time\;[s]}", latexstring("\$\\eta(t)\$ [kPa s]"); xlims=(0, end_obs_win))
                    Plots.savefig(plt_η, plotpath(win_exp_path,"η.pdf"))

                    plt_β = set_plot(fs, sz=(eta_beta_plt_width, eta_beta_plt_height), legend_column=3, left_margin = plt_lft_margin, right_margin = plt_right_margin, top_margin = plt_top_margin)
                    t_prev = 0.1
                    prev_β = 0.0
                    for ti::Int in 1:(size(data_ranges_, 1))
                        t = t_windows[ti]
                        Plots.vline!(plt_β, [t], color=:gray, linestyle=:dash, label=false)
                        data_range_ = data_ranges_[ti]
                        t_win = collect(range(start=t_prev, stop=t, step=t_steps))
                        if ti == 1
                            Plots.plot!(plt_β, t_win, est_βpList[data_range_], label=L"\beta_{\mathrm{est}}(t)", color=def_red)
                            Plots.plot!(plt_β, [], label=L"\beta_{\mathrm{pred}}(t)", color=def_blue, linestyle=:dash)
                        else
                            Plots.plot!(plt_β, t_win, est_βpList[data_range_], color=def_red, label=false)
                            Plots.plot!(plt_β, t_win, prev_β*ones(length(t_win)), label=false, color=def_blue, linestyle=:dash)
                        end
                        prev_β = est_βpList[data_range_[end]]
                        t_prev = t+t_steps
                    end
                    if data_type != "physical"
                        Plots.hline!(plt_β, β_gt, label=L"\beta_{\mathrm{gt}}", color=def_green)
                    end
                    _label!(plt_β, L"\mathrm{Time\;[s]}", latexstring("\$\\beta(t)\$ [MPa s m\$^{-1}\$]"); xlims=(0, end_obs_win))
                    Plots.savefig(plt_β, plotpath(win_exp_path,"β.pdf"))

                    h_plt = _fig(margins=:all, legend_column=3)
                    Plots.plot!(h_plt, t_full_h, est_h_list, label="Estimation", color=def_red)
                    pred_range_start = 0
                    for ti::Int in 1:(size(data_ranges_, 1))
                        range_ = collect(range(start=data_ranges_[ti][1], stop=(data_ranges_[ti][end]+1), step=1))
                        range_pred = if ti === 1
                                        collect(range(start=data_ranges_[ti][1], stop=(data_ranges_[ti][end]+1), step=1))
                                    else
                                        collect(range(start=pred_range_start, stop=(pred_range_start+size(data_ranges_[ti], 1)), step=1))
                                    end
                        t = t_windows[ti]
                        if ti == 1
                            Plots.plot!(h_plt, [], label="Prediction", linestyle=:dash, color=def_blue)
                        else
                            Plots.plot!(h_plt, t_full_h[range_], pred_h_list[range_pred], linestyle=:dash, color=def_blue, label=false)
                        end
                        Plots.vline!(h_plt, [t], color=:gray, linestyle=:dash, label=false)
                        pred_range_start = range_pred[end] + 1
                    end
                    if data_type != "physical"
                        Plots.plot!(h_plt, t_full_h, gt_h, label="Ground truth", color=def_green)
                    else
                        Plots.plot!(h_plt, t_full_h, gt_h, label="Ground truth", color=def_green)
                    end
                    _label!(h_plt, L"\mathrm{Time\;[s]}", L"\mathrm{Height\;[mm]}"; xlims=(0, end_obs_win), ylims=(minimum(vcat(gt_h, est_h_list))*0.8, maximum(vcat(gt_h, est_h_list))*1.2))
                    Plots.savefig(plotpath(win_exp_path,"h.pdf"))

                    h_normalized_plt = _fig(margins=:all, legend_column=2)
                    pred_range_start = 0
                    for ti::Int in 1:(size(data_ranges_, 1))
                        range_ = collect(range(start=data_ranges_[ti][1], stop=(data_ranges_[ti][end]+1), step=1))
                        range_pred = if ti === 1
                                        collect(range(start=data_ranges_[ti][1], stop=(data_ranges_[ti][end]+1), step=1))
                                    else
                                        collect(range(start=pred_range_start, stop=(pred_range_start+size(data_ranges_[ti], 1)), step=1))
                                    end
                        t = t_windows[ti]
                        if ti == 1
                            if data_type != "physical"
                                Plots.plot!(h_normalized_plt, t_full_h[range_], est_h_list[range_]./gt_h[range_], label="Estimation", color=def_red)
                                Plots.plot!(h_normalized_plt, [], label="Prediction", linestyle=:dash, color=def_blue)
                            else
                                Plots.plot!(h_normalized_plt, t_full_h[range_], est_h_list[range_]./gt_h[range_], label="Estimation", color=def_red)
                                Plots.plot!(h_normalized_plt, [], label="Prediction", linestyle=:dash, color=def_blue)
                            end
                        else
                            Plots.plot!(h_normalized_plt, t_full_h[range_], pred_h_list[range_pred]./gt_h[range_], linestyle=:dash, color=def_blue, label=false)
                            Plots.plot!(h_normalized_plt, t_full_h[range_], est_h_list[range_]./gt_h[range_], color=def_red, label=false)
                        end
                        Plots.vline!(h_normalized_plt, [t], color=:gray, linestyle=:dash, label=false)
                        pred_range_start = range_pred[end] + 1 
                    end
                    _label!(h_normalized_plt, L"\mathrm{Time\;[s]}", L"h/h_{\mathrm{gt}}"; xlims=(0, end_obs_win), ylims=(y_lims_h_norm))
                    Plots.savefig(plotpath(win_exp_path,"h_normalized.pdf"))

                    error_plt = _fig(margins=:all, legend_column=2)
                    Plots.plot!(error_plt, t_full_h, abs.(est_h_list-gt_h), label="Estimation", color=def_red)
                    for ti::Int in 1:(size(data_ranges_, 1))
                        range_ = collect(range(start=data_ranges_[ti][1], stop=(data_ranges_[ti][end]+1), step=1))
                        range_pred = if ti === 1
                                        collect(range(start=data_ranges_[ti][1], stop=(data_ranges_[ti][end]+1), step=1))
                                    else
                                        collect(range(start=pred_range_start, stop=(pred_range_start+size(data_ranges_[ti], 1)), step=1))
                                    end
                        t = t_windows[ti]
                        if ti == 1
                            Plots.plot!(error_plt, [], label="Prediction", linestyle=:dash, color=def_blue)
                        else
                            Plots.plot!(error_plt, t_full_h[range_], abs.(pred_h_list[range_pred]-gt_h[range_]), linestyle=:dash, color=def_blue, label=false)
                        end
                        Plots.vline!(error_plt, [t], color=:gray, linestyle=:dash, label=false)
                        pred_range_start = range_pred[end] + 1 
                    end
                    _label!(error_plt, L"\mathrm{Time\;[s]}", L"\mathrm{Height\;Error\;[mm]}"; xlims=(0, end_obs_win))
                    Plots.savefig(plotpath(win_exp_path,"h_est_error.pdf"))

                    rel_error_plt = _fig(margins=:all, legend_column=2)
                    Plots.plot!(rel_error_plt, t_full_h, abs.(est_h_list-gt_h)./gt_h*100, label="Estimation", color=def_red)
                    for ti::Int in 1:(size(data_ranges_, 1))
                        range_ = collect(range(start=data_ranges_[ti][1], stop=(data_ranges_[ti][end]+1), step=1))
                        range_pred = if ti === 1
                                        collect(range(start=data_ranges_[ti][1], stop=(data_ranges_[ti][end]+1), step=1))
                                    else
                                        collect(range(start=pred_range_start, stop=(pred_range_start+size(data_ranges_[ti], 1)), step=1))
                                    end
                        t = t_windows[ti]
                        if ti == 1
                            Plots.plot!(rel_error_plt, [], label="Prediction", linestyle=:dash, color=def_blue)
                        else
                            Plots.plot!(rel_error_plt, t_full_h[range_], abs.(pred_h_list[range_pred]-gt_h[range_])./gt_h[range_]*100, linestyle=:dash, color=def_blue, label=false)
                        end
                        Plots.vline!(rel_error_plt, [t], color=:gray, linestyle=:dash, label=false)
                        pred_range_start = range_pred[end] + 1 
                    end
                    _label!(rel_error_plt, L"\mathrm{Time\;[s]}", latexstring("Relative Error [\$\\%\$]"); xlims=(0, end_obs_win), ylims=(y_lims_rel_error))
                    Plots.savefig(plotpath(win_exp_path,"rel_error.pdf"))

                end
            end
        end
    end
    return
end

"""
    post_analysis_const(filepath_gt_, filepath, avoid_list)

Cross-experiment post-analysis and figure generation for `viscosity_type ==
"constant"` results (the counterpart of `replot`, which handles a
single experiment). Walks the nested `dir/elem_size_folder/simtime_folder/
noise_folder` result tree under `filepath` (top-level entries in `avoid_list`,
or `"post_analysis_global"`, are skipped), reads each leaf experiment's
already-computed `optimize`/`replot` outputs (η/β estimates, cost history,
height, contours), and overlays them onto shared comparison figures at four
nested levels:

- **global** (across all top-level `dir`s): cost convergence, normalized
  η/β, height/normalized-height/relative-error, and contour comparison,
  bucketed by simulation-time case (2/5/10/20/30 s).
- **per-`dir`** (across element sizes): mesh-convergence plots of cost, η, β,
  and height/relative-height-error.
- **per-`elem_size_folder`** (across simulation-time windows): the same set
  of comparison plots, plus cost-surface slice-direction plots.
- **per-`sim_time_folder`** (across noise levels): covariance ellipse and
  noisy height/error distribution plots.

# Arguments
- `filepath_gt_::String`: ground-truth data root each experiment is compared against.
- `filepath::String`: root of the result tree to walk.
- `avoid_list`: top-level `dir` names under `filepath` to skip.

# Returns
None. All outputs are PDFs written under `plot_path_global`/`plot_path_elems`/
`plot_path_sim_time`/`plot_path_noise`.
"""
function post_analysis_const(filepath_gt_::String, filepath::String, avoid_list; method::Union{Nothing,String}=nothing)

    # Validation metrics, normalized per experiment so they are comparable across runs made
    # at different ground-truth parameters.
    norm_est_runs    = Vector{Vector{Float64}}()   # h_est / h_gt
    rel_err_est_runs = Vector{Vector{Float64}}()   # |h_est - h_gt| / h_gt · 100
    stat_time = Float64[]
    dir_list = readdir(filepath)

    # figure for legend

    # plot for convergence per slip case
    plot_conv_2 = _fig(margins=:all, legend_column=3)
    Plots.plot!(plot_conv_2, [],  label=false)
    _label!(plot_conv_2, L"\mathrm{Iterations}", L"\mathrm{Cost\;[px]}")

    plot_conv_5 = _fig(margins=:all, legend_column=3)
    Plots.plot!(plot_conv_5, [],  label=false)
    _label!(plot_conv_5, L"\mathrm{Iterations}", L"\mathrm{Cost\;[px]}")

    plot_conv_10 = _fig(margins=:all, legend_column=3)
    Plots.plot!(plot_conv_10, [],  label=false)
    _label!(plot_conv_10, L"\mathrm{Iterations}", L"\mathrm{Cost\;[px]}")

    plot_conv_20 = _fig(margins=:all, legend_column=3)
    Plots.plot!(plot_conv_20, [],  label=false)
    _label!(plot_conv_20, L"\mathrm{Iterations}", L"\mathrm{Cost\;[px]}")

    plot_conv_30 = _fig(margins=:all, legend_column=3)
    Plots.plot!(plot_conv_30, [],  label=false)
    _label!(plot_conv_30, L"\mathrm{Iterations}", L"\mathrm{Cost\;[px]}")

    plot_conv_log_2 = _fig(margins=:all, legend_column=3)
    Plots.plot!(plot_conv_log_2, [],  label=false)
    _label!(plot_conv_log_2, L"\mathrm{Iterations}", L"\mathrm{Cost\;[px]}")

    plot_conv_log_5 = _fig(margins=:all, legend_column=3)
    Plots.plot!(plot_conv_log_5, [],  label=false)
    _label!(plot_conv_log_5, L"\mathrm{Iterations}", L"\mathrm{Cost\;[px]}")

    plot_conv_log_10 = _fig(margins=:all, legend_column=3)
    Plots.plot!(plot_conv_log_10, [],  label=false)
    _label!(plot_conv_log_10, L"\mathrm{Iterations}", L"\mathrm{Cost\;[px]}")

    plot_conv_log_20 = _fig(margins=:all, legend_column=3)
    Plots.plot!(plot_conv_log_20, [],  label=false)
    _label!(plot_conv_log_20, L"\mathrm{Iterations}", L"\mathrm{Cost\;[px]}")

    plot_conv_log_30 = _fig(margins=:all, legend_column=3)
    Plots.plot!(plot_conv_log_30, [],  label=false)
    _label!(plot_conv_log_30, L"\mathrm{Iterations}", L"\mathrm{Cost\;[px]}")

    η_norm_plot_2 = _fig(margins=:all, legend_column=3)
    Plots.hline!(η_norm_plot_2, [1.0],  linestyle=:dash, label=false, color=:black)
    _label!(η_norm_plot_2, L"\mathrm{Iterations}", L"\eta_{\mathrm{est}}/\eta_{\mathrm{gt}}")

    β_norm_plot_2 = _fig(margins=:all, legend_column=3)
    Plots.hline!(β_norm_plot_2, [1.0],  linestyle=:dash, label=false, color=:black)
    _label!(β_norm_plot_2, L"\mathrm{Iterations}", L"\beta_{\mathrm{est}}/\beta_{\mathrm{gt}}")

    η_norm_plot_5 = _fig(margins=:all, legend_column=3)
    Plots.hline!(η_norm_plot_5, [1.0],  linestyle=:dash, label=false, color=:black)
    _label!(η_norm_plot_5, L"\mathrm{Iterations}", L"\eta_{\mathrm{est}}/\eta_{\mathrm{gt}}")

    β_norm_plot_5 = _fig(margins=:all, legend_column=3)
    Plots.hline!(β_norm_plot_5, [1.0],  linestyle=:dash, label=false, color=:black, yscale=:log10)
    _label!(β_norm_plot_5, L"\mathrm{Iterations}", L"\beta_{\mathrm{est}}/\beta_{\mathrm{gt}}")

    η_norm_plot_20 = _fig(margins=:all, legend_column=3)
    Plots.hline!(η_norm_plot_20, [1.0],  linestyle=:dash, label=false, color=:black)
    _label!(η_norm_plot_20, L"\mathrm{Iterations}", L"\eta_{\mathrm{est}}/\eta_{\mathrm{gt}}")

    β_norm_plot_20 = _fig(margins=:all, legend_column=3)
    Plots.hline!(β_norm_plot_20, [1.0],  linestyle=:dash, label=false, color=:black)
    _label!(β_norm_plot_20, L"\mathrm{Iterations}", L"\beta_{\mathrm{est}}/\beta_{\mathrm{gt}}")

    η_norm_plot_30 = _fig(margins=:all, legend_column=3)
    Plots.hline!(η_norm_plot_30, [1.0],  linestyle=:dash, label=false, color=:black)
    _label!(η_norm_plot_30, L"\mathrm{Iterations}", L"\eta_{\mathrm{est}}/\eta_{\mathrm{gt}}")

    β_norm_plot_30 = _fig(margins=:all, legend_column=3)
    Plots.hline!(β_norm_plot_30, [1.0],  linestyle=:dash, label=false, color=:black)
    _label!(β_norm_plot_30, L"\mathrm{Iterations}", L"\beta_{\mathrm{est}}/\beta_{\mathrm{gt}}")

    h_glob_plot_2 = _fig(margins=:all, legend_column=3)
    Plots.plot([5,end_obs_win],[0.998,0.998], arrow=arrow(:closed, :both), color=:black, label=false)
    Plots.annotate!(h_glob_plot_2, 15, 0.997, ("Prediction",:black, :center,10,"computer modern"))
    _label!(h_glob_plot_2, L"\mathrm{Time\;[s]}", L"h_{\mathrm{est}}\;\mathrm{[mm]}"; xlims=(0, end_obs_win))

    h_glob_plot_5 = _fig(margins=:all, legend_column=3)
    _label!(h_glob_plot_5, L"\mathrm{Time\;[s]}", L"h_{\mathrm{est}}\;\mathrm{[mm]}"; xlims=(0, end_obs_win))

    h_norm_plot_2 = _fig(margins=:all, legend_column=3)
    Plots.hline!(h_norm_plot_2, [1.0],  left_margin=plt_lft_margin, linestyle=:dash, label=false, color=:black)
    Plots.plot!(h_norm_plot_2,[5,end_obs_win],[0.998,0.998], arrow=arrow(:closed, :both), color=:black, label=false)
    Plots.annotate!(h_norm_plot_2, 15, 0.997, ("Prediction",:black, :center,10,"computer modern"))
    _label!(h_norm_plot_2, L"\mathrm{Time\;[s]}", L"h_{\mathrm{est}}/h_{\mathrm{gt}}"; xlims=(0, end_obs_win), ylims=(y_lims_h_norm))

    h_norm_plot_5 = _fig(margins=:all, legend_column=3)
    Plots.hline!(h_norm_plot_5, [1.0],  left_margin=plt_lft_margin, linestyle=:dash, label=false, color=:black)

    Plots.plot!(h_norm_plot_5,[5,end_obs_win],[1.01,1.01], arrow=arrow(:closed, :both), color=:black, label=false)
    Plots.annotate!(h_norm_plot_5, 15, 1.0125, ("Prediction",:black, :center, 8,"computer modern"))
    Plots.vline!(h_norm_plot_5, [5.0], color=:black, linestyle=:dash, label=false)
    _label!(h_norm_plot_5, L"\mathrm{Time\;[s]}", L"h_{\mathrm{est}}/h_{\mathrm{gt}}"; xlims=(0, end_obs_win), ylims=(y_lims_h_norm))

    rel_height_error_glob_plot_2 = _fig(margins=:all, legend_column=3)
    _label!(rel_height_error_glob_plot_2, L"\mathrm{Time\;[s]}", latexstring("Relative Height Error [\$\\%\$]"); xlims=(0, end_obs_win), ylims=(y_lims_rel_error))

    rel_height_error_glob_plot_5 = _fig(margins=:all, legend_column=3)
    Plots.vline!(rel_height_error_glob_plot_5, [5.0], color=:black, linestyle=:dash, label=false)

    Plots.plot!(rel_height_error_glob_plot_5,[5,end_obs_win],[1.01,1.01], arrow=arrow(:closed, :both), color=:black, label=false)
    Plots.annotate!(rel_height_error_glob_plot_5, 15, 1.0175, ("Prediction",:black, :center, 8,"computer modern"))
    _label!(rel_height_error_glob_plot_5, L"\mathrm{Time\;[s]}", latexstring("Relative Height Error [\$\\%\$]"); xlims=(0, end_obs_win), ylims=(y_lims_rel_error))

    rel_height_error_glob_plot_20 = _fig(margins=:all, legend_column=3)
    Plots.vline!(rel_height_error_glob_plot_20, [20.0], color=:black, linestyle=:dash, label=false)
    _label!(rel_height_error_glob_plot_20, L"\mathrm{Time\;[s]}", latexstring("Relative Height Error [\$\\%\$]"); xlims=(0, end_obs_win), ylims=(y_lims_rel_error))

    rel_height_error_glob_plot_30 = _fig(margins=:all, legend_column=3)
    _label!(rel_height_error_glob_plot_30, L"\mathrm{Time\;[s]}", latexstring("Relative Height Error [\$\\%\$]"); xlims=(0, end_obs_win), ylims=(y_lims_rel_error))

    cont_y_lims = [400, 1200] 
    cont_x_lims = [1380, 1430]
    contour_plt = set_plot(fs, sz=(round(Int,plt_width*1.5), plt_height), legend_column=3)
    Plots.yflip!(true)
    Plots.xlims!(contour_plt, 480, 1520)
    _label!(contour_plt, L"x\;\mathrm{[px]}", L"y\;\mathrm{[px]}")

    cnt_plt_width = 179
    contour_plt_zoom = set_plot(fs, sz=(cnt_plt_width, plt_height), legend_column=3)
    Plots.xticks!(contour_plt_zoom, 1100:200:1520)
    Plots.yflip!(true)
    Plots.xlims!(contour_plt_zoom, cont_x_lims[1], cont_x_lims[2])
    _label!(contour_plt_zoom, L"x\;\mathrm{[px]}", L"y\;\mathrm{[px]}")

    contour_plt_zoom_05 = set_plot(fs, sz=(cnt_plt_width, plt_height), legend_column=3)
    Plots.xticks!(contour_plt_zoom_05, 1100:200:1520)
    Plots.yflip!(true)
    Plots.xlims!(contour_plt_zoom_05, cont_x_lims[1], cont_x_lims[2])
    _label!(contour_plt_zoom_05, L"x\;\mathrm{[px]}", L"y\;\mathrm{[px]}")

    contour_plt_zoom_10 = set_plot(fs, sz=(cnt_plt_width, plt_height), legend_column=3)
    Plots.xticks!(contour_plt_zoom_10, 1100:200:1520)
    Plots.yflip!(true)
    Plots.xlims!(contour_plt_zoom_10, cont_x_lims[1], cont_x_lims[2])
    _label!(contour_plt_zoom_10, L"x\;\mathrm{[px]}", L"y\;\mathrm{[px]}")

    contour_plt_zoom_15 = set_plot(fs, sz=(cnt_plt_width, plt_height), legend_column=3)
    Plots.xticks!(contour_plt_zoom_15, 1100:200:1520)
    Plots.yflip!(true)
    Plots.xlims!(contour_plt_zoom_15, cont_x_lims[1], cont_x_lims[2])
    _label!(contour_plt_zoom_15, L"x\;\mathrm{[px]}", L"y\;\mathrm{[px]}")

    contour_plt_zoom_20 = set_plot(fs, sz=(cnt_plt_width, plt_height), legend_column=3)
    Plots.xticks!(contour_plt_zoom_20, 1100:200:1520)
    Plots.yflip!(true)
    Plots.xlims!(contour_plt_zoom_20, cont_x_lims[1], cont_x_lims[2])
    _label!(contour_plt_zoom_20, L"x\;\mathrm{[px]}", L"y\;\mathrm{[px]}")

    cont_plt_legend = set_plot(fs, sz=(cnt_plt_width, plt_height), legend=:outertopright)

    # Create a thin horizontal legend strip for the outer loop (directory/viscosity slip case legend)
    slip_case_legend = set_plot(11, sz=(round(Int, plt_width*1.7), 50), legend=:bottom, legend_column=3, bottom_margin=-35mm, top_margin=2mm, left_margin=-25mm, right_margin=-15mm)
    Plots.plot!(slip_case_legend, [0, 1], [0, 0], label=false, color=:white, linewidth=0)  # Hidden base plot for axes
    Plots.xlims!(slip_case_legend, -0.2, 1.2)
    Plots.ylims!(slip_case_legend, -0.5, 0.5)

    global color_palette

    max_iter = 1
    for dir in dir_list

        if dir in avoid_list || dir == "post_analysis_global"
            @debug "Skipping directory: $dir"
            continue
        end
        filepath_dir = joinpath(filepath, dir)
        filepath_gt = joinpath(filepath_gt_, dir)

        printstyled("Processing directory: $(filepath_dir)\n", color=:green)
        sim_params = read_json(datapath(filepath_gt,"sim_params"))
        η_gt = sim_params["η"]
        β_gt = sim_params["β"]  
        sim_time = sim_params["simulation_time"]    
        t_steps = sim_params["time_steps"]
        gt_h = readdlm(datapath(filepath_gt,"h.csv"), ',', Float64)

        time = collect(Float64, range(start=0, stop=sim_time, step=t_steps))

        @debug "Ground truth η: $(η_gt[1])"

        Plots.plot!(slip_case_legend, [0, 0.2], [0, 0], label=latexstring("\$\\beta_{\\mathrm{gt}}:$(β_gt[1])\\,\\mathrm{MPa\\,s\\,\\mathrm{m^{-1}}}\$"), color=_exp_color(dir), linewidth=2, markerstrokewidth=0)

        groups = collect_experiment_groups(filepath_dir; method=method)
        elem_size_folders = _elem_size_folders(groups)

        elem_conv_plt = _fig(margins=:all, legend_column=3)
        _label!(elem_conv_plt, L"\mathrm{Iterations}", L"\mathrm{Cost\;[px]}")

        elem_conv_plt_log = _fig(margins=:all, legend_column=3)
        _label!(elem_conv_plt_log, L"\mathrm{Iterations}", L"\mathrm{Cost\;[px]}")

        elem_η_plt = _fig(margins=:all, legend_column=3)
        Plots.hline!(elem_η_plt, [η_gt[1]], label="Ground truth η",  left_margin=plt_lft_margin)
        _label!(elem_η_plt, L"\mathrm{Iterations}", L"\eta\;\mathrm{[kPa\, s]}")

        elem_β_plt = _fig(margins=:all, legend_column=3)
        Plots.hline!(elem_β_plt, [β_gt[1]], label="Ground truth β",  left_margin=plt_lft_margin)
        _label!(elem_β_plt, L"\mathrm{Iterations}", L"\beta\;\mathrm{[Pa\, s \, m]}")

        elem_η_norm_plt = _fig(margins=:all, legend_column=3)
        Plots.hline!(elem_η_norm_plt, [1.0], label=false,  left_margin=plt_lft_margin, linestyle=:dash, color=:black)
        _label!(elem_η_norm_plt, L"\mathrm{Iterations}", L"\eta_{\mathrm{est}}/\eta_{\mathrm{gt}}")

        elem_β_norm_plt = _fig(margins=:all, legend_column=3)
        Plots.hline!(elem_β_norm_plt, [1.0], label=false,  left_margin=plt_lft_margin, linestyle=:dash, color=:black)
        _label!(elem_β_norm_plt, L"\mathrm{Iterations}", L"\beta_{\mathrm{est}}/\beta_{\mathrm{gt}}")

        elem_rel_height_error_plt = _fig(margins=:all, legend_column=3)
        Plots.plot!(elem_rel_height_error_plt, [], left_margin=plt_lft_margin, label=false)
        _label!(elem_rel_height_error_plt, L"\mathrm{Time\;[s]}", latexstring("Relative Height Error [\$\\%\$]"); xlims=(0, end_obs_win))

        elem_height_plt = _fig(margins=:all, legend_column=3)
        Plots.plot!(elem_height_plt, [], left_margin=plt_lft_margin, label=false)
        _label!(elem_height_plt, L"\mathrm{Time\;[s]}", L"\mathrm{Height\;[mm]}"; xlims=(0, end_obs_win))

        for elem_size_folder_ in elem_size_folders
            if elem_size_folder_ == "post_analysis" || elem_size_folder_ == "Q2_16"
                continue
            end

            elem_size_folder = joinpath(filepath_dir, elem_size_folder_)
            printstyled("Processing element size folder: $(elem_size_folder)\n", color=:blue)  
            sim_time_folders = _sim_time_folders(groups, elem_size_folder_)

            sim_window_η_plt = _fig(margins=:none)
            Plots.hline!(sim_window_η_plt, [η_gt[1]], label="Ground truth η")
            _label!(sim_window_η_plt, L"\mathrm{Iterations}", L"\eta\;\mathrm{[kPa\, s]}")

            sim_window_β_plt = _fig(margins=:none)
            Plots.hline!(sim_window_β_plt, [β_gt[1]], label="Ground truth β")
            _label!(sim_window_β_plt, L"\mathrm{Iterations}", L"\beta\;\mathrm{[Pa\, s \, m]}")

            sim_window_η_norm_plt = _fig(margins=:none)
            Plots.hline!(sim_window_η_norm_plt, [1.0], label=false,  left_margin=plt_lft_margin, linestyle=:dash, color=:black)
            _label!(sim_window_η_norm_plt, L"\mathrm{Iterations}", L"\eta_{\mathrm{est}}/\eta_{\mathrm{gt}}")

            sim_window_β_norm_plt = _fig(margins=:none)
            Plots.hline!(sim_window_β_norm_plt, [1.0], label=false,  left_margin=plt_lft_margin, linestyle=:dash, color=:black)
            _label!(sim_window_β_norm_plt, L"\mathrm{Iterations}", L"\beta_{\mathrm{est}}/\beta_{\mathrm{gt}}")

            sim_window_rel_height_error_plt = _fig(margins=:none)
            Plots.plot!(sim_window_rel_height_error_plt, [], left_margin=plt_lft_margin, label=false)
            _label!(sim_window_rel_height_error_plt, L"\mathrm{Time\;[s]}", latexstring("Relative Height Error [\$\\%\$]"); xlims=(0, end_obs_win))

            sim_window_height_plt = _fig(margins=:none)
            Plots.plot!(sim_window_height_plt, [], left_margin=plt_lft_margin, label=false)
            _label!(sim_window_height_plt, L"\mathrm{Time\;[s]}", L"\mathrm{Height\;[mm]}"; xlims=(0, end_obs_win))

            plt_slices = set_plot(fs, sz=(350, 750))
            Plots.vline!(plt_slices, [0.0], color=:blue, linestyle=:dash, label="Minimum",  lw=1)
            _label!(plt_slices, L"\mathrm{Distance\;from\;minimum\;[px]}", L"\mathrm{Cost}"; ylims=(0, 50))

            for sim_time_folder_ in reverse(sort(sim_time_folders))

                if sim_time_folder_ == "post_analysis_time" || sim_time_folder_ == "Results" || sim_time_folder_ == "simtime_2.0" || sim_time_folder_ == "simtime_20.0" || sim_time_folder_ == "simtime_10.0" || sim_time_folder_ == "simtime_30.0"
                    continue
                end

                height_error_plt = _fig(margins=:none)
                Plots.plot!(height_error_plt, [], label=false)
                _label!(height_error_plt, L"\mathrm{Time\;[s]}", L"\mathrm{Height\;Error\;[mm]}"; xlims=(0, end_obs_win))

                rel_height_error_plt = _fig(margins=:none)
                Plots.plot!(rel_height_error_plt, [], label=false)
                _label!(rel_height_error_plt, L"\mathrm{Time\;[s]}", latexstring("Relative Height Error [\$\\%\$]"); xlims=(0, end_obs_win))

                noise_cols::Int = 2
                covarience_plt = _fig(margins=:all, legend_column=noise_cols)
                _label!(covarience_plt, L"\eta/\eta_{\mathrm{gt}}", L"\beta/\beta_{\mathrm{gt}}"; xlims=(0.95, 1.05), ylims=(0.0, 1.75))

                height_noise_plt = _fig(margins=:all, legend_column=noise_cols)
                _label!(height_noise_plt, L"\mathrm{Time\;[s]}", L"\mathrm{Height\;[mm]}"; xlims=(0, end_obs_win))

                height_error_noise_plt = _fig(margins=:all, legend_column=noise_cols)
                Plots.vline!(height_error_noise_plt, [5.0], color=:black, linestyle=:dash, label=false)
                _label!(height_error_noise_plt, L"\mathrm{Time\;[s]}", L"\mathrm{Height\;Error\;[mm]}"; xlims=(0, end_obs_win))

                rel_height_error_noise_plt = _fig(margins=:all, legend_column=noise_cols)
                Plots.vline!(rel_height_error_noise_plt, [5.0], color=:black, linestyle=:dash, label=false)
                _label!(rel_height_error_noise_plt, L"\mathrm{Time\;[s]}", latexstring("Relative Height Error [\$\\%\$]"); xlims=(0, end_obs_win))

                normalized_height_noise_plt = _fig(margins=:all, legend_column=noise_cols)
                Plots.vline!(normalized_height_noise_plt, [5.0], color=:black, linestyle=:dash, label=false)
                Plots.hline!(normalized_height_noise_plt, [1.0], linestyle=:dash, color=:black, label=false)
                Plots.ylabel!(normalized_height_noise_plt, L"h_{\mathrm{est}}/h_{\mathrm{gt}}")
                Plots.xlabel!(normalized_height_noise_plt, L"\mathrm{Time\;[s]}")
                Plots.xlims!(normalized_height_noise_plt, 0, end_obs_win)

                η_noise_norm_plt = _fig(margins=:all, legend_column=noise_cols)
                Plots.hline!(η_noise_norm_plt, [1.0], linestyle=:dash, color=:black, label=false)
                _label!(η_noise_norm_plt, L"\mathrm{Time\;[s]}", L"\eta_{\mathrm{est}}/\eta_{\mathrm{gt}}")

                sim_time_folder = joinpath(elem_size_folder, sim_time_folder_)
                noise_folders = _noise_folders(groups, elem_size_folder_, sim_time_folder_)

                printstyled("Processing simulation time folder: $(sim_time_folder)\n", color=:cyan)

                noise_iter = 0
                min_η_norm, max_η_norm, min_β_norm, max_β_norm = 0.0, 0.0, 0.0, 0.0
                for noise_folder_ in reverse(sort(noise_folders))

                    if noise_folder_ == "post_analysis_noise" || noise_folder_ == "Results" 
                        continue
                    end

                    group = _find_group(groups, elem_size_folder_, sim_time_folder_, noise_folder_)
                    group === nothing && continue
                    leaf = _single_leaf(group)
                    leaf === nothing && continue

                    exp_path = leaf.exp_path
                    _skip_incomplete(exp_path) && continue
                    printstyled("Processing experiment folder: $(exp_path)\n", color=:magenta)
                    exp_params = read_json(datapath(exp_path,"experiment_parameters"))

                    noise_level = exp_params["noise_level"]
                    ne = exp_params["ne_exp"]
                    sim_time = exp_params["sim_time_exp"]
                    data_type = exp_params["data_type"]
                    exp_path_n0 = replace(exp_path, "$noise_folder_" => "noise_0.0")
                    num_exp_points::Int = round(Int,sim_time/t_steps)

                    printstyled("Processing for noise level: $(noise_level): $(exp_path)\n", color=:yellow)
                    obs_border_pt_lst, sim_border_pt_lst, gt_Splinex, gt_Spliney, splinex, spliney = _get_borders(data_type, filepath_gt, exp_path_n0, num_exp_points; view_folder=leaf.view_folder)
                    if noise_level == 0.0

                        est_η = readdlm(datapath(exp_path,"η.csv"), ',', Float64)
                        est_β = readdlm(datapath(exp_path,"β.csv"), ',', Float64)
                        est_h = readdlm(datapath(exp_path,"est_h.csv"), ',', Float64)
                        iter = collect(range(start=1, stop=size(est_η,1), step=1))

                        max_iter = max(max_iter, size(iter,1))

                        if length(est_h) != length(gt_h)
                            @warn "Length mismatch between estimated height ($(length(est_h))) and ground truth height ($(length(gt_h))). Adjusting..."
                            min_length = min(length(est_h), length(gt_h))
                            est_h = est_h[1:min_length]
                            gt_h = gt_h[1:min_length]
                        end

                        height_error = abs.(est_h - gt_h)
                        rel_height_error = height_error ./ gt_h * 100.0

                        # normalize by ground-truth                        
                        est_η_norm = est_η./η_gt
                        est_β_norm = est_β./β_gt
                        h_norm = est_h ./ gt_h

                        cost_list = readdlm(datapath(exp_path,"cost_iter.csv"), ',', Float64)

                        Plots.plot!(sim_window_η_plt, iter, est_η, label=string("Window - $(sim_time)s"), marker=1, legend=:outerbottom, legend_column=2)
                        Plots.plot!(sim_window_β_plt, iter, est_β, label=string("Window - $(sim_time)s"), marker=1, legend=:outerbottom, legend_column=2)

                        time_h = copy(time)
                        n_time = min(length(time), length(rel_height_error))
                        if length(rel_height_error) != length(time) 
                            @warn "Length mismatch between height plot x-axis ($(length(rel_height_error))) and ground truth height ($(length(gt_h))). Adjusting..."
                            min_length = min(length(rel_height_error), length(time))
                            time_h = time[1:min_length]
                            rel_height_error = rel_height_error[1:min_length]
                            est_h = est_h[1:min_length]
                            gt_h = gt_h[1:min_length]
                        end
                        Plots.plot!(sim_window_rel_height_error_plt, time_h, rel_height_error, label=string("Window - $(sim_time)s"," - ne: ",ne), legend=:outerbottom, legend_column=2)
                        Plots.plot!(sim_window_height_plt, time_h, est_h, label=string("Window - $(sim_time)s"), legend=:outerbottom, legend_column=2)
                        Plots.plot!(sim_window_height_plt, time_h, gt_h, label="Ground truth", legend=:outerbottom, legend_column=2)

                        Plots.plot!(sim_window_η_norm_plt, iter, est_η_norm, label=string("Window - $(sim_time)s"), marker=1, legend=:outerbottom, legend_column=2)
                        Plots.plot!(sim_window_β_norm_plt, iter, est_β_norm, label=string("Window - $(sim_time)s"), marker=1, legend=:outerbottom, legend_column=2)

                        try
                            slice_data = read_json(datapath(exp_path,"slice_data"))
                            t_steep = Float64.(collect(slice_data["steep"]["t"]))
                            zs_steep = Float64.(collect(slice_data["steep"]["zs"]))
                            t_flat = Float64.(collect(slice_data["flat"]["t"]))
                            zs_flat = Float64.(collect(slice_data["flat"]["zs"]))

                            # 2D slices: cost vs distance along the two directions
                            if length(t_steep) > 0 && length(zs_steep) == length(t_steep)
                                Plots.plot!(plt_slices, t_steep, zs_steep, legend = :outerright, label = "Steepest direction; Window = $(sim_time)s")
                            else
                                @warn "Skipping steep slice plot: empty or mismatched lengths: $(length(t_steep)) vs $(length(zs_steep))"
                            end

                            if length(t_flat) > 0 && length(zs_flat) == length(t_flat)
                                Plots.plot!(plt_slices, t_flat,  zs_flat, legend = :outerright, label = "Flattest direction; Window = $(sim_time)s", legendfontsize=20)
                            else
                                @warn "Skipping flat slice plot: empty or mismatched lengths: $(length(t_flat)) vs $(length(zs_flat))"
                            end
                        catch e
                            @warn "Failed to read or plot slice data: $e"
                        end

                        if sim_time == 2.0
                            if ne == 6
                                Plots.plot!(η_norm_plot_2, iter, est_η_norm, label=latexstring("\$\\beta_{\\mathrm{gt}}:$(β_gt[1])\\,\\mathrm{MPa\\,s\\,\\mathrm{m^{-1}}}\$"), marker=1, color=_exp_color(dir))
                                Plots.xticks!(η_norm_plot_2, 0:2:(max_iter+1))

                                Plots.plot!(β_norm_plot_2, iter, est_β_norm, label=latexstring("\$\\beta_{\\mathrm{gt}}:$(β_gt[1])\\,\\mathrm{MPa\\,s\\,\\mathrm{m^{-1}}}\$"), marker=1, color=_exp_color(dir))
                                Plots.xticks!(β_norm_plot_2, 0:2:(max_iter+1))

                                Plots.plot!(plot_conv_2, iter, cost_list, label=latexstring("\$\\beta_{\\mathrm{gt}}:$(β_gt[1])\\,\\mathrm{MPa\\,s\\,\\mathrm{m^{-1}}}\$"), marker=1, color=_exp_color(dir))
                                Plots.xticks!(plot_conv_2, 0:2:(max_iter+1))

                                Plots.plot!(plot_conv_log_2, iter, cost_list, label=latexstring("\$\\beta_{\\mathrm{gt}}:$(β_gt[1])\\,\\mathrm{MPa\\,s\\,\\mathrm{m^{-1}}}\$"), marker=1, color=_exp_color(dir), yscale=:log10)
                                Plots.xticks!(plot_conv_log_2, 1:10:max_iter)

                                Plots.plot!(rel_height_error_glob_plot_2, time_h, rel_height_error, label=latexstring("\$\\beta_{\\mathrm{gt}}:$(β_gt[1])\\,\\mathrm{MPa\\,s\\,\\mathrm{m^{-1}}}\$"), color=_exp_color(dir))
                                Plots.plot!(h_norm_plot_2, h_norm, label=latexstring("\$\\beta_{\\mathrm{gt}}:$(β_gt[1])\\,\\mathrm{MPa\\,s\\,\\mathrm{m^{-1}}}\$"), color=_exp_color(dir))
                            end
                        elseif sim_time == 5.0
                            if ne == 6
                                Plots.plot!(η_norm_plot_5, iter, est_η_norm, label=latexstring("\$\\beta_{\\mathrm{gt}}:$(β_gt[1])\\,\\mathrm{MPa\\,s\\,\\mathrm{m^{-1}}}\$"), marker=1, color=_exp_color(dir))
                                Plots.xticks!(η_norm_plot_5, 0:2:(max_iter+1))

                                Plots.plot!(β_norm_plot_5, iter, est_β_norm, label=latexstring("\$\\beta_{\\mathrm{gt}}:$(β_gt[1])\\,\\mathrm{MPa\\,s\\,\\mathrm{m^{-1}}}\$"), marker=1, color=_exp_color(dir))
                                Plots.xticks!(β_norm_plot_5, 0:2:(max_iter+1))

                                Plots.plot!(plot_conv_5, iter, cost_list, label=latexstring("\$\\beta_{\\mathrm{gt}}:$(β_gt[1])\\,\\mathrm{MPa\\,s\\,\\mathrm{m^{-1}}}\$"), marker=1, color=_exp_color(dir))
                                Plots.xticks!(plot_conv_5, 0:2:(max_iter+1))

                                Plots.plot!(plot_conv_log_5, iter, cost_list, label=latexstring("\$\\beta_{\\mathrm{gt}}:$(β_gt[1])\\,\\mathrm{MPa\\,s\\,\\mathrm{m^{-1}}}\$"), marker=1, color=_exp_color(dir), yscale=:log10)
                                Plots.xticks!(plot_conv_log_5, 1:2:max_iter)

                                Plots.plot!(rel_height_error_glob_plot_5, time_h[1:51], rel_height_error[1:51], label=latexstring("\$\\beta_{\\mathrm{gt}}:$(β_gt[1])\\,\\mathrm{MPa\\,s\\,\\mathrm{m^{-1}}}\$"), color=_exp_color(dir))
                                Plots.plot!(rel_height_error_glob_plot_5, time_h, rel_height_error, label=false, color=_exp_color(dir), linestyle=:dash)

                                Plots.plot!(h_norm_plot_5, time_h[1:51], h_norm[1:51], label=latexstring("\$\\beta_{\\mathrm{gt}}:$(β_gt[1])\\,\\mathrm{MPa\\,s\\,\\mathrm{m^{-1}}}\$"), color=_exp_color(dir))
                                Plots.plot!(h_norm_plot_5, time_h, h_norm[1:n_time], label=false, color=_exp_color(dir), linestyle=:dash)

                                push!(norm_est_runs, vec(h_norm[1:n_time]))
                                push!(rel_err_est_runs, vec(rel_height_error[1:min(n_time, length(rel_height_error))]))
                                isempty(stat_time) && (stat_time = collect(time_h))

                                Plots.plot!(h_glob_plot_5, time_h[1:51], h_norm[1:51], label=latexstring("\$\\beta_{\\mathrm{gt}}:$(β_gt[1])\\,\\mathrm{MPa\\,s\\,\\mathrm{m^{-1}}}\$"), color=_exp_color(dir))
                                Plots.plot!(h_glob_plot_5, time_h, h_norm[1:n_time], label=false, color=_exp_color(dir), linestyle=:dash)
                                Plots.plot!(h_glob_plot_5, time_h, gt_h[1:n_time], label=false, style=:dash, color=_exp_color(dir))

                                Plots.plot!(contour_plt, gt_Splinex[end], gt_Spliney[end], label=latexstring("\$\\beta_{\\mathrm{gt}}:$(β_gt[1])\\,\\mathrm{MPa\\,s\\,\\mathrm{m^{-1}}}\$"), color=_exp_color(dir))
                                Plots.plot!(contour_plt_zoom, gt_Splinex[end], gt_Spliney[end], label=latexstring("\$\\beta_{\\mathrm{gt}}:$(β_gt[1])\\,\\mathrm{MPa\\,s\\,\\mathrm{m^{-1}}}\$"), color=_exp_color(dir))
                            end
                            @debug "Size mismatch: $(size(height_error, 1)) vs $(size(gt_h)) vs $(size(time_h)) vs $(length(rel_height_error))"
                            Plots.plot!(height_error_plt, time_h, height_error[1:n_time], label=latexstring("\$$(ne)\\times$(ne)\\times$(ne)\$"), marker=1, legend=:outerbottom)

                            Plots.plot!(rel_height_error_plt, time_h, rel_height_error[1:n_time], label=latexstring("\$$(ne)\\times$(ne)\\times$(ne)\$"), marker=1, legend=:outerbottom)

                            Plots.plot!(elem_conv_plt, iter, cost_list, label=latexstring("\$$(ne)\\times$(ne)\\times$(ne)\$"), marker=1, legend=:outerbottom)
                            Plots.xticks!(elem_conv_plt, 0:2:(max_iter+1))

                            Plots.plot!(elem_conv_plt_log, iter, cost_list, label=latexstring("\$$(ne)\\times$(ne)\\times$(ne)\$"), marker=1, legend=:outerbottom, yscale=:log10)
                            Plots.xticks!(elem_conv_plt_log, 0:2:(max_iter+1))

                            Plots.plot!(elem_η_plt, iter, est_η, label=latexstring("\$$(ne)\\times$(ne)\\times$(ne)\$"), marker=1, legend=:outerbottom)
                            Plots.xticks!(elem_η_plt, 0:2:(max_iter+1))

                            Plots.plot!(elem_β_plt, iter, est_β, label=latexstring("\$$(ne)\\times$(ne)\\times$(ne)\$"), marker=1, legend=:outerbottom)
                            Plots.xticks!(elem_β_plt, 0:2:(max_iter+1))

                            Plots.plot!(elem_rel_height_error_plt, time_h, rel_height_error[1:n_time], label=latexstring("\$$(ne)\\times$(ne)\\times$(ne)\$"), legend=:outerbottom)
                            Plots.plot!(elem_height_plt, time_h, est_h[1:n_time], label=latexstring("\$$(ne)\\times$(ne)\\times$(ne)\$"), legend=:outerbottom)
                            Plots.plot!(elem_height_plt, time_h, gt_h[1:n_time], label="Ground truth", legend=:outerbottom)

                            Plots.plot!(elem_η_norm_plt, iter, est_η_norm, label=latexstring("\$$(ne)\\times$(ne)\\times$(ne)\$"), marker=1, legend=:outerbottom)
                            Plots.xticks!(elem_η_norm_plt, 0:2:(max_iter+1))
                            Plots.plot!(elem_β_norm_plt, iter, est_β_norm, label=latexstring("\$$(ne)\\times$(ne)\\times$(ne)\$"), marker=1, legend=:outerbottom)
                            Plots.xticks!(elem_β_norm_plt, 0:2:(max_iter+1))
                        end
                    else
                        η_pred = readdlm(datapath(exp_path,"eta_est.csv"), ',', Float64) # estimated η values per sample
                        β_pred = readdlm(datapath(exp_path,"beta_est.csv"), ',', Float64) # estimated β values per sample
                        h_pred = readdlm(datapath(exp_path,"h_est.csv"), ',', Float64) # estimated height values per sample

                        n_obs_border_pt_lst, n_gt_Splinex, n_gt_Spliney = add_noise(obs_border_pt_lst, nFactor=noise_level)

                        η_norm = η_pred[:,1] ./ η_gt
                        β_norm = β_pred[:,1] ./ β_gt

                        if minimum(η_norm) < min_η_norm || min_η_norm == 0.0
                            min_η_norm = minimum(η_norm)
                        end
                        if maximum(η_norm) > max_η_norm || max_η_norm == 0.0
                            max_η_norm = maximum(η_norm)
                        end
                        if minimum(β_norm) < min_β_norm || min_β_norm == 0.0
                            min_β_norm = minimum(β_norm)
                        end
                        if maximum(β_norm) > max_β_norm || max_β_norm == 0.0
                            max_β_norm = maximum(β_norm)
                        end

                        plot_covariance!(covarience_plt, η_pred[:,1]./η_gt, β_pred[:,1]./β_gt, label=string(L"\sigma:\;",(round(noise_level,digits=2))," px  "), legend_column=noise_cols, color_ellipse=_palette(noise_iter + 1))

                        if noise_level == 0.5
                            Plots.plot!(contour_plt_zoom_05, n_gt_Splinex[end], n_gt_Spliney[end], label=latexstring("\$\\beta_{\\mathrm{gt}}:$(β_gt[1])\\,\\mathrm{MPa\\,s\\,\\mathrm{m^{-1}}}\$"), color=_exp_color(dir))
                            Plots.plot!(cont_plt_legend, [], [], label=latexstring("\$\\beta_{\\mathrm{gt}}:$(β_gt[1])\\,\\mathrm{MPa\\,s\\,\\mathrm{m^{-1}}}\$"), color=_exp_color(dir))
                        elseif noise_level == 1.0
                            Plots.plot!(contour_plt_zoom_10, n_gt_Splinex[end], n_gt_Spliney[end], label=latexstring("\$\\beta_{\\mathrm{gt}}:$(β_gt[1])\\,\\mathrm{MPa\\,s\\,\\mathrm{m^{-1}}}\$"), color=_exp_color(dir))
                        elseif noise_level == 1.5
                            Plots.plot!(contour_plt_zoom_15, n_gt_Splinex[end], n_gt_Spliney[end], label=latexstring("\$\\beta_{\\mathrm{gt}}:$(β_gt[1])\\,\\mathrm{MPa\\,s\\,\\mathrm{m^{-1}}}\$"), color=_exp_color(dir))
                        elseif noise_level == 2.0
                            Plots.plot!(contour_plt_zoom_20, n_gt_Splinex[end], n_gt_Spliney[end], label=latexstring("\$\\beta_{\\mathrm{gt}}:$(β_gt[1])\\,\\mathrm{MPa\\,s\\,\\mathrm{m^{-1}}}\$"), color=_exp_color(dir))
                        end

                        n_time = min(length(time), size(h_pred, 2), length(gt_h))
                        if n_time < length(time) || n_time < size(h_pred, 2) || n_time < length(gt_h)
                            @warn "Time and height vectors have mismatched lengths: time=$(length(time)), est_h=$(size(h_pred, 2)), gt_h=$(length(gt_h)). Truncating to $n_time samples for plotting."
                        end
                        height_error = abs.(h_pred[:,1:n_time]' .- gt_h[1:n_time])
                        normalized_height = h_pred[:,1:n_time]' ./ gt_h[1:n_time]
                        rel_height_error = abs.(h_pred[:,1:n_time]' .- gt_h[1:n_time]) ./ gt_h[1:n_time] * 100.0 # in percentage

                        time_h = copy(time)
                        if size(rel_height_error,1) != length(time) 
                            @warn "Length mismatch between height plot x-axis ($(size(rel_height_error,1))) and ground truth height ($(length(gt_h))). Adjusting..."
                            min_length = min(size(rel_height_error,1), length(time))
                            time_h = time[1:min_length] 
                            h_pred = h_pred[:,1:min_length]
                            gt_h = gt_h[1:min_length]
                        end
                        Plots.plot!(height_noise_plt, time_h[1:n_time], gt_h[1:n_time], style=:dash, label=false)
                        StatsPlots.errorline!(height_noise_plt, time_h[1:n_time], h_pred[:,1:n_time]', label=false, groupcolor=_palette(noise_iter + 1), linestyle=:dash)
                        StatsPlots.errorline!(height_noise_plt, time_h[1:51], h_pred'[1:51,:], label=string(L"\sigma:\;",(round(noise_level,digits=2))," px  "), groupcolor=_palette(noise_iter + 1))

                        StatsPlots.errorline!(height_error_noise_plt, time_h, height_error', label=false, groupcolor=_palette(noise_iter + 1), linestyle=:dash)
                        StatsPlots.errorline!(height_error_noise_plt, time_h[1:51], height_error[1:51,:]', label=string(L"\sigma:\;",(round(noise_level,digits=2))," px  "), color=_palette(noise_iter + 1))

                        StatsPlots.errorline!(rel_height_error_noise_plt, time_h, rel_height_error', label=false, groupcolor=_palette(noise_iter + 1), linestyle=:dash)
                        StatsPlots.errorline!(rel_height_error_noise_plt, time_h[1:51], rel_height_error[1:51,:]', label=string(L"\sigma:\;",(round(noise_level,digits=2))," px  "), groupcolor=_palette(noise_iter + 1))

                        StatsPlots.errorline!(normalized_height_noise_plt, time_h, normalized_height', label=false, groupcolor=_palette(noise_iter + 1), linestyle=:dash)
                        StatsPlots.errorline!(normalized_height_noise_plt, time_h[1:51], normalized_height[1:51,:]', label=string(L"\sigma:\;",(round(noise_level,digits=2))," px  "), groupcolor=_palette(noise_iter + 1))

                        noise_iter += 1
                    end
                end
                xlims = [min_η_norm, max_η_norm]
                ylims = [min_β_norm, max_β_norm]
                Plots.hline!(covarience_plt, [1], linestyle=:dash, label=false, color=:black)
                Plots.vline!(covarience_plt, [1], linestyle=:dash, label=false, color=:black)
                Plots.plot!(covarience_plt, bottom_margin=-20mm)

                plot_path_noise = _method_dir(joinpath(sim_time_folder,"post_analysis_noise","plots"), method)
                set_file(plot_path_noise)
                Plots.savefig(covarience_plt, joinpath(plot_path_noise,"covariance_$dir.pdf"))
                Plots.savefig(height_noise_plt, joinpath(plot_path_noise,"height_$dir.pdf"))
                Plots.savefig(height_error_noise_plt, joinpath(plot_path_noise,"height_error_$dir.pdf"))
                Plots.savefig(rel_height_error_noise_plt, joinpath(plot_path_noise,"relative_height_error_$dir.pdf"))
                Plots.savefig(normalized_height_noise_plt, joinpath(plot_path_noise,"normalized_height_$dir.pdf"))
                Plots.savefig(η_noise_norm_plt, joinpath(plot_path_noise,"eta_normalized_$dir.pdf"))

                @info "Saved plots to $plot_path_noise"
            end
            plot_path_sim_time = _method_dir(joinpath(elem_size_folder,"post_analysis_time","plots"), method)
            set_file(plot_path_sim_time)

            Plots.savefig(sim_window_η_plt, joinpath(plot_path_sim_time,"η_$dir.pdf"))
            Plots.savefig(sim_window_β_plt, joinpath(plot_path_sim_time,"β_$dir.pdf"))
            Plots.savefig(sim_window_rel_height_error_plt, joinpath(plot_path_sim_time,"relative_height_error_$dir.pdf"))
            Plots.savefig(sim_window_height_plt, joinpath(plot_path_sim_time,"height_$dir.pdf"))
            Plots.savefig(sim_window_η_norm_plt, joinpath(plot_path_sim_time,"η_normalized_$dir.pdf"))
            Plots.savefig(sim_window_β_norm_plt, joinpath(plot_path_sim_time,"β_normalized_$dir.pdf"))
            Plots.savefig(plt_slices, joinpath(plot_path_sim_time,"cost_slices_along_directions_$dir.pdf"))
            @info "Saved plots to $plot_path_sim_time"
        end

        plot_path_elems = _method_dir(joinpath(filepath_dir,"post_analysis","plots"), method)
        set_file(plot_path_elems)
        Plots.savefig(elem_conv_plt, joinpath(plot_path_elems,"conv_$dir.pdf"))
        Plots.savefig(elem_conv_plt_log, joinpath(plot_path_elems,"conv_log_$dir.pdf"))
        Plots.savefig(elem_η_plt, joinpath(plot_path_elems,"η_$dir.pdf"))
        Plots.savefig(elem_β_plt, joinpath(plot_path_elems,"β_$dir.pdf"))
        Plots.savefig(elem_η_norm_plt, joinpath(plot_path_elems,"η_normalized_$dir.pdf"))
        Plots.savefig(elem_β_norm_plt, joinpath(plot_path_elems,"β_normalized_$dir.pdf"))
        Plots.savefig(elem_rel_height_error_plt, joinpath(plot_path_elems,"relative_height_error_$dir.pdf"))
        Plots.savefig(elem_height_plt, joinpath(plot_path_elems,"height_$dir.pdf"))
        @info "Saved plots to $plot_path_elems"
    end

    plot_path_global = _method_dir(joinpath(filepath,"post_analysis_global","plots"), method)
    set_file(plot_path_global)

    @info "Saving plots to $plot_path_global"

    Plots.savefig(plot_conv_5, joinpath(plot_path_global,"conv_5.pdf"))
    Plots.savefig(plot_conv_log_5, joinpath(plot_path_global,"conv_log_5.pdf"))
    Plots.savefig(η_norm_plot_5, joinpath(plot_path_global,"η_normalized_5.pdf"))
    Plots.savefig(β_norm_plot_5, joinpath(plot_path_global,"β_normalized_5.pdf"))
    Plots.savefig(rel_height_error_glob_plot_5, joinpath(plot_path_global,"relative_height_error_5.pdf"))
    Plots.savefig(h_norm_plot_5, joinpath(plot_path_global,"h_normalized_5.pdf"))
    Plots.savefig(h_glob_plot_5, joinpath(plot_path_global,"height_comparison_5.pdf"))

    height_stats = filter(!isnothing, [
        normalized_replicate_stats(norm_est_runs,    "h_est/h_gt";         reference=1.0),
        normalized_replicate_stats(rel_err_est_runs, "rel. error est [%]"; reference=0.0, signed=false),
    ])
    if replicate_report(height_stats, "EXPERIMENT STATISTICS — normalized height",
                        joinpath(plot_path_global, "height_replicate_statistics.csv"))
        by_label = Dict(st.label => st for st in height_stats)
        _stat_fig(((norm_est_runs, by_label["h_est/h_gt"], def_red, :solid, "Estimation"),),
                  stat_time, joinpath(plot_path_global, "h_normalized_5_ci.pdf");
                  xlabel=L"\mathrm{Time\;[s]}", ylabel=L"h_{\mathrm{est}}/h_{\mathrm{gt}}",
                  xlims=(0, end_obs_win), hline=1.0)
        _stat_fig(((rel_err_est_runs, by_label["rel. error est [%]"], def_red, :solid, "Estimation"),),
                  stat_time, joinpath(plot_path_global, "relative_height_error_5_ci.pdf");
                  xlabel=L"\mathrm{Time\;[s]}",
                  ylabel=latexstring("Relative Height Error [\$\\%\$]"),
                  xlims=(0, end_obs_win))
    end

    Plots.savefig(contour_plt, joinpath(plot_path_global,"contour_comparison_5.pdf"))
    Plots.savefig(contour_plt_zoom, joinpath(plot_path_global,"contour_comparison_5_zoomed.pdf"))
    Plots.savefig(contour_plt_zoom_05, joinpath(plot_path_global,"contour_comparison_zoom_05.pdf"))
    Plots.savefig(contour_plt_zoom_10, joinpath(plot_path_global,"contour_comparison_zoom_10.pdf"))
    Plots.savefig(contour_plt_zoom_15, joinpath(plot_path_global,"contour_comparison_zoom_15.pdf"))
    Plots.savefig(contour_plt_zoom_20, joinpath(plot_path_global,"contour_comparison_zoom_20.pdf"))
    Plots.savefig(cont_plt_legend, joinpath(plot_path_global,"slip_legend_vertical.pdf"))
    Plots.savefig(slip_case_legend, joinpath(plot_path_global,"slip_case_legend.pdf"))
    @info "Saved plots to $plot_path_global"
end

"""
    post_analysis_bulk(filepath_gt_, filepath, avoid_list)

Cross-experiment post-analysis and figure generation for `viscosity_type ==
"bulk_viscosity"` results — the bulk-viscosity counterpart of
`post_analysis_const`. Walks the same nested `dir/elem_size_folder/
simtime_folder/window` result tree under `filepath` (skipping entries in
`avoid_list`, `"post_analysis_global"`, and a few hardcoded element-size
folders), and overlays each leaf's time-windowed η/β estimates and
height/relative-height-error against the ground truth in `filepath_gt_`, at
three nested levels:

- **global** (across all top-level `dir`s): normalized η/β and
  height/normalized-height/relative-error, bucketed by simulation-time case
  (5/10 s).
- **per-`dir`** (across element sizes): mesh-convergence plots of η, β, and
  height/relative-height-error.
- **per-`elem_size_folder`** (across simulation-time windows): the same
  comparison plots, plus cost-surface slice-direction plots and (per
  simulation-time folder) height-error plots aggregated across windows.

# Arguments
- `filepath_gt_::String`: ground-truth data root each experiment is compared against.
- `filepath::String`: root of the result tree to walk.
- `avoid_list`: top-level `dir` names under `filepath` to skip.

# Returns
None. All outputs are PDFs written under `plot_path_global`/
`plot_path_elems`/`plot_path_sim_time`/`plot_path_noise`.
"""
function post_analysis_bulk(filepath_gt_::String, filepath::String, avoid_list; method::Union{Nothing,String}=nothing)

    # Relative MOSD (Max/Mean Of Surface Depth), one series per experiment, for the
    # cross-experiment statistics at the end. Synthetic/simulated ground truth only: the
    # physical data has no ground-truth mesh, so there is no surface to normalize against
    # and `post_analysis_real` does not compute this.
    mosd_est_runs  = Vector{Vector{Float64}}()   # MOSD_est  / MOSD_gt
    mosd_pred_runs = Vector{Vector{Float64}}()   # MOSD_pred / MOSD_gt
    mosd_pred_segs = Vector{Vector{Tuple{UnitRange{Int},Vector{Float64}}}}()  # per-window, for drawing
    mosd_time = Float64[]

    # Validation metrics, normalized per experiment so they are comparable across runs made
    # at different η/β.
    norm_est_runs    = Vector{Vector{Float64}}()   # h_est / h_gt
    rel_err_est_runs = Vector{Vector{Float64}}()   # |h_est - h_gt| / h_gt · 100
    norm_pred_runs    = Vector{Vector{Float64}}()   # h_pred / h_gt
    rel_err_pred_runs = Vector{Vector{Float64}}()
    norm_pred_segs    = Vector{Vector{Tuple{UnitRange{Int},Vector{Float64}}}}()
    rel_err_pred_segs = Vector{Vector{Tuple{UnitRange{Int},Vector{Float64}}}}()
    stat_time = Float64[]
    stat_windows = Float64[]                        # window boundaries, for the vlines

    dir_list = readdir(filepath)

    η_norm_plot_5 = _fig(legend_column=3)
    Plots.hline!(η_norm_plot_5, [1.0],  linestyle=:dash, label=false, color=:black)
    _label!(η_norm_plot_5, L"\mathrm{Iterations}", L"\eta_{\mathrm{est}}/\eta^{\mathrm{avg}}_{\mathrm{gt}}(t)")

    β_norm_plot_5 = _fig(legend_column=3)
    Plots.hline!(β_norm_plot_5, [1.0],  linestyle=:dash, label=false, color=:black)
    _label!(β_norm_plot_5, L"\mathrm{Iterations}", L"\beta_{\mathrm{est}}/\beta_{\mathrm{gt}}")

    h_plot_5 = _fig(legend_column=3)
    _label!(h_plot_5, L"\mathrm{Time\;[s]}", L"h_{\mathrm{est}}\;\mathrm{[mm]}"; xlims=(0, end_obs_win))

    h_norm_plot_5 = _fig(legend_column=3)
    Plots.hline!(h_norm_plot_5, [1.0],  linestyle=:dash, label=false, color=:black)
    _label!(h_norm_plot_5, L"\mathrm{Time\;[s]}", L"h_{\mathrm{est}}/h_{\mathrm{gt}}"; xlims=(0, end_obs_win), ylims=(y_lims_h_norm))

    rel_height_error_glob_plot_5 = _fig(legend_column=3)
    _label!(rel_height_error_glob_plot_5, L"\mathrm{Time\;[s]}", latexstring("Relative Height Error [\$\\%\$]"); xlims=(0, end_obs_win), ylims=(y_lims_rel_error))

    η_y_max = 0
    β_y_max = 0
    for dir in dir_list
        if dir in avoid_list || dir == "post_analysis_global"
            @debug "Skipping directory: $dir"
            continue
        end
        filepath_dir = joinpath(filepath, dir)
        filepath_gt = joinpath(filepath_gt_, dir)

        printstyled("Processing directory: $(filepath_dir)\n", color=:green)
        sim_params = read_json(datapath(filepath_gt,"sim_params"))

        η_gt = sim_params["η"]
        β_gt = sim_params["β"]  
        gt_h_ = readdlm(datapath(filepath_gt,"h.csv"), ',', Float64)

        sim_time = sim_params["simulation_time"]    
        t_steps = sim_params["time_steps"]

        time = collect(Float64, range(start=0, stop=sim_time, step=t_steps))

        groups = collect_experiment_groups(filepath_dir; method=method)
        elem_size_folders = _elem_size_folders(groups)

        elem_η_plt = _fig(margins=:none)
        Plots.hline!(elem_η_plt, [η_gt[1]], label="Ground truth η",  left_margin=plt_lft_margin)
        _label!(elem_η_plt, L"\mathrm{Iterations}", L"\eta\;\mathrm{[kPa\, s]}")

        elem_β_plt = _fig(margins=:none)
        Plots.hline!(elem_β_plt, [β_gt[1]], label="Ground truth β",  left_margin=plt_lft_margin)
        _label!(elem_β_plt, L"\mathrm{Iterations}", L"\beta\;\mathrm{[Pa\, s \, m]}")

        elem_η_norm_plt = _fig(margins=:none)
        Plots.hline!(elem_η_norm_plt, [1.0], label=false,  left_margin=plt_lft_margin, linestyle=:dash)
        _label!(elem_η_norm_plt, L"\mathrm{Iterations}", L"\eta_{\mathrm{est}}/\eta_{\mathrm{gt}}")

        elem_β_norm_plt = _fig(margins=:none)
        Plots.hline!(elem_β_norm_plt, [1.0], label=false,  left_margin=plt_lft_margin, linestyle=:dash)
        _label!(elem_β_norm_plt, L"\mathrm{Iterations}", L"\beta_{\mathrm{est}}/\beta_{\mathrm{gt}}")

        elem_rel_height_error_plt = _fig(margins=:none)
        Plots.plot!(elem_rel_height_error_plt, [], left_margin=plt_lft_margin, label=false)
        _label!(elem_rel_height_error_plt, L"\mathrm{Time\;[s]}", L"\mathrm{Relative\;Height\;Error}"; xlims=(0, end_obs_win))

        elem_height_plt = _fig(margins=:none)
        Plots.plot!(elem_height_plt, [], left_margin=plt_lft_margin, label=false)
        _label!(elem_height_plt, L"\mathrm{Time\;[s]}", L"\mathrm{Height\;[mm]}"; xlims=(0, end_obs_win))

        for elem_size_folder_ in elem_size_folders
            if elem_size_folder_ == "post_analysis" || elem_size_folder_ == "Q2_2" || elem_size_folder_ == "Q2_4" || elem_size_folder_ == "Q2_8"
                continue
            end

            elem_size_folder = joinpath(filepath_dir, elem_size_folder_)
            printstyled("Processing element size folder: $(elem_size_folder)\n", color=:blue)  
            sim_time_folders = _sim_time_folders(groups, elem_size_folder_)

            sim_window_η_plt = _fig(margins=:none)
            Plots.hline!(sim_window_η_plt, [η_gt[1]], label="Ground truth η")
            _label!(sim_window_η_plt, L"\mathrm{Iterations}", latexstring("\$\\eta\$ [kPa s]"))

            sim_window_β_plt = _fig(margins=:none)
            Plots.hline!(sim_window_β_plt, [β_gt[1]], label="Ground truth β")
            _label!(sim_window_β_plt, L"\mathrm{Iterations}", latexstring("\$\\beta\$ [MPa s m\$^{-1}\$]"))

            sim_window_η_norm_plt = _fig(margins=:none)
            Plots.hline!(sim_window_η_norm_plt, [1.0], label=false,  left_margin=plt_lft_margin, linestyle=:dash)
            _label!(sim_window_η_norm_plt, L"\mathrm{Iterations}", L"\eta_{\mathrm{est}}/\eta_{\mathrm{gt}}")

            sim_window_β_norm_plt = _fig(margins=:none)
            Plots.hline!(sim_window_β_norm_plt, [1.0], label=false,  left_margin=plt_lft_margin, linestyle=:dash)
            _label!(sim_window_β_norm_plt, L"\mathrm{Iterations}", L"\beta_{\mathrm{est}}/\beta_{\mathrm{gt}}")

            sim_window_rel_height_error_plt = _fig(margins=:none)
            Plots.plot!(sim_window_rel_height_error_plt, [], left_margin=plt_lft_margin, label=false)
            _label!(sim_window_rel_height_error_plt, L"\mathrm{Time\;[s]}", L"\mathrm{Relative\;Height\;Error}"; xlims=(0, end_obs_win))

            sim_window_height_plt = _fig(margins=:none)
            Plots.plot!(sim_window_height_plt, [], left_margin=plt_lft_margin, label=false)
            _label!(sim_window_height_plt, L"\mathrm{Time\;[s]}", L"\mathrm{Height\;[mm]}"; xlims=(0, end_obs_win))

            plt_slices = set_plot(fs, sz=(350, 750))
            Plots.vline!(plt_slices, [0.0], color=:blue, linestyle=:dash, label="Minimum",  lw=1)
            _label!(plt_slices, L"\mathrm{Distance\;from\;minimum\;[px]}", L"\mathrm{Cost}"; ylims=(0, 50))

            for sim_time_folder_ in sim_time_folders

                if sim_time_folder_ == "post_analysis_time" || sim_time_folder_ == "Results" || sim_time_folder_ == "simtime_2.0" || sim_time_folder_ == "simtime_30.0"
                    continue
                end

                height_error_plt = _fig(margins=:none)
                Plots.plot!(height_error_plt, [], label=false)
                _label!(height_error_plt, L"\mathrm{Time\;[s]}", L"\mathrm{Height\;Error\;[mm]}"; xlims=(0, end_obs_win))

                rel_height_error_plt = _fig(margins=:none)
                Plots.plot!(rel_height_error_plt, [], label=false)
                _label!(rel_height_error_plt, L"\mathrm{Time\;[s]}", L"\mathrm{Height\;Error\;[mm]}"; xlims=(0, end_obs_win))

                group = _find_group(groups, elem_size_folder_, sim_time_folder_, "noise_0.0")
                if group === nothing
                    @debug "Skipping $(joinpath(elem_size_folder_, sim_time_folder_)): no noise_0.0 group found"
                    continue
                end
                leaf = _single_leaf(group)
                leaf === nothing && continue

                sim_time_folder = group.step_path

                printstyled("Processing simulation time folder: $(sim_time_folder)\n", color=:cyan)

                window_dirs = readdir(leaf.exp_path)
                for window_dir in window_dirs
                    if window_dir == "Results" || window_dir == "post_analysis_window" || window_dir == "single_window" || window_dir == "post_analysis_noise"
                        @debug "Skipping directory: $window_dir"
                        continue
                    end
                    win_exp_path = joinpath(leaf.exp_path, window_dir)
                    _skip_incomplete(win_exp_path) && continue

                    @debug "Processing window: $win_exp_path"
                    exp_params = read_json(datapath(win_exp_path,"experiment_parameters"))
                    sim_time_exp = exp_params["sim_time_exp"]
                    data_type = exp_params["data_type"]
                    noise_level = exp_params["noise_level"]
                    ne = exp_params["ne_exp"]

                    if data_type  == "synthetic"
                        ObsDataList, splinexObs, splineyObs = read_csv(datapath(filepath_gt,"img_data",leaf.view_folder,"contour_data"))
                        @info "Data type $data_type Reading synthetic contour data of $(length(ObsDataList)) time steps"
                    elseif data_type == "simulated"
                        ObsDataList, splinexObs, splineyObs = read_csv(datapath(filepath_gt,"sim_data",leaf.view_folder,"contour_data"))
                        @info "Data type $data_type Reading simulated contour data from $(datapath(filepath_gt,"sim_data",leaf.view_folder,"contour_data")) of $(length(ObsDataList)) time steps"
                    else
                        error("Unknown data type: $data_type")
                    end

                    obs_border_pt_lst, nSplinex, nSpliney, pd = add_noise(ObsDataList, nFactor=0.0)

                    local η_gt::Vector{Float64} = float.(sim_params["η"])
                    local β_gt  = sim_params["β"]

                    est_ηpList = readdlm(datapath(win_exp_path,"est_η.csv"), ',', Float64)
                    est_βpList = readdlm(datapath(win_exp_path,"est_β.csv"), ',', Float64)
                    avg_ηList = readdlm(datapath(win_exp_path,"avg_η.csv"), ',', Float64)
                    est_h_list = readdlm(datapath(win_exp_path,"est_h.csv"), ',', Float64)
                    pred_h_list = isfile(datapath(win_exp_path,"pred_h.csv")) ?
                                  readdlm(datapath(win_exp_path,"pred_h.csv"), ',', Float64) : nothing

                    data_ranges_ = get_time_windows(datapath(win_exp_path,"window_data","data_ranges.csv"))
                    t_windows = readdlm(datapath(win_exp_path,"window_data","t_windows.csv"),',',Float64)
                    time_windows = readdlm(datapath(win_exp_path,"window_data","time_windows.csv"),',',Float64)

                    sim_border_pt_lst, splinex, spliney = read_csv(datapath(win_exp_path,"sim_data","view_1","2D_border_points_est"))

                    @debug "Time windows: $(time_windows)"
                    obs_time = sum(time_windows)

                    if obs_time != sim_time
                        @warn "Observation time frame $obs_time is less than preset ground truth time frame $sim_time, updating time frame"
                        min_time = min(obs_time, sim_time)
                        sim_time = min_time
                        obs_time = min_time
                    end

                    if obs_time < sim_time_exp
                        @warn "Observation time frame $obs_time is less than experimental simulation time frame $sim_time_exp, switching to observation time frame"
                        sim_time_exp = obs_time
                    end

                    data_point_len = round(Int, obs_time/t_steps)
                    η_gt = η_gt[1:data_point_len]
                    gt_h = gt_h_[1:(data_point_len+1)]
                    sim_border_pt_lst = sim_border_pt_lst[1:data_point_len+1, :]
                    obs_border_pt_lst = obs_border_pt_lst[1:data_point_len+1, :]
                    est_h = est_h_list[1:data_point_len+1, :]
                    # The ground truth can run longer than the fitted window — the Carreau runs do —
                    # so the time axis has to be cut to the data, not to the ground truth's duration.
                    time = time[1:min(length(time), length(gt_h))]

                    height_error = abs.(est_h - gt_h)
                    rel_height_error = height_error ./ gt_h * 100.0

                    # Relative MOSD, computed the same way `replot` plots it. Guarded on the
                    # ground-truth surface existing, so a partially written experiment is
                    # skipped rather than aborting the whole post-analysis.
                    if data_type != "physical" && isdir(datapath(filepath_gt,"sim_data","surface_nodes"))
                        gt_surf, _, _   = read_csv(datapath(filepath_gt,"sim_data","surface_nodes"))
                        est_surf, _, _  = read_csv(datapath(win_exp_path,"sim_data","3D_surface_points_est"))
                        pred_surf, _, _ = read_csv(datapath(win_exp_path,"sim_data","3D_surface_points_pred"))

                        # Carreau ground truth stores h as an Int; get_surface_mosd wants Float64.
                        h_gt      = float(sim_params["h"])
                        pose_gt   = _read_obj_pose(sim_params)
                        gt_surf   = gt_surf[1:(data_point_len+1)]
                        mosd_gt   = get_surface_mosd(gt_surf, pose_gt, h_gt)
                        mosd_est  = get_surface_mosd(est_surf, pose_gt, h_gt)
                        mosd_pred = get_surface_mosd(pred_surf, pose_gt, h_gt)

                        # `mosd_est` is a time series, but `mosd_pred` is the per-window
                        # predictions concatenated — one sample longer per window boundary —
                        # so it has to be mapped onto the time axis before dividing by a
                        # time-indexed ground truth.
                        ne_ = min(length(mosd_est), length(mosd_gt))
                        push!(mosd_est_runs, mosd_est[1:ne_] ./ mosd_gt[1:ne_])
                        push!(mosd_pred_runs,
                              _align_windowed(mosd_pred, data_ranges_, ne_) ./ mosd_gt[1:ne_])
                        push!(mosd_pred_segs,
                              [(r, v ./ mosd_gt[r]) for (r, v) in _windowed_series(mosd_pred, data_ranges_, ne_)])
                        isempty(mosd_time) && (mosd_time = collect(time[1:min(ne_, length(time))]))
                    end

                    # normalize by ground-truth                      
                    est_η_norm = est_ηpList./avg_ηList
                    est_β_norm = est_βpList./β_gt

                    Plots.plot!(sim_window_η_plt, est_ηpList, label=string("Window - $(sim_time)s"), marker=1, legend=:outerbottom, legend_column=2)
                    Plots.plot!(sim_window_β_plt, est_βpList, label=string("Window - $(sim_time)s"), marker=1, legend=:outerbottom, legend_column=2)

                    Plots.plot!(sim_window_rel_height_error_plt, time, rel_height_error, label=string("Window - $(sim_time)s"," - ne: ",ne), legend=:outerbottom, legend_column=2)
                    Plots.plot!(sim_window_height_plt, time, est_h, label=string("Window - $(sim_time)s"), legend=:outerbottom, legend_column=2)
                    Plots.plot!(sim_window_height_plt, time, gt_h, label="Ground truth", legend=:outerbottom, legend_column=2)

                    Plots.plot!(sim_window_η_norm_plt, est_η_norm, label=string("Window - $(sim_time)s"), marker=1, legend=:outerbottom, legend_column=2)
                    Plots.plot!(sim_window_β_norm_plt, est_β_norm, label=string("Window - $(sim_time)s"), marker=1, legend=:outerbottom, legend_column=2)

                    if sim_time_exp == 5.0
                        if ne == 6
                            η_y_max = max(η_y_max, maximum(est_η_norm))
                            β_y_max = max(β_y_max, maximum(est_β_norm))
                            t_prev = 0.1
                            for ti::Int in 1:(size(data_ranges_, 1)-1)
                                t = t_windows[ti]
                                data_range_ = data_ranges_[ti]
                                t_win = collect(range(start=t_prev, stop=t, step=t_steps))

                                Plots.vline!(η_norm_plot_5, [t], color=:gray, linestyle=:dash, label=false)
                                Plots.vline!(β_norm_plot_5, [t], color=:gray, linestyle=:dash, label=false)

                                if ti == 1
                                    Plots.plot!(η_norm_plot_5, t_win, est_η_norm[data_range_], label=latexstring("\$\\beta_{\\mathrm{gt}}:$(β_gt[1])\\,\\mathrm{MPa\\,s\\,\\mathrm{m^{-1}}}\$"), color=_exp_color(dir))
                                    Plots.plot!(β_norm_plot_5, t_win, est_β_norm[data_range_], label=latexstring("\$\\beta_{\\mathrm{gt}}:$(β_gt[1])\\,\\mathrm{MPa\\,s\\,\\mathrm{m^{-1}}}\$"), color=_exp_color(dir))
                                else
                                    Plots.plot!(η_norm_plot_5, t_win, est_η_norm[data_range_], label=false, color=_exp_color(dir))
                                    Plots.plot!(β_norm_plot_5, t_win, est_β_norm[data_range_], label=false, color=_exp_color(dir))
                                end
                                t_prev = t+t_steps
                            end
                            Plots.ylims!(η_norm_plot_5, 0, max(η_y_max*1.1, 2))
                            Plots.ylims!(β_norm_plot_5, 0, max(β_y_max*1.1, 2))

                            Plots.plot!(h_plot_5, time, est_h, label=latexstring("\$\\beta_{\\mathrm{gt}}:$(β_gt[1])\\,\\mathrm{MPa\\,s\\,\\mathrm{m^{-1}}}\$"), color=_exp_color(dir))
                            Plots.plot!(h_plot_5, time, gt_h, label=false, linestyle=:dash, color=_exp_color(dir))

                            t = collect(range(start=0, stop=sim_time, step=t_steps))
                            Plots.plot!(rel_height_error_glob_plot_5, t, rel_height_error, label=latexstring("\$\\beta_{\\mathrm{gt}}:$(β_gt[1])\\,\\mathrm{MPa\\,s\\,\\mathrm{m^{-1}}}\$"), color=_exp_color(dir))
                            Plots.plot!(h_norm_plot_5, t, est_h./gt_h, label=latexstring("\$\\beta_{\\mathrm{gt}}:$(β_gt[1])\\,\\mathrm{MPa\\,s\\,\\mathrm{m^{-1}}}\$"), color=_exp_color(dir))

                            push!(norm_est_runs, vec(est_h ./ gt_h))
                            push!(rel_err_est_runs, vec(rel_height_error))
                            if !isnothing(pred_h_list)
                                gtp = vec(gt_h)
                                ph  = _align_windowed(vec(pred_h_list), data_ranges_, length(gtp))
                                push!(norm_pred_runs, ph ./ gtp)
                                push!(rel_err_pred_runs, abs.(ph .- gtp) ./ gtp .* 100.0)
                                segs = _windowed_series(vec(pred_h_list), data_ranges_, length(gtp))
                                push!(norm_pred_segs, [(r, v ./ gtp[r]) for (r, v) in segs])
                                push!(rel_err_pred_segs, [(r, abs.(v .- gtp[r]) ./ gtp[r] .* 100.0) for (r, v) in segs])
                            end
                            isempty(stat_time) && (stat_time = collect(t))
                            isempty(stat_windows) && (stat_windows = vec(Float64.(t_windows)))
                            if dir == "1"
                                for ti::Int in 1:(size(data_ranges_, 1)-1)
                                    t_win = t_windows[ti]
                                    data_range_ = data_ranges_[ti]

                                    Plots.vline!(h_norm_plot_5, [t_win], color=:gray, linestyle=:dash, label=false)
                                    Plots.vline!(rel_height_error_glob_plot_5, [t_win], color=:gray, linestyle=:dash, label=false)
                                end
                            end

                        end

                    end
                end
                plot_path_noise = _method_dir(joinpath(sim_time_folder,"post_analysis_noise","plots"), method)
                set_file(plot_path_noise)
                @info "Saving plots to $plot_path_noise"
                Plots.savefig(height_error_plt, joinpath(plot_path_noise,"height_error.pdf"))
                Plots.savefig(rel_height_error_plt, joinpath(plot_path_noise,"relative_height_error.pdf"))
                @info "Saved plots to $plot_path_noise"
            end
            plot_path_sim_time = _method_dir(joinpath(elem_size_folder,"post_analysis_time","plots"), method)
            set_file(plot_path_sim_time)

            Plots.savefig(sim_window_η_plt, joinpath(plot_path_sim_time,"η.pdf"))
            Plots.savefig(sim_window_β_plt, joinpath(plot_path_sim_time,"β.pdf"))
            Plots.savefig(sim_window_rel_height_error_plt, joinpath(plot_path_sim_time,"relative_height_error.pdf"))
            Plots.savefig(sim_window_height_plt, joinpath(plot_path_sim_time,"height.pdf"))
            Plots.savefig(sim_window_η_norm_plt, joinpath(plot_path_sim_time,"η_normalized.pdf"))
            Plots.savefig(sim_window_β_norm_plt, joinpath(plot_path_sim_time,"β_normalized.pdf"))
            Plots.savefig(plt_slices, joinpath(plot_path_sim_time,"cost_slices_along_directions.pdf"))
            @info "Saved plots to $plot_path_sim_time"
        end

        plot_path_elems = _method_dir(joinpath(filepath_dir,"post_analysis","plots"), method)
        set_file(plot_path_elems)

        Plots.savefig(elem_η_plt, joinpath(plot_path_elems,"η.pdf"))
        Plots.savefig(elem_β_plt, joinpath(plot_path_elems,"β.pdf"))
        Plots.savefig(elem_η_norm_plt, joinpath(plot_path_elems,"η_normalized.pdf"))
        Plots.savefig(elem_β_norm_plt, joinpath(plot_path_elems,"β_normalized.pdf"))
        Plots.savefig(elem_rel_height_error_plt, joinpath(plot_path_elems,"relative_height_error.pdf"))
        Plots.savefig(elem_height_plt, joinpath(plot_path_elems,"height.pdf"))
        @info "Saved plots to $plot_path_elems"
    end
    plot_path_global = _method_dir(joinpath(filepath,"post_analysis_global","plots"), method)
    set_file(plot_path_global)

    Plots.savefig(η_norm_plot_5, joinpath(plot_path_global,"η_normalized_5.pdf"))
    Plots.savefig(β_norm_plot_5, joinpath(plot_path_global,"β_normalized_5.pdf"))
    Plots.savefig(rel_height_error_glob_plot_5, joinpath(plot_path_global,"relative_height_error_5.pdf"))
    Plots.savefig(h_plot_5, joinpath(plot_path_global,"height_comparison_5.pdf"))
    Plots.savefig(h_norm_plot_5, joinpath(plot_path_global,"height_normalized_5.pdf"))

    @info "Saved plots to $plot_path_global"

    # ---- cross-experiment statistics on the relative MOSD --------------------------------
    # Relative MOSD is normalized per experiment against its own ground truth, so it is
    # comparable across experiments that were run at different η/β — the spread describes how
    # well the surface depth is recovered, not how the ground truths differ.
    height_stats = filter(!isnothing, [
        normalized_replicate_stats(norm_pred_runs,    "h_pred/h_gt";         reference=1.0),
        normalized_replicate_stats(norm_est_runs,    "h_est/h_gt";          reference=1.0),
        normalized_replicate_stats(rel_err_est_runs, "rel. error est [%]";  reference=0.0, signed=false),
        normalized_replicate_stats(rel_err_pred_runs, "rel. error pred [%]"; reference=0.0, signed=false),
    ])
    if replicate_report(height_stats, "EXPERIMENT STATISTICS — normalized height",
                        joinpath(plot_path_global, "height_replicate_statistics.csv"))
        by_label = Dict(st.label => st for st in height_stats)
        _stat_fig(((norm_est_runs, by_label["h_est/h_gt"], def_red, :solid, "Estimation"),
                   (norm_pred_segs, get(by_label, "h_pred/h_gt", nothing), def_blue, :dash, "Prediction")),
                  stat_time, joinpath(plot_path_global, "height_normalized_5_ci.pdf");
                  xlabel=L"\mathrm{Time\;[s]}", ylabel=L"h_{\mathrm{est}}/h_{\mathrm{gt}}",
                  xlims=(0, end_obs_win), hline=1.0, vlines=stat_windows)
        _stat_fig(((rel_err_est_runs, by_label["rel. error est [%]"], def_red, :solid, "Estimation"),
                   (rel_err_pred_segs, get(by_label, "rel. error pred [%]", nothing), def_blue, :dash, "Prediction")),
                  stat_time, joinpath(plot_path_global, "relative_height_error_5_ci.pdf");
                  xlabel=L"\mathrm{Time\;[s]}",
                  ylabel=latexstring("Relative Height Error [\$\\%\$]"),
                  xlims=(0, end_obs_win), vlines=stat_windows)
    end

    mosd_stats = filter(!isnothing, [
        normalized_replicate_stats(mosd_est_runs,  "MOSD_est/MOSD_gt";  reference=1.0),
        normalized_replicate_stats(mosd_pred_runs, "MOSD_pred/MOSD_gt"; reference=1.0),
    ])

    if !replicate_report(mosd_stats, "EXPERIMENT STATISTICS — relative MOSD",
                         joinpath(plot_path_global, "mosd_statistics.csv"))
        @info "No ground-truth surface found — skipping MOSD statistics"
    else
        # Same combined figure as the height metrics: estimation solid red, prediction
        # dashed blue, both with their runs underneath.
        _stat_fig(((mosd_est_runs,  mosd_stats[1], def_red,  :solid, "Estimation"),
                   (mosd_pred_segs, length(mosd_stats) > 1 ? mosd_stats[2] : nothing,
                    def_blue, :dash, "Prediction")),
                  mosd_time, joinpath(plot_path_global, "mosd_relative_ci.pdf");
                  xlabel=L"\mathrm{Time\;[s]}", ylabel=L"\mathrm{Relative\;MOSD}",
                  hline=1.0, vlines=stat_windows)
    end
end

"""
    post_analysis_real(filepath_gt_, filepath, avoid_list)

Cross-experiment post-analysis and figure generation for physical (real,
measured) experiments — the counterpart of `post_analysis_bulk` for
`data_type == "physical"`, where no ground-truth η/β parameters exist, only a
measured height `h_m(t)`. Walks the same nested `dir/elem_size_folder/
simtime_folder/window` result tree under `filepath` (skipping entries in
`avoid_list`, `"post_analysis_global"`, a few hardcoded element-size folders,
and `simtime_2.0`), and overlays each leaf's time-windowed η(t)/β(t)
estimates and height/relative-height-error against
the measured height in `filepath_gt_`, at three nested levels: global (across
all top-level `dir`s, bucketed by simulation-time case), per-`dir` (across
element sizes), and per-`elem_size_folder` (across simulation-time windows,
plus per-window height-error plots).

# Arguments
- `filepath_gt_::String`: ground-truth (measured height) data root each experiment
  is compared against.
- `filepath::String`: root of the result tree to walk.
- `avoid_list`: top-level `dir` names under `filepath` to skip.

# Returns
None. All outputs are PDFs written under `plot_path_global`/
`plot_path_elems`/`plot_path_sim_time`.
"""
function post_analysis_real(filepath_gt_::String, filepath::String, avoid_list; method::Union{Nothing,String}=nothing)
    dir_list = readdir(filepath)

    # Normalized height series, one entry per replicate, collected for the cross-run
    # statistics at the end. Normalized exactly as the overlay plots below are, so the
    # numbers and the figures describe the same quantity.
    norm_est_runs  = Vector{Vector{Float64}}()   # h_est  / h_m
    norm_pred_runs = Vector{Vector{Float64}}()   # h_pred / h_m
    rel_err_est_runs  = Vector{Vector{Float64}}()   # |h_est  - h_m| / h_m · 100
    rel_err_pred_runs = Vector{Vector{Float64}}()
    cpd_runs = Vector{Vector{Float64}}()            # closest-point distance error [px]
    stat_time = Float64[]                           # shared time grid for the CI bands
    stat_windows = Float64[]                        # window boundaries, for the vlines
    norm_pred_segs    = Vector{Vector{Tuple{UnitRange{Int},Vector{Float64}}}}()
    rel_err_pred_segs = Vector{Vector{Tuple{UnitRange{Int},Vector{Float64}}}}()

    η_plot_5 = _fig(legend_column=3)
    _label!(η_plot_5, L"\mathrm{Time\;[s]}", latexstring("\$\\eta_{\\mathrm{est}}(t)\$ [kPa s]"))

    β_plot_5 = _fig(legend_column=3)
    _label!(β_plot_5, L"\mathrm{Time\;[s]}", latexstring("\$\\beta_{\\mathrm{est}}(t)\$ [MPa s m\$^{-1}\$]"))

    gt_h_plot = _fig(legend_column=3)
    _label!(gt_h_plot, L"\mathrm{Time\;[s]}", L"h_{\mathrm{m}}\;\mathrm{[mm]}"; xlims=(0, end_obs_win))

    h_plot_5 = _fig(legend_column=3)
    _label!(h_plot_5, L"\mathrm{Time\;[s]}", L"h\;\mathrm{[mm]}"; xlims=(0, end_obs_win))

    h_plot_est_5 = _fig(legend_column=3)
    _label!(h_plot_est_5, L"\mathrm{Time\;[s]}", L"h_{\mathrm{est}}\;\mathrm{[mm]}"; xlims=(0, end_obs_win))

    h_norm_plot_5 = _fig(legend_column=3)
    Plots.hline!(h_norm_plot_5, [1.0],  linestyle=:dash, label=false, color=:black)
    _label!(h_norm_plot_5, L"\mathrm{Time\;[s]}", L"h_{\mathrm{est}}/h_{\mathrm{m}}"; xlims=(0, end_obs_win), ylims=(y_lims_h_norm))

    rel_height_error_glob_plot_5 = _fig(legend_column=3)
    _label!(rel_height_error_glob_plot_5, L"\mathrm{Time\;[s]}", latexstring("Relative Height Error [\$\\%\$]"); ylims=(y_lims_rel_error))

    for dir in dir_list
        if dir in avoid_list || dir == "post_analysis_global"
            @debug "Skipping directory: $dir"
            continue
        end
        filepath_dir = joinpath(filepath, dir)
        filepath_gt = joinpath(filepath_gt_, dir)

        printstyled("Processing directory: $(filepath_dir)\n", color=:green)
        sim_params = read_json(datapath(filepath_gt,"sim_params"))
        sim_time = sim_params["simulation_time"]    
        t_steps = sim_params["time_steps"]
        gt_h_ = readdlm(datapath(filepath_gt,"h.csv"), ',', Float64)
        time = collect(Float64, range(start=0, stop=sim_time, step=t_steps))

        groups = collect_experiment_groups(filepath_dir; method=method)
        elem_size_folders = _elem_size_folders(groups)

        for elem_size_folder_ in elem_size_folders
            if elem_size_folder_ == "post_analysis" || elem_size_folder_ == "Q2_2" || elem_size_folder_ == "Q2_4" || elem_size_folder_ == "Q2_8"
                continue
            end

            elem_size_folder = joinpath(filepath_dir, elem_size_folder_)
            printstyled("Processing element size folder: $(elem_size_folder)\n", color=:blue)  
            sim_time_folders = _sim_time_folders(groups, elem_size_folder_)

            sim_window_η_plt = _fig(margins=:none)
            _label!(sim_window_η_plt, L"\mathrm{Iterations}", L"\eta\;\mathrm{[kPa\, s]}")

            sim_window_β_plt = _fig(margins=:none)
            _label!(sim_window_β_plt, L"\mathrm{Iterations}", L"\beta\;\mathrm{[Pa\, s \, m]}")

            sim_window_rel_height_error_plt = _fig(margins=:none)
            Plots.plot!(sim_window_rel_height_error_plt, [], left_margin=plt_lft_margin, label=false)
            _label!(sim_window_rel_height_error_plt, L"\mathrm{Time\;[s]}", L"\mathrm{Relative\;Height\;Error}"; xlims=(0, end_obs_win))

            sim_window_height_plt = _fig(margins=:none)
            Plots.plot!(sim_window_height_plt, [], left_margin=plt_lft_margin, label=false)
            _label!(sim_window_height_plt, L"\mathrm{Time\;[s]}", L"\mathrm{Height\;[mm]}"; xlims=(0, end_obs_win))

            for sim_time_folder_ in sim_time_folders

                if sim_time_folder_ == "post_analysis_time" || sim_time_folder_ == "Results" || sim_time_folder_ == "simtime_2.0" #|| sim_time_folder_ == "simtime_10.0"
                    continue
                end

                height_error_plt = _fig(margins=:none)
                Plots.plot!(height_error_plt, [], label=false)
                _label!(height_error_plt, L"\mathrm{Time\;[s]}", L"\mathrm{Height\;Error\;[mm]}"; xlims=(0, end_obs_win))

                rel_height_error_plt = _fig(margins=:none)
                Plots.plot!(rel_height_error_plt, [], label=false)
                _label!(rel_height_error_plt, L"\mathrm{Time\;[s]}", L"\mathrm{Height\;Error\;[mm]}"; xlims=(0, end_obs_win))

                group = _find_group(groups, elem_size_folder_, sim_time_folder_, "noise_0.0")
                if group === nothing
                    @debug "Skipping $(joinpath(elem_size_folder_, sim_time_folder_)): no noise_0.0 group found"
                    continue
                end
                leaf = _single_leaf(group)
                leaf === nothing && continue

                sim_time_folder = group.step_path

                printstyled("Processing simulation time folder: $(sim_time_folder)\n", color=:cyan)

                window_dirs = readdir(leaf.exp_path)
                for window_dir in window_dirs
                    if window_dir == "Results" || window_dir == "post_analysis_window" || window_dir == "single_window" || window_dir == "post_analysis_noise"
                        @debug "Skipping directory: $window_dir"
                        continue
                    end
                    win_exp_path = joinpath(leaf.exp_path, window_dir)
                    _skip_incomplete(win_exp_path) && continue

                    @debug "Processing window: $win_exp_path"
                    exp_params = read_json(datapath(win_exp_path,"experiment_parameters"))
                    sim_time_exp = exp_params["sim_time_exp"]
                    data_type = exp_params["data_type"]
                    ne = exp_params["ne_exp"]

                    ObsDataList, splinexObs, splineyObs = read_csv(datapath(filepath_gt,"img_data","contour_data"))
                    # The video starts before the press engages, so the contour series has to be
                    # indexed by `compressed_frames` exactly as `optimize` does it. Taking raw
                    # frame k against simulation step k lands ~100 frames early and scores a fit
                    # that never happened.
                    _meta = read_json(datapath(filepath_gt,"video_metadata"))
                    ObsDataList = ObsDataList[_meta["compressed_frames"]]
                    obs_border_pt_lst, nSplinex, nSpliney, pd = add_noise(ObsDataList, nFactor=0.0)

                    est_ηpList = readdlm(datapath(win_exp_path,"est_η.csv"), ',', Float64)
                    est_βpList = readdlm(datapath(win_exp_path,"est_β.csv"), ',', Float64)
                    est_h_list = readdlm(datapath(win_exp_path,"est_h.csv"), ',', Float64)
                    pred_h_list = isfile(datapath(win_exp_path,"pred_h.csv")) ?
                                  readdlm(datapath(win_exp_path,"pred_h.csv"), ',', Float64) : nothing

                    data_ranges_ = get_time_windows(datapath(win_exp_path,"window_data","data_ranges.csv"))
                    t_windows = readdlm(datapath(win_exp_path,"window_data","t_windows.csv"),',',Float64)
                    time_windows = readdlm(datapath(win_exp_path,"window_data","time_windows.csv"),',',Float64)

                    sim_border_pt_lst, splinex, spliney = read_csv(datapath(win_exp_path,"sim_data","view_1","2D_border_points_est"))

                    @debug "Time windows: $(time_windows)"
                    obs_time = sum(time_windows)

                    if obs_time != sim_time
                        @warn "Observation time frame $obs_time is less than preset ground truth time frame $sim_time, updating time frame"
                        min_time = min(obs_time, sim_time)
                        sim_time = min_time
                        obs_time = min_time
                    end

                    if obs_time < sim_time_exp
                        @warn "Observation time frame $obs_time is less than experimental simulation time frame $sim_time_exp, switching to observation time frame"
                        sim_time_exp = obs_time
                    end

                    data_point_len = round(Int, obs_time/t_steps)
                    gt_h = gt_h_[1:(data_point_len+1)]
                    sim_border_pt_lst = sim_border_pt_lst[1:data_point_len+1, :]
                    obs_border_pt_lst = obs_border_pt_lst[1:data_point_len+1, :]
                    est_h = est_h_list[1:data_point_len+1, :]
                    time = collect(range(start=0, stop=sim_time, step=t_steps))

                    height_error = abs.(est_h - gt_h)
                    rel_height_error = height_error ./ gt_h * 100.0

                    Plots.plot!(sim_window_η_plt, est_ηpList, label=string("Window - $(sim_time)s"), marker=1, legend=:outerbottom, legend_column=2)
                    Plots.plot!(sim_window_β_plt, est_βpList, label=string("Window - $(sim_time)s"), marker=1, legend=:outerbottom, legend_column=2)

                    Plots.plot!(sim_window_rel_height_error_plt, time, rel_height_error, label=string("Window - $(sim_time)s"," - ne: ",ne), legend=:outerbottom, legend_column=2)
                    Plots.plot!(sim_window_height_plt, time[1:length(est_h)], est_h, label=string("Window - $(sim_time)s"), legend=:outerbottom, legend_column=2)
                    Plots.plot!(sim_window_height_plt, time, gt_h, label="Ground truth", legend=:outerbottom, legend_column=2)

                    if sim_time_exp == 0.5
                        if ne == 6
                            t_prev = 0.1
                            for ti::Int in 1:(size(data_ranges_, 1)-1)
                                t = t_windows[ti]
                                data_range_ = data_ranges_[ti]
                                t_win = collect(range(start=t_prev, stop=t, step=t_steps))

                                Plots.vline!(η_plot_5, [t], color=:gray, linestyle=:dash, label=false)
                                Plots.vline!(β_plot_5, [t], color=:gray, linestyle=:dash, label=false)

                                if ti == 1
                                    Plots.plot!(η_plot_5, t_win, est_ηpList[data_range_], label=string(L"\mathrm{Exp}:\;",dir), color=_exp_color(dir))
                                    Plots.plot!(β_plot_5, t_win, est_βpList[data_range_], label=string(L"\mathrm{Exp}:\;",dir), color=_exp_color(dir))
                                else
                                    Plots.plot!(η_plot_5, t_win, est_ηpList[data_range_], label=false, color=_exp_color(dir))
                                    Plots.plot!(β_plot_5, t_win, est_βpList[data_range_], label=false, color=_exp_color(dir))
                                end
                                t_prev = t+t_steps
                            end

                            Plots.plot!(rel_height_error_glob_plot_5, time, rel_height_error, label=string(L"\mathrm{Exp}:\;",dir), color=_exp_color(dir))
                            Plots.plot!(h_norm_plot_5, time, est_h./gt_h, label=string(L"\mathrm{Exp}:\;",dir), color=_exp_color(dir))

                            # Same normalization as the two plots above, kept for the
                            # cross-replicate statistics written at the end of this function.
                            push!(norm_est_runs, vec(est_h ./ gt_h))
                            push!(rel_err_est_runs, vec(rel_height_error))
                            push!(cpd_runs, vec(contour_cost(sim_border_pt_lst, obs_border_pt_lst)[1]))
                            isempty(stat_time) && (stat_time = collect(time))
                            isempty(stat_windows) && (stat_windows = vec(Float64.(t_windows)))
                            if !isnothing(pred_h_list)
                                # Per-window segments, so map them onto the time axis before
                                # dividing by a time-indexed height.
                                gtp = vec(gt_h)
                                ph  = _align_windowed(vec(pred_h_list), data_ranges_, length(gtp))
                                push!(norm_pred_runs, ph ./ gtp)
                                push!(rel_err_pred_runs, abs.(ph .- gtp) ./ gtp .* 100.0)
                                # Per-window values for drawing; the aligned vectors above
                                # remain what the statistics are computed from.
                                segs = _windowed_series(vec(pred_h_list), data_ranges_, length(gtp))
                                push!(norm_pred_segs, [(r, v ./ gtp[r]) for (r, v) in segs])
                                push!(rel_err_pred_segs,
                                      [(r, abs.(v .- gtp[r]) ./ gtp[r] .* 100.0) for (r, v) in segs])
                            end
                            Plots.plot!(h_plot_est_5, time, est_h, label=string(L"\mathrm{Exp}:\;",dir), color=_exp_color(dir))
                            Plots.plot!(gt_h_plot, time, gt_h, label=string(L"\mathrm{Exp}:\;",dir), color=_exp_color(dir))
                            Plots.plot!(h_plot_5, time, gt_h, label=string(L"\mathrm{Exp}:\;",dir), color=_exp_color(dir), linestyle=:dash)
                            Plots.plot!(h_plot_5, time, est_h, label=false, color=_exp_color(dir))
                            if dir == "1"
                                for ti::Int in 1:(size(data_ranges_, 1))
                                    t_win = t_windows[ti]
                                    data_range_ = data_ranges_[ti]
                                    Plots.vline!(h_norm_plot_5, [t_win], color=:gray, linestyle=:dash, label=false)
                                    Plots.vline!(rel_height_error_glob_plot_5, [t_win], color=:gray, linestyle=:dash, label=false)
                                    Plots.vline!(h_plot_est_5, [t_win], color=:gray, linestyle=:dash, label=false)
                                    Plots.vline!(h_plot_5, [t_win], color=:gray, linestyle=:dash, label=false)
                                end
                            end

                        end

                        Plots.plot!(height_error_plt, time, height_error, label=string("Number of elements: ",ne), marker=1, legend=:outerbottom, legend_column=2)
                        Plots.plot!(rel_height_error_plt, time, rel_height_error, label=string("Number of elements: ",ne), marker=1, legend=:outerbottom, legend_column=2)

                    end
                end
            end
            plot_path_sim_time = _method_dir(joinpath(elem_size_folder,"post_analysis_time","plots"), method)
            set_file(plot_path_sim_time)

            Plots.savefig(sim_window_η_plt, joinpath(plot_path_sim_time,"η.pdf"))
            Plots.savefig(sim_window_β_plt, joinpath(plot_path_sim_time,"β.pdf"))
            Plots.savefig(sim_window_height_plt, joinpath(plot_path_sim_time,"height.pdf"))
            @info "Saved plots to $plot_path_sim_time"
        end
    end
    plot_path_global = _method_dir(joinpath(filepath,"post_analysis_global","plots"), method)
    set_file(plot_path_global)

    Plots.savefig(η_plot_5, joinpath(plot_path_global,"η_plot_5.pdf"))
    Plots.savefig(β_plot_5, joinpath(plot_path_global,"β_plot_5.pdf"))
    Plots.savefig(h_plot_est_5, joinpath(plot_path_global,"height_5.pdf"))
    Plots.savefig(h_plot_5, joinpath(plot_path_global,"height_comparison_5.pdf"))
    Plots.savefig(rel_height_error_glob_plot_5, joinpath(plot_path_global,"relative_height_error_5.pdf"))
    Plots.savefig(h_norm_plot_5, joinpath(plot_path_global,"height_normalized_5.pdf"))
    Plots.savefig(gt_h_plot, joinpath(plot_path_global,"gt_height.pdf"))
    @info "Saved plots to $plot_path_global"

    # ---- cross-replicate statistics on the normalized heights ---------------------------
    # Height is the only quantity these replicates share. η and β belong to the individual
    # specimen and legitimately differ from run to run, so their spread would describe the
    # specimens rather than the method; normalizing by each run's own h_m removes the
    # specimen-to-specimen difference and leaves the agreement between model and measurement.
    stats = filter(!isnothing, [
        normalized_replicate_stats(norm_est_runs,     "h_est/h_m";           reference=1.0),
        normalized_replicate_stats(norm_pred_runs,    "h_pred/h_m";          reference=1.0),
        normalized_replicate_stats(rel_err_est_runs,  "rel. error est [%]";  reference=0.0, signed=false),
        normalized_replicate_stats(rel_err_pred_runs, "rel. error pred [%]"; reference=0.0, signed=false),
        normalized_replicate_stats(cpd_runs,          "closest-pt dist [px]"; reference=0.0, signed=false),
    ])

    if !replicate_report(stats, "REPLICATE STATISTICS — normalized height",
                         joinpath(plot_path_global, "height_replicate_statistics.csv"))
        @warn "No normalized height series collected — skipping replicate statistics"
    else
        by_label = Dict(st.label => st for st in stats)
        out(name) = joinpath(plot_path_global, name)

        # Estimation and prediction share a figure, as they do per experiment.
        _stat_fig(((norm_est_runs,  by_label["h_est/h_m"],  def_red,  :solid, "Estimation"),
                   (norm_pred_segs, get(by_label, "h_pred/h_m", nothing), def_blue, :dash,
                    "Prediction")),
                  stat_time, out("height_normalized_5_ci.pdf");
                  xlabel=L"\mathrm{Time\;[s]}", ylabel=L"h/h_{\mathrm{m}}",
                  xlims=(0, end_obs_win), hline=1.0, vlines=stat_windows)

        _stat_fig(((rel_err_est_runs,  by_label["rel. error est [%]"],  def_red,  :solid, "Estimation"),
                   (rel_err_pred_segs, get(by_label, "rel. error pred [%]", nothing), def_blue, :dash,
                    "Prediction")),
                  stat_time, out("relative_height_error_5_ci.pdf");
                  xlabel=L"\mathrm{Time\;[s]}",
                  ylabel=latexstring("Relative Height Error [\$\\%\$]"),
                  xlims=(0, end_obs_win), vlines=stat_windows)

        if !isempty(cpd_runs)
            _stat_fig(((cpd_runs, by_label["closest-pt dist [px]"], def_red, :solid, "Estimation"),),
                      stat_time, out("closest_point_distance_error_ci.pdf");
                      xlabel=L"\mathrm{Time\;[s]}",
                      ylabel=L"\mathrm{Closest\;Point\;Distance\;[px]}",
                      xlims=(0, end_obs_win), vlines=stat_windows)
        end

    end
end

"""
    optimize_sim(use_parallel=true)

Driver for the mesh-convergence study: builds one `optimize` call
(`exp_params`) per combination of mesh element count (`nz_list`), noise level,
simulation time, and time step for `data_type == "simulated"`, constant
viscosity, reading ground truth from `ground_truth/sim_data/Stokes` and
writing results under `experiments/sim_data/convergence_analysis/
stokes_convergence`. Dispatches the resulting parameter list to
`run_parallel_tasks` when `use_parallel`, otherwise runs `optimize`
sequentially, logging (not raising) any per-task failure via
`_handle_worker_error`.

# Arguments
- `use_parallel::Bool`: dispatch to `run_parallel_tasks` instead of running
  sequentially (default: `true`).

# Returns
None.
"""
function optimize_sim(use_parallel::Bool=true)

    element_shape_u::Symbol = :Hex
    basis_order_u::Int = 2
    element_shape_p::Symbol = :Hex
    basis_order_p::Int = 1
    element_shape_x::Symbol = :Hex
    basis_order_x::Int = 2
    nz_list = Union{Int,Float64}[6]
    mode::Symbol = :exp  # :exp or :conv_exp_mesh
    cost_function_list = [:closest_point] # :chamfer or :closest_point

    if mode == :conv_exp_mesh
        nz_list = Union{Int,Float64}[6]
    end
    dt_list = [0.1] 
    control = "force" # "force" or "velocity"
    viscosity_type_list = ["constant"]
    window = "multi_window"
    camera_matrix::AbstractArray = [[2.39642674e+03, 0.0, 1.00429248e+03] [0.0, 2.40565353e+03, 7.57028161e+02] [0.0, 0.0, 1.0]]'
    filepath_res::String = ""
    param_list = Vector{Dict}(undef, 0)
    geometry::Symbol = :cylinder # :cylinder or :cube

    # Optimization method and its `fit_model` keyword arguments. `opt_method` also
    # becomes a path segment under `dt_*`, so runs with different methods never
    # overwrite each other. Keys in `opt_kwargs` are passed straight through to
    # `fit_model`, so anything it accepts works here (see `_opt_kwargs`).
    opt_method::Symbol = :gn_tikhonov                # :gn, :lm or :gn_tikhonov
    opt_kwargs = Dict{String,Any}()         # e.g. "λ_scale" => 0.01, "line_search_method" => :armijo
    opt_kwargs["λ_scale"] = [0.0001, 0.00002]

    avoid_dirs = ["post_analysis_global","1","2"]
    for viscosity_type in viscosity_type_list
        _filepath_gt = joinpath(resolve_data_path("ground_truth/sim_data/Stokes"), control, viscosity_type, "Hex2_16", string(geometry))
        if !isdir(_filepath_gt)
            @warn "Ground truth directory not found, skipping: $_filepath_gt"
            continue
        end
        dir_list = filter(d -> !startswith(d, ".") && isdir(joinpath(_filepath_gt, d)), readdir(_filepath_gt))
        for dir in dir_list
            if dir in avoid_dirs
                @debug "Skipping dir $dir"
                continue
            end
            filepath_gt = if mode === :conv_exp_mesh
                            joinpath(resolve_data_path("ground_truth/sim_data/Stokes/force/constant/Hex2_16/cylinder/4/"))
                          else    
                            joinpath(_filepath_gt, dir)
                          end
            @debug "Ground truth directory: $filepath_gt"
            for nz in nz_list
                if nz == 6 && mode == :conv_exp_mesh
                    dt_list = [0.1, 0.2, 0.5, 0.7, 0.8, 1.0]
                else
                    dt_list = [0.1]
                end
                if nz == 6 && viscosity_type == "constant"
                    noise_level_list = [0.0, 0.5, 1.0, 1.5, 2.0]
                else
                    noise_level_list = [0.0]
                end
                for noise_level in noise_level_list 
                    sim_time_exp_list = [5.0]  # simulation time in seconds
                    @debug "Simulation time experiments to run: $sim_time_exp_list"
                    for sim_time_exp::Float16 in sim_time_exp_list
                        if viscosity_type == "constant"
                            window = ""
                        end
                        @info "Running optimization with ne = $nz and simulation time = $sim_time_exp with noise level = $noise_level"
                        for cost_function in cost_function_list
                            for dt in dt_list
                                filepath_res = if mode === :conv_exp_mesh
                                                    # `window` must stay last: optimize splices view_$i before the basename.
                                                    joinpath(resolve_data_path(joinpath("experiments","sim_data","convergence_analysis","stokes_convergence")), "experiment_mesh_conv", string(geometry), dir, "$(element_shape_x)$(basis_order_x)_$(nz)", "dt_$(dt)", _run_dir(opt_method, cost_function), window)
                                                else
                                                    joinpath(resolve_data_path(joinpath("experiments","sim_data","optimization","Stokes")), control, viscosity_type, "Hex2_16", string(geometry), dir, "$(element_shape_x)$(basis_order_x)_$(nz)", "simtime_$(sim_time_exp)", "noise_$(noise_level)", "dt_$(dt)", _run_dir(opt_method, cost_function), window)
                                                end
                    
                                                @info "Running optimization with element_shape_x = $element_shape_x, basis_order_x = $basis_order_x, ne = $nz, dt = $dt"

                                exp_params = Dict(
                                    "element_shape_u" => element_shape_u, "basis_order_u" => basis_order_u,
                                    "element_shape_p" => element_shape_p, "basis_order_p" => basis_order_p,
                                    "element_shape_x" => element_shape_x, "basis_order_x" => basis_order_x,
                                    "ne_exp" => nz, "sim_time_exp" => sim_time_exp, "dt" => dt,
                                    "filepath_res" => filepath_res, "filepath_gt" => filepath_gt,
                                    "control" => control, "data_type" => "simulated",
                                    "camera_matrix" => camera_matrix, "WRITE_GT" => false,
                                    "noise_level" => noise_level, "mode" => window,
                                    "opt_method" => opt_method, "opt_kwargs" => opt_kwargs,
                                    "cost_function" => cost_function)

                                push!(param_list, exp_params)
                            end
                        end
                    end
                end
            end
        end
        if use_parallel
            run_parallel_tasks(param_list, optimize; max_workers=15, memory_per_task_mb=8192.0)
        else
            for (i, params) in enumerate(param_list)
                @info "Sequential execution: calling write_gt_data for index $i / $(length(param_list))"
                try
                    optimize(params)
                    @info "Completed write_gt_data for index $i"
                catch err
                    _handle_worker_error(err, i, params)
                end
            end
            @debug "All experiments completed."
        end
    end
end

"""
    optimize_syn(use_parallel=true)

Driver for synthetic (rendered-image) bulk-viscosity experiments: builds one
`optimize` call per combination of mesh element count (`ne_list`),
noise level, simulation time, and time step for `data_type == "synthetic"`,
`viscosity_type == "bulk_viscosity"`, reading ground truth from
`ground_truth/sim_data/Stokes` and writing results under
`experiments/syn_data/optimization/Stokes`. Dispatches the resulting
parameter list to `run_parallel_tasks` when `use_parallel`, otherwise runs
`optimize` sequentially, logging (not raising) any per-task failure via
`_handle_worker_error`.

# Arguments
- `use_parallel::Bool`: dispatch to `run_parallel_tasks` instead of running
  sequentially (default: `true`).

# Returns
None.
"""
function optimize_syn(use_parallel::Bool=true)

    element_shape_u::Symbol = :Hex
    basis_order_u::Int = 2
    element_shape_p::Symbol = :Hex
    basis_order_p::Int = 1
    element_shape_x::Symbol = :Hex
    basis_order_x::Int = 2

    ne_list = [6]
    noise_level_list = [0.0] # No noise for synthetic data
    control = "force" # "force" or "velocity"
    viscosity_type_list = ["bulk_viscosity"]
    window = "multi_window" 
    camera_matrix::AbstractArray = [[2.39642674e+03, 0.0, 1.00429248e+03] [0.0, 2.40565353e+03, 7.57028161e+02] [0.0, 0.0, 1.0]]'
    filepath_res::String = ""
    param_list = Vector{Dict}(undef, 0)

    # Optimization method and its `fit_model` keyword arguments. `opt_method` also
    # becomes a path segment under `dt_*`, so runs with different methods never
    # overwrite each other. Keys in `opt_kwargs` are passed straight through to
    # `fit_model`, so anything it accepts works here (see `_opt_kwargs`).
    cost_function_list = [:closest_point]   # :chamfer or :closest_point
    opt_method::Symbol = :gn                # :gn, :lm or :gn_tikhonov
    opt_kwargs = Dict{String,Any}()         # e.g. "λ_scale" => 0.01, "line_search_method" => :armijo
    geometry_list::Vector{Symbol} = [:cylinder, :cube]
    avoid_dirs = ["post_analysis_global"]
    dt_list = [0.1]

    for geometry::Symbol in geometry_list
        for viscosity_type in viscosity_type_list
            _filepath_gt = joinpath(resolve_data_path("ground_truth/sim_data/Stokes"), control, viscosity_type, "Hex2_16", string(geometry))
            if !isdir(_filepath_gt)
                @warn "Ground truth directory not found, skipping: $_filepath_gt"
                continue
            end
            dir_list = filter(d -> !startswith(d, ".") && isdir(joinpath(_filepath_gt, d)), readdir(_filepath_gt))
            for dir in dir_list
                if dir in avoid_dirs
                    @debug "Skipping dir $dir"
                    continue
                end
                filepath_gt = joinpath(_filepath_gt, dir)
                for ne in ne_list
                    for noise_level in noise_level_list
                        sim_time_exp_list = [5.0]
                        @debug "Simulation time experiments to run: $sim_time_exp_list"
                        for sim_time_exp::Float16 in sim_time_exp_list
                            if viscosity_type == "constant"
                                window = ""
                            end
                            @info "Running optimization with ne = $ne and simulation time = $sim_time_exp with noise level = $noise_level"
                            for cost_function in cost_function_list
                                for dt in dt_list
                                    filepath_res = joinpath(resolve_data_path(joinpath("experiments","syn_data","optimization","Stokes")), control, viscosity_type, "Hex2_16", string(geometry), dir, "$(element_shape_x)$(basis_order_x)_$(ne)", "simtime_$(sim_time_exp)", "noise_$(noise_level)", "dt_$(dt)", _run_dir(opt_method, cost_function), window)
                                    @info "Running optimization with element_shape_x = $element_shape_x, basis_order_x = $basis_order_x with $ne elements"
    
                                    exp_params = Dict(
                                        "element_shape_u" => element_shape_u, "basis_order_u" => basis_order_u,
                                        "element_shape_p" => element_shape_p, "basis_order_p" => basis_order_p,
                                        "element_shape_x" => element_shape_x, "basis_order_x" => basis_order_x,
                                        "ne_exp" => ne, "sim_time_exp" => sim_time_exp,
                                        "dt" => dt,
                                        "filepath_res" => filepath_res, "filepath_gt" => filepath_gt, "control" => control, "data_type" => "synthetic", "camera_matrix" => camera_matrix, "WRITE_GT" => false,
                                        "noise_level"=>noise_level, "mode"=>window,
                                        "opt_method" => opt_method, "opt_kwargs" => opt_kwargs,
                                        "cost_function" => cost_function)
    
                                    push!(param_list, exp_params)
                                end
                            end
                        end
                    end
                end
            end
        end
        if use_parallel
            run_parallel_tasks(param_list, optimize; max_workers=15, memory_per_task_mb=8192.0)
        else
            for (i, params) in enumerate(param_list)
                @info "Sequential execution: calling write_gt_data for index $i / $(length(param_list))"
                try
                    optimize(params)
                    @info "Completed write_gt_data for index $i"
                catch err
                    _handle_worker_error(err, i, params)
                end
            end
            @debug "All experiments completed."
        end
    end
end

"""
    optimize_real(use_parallel=true)

Driver for physically-measured (or Carreau-model synthetic) bulk-viscosity
experiments: builds one `optimize` call per ground-truth directory
and mesh element count (`ne_list`). When `model_type == "carreau"`, reads
ground truth from `ground_truth/sim_data/Carreau` and sets `data_type =
"synthetic"`; otherwise reads from `ground_truth/physical_data` and sets
`data_type = "physical"` (supplying the applied force `F_ext`, cylinder
radius `r`, and height `h`). Dispatches the resulting parameter list to
`run_parallel_tasks` when `use_parallel`, otherwise runs `optimize`
sequentially, logging (not raising) any per-task failure via
`_handle_worker_error`.

# Arguments
- `use_parallel::Bool`: dispatch to `run_parallel_tasks` instead of running
  sequentially (default: `true`).

# Returns
None.
"""
function optimize_real(use_parallel::Bool=true)
    F_ext::Float64 = 9.812*1e3 # force applied to the cylinder in N
    sim_time_exp_list::Vector{Float64} = [0.5] # simulation time in seconds

    element_shape_u::Symbol = :Hex
    basis_order_u::Int = 2
    element_shape_p::Symbol = :Hex
    basis_order_p::Int = 1
    element_shape_x::Symbol = :Hex
    basis_order_x::Int = 2

    ne_list = [6]
    control = "force" # "force" or "velocity"
    η_start = 50.0
    β_start = 50.0
    viscosity_type = "bulk_viscosity"

    r::Float64 = 25.0  # radius of the cylinder in mm
    h::Float64 = 40.0  # height of the cylinder in mm
    camera_matrix::AbstractArray = [[2.39642674e+03, 0.0, 1.00429248e+03] [0.0, 2.40565353e+03, 7.57028161e+02] [0.0, 0.0, 1.0]]'

    model_type = "carreau" # "carreau" or "Stokes"
    window = "multi_window" # "multi_window" or "single_window"
    filepath_res::String = ""
    param_list = Vector{Dict}(undef, 0)
    avoid_dirs = ["post_analysis_global","1","2","3","5","4","7","8","9"] 

    # Optimization method and its `fit_model` keyword arguments. `opt_method` also
    # becomes a path segment under `dt_*`, so runs with different methods never
    # overwrite each other. Keys in `opt_kwargs` are passed straight through to
    # `fit_model`, so anything it accepts works here (see `_opt_kwargs`).
    cost_function_list = [:closest_point]   # :chamfer or :closest_point
    opt_method::Symbol = :gn                # :gn, :lm or :gn_tikhonov
    opt_kwargs = Dict{String,Any}()         # e.g. "λ_scale" => 0.01, "line_search_method" => :armijo

    dt_carreau::Float64 = 0.1 # must match `time_steps` in the Carreau ground-truth sim_params
    dt_physical::Float64 = 0.1 # must match `time_steps` in the physical ground-truth sim_params

    if model_type == "carreau" && viscosity_type == "bulk_viscosity"
        _filepath_gt = joinpath(resolve_data_path("ground_truth/sim_data/Carreau"))
        sim_time_exp_list = [5.0]
        F_ext = 1e4 # force applied to the cylinder in N
    else
        _filepath_gt = joinpath(resolve_data_path("ground_truth/physical_data"))
        sim_time_exp_list = [0.5] # simulation time in seconds
        F_ext = 9.812*1e3 # force applied to the cylinder in N
        h = 38.5  # height of the cylinder in mm
    end

    # Extract base paths to variables using resolve_data_path for configuration-aware resolution
    base_experiments_path = resolve_data_path("experiments")
    base_gt_path = resolve_data_path("ground_truth")

    if !isdir(_filepath_gt)
        @warn "Ground truth directory not found, skipping optimize_real: $_filepath_gt"
        return
    end
    dir_list = filter(d -> !startswith(d, ".") && isdir(joinpath(_filepath_gt, d)), readdir(_filepath_gt))
    for dir in dir_list
        if dir in avoid_dirs
            @debug "Skipping dir $dir"
            continue
        end
        @info "Processing ground truth directory: $dir for $viscosity_type viscosity ..."
        filepath_gt = joinpath(_filepath_gt, dir)
        for ne in ne_list
            @debug "Simulation time experiments to run: $sim_time_exp_list"
            for sim_time_exp in sim_time_exp_list
                @debug "Simulation time: $sim_time_exp seconds"
                for cost_function in cost_function_list
                    for ne in ne_list
                        if model_type == "carreau" && viscosity_type == "bulk_viscosity"
                            filepath_res = joinpath(base_experiments_path, "syn_data", "optimization", "Carreau", dir, "$(element_shape_x)$(basis_order_x)_$(ne)", "simtime_$(sim_time_exp)", "noise_0.0", "dt_$(dt_carreau)", _run_dir(opt_method, cost_function), window)
                            exp_params = Dict(
                                "element_shape_u" => element_shape_u, "basis_order_u" => basis_order_u,
                                "element_shape_p" => element_shape_p, "basis_order_p" => basis_order_p,
                                "element_shape_x" => element_shape_x, "basis_order_x" => basis_order_x,
                                "ne_exp" => ne, "sim_time_exp" => sim_time_exp,
                                "dt" => dt_carreau,
                                "η_start" => η_start, "β_start" => β_start,
                                "filepath_res" => filepath_res, "filepath_gt" => filepath_gt,
                                "control" => control, "data_type" => "synthetic",
                                "camera_matrix" => camera_matrix, "WRITE_GT" => false, "noise_level" => 0.0, "mode" => window,
                                "opt_method" => opt_method, "opt_kwargs" => opt_kwargs,
                                "cost_function" => cost_function)
                        else
                            filepath_res = joinpath(base_experiments_path, "physical_data", "optimization", dir, "$(element_shape_x)$(basis_order_x)_$(ne)", "simtime_$(sim_time_exp)", "noise_0.0", "dt_$(dt_physical)", _run_dir(opt_method, cost_function), window)
                            exp_params = Dict(
                                "element_shape_u" => element_shape_u, "basis_order_u" => basis_order_u,
                                "element_shape_p" => element_shape_p, "basis_order_p" => basis_order_p,
                                "element_shape_x" => element_shape_x, "basis_order_x" => basis_order_x,
                                "ne_exp" => ne, "sim_time_exp" => sim_time_exp,
                                "dt" => dt_physical,
                                "η_start" => η_start, "β_start" => β_start,
                                "filepath_res" => filepath_res, "filepath_gt" => filepath_gt,
                                "control" => control, "viscosity_type" => viscosity_type,
                                "data_type" => "physical", "r" => r, "h" => h,
                                "camera_matrix" => camera_matrix, "F_ext" => F_ext, "mode" => window,
                                "opt_method" => opt_method, "opt_kwargs" => opt_kwargs,
                                "cost_function" => cost_function)
                        end
                        @info "Running optimization with element_shape_x = $element_shape_x, basis_order_x = $basis_order_x with $ne elements"
    
                        push!(param_list, exp_params)
                    end
                end
            end
        end
    end
    if use_parallel
        run_parallel_tasks(param_list, optimize; max_workers=5, memory_per_task_mb=8192.0)
    else
        for (i, params) in enumerate(param_list)
            @info "Sequential execution: calling write_gt_data for index $i / $(length(param_list))"
            try
                optimize(params)
                @info "Completed write_gt_data for index $i"
            catch err
                _handle_worker_error(err, i, params)
            end
        end
        @debug "All experiments completed."
    end
end

"""
    plot_results()

Top-level driver that regenerates plots for every result tree under
`experiments/` in the configured data directory (see `resolve_data_path`),
across the configured `data_type_list` (`"simulated"`/`"synthetic"`/`"physical"`)
and their applicable `model_type` (`"Stokes"`/`"carreau"`) and `viscosity_type`
combinations. For each ground-truth directory (skipping `avoid_dirs`), applies
`set_plot_config` for that data/viscosity type and calls `replot` on the matching
experiment directory, then runs the `post_analysis_*` cross-experiment summary
matching the viscosity/data type.

# Returns
None. All outputs are the side effects of `replot`.
"""
function plot_results()
    control::String = "force"
    viscosity_type_list = [] # "constant" or "bulk_viscosity"
    model_type = [] # "carreau" or "Stokes"
    avoid_dirs = ["post_analysis_global"] # directories to skip in post-analysis and plotting
    data_type_list = ["simulated","synthetic","physical"] # ["simulated", "synthetic", "physical"]
    base_path = ""
    geometry::Symbol = :cylinder # :cylinder or :cube

    for data_type in data_type_list
        avoid_dirs = data_type == "physical" ? ["post_analysis_global"] : ["post_analysis_global", "5"]
        if data_type == "physical"
            model_type = ["Stokes"] # for physical data, we only have Stokes model results for now
        elseif data_type == "synthetic"
            model_type = ["carreau", "Stokes"]
        elseif data_type == "simulated"
            model_type = ["Stokes"]
        end

        base_experiments_path = resolve_data_path("experiments")
        base_gt_path = resolve_data_path("ground_truth")
        for model_type::String in model_type
            if data_type == "synthetic" 
                if model_type == "carreau"
                    base_path = joinpath(base_experiments_path, "syn_data", "optimization", "Carreau")
                    viscosity_type_list = ["bulk_viscosity"] # for Carreau model, we only have bulk viscosity results for now
                else
                    base_path = joinpath(base_experiments_path, "syn_data", "optimization", "Stokes")
                    viscosity_type_list = ["bulk_viscosity"]
                end
            elseif data_type == "simulated"
                base_path = joinpath(base_experiments_path, "sim_data", "optimization", "Stokes")
                viscosity_type_list = ["constant"]
            elseif data_type == "physical"
                base_path = joinpath(base_experiments_path, "physical_data", "optimization")
                viscosity_type_list = ["bulk_viscosity"]
            end
            for viscosity_type::String in viscosity_type_list
                set_plot_config(data_type, viscosity_type)
                if data_type == "physical"
                    filepath_gt = joinpath(base_gt_path, "physical_data")
                    filepath_res = base_path
                else
                    if model_type == "carreau" && viscosity_type == "bulk_viscosity"
                        filepath_gt = joinpath(base_gt_path, "sim_data", "Carreau")
                        filepath_res = base_path
                    else
                        filepath_gt = joinpath(base_gt_path, "sim_data", "Stokes", control, viscosity_type, "Hex2_16" , string(geometry))
                        filepath_res = joinpath(base_path, control, viscosity_type, "Hex2_16" , string(geometry))
                    end
                end

                if !isdir(filepath_gt)
                    @warn "Ground truth directory not found, skipping: $filepath_gt"
                    continue
                end
                dirs = filter(d -> !startswith(d, ".") && isdir(joinpath(filepath_gt, d)), readdir(filepath_gt))

                # One pass per optimizer, so `gn` and `gn_tikhonov` results are never mixed
                # into the same figures or overwritten by each other. Trees written before
                # the method level existed have `view_*` straight under `dt_*`; those yield
                # no method folder and are processed once with `nothing`, keeping the flat
                # output layout they already have.
                for opt_method in _result_methods(filepath_res, dirs, avoid_dirs)
                    @info "Post-processing" data_type model_type viscosity_type method=something(opt_method, "(none)")
                    for dir::String in dirs
                        if dir in avoid_dirs
                            @debug "Skipping dir $dir"
                            continue
                        end
                        filepath_gt_dir = joinpath(filepath_gt, dir)
                        filepath_res_dir = joinpath(filepath_res, dir)
                        @debug "Processing ground truth directory: $filepath_gt_dir for $viscosity_type viscosity ..."
                        replot(filepath_res_dir, filepath_gt_dir; method=opt_method)
                    end
                    if viscosity_type == "constant"
                        post_analysis_const(filepath_gt, filepath_res, avoid_dirs; method=opt_method)
                    elseif viscosity_type == "bulk_viscosity" && data_type != "physical"
                        post_analysis_bulk(filepath_gt, filepath_res, avoid_dirs; method=opt_method)
                    elseif data_type == "physical"
                        post_analysis_real(filepath_gt, filepath_res, avoid_dirs; method=opt_method)
                    end
                end
            end
        end
    end
end

"""
    main()

Run the full pipeline: the simulated, synthetic, and physical optimization
batches followed by result plotting. Invoked automatically when this file is
run as a script (`julia test_opt_stokes.jl`), but not when it is `include`d
from the REPL or another script.
"""
function main()
    optimize_sim(false)
    # optimize_syn(false)
    # optimize_real(false)
    # plot_results()
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

