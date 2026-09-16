# Figures and numbers for the paper, in two parts:
#
#   1. Contour-extraction (segmentation) error, constant-viscosity experiments.
#   2. The Tikhonov λ sweep on the synthetic bulk-viscosity experiments.
#
# Part 2 fits models, so this file includes the pipeline rather than duplicating it; the
# colour constants and the `datapath`/`plotpath` helpers come from there too.
#
# --- part 1 ---------------------------------------------------------------------------
# Each constant-viscosity ground-truth run carries two contour sets per frame:
#   data/sim_data/<view>/contour_data — the FEM surface projected to the image plane
#   data/img_data/<view>/contour_data — segmented from the rendered image
# The projected contour is exact by construction, so the distance between the two is
# the error the segmentation injects before the optimizer ever sees the data. Readers
# and metrics are the ones `test_opt_stokes.jl` uses (`read_csv`, `compare_pt_clouds`),
# so these numbers sit on the same px scale as the contour costs reported there.

using smearFEM
using NearestNeighbors
using StatsPlots
using Statistics
using Printf
using DelimitedFiles

# Brings `optimize` (part 2 fits models), the plotting helpers `_fig`/`_label!`, the
# `datapath`/`plotpath` path helpers, the `def_*` colours, `_align_windowed` and `_run_dir`.
# `main()` there is guarded by a PROGRAM_FILE check, so including it runs nothing.
include(joinpath(@__DIR__, "model_optimization", "test_opt_stokes.jl"))

const GT_CONST  = "ground_truth/sim_data/Stokes/force/constant/Hex2_16/cylinder"
const RES_CONST = "experiments/sim_data/optimization/Stokes/force/constant/Hex2_16/cylinder"
const GT_BULK   = "ground_truth/sim_data/Stokes/force/bulk_viscosity/Hex2_16/cylinder"
const RES_BULK  = "experiments/sim_data/optimization/Stokes/force/bulk_viscosity/Hex2_16/cylinder"

"""
    _frame_window(gt_pts, seg_pts, skip_frames, max_frames, filepath_gt) -> UnitRange

The frames both contour sets share, after dropping `skip_frames` leading ones and keeping at
most `max_frames` of the rest. Warns when a run is too short to fill the window, since that
run then covers less time than the others it is averaged with.
"""
function _frame_window(gt_pts, seg_pts, skip_frames::Int, max_frames::Union{Nothing,Int},
                       filepath_gt::String)
    nf = min(length(gt_pts), length(seg_pts))
    last = isnothing(max_frames) ? nf : min(nf, skip_frames + max_frames)
    frames = (1 + skip_frames):last
    isempty(frames) && throw(ArgumentError(
        "No frames left in $filepath_gt: $nf shared frames, skip_frames=$skip_frames"))
    if !isnothing(max_frames) && length(frames) < max_frames
        @warn "Only $(length(frames)) of the $max_frames requested frames in $filepath_gt " *
              "(projected $(length(gt_pts)), segmented $(length(seg_pts)))"
    end
    return frames
end

"""
    _segmented_runs(filepath_gt, view_folder, avoid_dirs) -> Vector{String}

Run folders under `filepath_gt` that carry a segmented contour for `view_folder`.
"""
function _segmented_runs(filepath_gt::String, view_folder::String, avoid_dirs)
    runs = filter(d -> !startswith(d, ".") && !(d in avoid_dirs) &&
                       isdir(datapath(joinpath(filepath_gt, d), "img_data", view_folder, "contour_data")),
                  readdir(filepath_gt))
    isempty(runs) && throw(ArgumentError("No runs with segmented contours under $filepath_gt"))
    return runs
end

"""
    contour_extraction_error(filepath_gt; view_folder="view_1", squared_chamfer=false, skip_frames=0, max_frames=nothing) -> Dict

Per-frame distance between the segmented contour (`img_data`) and the projected
ground-truth contour (`sim_data`) of one experiment.

Frames are index-aligned — segmented frame `k` is the render of simulation step `k` —
and truncated to the shorter of the two sequences, since a render can stop before the
simulation does.

# Arguments
- `filepath_gt::String`: ground-truth run directory (the one holding `data/`).

# Keyword Arguments
- `view_folder::String`: camera-view subdirectory (default: `"view_1"`).
- `squared_chamfer::Bool`: report the squared Chamfer distance, comparable to the
  optimizer's `ChamferCost`, instead of the px-unit form (default: `false`).
- `skip_frames::Int`: drop this many leading frames (default: `0`). Frame 1 is an
  outlier in the constant-viscosity data — see [`contour_extraction_error_const`](@ref).
- `max_frames::Union{Nothing,Int}`: keep at most this many frames after the skip
  (default: `nothing`, all of them).

# Returns
- `Dict` with `"time"` (s), and the per-frame series `"hausdorff"`, `"chamfer"`,
  `"closest_pt"` and `"directed_hausdorff"` (segmented → projected), in px — px² for the
  Chamfer entry when `squared_chamfer`.
"""
function contour_extraction_error(filepath_gt::String; view_folder::String="view_1",
                                  squared_chamfer::Bool=false, skip_frames::Int=0,
                                  max_frames::Union{Nothing,Int}=nothing)

    gt_pts, _, _  = read_csv(datapath(filepath_gt, "sim_data", view_folder, "contour_data"))
    seg_pts, _, _ = read_csv(datapath(filepath_gt, "img_data", view_folder, "contour_data"))

    frames = _frame_window(gt_pts, seg_pts, skip_frames, max_frames, filepath_gt)
    h, c, cp, dh = compare_pt_clouds(seg_pts[frames], gt_pts[frames]; squared_chamfer=squared_chamfer)
    frame_rate = read_json(datapath(filepath_gt, "video_metadata"))["frame_rate"]

    return Dict{String,Any}("time" => collect(frames .- 1) ./ frame_rate,
                            "hausdorff" => h, "chamfer" => c, "closest_pt" => cp,
                            "directed_hausdorff" => dh)
end

"""
    contour_extraction_error_const(filepath_gt; kwargs...) -> Dict{String,Dict}

Run [`contour_extraction_error`](@ref) over every constant-viscosity replicate, print a
per-run summary table, and write the series and figures under `post_analysis_global`.

# Arguments
- `filepath_gt::String`: directory holding the ground-truth run folders
  (default: the constant-viscosity cylinder tree).

# Keyword Arguments
- `view_folder::String`: camera-view subdirectory (default: `"view_1"`).
- `filepath_res::String`: results root the outputs are written under
  (default: the matching `experiments/sim_data` tree).
- `avoid_dirs`: run folders to skip (default: `["post_analysis_global"]`).
- `max_frames::Int`: frames to keep after the skip (default: `200`, i.e. the 20 s
  observation window at the 10 Hz frame rate; every run is long enough to fill it).
- `skip_frames::Int`: leading frames to drop (default: `1`). Frame 1 (t = 0, the
  undeformed cylinder) is an outlier in every replicate: its projected contour carries
  arcs that run *inside* the silhouette the segmentation sees, so the Hausdorff distance
  is ~340 px there against ~10 px for every later frame — a defect of the projection at
  t = 0, not of the segmentation. `replot` likewise scores contours from the second frame
  on. Set `0` to keep it.
- `write_results::Bool`, `plot_results::Bool`: write CSVs / figures (default: `true`).

# Returns
- `Dict` mapping run name to that run's `contour_extraction_error` result.
"""
function contour_extraction_error_const(filepath_gt::String=resolve_data_path(GT_CONST);
                                        view_folder::String="view_1",
                                        filepath_res::String=resolve_data_path(RES_CONST),
                                        avoid_dirs=["post_analysis_global"], skip_frames::Int=1,
                                        max_frames::Int=200,
                                        write_results::Bool=true, plot_results::Bool=true)

    runs = _segmented_runs(filepath_gt, view_folder, avoid_dirs)

    errors = Dict{String,Dict}()
    # Means are over frames; `hausdorff_max` and the time it occurs at expose any single
    # frame the segmentation lost, which the means alone would hide.
    @printf("%-6s %7s %11s %11s %11s %11s %11s %9s\n",
            "run", "frames", "closest_pt", "chamfer", "dir_hausd", "hausdorff", "hausd_max",
            "t@max[s]")
    for run in runs
        e = contour_extraction_error(joinpath(filepath_gt, run); view_folder=view_folder,
                                     skip_frames=skip_frames, max_frames=max_frames)
        errors[run] = e
        @printf("%-6s %7d %11.3f %11.3f %11.3f %11.3f %11.3f %9.1f\n", run, length(e["time"]),
                mean(e["closest_pt"]), mean(e["chamfer"]), mean(e["directed_hausdorff"]),
                mean(e["hausdorff"]), maximum(e["hausdorff"]),
                e["time"][argmax(e["hausdorff"])])
    end

    pooled = Dict(k => reduce(vcat, errors[r][k] for r in runs)
                  for k in ("closest_pt", "chamfer", "hausdorff", "directed_hausdorff"))
    @printf("%-6s %7d %11.3f %11.3f %11.3f %11.3f %11.3f\n", "all", length(pooled["hausdorff"]),
            mean(pooled["closest_pt"]), mean(pooled["chamfer"]),
            mean(pooled["directed_hausdorff"]),
            mean(pooled["hausdorff"]), maximum(pooled["hausdorff"]))
    tspan = extrema(errors[runs[1]]["time"])
    @info "Segmentation error over $(length(runs)) runs, t = $(tspan[1])–$(tspan[2]) s [px]: " *
          "closest-point $(round(mean(pooled["closest_pt"]), sigdigits=4)), " *
          "Hausdorff $(round(mean(pooled["hausdorff"]), sigdigits=4))"

    out = joinpath(filepath_res, "post_analysis_global")
    set_file(out)
    if write_results
        for run in runs, k in ("time", "closest_pt", "chamfer", "hausdorff", "directed_hausdorff")
            write_csv(datapath(out, "contour_extraction_error", run, k), errors[run][k])
        end
        # One row per run: run, mean closest-point, mean chamfer, mean and max Hausdorff.
        write_csv(datapath(out, "contour_extraction_error", "summary"),
                  reduce(vcat, [r mean(errors[r]["closest_pt"]) mean(errors[r]["chamfer"]) mean(errors[r]["hausdorff"]) maximum(errors[r]["hausdorff"])] for r in runs))
    end

    if plot_results
        # All four are px here (the Chamfer is requested unsquared), so one figure can carry
        # them honestly: the run-averaged series of each, on a shared axis.
        cmp_plt = set_plot(12; sz=(480, 320), legend_column=2)
        Plots.plot!(cmp_plt, [], label=false)
        tgrid = errors[runs[1]]["time"]
        # The two Hausdorff curves nearly coincide, so the directed one is dashed: colour
        # alone cannot separate lines that sit within a few percent of each other.
        for (k, lbl, c, ls) in (("closest_pt", "closest point (RMSE)", def_blue, :solid),
                                ("chamfer", "Chamfer", def_green, :solid),
                                ("directed_hausdorff", "directed Hausdorff", def_orange, :dash),
                                ("hausdorff", "Hausdorff, symmetric", def_red, :solid))
            n = minimum(length(errors[r][k]) for r in runs)
            Plots.plot!(cmp_plt, tgrid[1:n],
                        mean(hcat((errors[r][k][1:n] for r in runs)...); dims=2),
                        label=lbl, color=c, linestyle=ls)
        end
        Plots.xlabel!(cmp_plt, L"\mathrm{Time\;[s]}")
        Plots.ylabel!(cmp_plt, L"\mathrm{Distance\;[px]}")
        Plots.ylims!(cmp_plt, 0, 1.1 * maximum(maximum(errors[r]["hausdorff"]) for r in runs))
        Plots.savefig(cmp_plt, plotpath(out, "contour_extraction_metrics.pdf"))

        for (k, ylabel) in (("closest_pt", L"\mathrm{Closest\;point\;distance\;[px]}"),
                            ("directed_hausdorff", L"\mathrm{Directed\;Hausdorff\;[px]}"),
                            ("hausdorff",  L"\mathrm{Hausdorff\;distance\;[px]}"))
            plt = set_plot(12; sz=(480, 320), legend_column=3)
            Plots.plot!(plt, [], label=false)
            for run in runs
                Plots.plot!(plt, errors[run]["time"], errors[run][k], label="run $run")
            end
            Plots.xlabel!(plt, L"\mathrm{Time\;[s]}")
            Plots.ylabel!(plt, ylabel)
            # Anchor at zero so the curves read as absolute px error, not as fluctuation.
            Plots.ylims!(plt, 0, 1.1 * maximum(maximum(errors[r][k]) for r in runs))
            Plots.savefig(plt, plotpath(out, "contour_extraction_$(k).pdf"))
        end
        @info "Wrote figures to $(plotpath(out))"
    end

    return errors
end


# ---------------------------------------------------------------------------------------
# Signed bias: does the segmentation sit systematically outside the true silhouette?
#
# The distances above are magnitudes, so they cannot separate a mask that is uniformly too
# large from one that is noisy but centred — only the first biases a parameter fit. Every
# quantity below is therefore signed, outward-positive, and each is aggregated to one number
# per run before any test: frames within a run are strongly autocorrelated, so the runs, not
# the frames within them, are the independent units.
# ---------------------------------------------------------------------------------------

# read_csv hands back points as (2, n); the geometry below wants one point per row.
_rows(p::AbstractArray) = size(p, 1) < size(p, 2) ? permutedims(p) : p

# Shoelace area of a closed polygon given by its ordered vertices.
function _polygon_area(p::AbstractMatrix)
    n = size(p, 1)
    return abs(sum(p[i,1]*p[mod1(i+1,n),2] - p[mod1(i+1,n),1]*p[i,2] for i in 1:n)) / 2
end

# Which side of the closed curve a point falls on, from the cross product of the local
# tangent at vertex `i` with the offset to the point. The sign on its own is meaningless —
# it flips with the traversal direction of the polygon — so callers calibrate it against a
# point known to lie inside.
function _side(p::AbstractMatrix, i::Int, qx::Real, qy::Real)
    n = size(p, 1)
    tx = p[mod1(i+1,n),1] - p[mod1(i-1,n),1]
    ty = p[mod1(i+1,n),2] - p[mod1(i-1,n),2]
    return sign(tx * (qy - p[i,2]) - ty * (qx - p[i,1]))
end

"""
    signed_offsets(q_pts, poly) -> Vector{Float64}

Signed nearest-neighbour distance from each point of `q_pts` to the closed curve `poly`,
positive outside it.

The sign convention is fixed per call by testing the polygon's centroid, which lies inside
these silhouettes (a squeezed cylinder projects to a barrel, never to a shape that excludes
its own centroid), rather than by assuming a traversal direction the writers do not
guarantee.

# Arguments
- `q_pts::AbstractMatrix`, `poly::AbstractMatrix`: query points and the closed curve, one
  point per row.

# Returns
- `Vector{Float64}`: one signed distance per query point, in the coordinate units (px).
"""
function signed_offsets(q_pts::AbstractMatrix, poly::AbstractMatrix)
    tree = KDTree(permutedims(poly))
    idx, d = nn(tree, permutedims(q_pts))
    c = vec(mean(poly; dims=1))
    inside = _side(poly, nn(tree, c)[1], c[1], c[2])
    return [_side(poly, idx[k], q_pts[k,1], q_pts[k,2]) == inside ? -d[k] : d[k]
            for k in eachindex(d)]
end

# The signed quantities, grouped for reporting, each with the value it takes when the
# segmentation is unbiased. `registration` is a reparametrization of `edges`: two opposite
# edges moving together are the mask sitting in the wrong place, moving apart are the mask
# being the wrong size, and only the second is something the segmentation itself did.
const BIAS_GROUPS = (
    "silhouette"   => (("normal_offset_px", 0.0), ("area_ratio", 1.0),
                       ("width_ratio", 1.0), ("height_ratio", 1.0)),
    "registration" => (("shift_x_px", 0.0), ("shift_y_px", 0.0),
                       ("dilation_x_px", 0.0), ("dilation_y_px", 0.0)),
    "edges"        => (("left_edge_px", 0.0), ("right_edge_px", 0.0),
                       ("top_edge_px", 0.0), ("bottom_edge_px", 0.0)),
)
const BIAS_QUANTITIES = reduce(vcat, [collect(g) for (_, g) in BIAS_GROUPS])

# What the paper reports for segmentation error: closest-point RMSE (magnitude) plus area
# ratio (larger/smaller), reusing quantities `contour_bias` already computes rather than
# introducing the signed-offset/registration decomposition as a separate metric. That
# decomposition stays available via `contour_bias` directly for diagnosis.
const PAPER_BIAS_QUANTITIES = (("closest_pt_px", 0.0), ("area_ratio", 1.0))

"""
    contour_bias(filepath_gt; view_folder="view_1", skip_frames=0, max_frames=nothing) -> Dict

Per-frame signed discrepancies between the segmented contour and the projected one, all
positive when the segmentation lies outside the truth.

`normal_offset_px` is the mean signed distance of the segmented points to the projected
curve — how much the mask is uniformly too large. The extent terms decompose that per image
edge, and `shift_*`/`dilation_*` recombine the opposite edges into the two things that
cannot be told apart from a single edge: a mask in the wrong place and a mask of the wrong
size.

# Arguments
- `filepath_gt::String`: ground-truth run directory (the one holding `data/`).

# Keyword Arguments
- `view_folder::String`: camera-view subdirectory (default: `"view_1"`).
- `skip_frames::Int`: drop this many leading frames (default: `0`); see
  [`contour_extraction_error_const`](@ref) for why frame 1 is dropped.
- `max_frames::Union{Nothing,Int}`: keep at most this many frames after the skip
  (default: `nothing`, all of them).

# Returns
- `Dict` with `"time"` (s), one per-frame series for every quantity in `BIAS_GROUPS`, and
  `"offset_rms_px"`, the RMS of the same signed offsets whose mean is `normal_offset_px`.
"""
function contour_bias(filepath_gt::String; view_folder::String="view_1", skip_frames::Int=0,
                      max_frames::Union{Nothing,Int}=nothing,
                      frame_rate::Union{Nothing,Real}=nothing)

    gt_pts, _, _  = read_csv(datapath(filepath_gt, "sim_data", view_folder, "contour_data"))
    seg_pts, _, _ = read_csv(datapath(filepath_gt, "img_data", view_folder, "contour_data"))

    frames = _frame_window(gt_pts, seg_pts, skip_frames, max_frames, filepath_gt)

    b = Dict(k => Float64[] for (k, _) in BIAS_QUANTITIES)
    # Not in BIAS_GROUPS: a spread is non-negative, so testing its mean against zero says
    # nothing. Kept because it is what the mean offset has to be read against.
    b["offset_rms_px"] = Float64[]
    for f in frames
        gt, sg = _rows(gt_pts[f]), _rows(seg_pts[f])
        gx, gy = extrema(gt[:,1]), extrema(gt[:,2])
        sx, sy = extrema(sg[:,1]), extrema(sg[:,2])

        offs = signed_offsets(sg, gt)
        push!(b["normal_offset_px"], mean(offs))
        push!(b["offset_rms_px"], sqrt(mean(offs .^ 2)))
        push!(b["area_ratio"],   _polygon_area(sg) / _polygon_area(gt))
        push!(b["width_ratio"],  (sx[2] - sx[1]) / (gx[2] - gx[1]))
        push!(b["height_ratio"], (sy[2] - sy[1]) / (gy[2] - gy[1]))
        # Outward-positive: y runs down the image, so the top edge is the smaller y.
        push!(b["left_edge_px"],   gx[1] - sx[1])
        push!(b["right_edge_px"],  sx[2] - gx[2])
        push!(b["top_edge_px"],    gy[1] - sy[1])
        push!(b["bottom_edge_px"], sy[2] - gy[2])
        # Common mode of the two opposite edges (the mask sitting off-centre, +x right and
        # +y down) and their differential (the mask being too large, outward-positive).
        push!(b["shift_x_px"],    ((sx[1] - gx[1]) + (sx[2] - gx[2])) / 2)
        push!(b["shift_y_px"],    ((sy[1] - gy[1]) + (sy[2] - gy[2])) / 2)
        push!(b["dilation_x_px"], ((gx[1] - sx[1]) + (sx[2] - gx[2])) / 2)
        push!(b["dilation_y_px"], ((gy[1] - sy[1]) + (sy[2] - gy[2])) / 2)
    end

    metadata = datapath(filepath_gt, "video_metadata")
    has_metadata = isfile(metadata) || isfile("$metadata.json")
    rate = if has_metadata
        read_json(metadata)["frame_rate"]
    elseif isnothing(frame_rate)
        error("No video_metadata in $filepath_gt; provide frame_rate")
    else
        frame_rate
    end
    b["time"] = collect(frames .- 1) ./ rate
    return b
end

# Mean ± SD across runs for one quantity, on the frames all runs share.
function _ribbon!(plt, series::Vector{Vector{Float64}}, tgrid; color, label)
    k = min(length.(series)...)
    k = min(k, length(tgrid))
    values = hcat((v[1:k] for v in series)...)
    StatsPlots.errorline!(plt, tgrid[1:k], values; centertype=:mean, errortype=:std,
                          errorstyle=:ribbon, groupcolor=color, fillalpha=0.30, label=label)
end

"""
    contour_bias_case(filepath_gt; kwargs...) -> (bias, stats)

Test the segmentation error across one viscosity case.

Reports `closest_pt_px` (the seg→gt nearest-neighbour RMSE — the same quantity
[`contour_extraction_error`](@ref) calls `"closest_pt"`) and `area_ratio`, each reduced to
one mean per run and tested against its unbiased value with a two-sided one-sample t-test
over the runs ([`normalized_replicate_stats`](@ref)). The finer signed-offset/registration
decomposition in [`contour_bias`](@ref) is still computed (and written to CSV per run) but
not part of this summary — see [`contour_bias`](@ref) directly to inspect it.

# Arguments / Keyword Arguments
As [`contour_extraction_error_const`](@ref).

# Returns
- `bias::Dict`: run name → that run's [`contour_bias`](@ref) result.
- `stats::Dict`: quantity name → its [`normalized_replicate_stats`](@ref) entry.
"""
function contour_bias_case(filepath_gt::String;
                            view_folder::String="view_1",
                            filepath_res::String=resolve_data_path(RES_CONST),
                            avoid_dirs=["post_analysis_global"], skip_frames::Int=1,
                            max_frames::Int=200,
                            write_results::Bool=true, plot_results::Bool=true,
                            case_label::String="constant viscosity",
                            frame_rate::Union{Nothing,Real}=nothing)

    runs = _segmented_runs(filepath_gt, view_folder, avoid_dirs)
    bias = Dict(r => contour_bias(joinpath(filepath_gt, r); view_folder=view_folder,
                                  skip_frames=skip_frames, max_frames=max_frames,
                                  frame_rate=frame_rate) for r in runs)

    per_run = Dict("closest_pt_px" => [bias[r]["offset_rms_px"] for r in runs],
                   "area_ratio"    => [bias[r]["area_ratio"] for r in runs])
    stats = Dict(k => normalized_replicate_stats(per_run[k], k; reference=ref)
                 for (k, ref) in PAPER_BIAS_QUANTITIES)
    print_replicate_table([stats[k] for (k, _) in PAPER_BIAS_QUANTITIES],
                          "Segmentation error — closest point & area, $case_label ($view_folder)")

    # What the bias is worth: within a frame the signed offsets have RMS² = mean² + SD², so
    # the squared residual the optimizer minimizes splits into a part the systematic offset
    # accounts for and a part it does not. Both are pooled per run first, then over runs.
    rms  = [sqrt(mean(bias[r]["offset_rms_px"] .^ 2)) for r in runs]
    syst = [abs(mean(bias[r]["normal_offset_px"])) for r in runs]
    rand_ = sqrt.(max.(rms .^ 2 .- syst .^ 2, 0.0))
    @printf("\nResidual split over %d runs [px]: total RMS %.3f = systematic %.3f + random %.3f (in quadrature)\n",
            length(runs), mean(rms), mean(syst), mean(rand_))
    @printf("  the systematic offset accounts for %.1f%% of the mean squared residual\n",
            100 * mean(syst .^ 2 ./ rms .^ 2))

    out = joinpath(filepath_res, "post_analysis_global")
    if write_results
        for r in runs, k in keys(bias[r])
            write_csv(datapath(out, "contour_bias", r, k), bias[r][k])
        end
        set_file(datapath(out, "contour_bias"))
        ordered = [stats[k] for (k, _) in PAPER_BIAS_QUANTITIES]
        write_replicate_stats(datapath(out, "contour_bias", "bias_stats.csv"), ordered)
        # Per-run means alongside the summary: at n = 6 every observation has to be visible.
        write_csv(datapath(out, "contour_bias", "per_run_means"),
                  reduce(vcat, [st.label reshape(st.per_run_mean, 1, :)] for st in ordered))
    end

    if plot_results
        set_file(plotpath(out))
        tgrid = bias[runs[1]]["time"]

        plt = set_plot(12; sz=(480, 320))
        for r in runs
            Plots.plot!(plt, bias[r]["time"], bias[r]["offset_rms_px"], lw=0.45, label=false,
                        color=:gray, linealpha=0.45, linestyle=:solid)
        end
        _ribbon!(plt, [bias[r]["offset_rms_px"] for r in runs], tgrid; color=def_blue,
                 label=L"\mathrm{mean}\;\pm\;1\,\mathrm{SD}")
        Plots.xlabel!(plt, L"\mathrm{Time\;[s]}")
        Plots.ylabel!(plt, L"\mathrm{Closest\;Point\;Distance\;[px]}")
        Plots.savefig(plt, plotpath(out, "contour_bias_closest_pt.pdf"))

        plt = set_plot(12; sz=(480, 320))
        for r in runs
            Plots.plot!(plt, bias[r]["time"], bias[r]["area_ratio"], lw=0.45, label=false,
                        color=:gray, linealpha=0.45, linestyle=:solid)
        end
        Plots.hline!(plt, [1.0], linestyle=:dash, color=:black, label=false)
        _ribbon!(plt, [bias[r]["area_ratio"] for r in runs], tgrid; color=def_orange,
                 label=L"\mathrm{mean}\;\pm\;1\,\mathrm{SD}")
        Plots.xlabel!(plt, L"\mathrm{Time\;[s]}")
        Plots.ylabel!(plt, L"\mathrm{Area\;ratio}")
        Plots.savefig(plt, plotpath(out, "contour_bias_area_ratio.pdf"))
        @info "Wrote bias figures to $(plotpath(out))"
    end

    return bias, stats
end

function contour_bias_const(filepath_gt::String=resolve_data_path(GT_CONST); kwargs...)
    return contour_bias_case(filepath_gt; kwargs...)
end

function contour_bias_bulk(; filepath_gt::String=resolve_data_path(GT_BULK),
                                                     filepath_res::String=resolve_data_path(RES_BULK),
                                                     frame_rate::Real=10.0, kwargs...)
    return contour_bias_case(filepath_gt; filepath_res=filepath_res,
                                                         case_label="bulk viscosity", frame_rate=frame_rate, kwargs...)
end

if abspath(PROGRAM_FILE) == @__FILE__
    contour_extraction_error_const()
    contour_bias_const()
    contour_bias_bulk()
end

# ---------------------------------------------------------------------------------------
# Part 2 — Tikhonov λ sweep on the synthetic bulk-viscosity experiments.
#
# η and β trade off against each other: the data misfit is nearly flat above the no-slip
# threshold (β ≳ 200), so β wanders over orders of magnitude between fitting windows while
# η absorbs the difference. The penalty ties β and leaves η free — `λ_scale = [0, λ_β]`,
# `θ_p = 0`, so `R(β) = ½·λ_β·β²/β_ref`, a soft magnitude bound rather than a physical
# prior.
#
# `Γ` is fixed rather than the optimizer's default `I ./ |θ|`, for two reasons: the default
# is `Inf` for a window starting from β = 0 (which happens), and normalizing by each
# window's own starting β would make the same λ mean a different strength in every window,
# which no sweep could be read through.
#
# The optimizer forms `ΓᵀΛΓ`, so Γ enters the penalty squared: `Γ = diag(1, 1/√β_ref)`
# gives `Γ² = diag(1, 1/β_ref)` and hence `R(β) = ½·λ_β·β²/β_ref`. The penalty therefore
# grows linearly in β_ref rather than quadratically, so at a given λ_β and β it is β_ref
# times stronger than the earlier `Γ = diag(1, 1/β_ref)` parametrization; the same penalty
# is reached at `λ_new = λ_old/β_ref`.
#
# `LAM_VALUES` is expressed in the new parametrization. The top three reproduce the earlier
# sweep's 0.7 / 2.2 / 7.0, and the bottom three extend into the gap below them, where the
# β = 1000 run was still unbiased. The `lam0p7`/`lam2p2`/`lam7p0` arms on disk are from the
# old Γ and are not comparable to anything fitted after this change.
# ---------------------------------------------------------------------------------------

const LAM_ROOT = "experiments/syn_data/lambda_sweep"     # summary CSV and figures land here
const LAM_SRC  = "experiments/syn_data/optimization/Stokes/force/bulk_viscosity/Hex2_16/cylinder"
const LAM_DEST = "experiments/syn_data/lambda_sweep/Stokes/force/bulk_viscosity/Hex2_16/cylinder"
const LAM_GT   = "ground_truth/sim_data/Stokes/force/bulk_viscosity/Hex2_16/cylinder"
const LAM_TAIL = "Hex2_6/simtime_5.0/noise_0.0/dt_0.1"
#                  old-equivalent λ:  0.02    0.07    0.2    0.7     2.2     7.0
const LAM_VALUES = [1.0e-5, 3.5e-5, 1.0e-4, 3.5e-4, 1.1e-3, 3.5e-3]
const BETA_REF = 2000.0          # penalty knee, between the true β range and the runaways

# Figure geometry for LaTeX inclusion at ½ and ⅓ of the text width, taken verbatim from the
# commented `PLOT_CONFIG` variants in `test_opt_stokes.jl` rather than derived — the margins
# differ between the two, not just the width, and a ⅓-width figure is *taller* than a ½-width
# one because it has less room for the legend.
const LAM_FIG_HALF  = (sz = (480, 320), left =  1pt, right =  5pt, top = 1pt)
const LAM_FIG_THIRD = (sz = (330, 360), left = -6pt, right = 10pt, top = 0pt)

"""
    _lam_fig(g; legend_column) -> Plots.Plot

A figure at one of the two paper widths, with that width's own margins.

`_fig` selects margins from the globals, which are set per data type; these two geometries are
per *figure size*, so they are passed explicitly.
"""
_lam_fig(g; legend_column::Int=1) =
    set_plot(fs; sz=g.sz, legend_column=legend_column,
             legend=:outerbottom, left_margin=g.left, right_margin=g.right, top_margin=g.top)

function _lam_legend(λs, outdir, run; suffix::String, gt_label=nothing, width::Int=480)
    legend = set_plot(11, sz=(round(Int, width * 1.7), 50), legend=:bottom,
                      legend_column=3,
                      bottom_margin=-35mm, top_margin=2mm,
                      left_margin=-25mm, right_margin=-15mm)
    Plots.plot!(legend, [0, 1], [0, 0], label=false, color=:white, linewidth=0)
    cmap = _lam_cmap(λs)
    if !isnothing(gt_label)
        Plots.plot!(legend, [0, 1], [0, 0], label=gt_label, color=def_green,
                    linewidth=2, markerstrokewidth=0)
    end
    for (λ, _) in _lam_series_order([(λ, nothing) for λ in λs])
        Plots.plot!(legend, [0, 1], [0, 0], label=_lam_label(λ),
                    color=(λ == 0 ? def_blue : cmap[float(λ)]), linewidth=2,
                    markerstrokewidth=0)
    end
    Plots.xlims!(legend, -0.2, 1.2)
    Plots.ylims!(legend, -0.5, 0.5)
    filename = isempty(suffix) ? "lambda_legend_run$(run).pdf" :
               "lambda_$(suffix)_legend_run$(run).pdf"
    Plots.savefig(legend, joinpath(outdir, filename))
end

_lam_tag(λ) = "lam$(replace(string(λ), "." => "p"))"

"""
    _lam_arm(run, λ) -> String

`data/` directory of one arm. λ = 0 is the unregularized baseline, which lives in the
production tree under `gn/`; every other λ lives in the sweep tree.
"""
_lam_arm(run::AbstractString, λ::Real) =
    λ == 0 ? datapath(joinpath(resolve_data_path(LAM_SRC), run, LAM_TAIL, "gn", "view_1", "multi_window")) :
             datapath(joinpath(resolve_data_path(LAM_DEST), run, LAM_TAIL,
                               "gn_tikhonov_$(_lam_tag(λ))", "view_1", "multi_window"))

"""
    run_lambda_sweep(; λs=LAM_VALUES, runs=nothing, β_ref=BETA_REF, force=false)

Fit every (run, λ) pair with `:gn_tikhonov`, cloning each run's existing `gn` parameters so
only the method differs.

Resumable: a fit whose `window_data/cost_windows.csv` exists is finished (that file is
written after the window loop) and is skipped unless `force`. One fit takes ~50 min, so the
full grid is an overnight job — run it detached, or it dies with the shell.

# Keyword Arguments
- `λs`: β penalty strengths to fit (default: `LAM_VALUES`).
- `runs`: run folders (default: every numeric folder in the ground-truth tree).
- `β_ref::Float64`: penalty knee (default: `BETA_REF`).
- `force::Bool`: refit even when a result is already complete (default: `false`).

# Returns
- `NamedTuple` with counts `fitted`, `skipped`, `failed`.
"""
function run_lambda_sweep(; λs=LAM_VALUES, runs=nothing, β_ref::Float64=BETA_REF,
                          force::Bool=false)
    gt_root = resolve_data_path(LAM_GT)
    runs = something(runs, sort(filter(d -> occursin(r"^[0-9]+$", d), readdir(gt_root))))
    fitted = skipped = failed = 0

    for run in runs, λβ in λs
        arm = _lam_arm(run, λβ)
        if !force && isfile(joinpath(arm, "window_data", "cost_windows.csv"))
            @info "Skipping run $run λ=$λβ (already complete)"
            skipped += 1
            continue
        end
        ep = read_json(joinpath(_lam_arm(run, 0), "experiment_parameters"))
        p = Dict{String,Any}(string(k) => v for (k, v) in pairs(ep))
        # Stored paths are relative to the data dir, so they need resolving before reuse.
        p["filepath_gt"]   = resolve_data_path(String(p["filepath_gt"]))
        p["filepath_res"]  = dirname(dirname(dirname(arm)))   # …/<method>/multi_window → …/<method>
        p["filepath_res"]  = joinpath(p["filepath_res"], "multi_window")
        p["opt_method"]    = :gn_tikhonov
        p["opt_kwargs"]    = Dict{String,Any}("λ_scale" => [0.0, float(λβ)],
                                              "θ_p"     => [0.0, 0.0],
                                              "Γ"       => [[1.0, 0.0], [0.0, 1/sqrt(β_ref)]])
        p["cost_function"] = :closest_point

        @info "Fitting run $run with λ_β = $λβ"
        t0 = time()
        try
            optimize(p)
            fitted += 1
            @info "run $run λ=$λβ done in $(round((time()-t0)/60, digits=1)) min"
        catch err
            failed += 1
            @warn "run $run λ=$λβ failed" exception=(err, catch_backtrace())
        end
    end
    @info "λ sweep complete" fitted skipped failed
    return (fitted=fitted, skipped=skipped, failed=failed)
end

"""
    read_lambda_arm(run, λ) -> NamedTuple

One arm's per-window estimates against the ground truth.

β is constant per window and stored per time step, so the returned `β`/`η` are the full
step series (a staircase) and `wβ`/`wη` their per-window values. `η_gt` is time-varying by
construction — these are bulk-viscosity experiments — so it is averaged within each window
before being compared.

# Returns
- `NamedTuple` with `time`, `β`, `η` (per step), `wβ`, `wη`, `wη_gt`, `ηerr` (%, per
  window), `β_gt`, `t_windows`, `rms_η`, `spread` (sd of log₁₀ β across windows) and
  `pred_err` (mean one-window-ahead height error, %).
"""
function read_lambda_arm(run::AbstractString, λ::Real)
    p = _lam_arm(run, λ)
    isdir(p) || return nothing
    sp = read_json(joinpath(resolve_data_path(LAM_GT), run, "data", "sim_params"))
    η_gt_full = Array(float.(sp["η"]))
    β_gt = Array(float.(sp["β"]))[1]

    β = vec(readdlm(joinpath(p, "est_β.csv"), ',', Float64))
    η = vec(readdlm(joinpath(p, "est_η.csv"), ',', Float64))
    dr = get_time_windows(joinpath(p, "window_data", "data_ranges.csv"))
    tw = vec(readdlm(joinpath(p, "window_data", "t_windows.csv"), ',', Float64))
    nw = size(dr, 1)

    wβ = [β[first(dr[i])] for i in 1:nw]
    wη = [η[first(dr[i])] for i in 1:nw]
    wη_gt = [mean(η_gt_full[clamp.(dr[i], 1, length(η_gt_full))]) for i in 1:nw]
    ηerr = 100 .* (wη .- wη_gt) ./ wη_gt

    gt_h = vec(readdlm(joinpath(p, "gt_h.csv"), ',', Float64))
    pred_err = NaN
    if isfile(joinpath(p, "pred_h.csv"))
        # Per-window predictions concatenated, so they must be mapped onto the time axis
        # before being compared against a time-indexed height.
        ph = _align_windowed(vec(readdlm(joinpath(p, "pred_h.csv"), ',', Float64)), dr, length(gt_h))
        ok = .!isnan.(ph)
        pred_err = 100 * mean(abs.(ph[ok] .- gt_h[ok]) ./ gt_h[ok])
    end

    est_h = vec(readdlm(joinpath(p, "est_h.csv"), ',', Float64))
    nh = min(length(est_h), length(gt_h))

    dt = length(β) > 1 ? tw[1] / length(dr[1]) : 0.1
    return (time=collect(0:length(β)-1) .* dt, β=β, η=η, wβ=wβ, wη=wη, wη_gt=wη_gt,
            ηerr=ηerr, β_gt=β_gt, t_windows=tw, data_ranges=dr, t_steps=dt,
            time_h=collect(0:nh-1) .* dt, h_rel=est_h[1:nh] ./ gt_h[1:nh],
            gt_h=gt_h[1:nh], pred_h=vec(readdlm(joinpath(p, "pred_h.csv"), ',', Float64)),
            rms_η=sqrt(mean(ηerr .^ 2)),
            spread=std(log10.(max.(wβ, 1e-3))), pred_err=pred_err)
end

"""
    lambda_sweep_table(; λs, runs, write_results=true) -> Dict

Print the sweep as a table and write `sweep_results.csv` beside the results.

Reports β per window against the constant `β_gt`, η per window against the window-averaged
`η_gt`, and the aggregate the choice of λ actually turns on: the rms η error and the
one-window-ahead prediction error, pooled over runs.

Note β is *not* scored as a percentage error. Above the no-slip threshold the profile is
insensitive to it — a 145× swing moves the data misfit by ~5% — so β = 500 and β = 5000 fit
equally well and a percentage would be meaningless. What matters there is the regime and
the spread across windows.
"""
function lambda_sweep_table(; λs=vcat(0.0, LAM_VALUES), runs=nothing, write_results::Bool=true)
    gt_root = resolve_data_path(LAM_GT)
    runs = something(runs, sort(filter(d -> occursin(r"^[0-9]+$", d), readdir(gt_root))))
    arms = Dict{Tuple{String,Float64},Any}()
    rows = Any[["run" "beta_gt" "lambda" "window" "eta_gt" "eta_est" "eta_err_pct" "beta_est"]]

    for run in runs
        r0 = read_lambda_arm(run, 0.0)
        r0 === nothing && continue
        @printf("\n=== run %s   β_gt = %-8g ===\n", run, r0.β_gt)
        @printf("%-12s %-6s %s\n", "", "", join([@sprintf("%9s", "win $i") for i in 1:length(r0.wβ)], ""))
        @printf("%-12s %-6s %s\n", "η_gt", "", join([@sprintf("%9.1f", x) for x in r0.wη_gt], ""))
        for λ in λs
            a = read_lambda_arm(run, λ)
            a === nothing && continue
            arms[(run, float(λ))] = a
            tag = λ == 0 ? "λ=0 (base)" : "λ=$λ"
            @printf("%-12s %-6s %s\n", tag, "β",  join([@sprintf("%9.1f", x) for x in a.wβ], ""))
            @printf("%-12s %-6s %s\n", "",  "η",  join([@sprintf("%9.1f", x) for x in a.wη], ""))
            @printf("%-12s %-6s %s\n", "",  "err%", join([@sprintf("%9.1f", x) for x in a.ηerr], ""))
            for i in eachindex(a.wβ)
                push!(rows, [run a.β_gt λ i round(a.wη_gt[i], digits=3) round(a.wη[i], digits=3) round(a.ηerr[i], digits=3) round(a.wβ[i], digits=3)])
            end
        end
    end

    @printf("\n=== aggregate over %d runs ===\n%-12s %10s %10s %10s %12s\n",
            length(runs), "arm", "rms_η%", "worst run", "max β", "pred_err%")
    for λ in λs
        as = [arms[(r, float(λ))] for r in runs if haskey(arms, (r, float(λ)))]
        isempty(as) && continue
        pe = filter(isfinite, [a.pred_err for a in as])
        @printf("%-12s %10.2f %10.2f %10.0f %12.3f\n", λ == 0 ? "λ=0 (base)" : "λ=$λ",
                sqrt(mean([a.rms_η for a in as] .^ 2)), maximum(a.rms_η for a in as),
                maximum(maximum(a.wβ) for a in as), isempty(pe) ? NaN : mean(pe))
    end

    if write_results
        out = joinpath(resolve_data_path(LAM_ROOT), "sweep_results.csv")
        set_file(dirname(out))
        writedlm(out, vcat(rows...), ',')
        @info "Wrote $out"
    end
    return arms
end

"""
    lambda_mosd(run, λ) -> NamedTuple

Maximum obscured surface depth of one λ arm: estimation, prediction and ground truth.

Computed on demand rather than in [`read_lambda_arm`](@ref): it projects every surface frame
into the camera frame, so it costs more than the rest of the arm put together and only the
MOSD figure needs it. All three surfaces are read from the arm itself, which keeps its own
copy of the ground truth, so no trimming against the ground-truth tree is required.

# Returns
- `NamedTuple` with `time`, `est` (already divided by the ground truth), `pred` and `gt`
  (both raw, since the prediction is indexed per window and divided segment by segment), and
  `data_ranges`; `nothing` if the arm was never fitted.
"""
function lambda_mosd(run::AbstractString, λ::Real)
    p = _lam_arm(run, λ)
    isdir(p) || return nothing
    sp = read_json(joinpath(resolve_data_path(LAM_GT), run, "data", "sim_params"))
    pose = _obj_pose_for(sp)
    hgt = Float64(sp["h"])

    est, _, _  = read_csv(joinpath(p, "sim_data", "3D_surface_points_est"))
    gt,  _, _  = read_csv(joinpath(p, "sim_data", "3D_surface_points_gt"))
    pred, _, _ = read_csv(joinpath(p, "sim_data", "3D_surface_points_pred"))
    n = min(length(est), length(gt))

    dr = get_time_windows(joinpath(p, "window_data", "data_ranges.csv"))
    tw = vec(readdlm(joinpath(p, "window_data", "t_windows.csv"), ',', Float64))
    dt = tw[1] / length(dr[1])

    m_est = get_surface_mosd(est[1:n], pose, hgt)
    m_gt  = get_surface_mosd(gt[1:n], pose, hgt)
    m_pr  = get_surface_mosd(pred, pose, hgt)

    return (time=collect(0:n-1) .* dt, est=m_est ./ m_gt, pred=m_pr, gt=m_gt, data_ranges=dr)
end

"""
    _lam_pred_segments(time, dr, pred, ref) -> Vector{Tuple{Vector,Vector}}

Split a per-window prediction into the `(t, value)` segments `replot` draws.

Each window's prediction is launched from that window's parameters and runs one sample past
it, so the segments are discontinuous at the boundaries: drawing them as a single polyline
would join the end of one forecast to the start of the next, a line the model never produced.
`ref` is divided out index by index, so the caller controls whether the result is a ratio to
the ground truth or to a fixed residual.

Layout detection is delegated to `_pred_range`, so a single-window arm — whose prediction is
already indexed by time — degrades to one segment covering the record.
"""
function _lam_pred_segments(time::AbstractVector, dr, pred::AbstractVector, ref::AbstractVector)
    segs = Tuple{Vector{Float64},Vector{Float64}}[]
    windowed = length(pred) > length(ref)
    start = 0
    for ti in 1:length(dr)
        range_ = collect(dr[ti][1]:(dr[ti][end] + 1))
        rp = _pred_range(ti, dr, start, range_, windowed)
        n = min(length(range_), length(rp), length(time), length(ref))
        n > 0 || continue
        range_, rp = range_[1:n], rp[1:n]
        push!(segs, (collect(time[range_]), pred[rp] ./ ref[range_]))
        start = rp[end] + 1
    end
    return segs
end

"""
    _window_series!(plt, a, values, label, color; pred=true)

Draw a per-window quantity the way `replot` draws η and β for bulk-viscosity fits: one solid
segment per window for the estimate, and — from the second window on — the previous window's
value held across the current one, dashed, as its prediction.

The estimates are piecewise constant by construction, one value per window held across it, so
a continuous line would draw a ramp across each boundary that no fit ever produced. Only the
first segment carries the label, so each λ contributes one legend entry.
"""
function _window_series!(plt, a, values::AbstractVector, label, color; pred::Bool=true)
    t_prev = a.t_steps
    prev = 0.0
    for ti in 1:length(a.t_windows)
        dr = a.data_ranges[ti]
        t_win = collect(range(start=t_prev, stop=a.t_windows[ti], step=a.t_steps))
        n = min(length(t_win), length(dr))
        Plots.plot!(plt, t_win[1:n], values[dr[1:n]], label=(ti == 1 ? label : false),
                    color=color)
        if pred && ti > 1
            Plots.plot!(plt, t_win[1:n], fill(prev, n), label=false, color=color,
                        linestyle=:dash)
        end
        prev = values[dr[end]]
        t_prev = a.t_windows[ti] + a.t_steps
    end
end

# Colour distinguishes the λ arms; linestyle keeps `replot`'s meaning, solid for the
# estimation and dashed for the prediction. def_green is reserved for the ground truth, and
# def_red and def_orange are close enough in hue that they are never placed adjacent here.
const LAM_LINE_COLORS = [def_blue, def_red, def_orange, color_palette[5], color_palette[8]]

"""
    _lam_cmap(λs) -> Dict{Float64,Any}

Map each λ in the order given onto a line colour, so an arm keeps its colour across every
figure in a set regardless of which arms are drawn.
"""
_lam_cmap(λs) = Dict(float(λ) => LAM_LINE_COLORS[mod1(i, length(LAM_LINE_COLORS))]
                     for (i, λ) in enumerate(λs))

_lam_series_order(series) = sort(series; by=x -> (x[1] == 0 ? 0 : 1, -float(x[1])))

"""
    _lam_paper(λ) -> Float64

The λ of the manuscript's formulation, given the λ this code takes.

The two are not the same number. `_fit_model_GN_tikhonov` forms `R(θ) = ½·θᵀΓᵀΛΓθ` with
`Γ = diag(1, 1/√β_ref)`, so its penalty is `½·λ·β²/β_ref`. The manuscript writes the penalty
as `½‖Λθ‖²` with `Λ = diag(0, λ_β)`, i.e. `½·λ_β²·β²`. They agree at `λ_β = √(λ/β_ref)`, so
every figure and table facing the paper must report the root, not the value passed to the fit.
"""
_lam_paper(λ) = sqrt(float(λ) / BETA_REF)

"""
    _lam_label(λ) -> LaTeXString

Legend entry for one λ arm, in the manuscript's units. See [`_lam_paper`](@ref).
"""
function _lam_label(λ)
    λ == 0 && return L"\lambda_\beta = 0"
    v = _lam_paper(λ)
    e = floor(Int, log10(v))
    m = round(Int, v / 10.0^e)
    return latexstring("\\lambda_\\beta = $m \\times 10^{$e}")
end

"""
    plot_lambda_sweep(; λs, runs, outdir) -> String

Per-run figures for the λ sweep, plus the trade-off figure the choice of λ turns on.

Four figures per run, following `replot`'s conventions for a bulk-viscosity fit so they read
against the per-experiment ones: β and η drawn window by window, and the height and MOSD
ratios drawn against the ground truth. Each carries one curve per λ, estimation solid and
prediction dashed — colour is spent on λ, so that convention belongs in the caption — with
the prediction segmented per window exactly as `replot` segments it.

# Returns
- `String`: the directory the figures were written to.
"""
function plot_lambda_sweep(; λs=vcat(0.0, LAM_VALUES), runs=nothing,
                           outdir::String=plotpath(resolve_data_path(LAM_ROOT)))
    gt_root = resolve_data_path(LAM_GT)
    runs = something(runs, sort(filter(d -> occursin(r"^[0-9]+$", d), readdir(gt_root))))
    set_file(outdir)
    cmap = _lam_cmap(λs)

    rms, pred, ok_λ = Float64[], Float64[], Float64[]
    for λ in λs
        as = [read_lambda_arm(r, λ) for r in runs]
        as = filter(!isnothing, as)
        isempty(as) && continue
        push!(ok_λ, float(λ))
        push!(rms, sqrt(mean([a.rms_η for a in as] .^ 2)))
        pe = filter(isfinite, [a.pred_err for a in as])
        push!(pred, isempty(pe) ? NaN : mean(pe))
    end

    for run in runs
        base = read_lambda_arm(run, 0.0)
        base === nothing && continue
        series = [(λ, read_lambda_arm(run, λ)) for λ in λs]
        series = _lam_series_order([(λ, a) for (λ, a) in series if a !== nothing])
        _lam_legend(λs, outdir, run; suffix="beta", gt_label=L"\beta_{\mathrm{gt}}")
        _lam_legend(λs, outdir, run; suffix="eta", gt_label=L"\mathrm{gt}")
        _lam_legend(λs, outdir, run; suffix="", width=330)
        # `replot` cuts every per-experiment figure at `end_obs_win`; match it so the
        # sweep figures can be laid beside them.
        t_end = end_obs_win

        plt_β = _lam_fig(LAM_FIG_HALF, legend_column=6)
        Plots.plot!(plt_β, [], label=false)
        for t in base.t_windows
            Plots.vline!(plt_β, [t], color=:gray, linestyle=:dash, label=false)
        end
        Plots.hline!(plt_β, [base.β_gt], color=def_green, label=L"\beta_{\mathrm{gt}}")
        for (λ, a) in series
            _window_series!(plt_β, a, a.β, _lam_label(λ), cmap[float(λ)]; pred=false)
        end
        # A collapsed window sits at exactly zero, which a log axis cannot render.
        if all(all(a.β .> 0) for (_, a) in series) && base.β_gt > 0
            Plots.plot!(plt_β, yscale=:log10)
        end
        _label!(plt_β, L"\mathrm{Time\;[s]}", latexstring("\$\\beta(t)\$ [MPa s m\$^{-1}\$]");
                xlims=(0, t_end))
        Plots.savefig(plt_β, joinpath(outdir, "lambda_beta_run$(run).pdf"))

        plt_η = _lam_fig(LAM_FIG_HALF, legend_column=6)
        Plots.plot!(plt_η, [], label=false)
        for t in base.t_windows
            Plots.vline!(plt_η, [t], color=:gray, linestyle=:dash, label=false)
        end
        η_gt_full = Array(float.(read_json(joinpath(gt_root, run, "data", "sim_params"))["η"]))
        n = min(length(η_gt_full), length(base.time))
        Plots.plot!(plt_η, base.time[1:n], η_gt_full[1:n], color=def_green,
                label=L"\mathrm{gt}")
        for (λ, a) in series
            _window_series!(plt_η, a, a.η, _lam_label(λ), cmap[float(λ)]; pred=false)
        end
        _label!(plt_η, L"\mathrm{Time\;[s]}", latexstring("\$\\eta(t)\$ [kPa s]");
                xlims=(0, t_end))
        Plots.savefig(plt_η, joinpath(outdir, "lambda_eta_run$(run).pdf"))

        plt_h = _lam_fig(LAM_FIG_THIRD, legend_column=5)
        Plots.plot!(plt_h, [], label=false)
        for t in base.t_windows
            Plots.vline!(plt_h, [t], color=:gray, linestyle=:dash, label=false)
        end
        Plots.hline!(plt_h, [1.0], color=:black, linestyle=:dot, label=false)
        for (λ, a) in series
            for (t, v) in _lam_pred_segments(a.time_h, a.data_ranges, a.pred_h, a.gt_h)
                Plots.plot!(plt_h, t, v, label=false, color=cmap[float(λ)], linestyle=:dash)
            end
            Plots.plot!(plt_h, a.time_h, a.h_rel, label=_lam_label(λ), color=cmap[float(λ)])
        end
        _label!(plt_h, L"\mathrm{Time\;[s]}", L"h/h_{\mathrm{gt}}";
                xlims=(0, t_end), ylims=y_lims_h_norm)
        Plots.savefig(plt_h, joinpath(outdir, "lambda_height_run$(run).pdf"))

        mosds = [(λ, lambda_mosd(run, λ)) for (λ, _) in series]
        mosds = [(λ, m) for (λ, m) in mosds if m !== nothing]
        if !isempty(mosds)
            plt_m = _lam_fig(LAM_FIG_THIRD, legend_column=5)
            Plots.plot!(plt_m, [], label=false)
            for t in base.t_windows
                Plots.vline!(plt_m, [t], color=:gray, linestyle=:dash, label=false)
            end
            Plots.hline!(plt_m, [1.0], color=:black, linestyle=:dot, label=false)
            for (λ, m) in mosds
                for (t, v) in _lam_pred_segments(m.time, m.data_ranges, m.pred, m.gt)
                    Plots.plot!(plt_m, t, v, label=false, color=cmap[float(λ)], linestyle=:dash)
                end
                Plots.plot!(plt_m, m.time, m.est, label=_lam_label(λ), color=cmap[float(λ)])
            end
            _label!(plt_m, L"\mathrm{Time\;[s]}", L"\mathrm{Relative\;MOSD}";
                    xlims=(0, t_end), ylims=y_lims_h_norm)
            Plots.savefig(plt_m, joinpath(outdir, "lambda_mosd_run$(run).pdf"))
        end
    end

    # The trade-off the choice of λ turns on: η accuracy against predictive quality.
    plt = _fig(margins=:all, legend_column=2)
    Plots.plot!(plt, [], label=false)
    Plots.plot!(plt, ok_λ, rms, marker=:circle, color=def_red, label=L"\mathrm{rms}\;\eta\;\mathrm{error}\;[\%]")
    _label!(plt, L"\lambda_\beta", L"\mathrm{rms}\;\eta\;\mathrm{error}\;[\%]")
    Plots.plot!(Plots.twinx(), ok_λ, pred, marker=:square, color=def_blue, linestyle=:dash,
                label=L"\mathrm{prediction\;error}\;[\%]", ylabel=L"\mathrm{prediction\;error}\;[\%]",
                legend=:topright)
    Plots.savefig(plt, joinpath(outdir, "lambda_tradeoff.pdf"))

    @info "Wrote λ sweep figures to $outdir"
    return outdir
end

"""
    lambda_contour_error(run, λ) -> NamedTuple

Per-frame contour error of one λ arm, for the estimated and the predicted contours.

Observations come through `_get_borders` and the contours from the arm itself, so both series
are by construction the ones the per-experiment figures show for that fit. The prediction is
scored window by window against the matching observations, as `replot` scores it.

# Returns
- `NamedTuple` with `time`, `d` (estimation, raw), `rel` (`d` over that fit's residual at the
  first Gauss--Newton iteration of the first window), `pred_segments` as `(t, d)` pairs, and
  `t_windows`; `nothing` if the arm was never fitted.
"""
function lambda_contour_error(run::AbstractString, λ::Real)
    p = _lam_arm(run, λ)
    isdir(p) || return nothing
    win = dirname(p)

    ep = read_json(joinpath(p, "experiment_parameters"))
    data_type = ep["data_type"]
    filepath_gt = resolve_data_path(String(ep["filepath_gt"]))

    dr = get_time_windows(joinpath(p, "window_data", "data_ranges.csv"))
    tw = vec(readdlm(joinpath(p, "window_data", "t_windows.csv"), ',', Float64))
    t_steps = tw[1] / length(dr[1])
    n = dr[end][end] + 1

    obs, sim, _, _, _, _ = _get_borders(data_type, filepath_gt, win, n)
    d_est, _ = contour_cost(sim, obs)
    cost_init = readdlm(joinpath(p, "window_data", "cost_windows.csv"), ',', '\n')[1, 1]
    time = collect(range(start=0, stop=(n - 1) * t_steps, step=t_steps))

    pred, _, _ = read_csv(joinpath(p, "sim_data", "view_1", "2D_border_points_pred"))
    windowed = length(pred) > n
    segs = Tuple{Vector{Float64},Vector{Float64}}[]
    start = 0
    for ti in 1:length(dr)
        range_ = collect(dr[ti][1]:(dr[ti][end] + 1))
        rp = _pred_range(ti, dr, start, range_, windowed)
        m = min(length(range_), length(rp), length(obs), length(pred))
        m > 0 || continue
        range_, rp = range_[1:m], rp[1:m]
        push!(segs, (time[range_], contour_cost(pred[rp], obs[range_])[1]))
        start = rp[end] + 1
    end

    return (time=time, d=d_est, rel=d_est ./ cost_init, cost_init=cost_init,
            pred_segments=segs, t_windows=tw)
end

"""
    plot_lambda_contour_error(; λs, runs, outdir) -> String

One figure per run: the contour error of every λ arm, estimation solid and prediction dashed.

Normalised by each fit's own residual at the first Gauss--Newton iteration of the first
window — `replot`'s convention for this quantity, and the one the tables report — so the
curves sit on the same scale as the per-experiment figures rather than on a scale defined by
one of the arms. λ = 0 is drawn as an ordinary curve, so its own contour fit is visible
instead of being flattened to the unit reference. The dashed line at unity is the
initialisation: a curve above it fits the contours worse than the starting guess did.

# Returns
- `String`: the directory the figures were written to.
"""
function plot_lambda_contour_error(; λs=vcat(0.0, LAM_VALUES), runs=nothing,
                                   outdir::String=plotpath(resolve_data_path(LAM_ROOT)))
    gt_root = resolve_data_path(LAM_GT)
    runs = something(runs, sort(filter(d -> occursin(r"^[0-9]+$", d), readdir(gt_root))))
    set_file(outdir)
    cmap = _lam_cmap(λs)

    for run in runs
        arms = [(λ, lambda_contour_error(run, λ)) for λ in λs]
        arms = _lam_series_order([(λ, a) for (λ, a) in arms if a !== nothing])
        isempty(arms) && continue

        plt = _lam_fig(LAM_FIG_THIRD, legend_column=5)
        Plots.plot!(plt, [], label=false)
        for t in arms[1][2].t_windows
            Plots.vline!(plt, [t], color=:gray, linestyle=:dash, label=false)
        end
        Plots.hline!(plt, [1.0], color=:black, linestyle=:dot, label=false)

        hi = 0.0
        for (λ, a) in arms
            for (t, v) in a.pred_segments
                r = v ./ a.cost_init
                hi = max(hi, maximum(r))
                Plots.plot!(plt, t, r, label=false, color=cmap[float(λ)], linestyle=:dash)
            end
            hi = max(hi, maximum(a.rel))
            Plots.plot!(plt, a.time, a.rel, label=_lam_label(λ), color=cmap[float(λ)])
        end

        _label!(plt, L"\mathrm{Time\;[s]}", L"\mathrm{Relative\;Cost}";
                xlims=(0, end_obs_win), ylims=(0, 1.1 * hi))
        Plots.savefig(plt, joinpath(outdir, "lambda_contour_error_run$(run).pdf"))
    end
    @info "Wrote λ contour-error figures to $outdir"
    return outdir
end

# ---------------------------------------------------------------------------------------
# Part 3 — sampling bias of the contour metrics
#
# The optimizer scores a coarse simulated contour against a dense observed one, and a
# nearest-neighbour distance measured to a sparse set of points overstates the distance to
# the curve those points lie on. That overstatement is a property of the sampling, not of
# the fit, and it is what has to be quantified before either cost can be trusted.
#
# The mesh-convergence tree gives a clean instrument: the same simulation discretized at
# ne = 2…16, whose `contour_data` is the border polyline subdivided ten times per segment,
# so the point count is tied to the mesh (501 points at ne = 6, 1961 at ne = 16).
# Comparing ne = 6 against ne = 16 as stored therefore mixes two effects — the ne = 6
# solution genuinely sits elsewhere, *and* it is described by four times fewer points. Only
# the second is metric bias.
#
# The two are separated by measuring the same pair twice: once with both sides resampled
# onto a common dense arclength grid, which is the sampling-free geometric difference, and
# once as stored, which is what the optimizer sees. The gap between them is the bias.
# ---------------------------------------------------------------------------------------

const MESH_CONV = "ground_truth/sim_data/Stokes/force/constant/Hex_2/convergence_analysis/mesh_convergence_analysis"
const MESH_STRIDE = 40      # 2001 frames at dt = 0.0025 s → 51 frames on the fits' dt = 0.1 s grid
const MESH_DENSITIES = [501, 1000, 2000, 4000, 8000]
const MESH_CONVERGED = 16000  # density at which both sides count as sampling-free
const MESH_BIAS_OUT = "experiments/sim_data/mesh_metric_bias"   # simulated ground truth, not syn
const MESH_OPT_ROOT = "experiments/sim_data/optimization/Stokes/force/constant/Hex2_16/cylinder"

"""
    read_mesh_contours(ne; stride, kind) -> Vector

Contour frames of one mesh from the convergence tree, subsampled onto the fits' time grid.

`kind` selects `"contour_data"` (the interpolated contour, ten points per border segment) or
`"2D_border_points"` (the raw projected mesh nodes).
"""
read_mesh_contours(ne::Int; stride::Int=MESH_STRIDE, kind::String="contour_data") =
    read_csv(joinpath(resolve_data_path(MESH_CONV), "mesh_sz_$ne", "data",
                      "sim_data", "view_1", kind))[1][1:stride:end]

"""
    mesh_contour_bias(; ne_test, ne_ref, stride, densities, n_conv, write_results) -> Dict

Sampling bias of the closest-point and Chamfer metrics, measured on the mesh-convergence
contours.

Three quantities per frame, all in px²:

- `*_geom`: both contours resampled to `n_conv` points, so the sampling contribution is
  negligible. This is the genuine discretization error of `ne_test` against `ne_ref`.
- `*_raw`: the contours exactly as stored, which is what the optimizer's cost sees.
- `*_up[N]`: the coarse side resampled to `N` points with the reference left as stored — the
  operational fix, since in a fit only the simulated side can be resampled.

`*_raw - *_geom` is the bias, and `*_up[N] - *_geom` shows how much of it upsampling removes.

# Returns
- `Dict` with the per-frame series and the density sweep; also written to
  `<LAM_ROOT>/../mesh_metric_bias/` when `write_results`.
"""
function mesh_contour_bias(; ne_test::Int=6, ne_ref::Int=16, stride::Int=MESH_STRIDE,
                           densities=MESH_DENSITIES, n_conv::Int=MESH_CONVERGED,
                           kind::String="contour_data", kind_ref::String=kind,
                           write_results::Bool=true)
    test = read_mesh_contours(ne_test; stride=stride, kind=kind)
    ref  = read_mesh_contours(ne_ref;  stride=stride, kind=kind_ref)
    n = min(length(test), length(ref))
    test, ref = test[1:n], ref[1:n]
    @info "mesh bias: ne=$ne_test ($(_npts(test[1])) pts) vs ne=$ne_ref ($(_npts(ref[1])) pts), $n frames"

    # Phase 1 — both sides dense, so what is left is geometry alone.
    tc = [upsample_contour(c, n_conv) for c in test]
    rc = [upsample_contour(c, n_conv) for c in ref]
    cp_geom = [closest_point_distance_kdtree(a, b)^2 for (a, b) in zip(tc, rc)]
    ch_geom = [chamfer_sq_distance_kdtree(a, b) for (a, b) in zip(tc, rc)]

    # Phase 2 — as stored.
    cp_raw = [closest_point_distance_kdtree(a, b)^2 for (a, b) in zip(test, ref)]
    ch_raw = [chamfer_sq_distance_kdtree(a, b) for (a, b) in zip(test, ref)]

    # Phase 3 — resample the coarse side only; the reference stays as the fit would find it.
    cp_up = Dict{Int,Vector{Float64}}()
    ch_up = Dict{Int,Vector{Float64}}()
    for N in densities
        u = [upsample_contour(c, N) for c in test]
        cp_up[N] = [closest_point_distance_kdtree(a, b)^2 for (a, b) in zip(u, ref)]
        ch_up[N] = [chamfer_sq_distance_kdtree(a, b) for (a, b) in zip(u, ref)]
        @info "  N=$N: closest point $(round(mean(cp_up[N]), sigdigits=4)) px², chamfer $(round(mean(ch_up[N]), sigdigits=4)) px²"
    end

    out = Dict{String,Any}("frames" => n, "ne_test" => ne_test, "ne_ref" => ne_ref,
                           "n_pts_test" => _npts(test[1]), "n_pts_ref" => _npts(ref[1]),
                           "densities" => collect(densities), "n_conv" => n_conv,
                           "cp_geom" => cp_geom, "ch_geom" => ch_geom,
                           "cp_raw" => cp_raw, "ch_raw" => ch_raw,
                           "cp_up" => cp_up, "ch_up" => ch_up)

    if write_results
        dir = joinpath(resolve_data_path(MESH_BIAS_OUT),
                       kind == kind_ref ? kind : "$(kind)_vs_$(kind_ref)")
        set_file(dir)
        writedlm(joinpath(dir, "per_frame_ne$(ne_test)_vs_ne$(ne_ref).csv"),
                 hcat(cp_geom, cp_raw, ch_geom, ch_raw), ',')
        writedlm(joinpath(dir, "density_sweep_ne$(ne_test)_vs_ne$(ne_ref).csv"),
                 hcat(collect(densities),
                      [mean(cp_up[N]) for N in densities],
                      [mean(ch_up[N]) for N in densities]), ',')
        @info "Wrote mesh metric-bias results to $dir"
    end
    return out
end

"""
    _decimate(c, n) -> Matrix

Keep `n` evenly spaced points of a contour, preserving its stored orientation.

Used to simulate having only a coarse mesh's worth of vertices *on a curve that is otherwise
known exactly*, which is what separates interpolation error from discretization error.
"""
function _decimate(c::AbstractArray, n::Int)
    flip = size(c, 1) != 2
    p = flip ? permutedims(c) : c
    idx = unique(round.(Int, range(1, size(p, 2), length=n)))
    out = p[:, idx]
    return flip ? permutedims(out) : out
end

"""
    interpolation_error(; ne_ref, n_verts, n_dense, stride, method, write_results) -> Dict

Error the upsampling interpolation itself introduces, isolated from FEM discretization.

Takes the finest mesh's contour as truth, decimates it to `n_verts` points, splines it back
to `n_dense`, and scores the reconstruction against the truth. Because both sides come from
the *same* solution, no discretization error enters: the residual is the interpolation error
alone, at the vertex count a simulated contour actually carries.

This is the quantity that licenses using an upsampled contour as the unbiased comparator. It
is only meaningful against the bias it is meant to remove — `mesh_contour_bias`'s
`*_raw - *_geom`, +0.43 px² for the closest-point cost — and upsampling is only worth doing
if it lands well below that.

# Returns
- `Dict` mapping each vertex count to per-frame closest-point and Chamfer errors, in px².
"""
function interpolation_error(; ne_ref::Int=16, n_verts=[18, 50, 196], n_dense::Int=16000,
                             stride::Int=MESH_STRIDE, method::Symbol=:spline,
                             write_results::Bool=true)
    ref = read_mesh_contours(ne_ref; stride=stride)
    truth = [upsample_contour(c, n_dense; method=method) for c in ref]
    @info "interpolation error: truth = ne=$ne_ref ($(_npts(ref[1])) pts) at $n_dense, $(length(ref)) frames"

    cp = Dict{Int,Vector{Float64}}()
    ch = Dict{Int,Vector{Float64}}()
    for nv in n_verts
        rec = [upsample_contour(_decimate(c, nv), n_dense; method=method) for c in ref]
        cp[nv] = [closest_point_distance_kdtree(a, b)^2 for (a, b) in zip(rec, truth)]
        ch[nv] = [chamfer_sq_distance_kdtree(a, b) for (a, b) in zip(rec, truth)]
        @info "  $nv vertices: closest point $(round(mean(cp[nv]), sigdigits=4)) px², chamfer $(round(mean(ch[nv]), sigdigits=4)) px²"
    end

    out = Dict{String,Any}("n_verts" => collect(n_verts), "n_dense" => n_dense,
                           "ne_ref" => ne_ref, "method" => String(method),
                           "cp" => cp, "ch" => ch)
    if write_results
        dir = resolve_data_path(MESH_BIAS_OUT)
        set_file(dir)
        writedlm(joinpath(dir, "interpolation_error_$(method).csv"),
                 hcat(collect(n_verts), [mean(cp[n]) for n in n_verts],
                      [mean(ch[n]) for n in n_verts]), ',')
        @info "Wrote interpolation-error results to $dir"
    end
    return out
end

"""
    _perimeter(c) -> Float64

Arclength of a contour as stored, whichever way round its axes are.
"""
function _perimeter(c::AbstractArray)
    p = size(c, 1) == 2 ? c : permutedims(c)
    return sum(sqrt((p[1, i+1] - p[1, i])^2 + (p[2, i+1] - p[2, i])^2) for i in 1:size(p, 2)-1)
end

"""
    obs_density_bias(; ne_test, ne_ref, obs_counts, stride, n_dense) -> Dict

Closest-point bias as a function of the **observation's** point spacing.

`ClosestPointCost` matches each simulated point to its nearest observed point, so the
distance floor is a property of the observation, not the simulation. Decimating the observed
contour and holding the simulated one fixed at its raw mesh nodes measures that directly,
rather than asserting it: the bias should scale with the observed spacing and be independent
of how the simulated side is sampled.

# Returns
- `Dict` with `obs_counts`, `spacing` (px between observed points), `cost`, `bias`
  (cost minus the both-dense geometric value) and `perimeter`.
"""
function obs_density_bias(; ne_test::Int=6, ne_ref::Int=16,
                          obs_counts=[50, 100, 200, 500, 1000, 1961],
                          stride::Int=MESH_STRIDE, n_dense::Int=16000)
    sim = read_mesh_contours(ne_test; kind="2D_border_points", stride=stride)
    obs = read_mesh_contours(ne_ref;  kind="contour_data", stride=stride)
    n = min(length(sim), length(obs))
    sim, obs = sim[1:n], obs[1:n]

    sd = [upsample_contour(c, n_dense) for c in sim]
    od = [upsample_contour(c, n_dense) for c in obs]
    geom = mean([closest_point_distance_kdtree(a, b)^2 for (a, b) in zip(sd, od)])
    per = mean(_perimeter.(obs))

    spacing, cost = Float64[], Float64[]
    for m in obs_counts
        o = [_decimate(c, m) for c in obs]
        push!(cost, mean([closest_point_distance_kdtree(a, b)^2 for (a, b) in zip(sim, o)]))
        push!(spacing, per / (m - 1))
        @info "  obs $m pts ($(round(per/(m-1), digits=2)) px spacing): closest point $(round(cost[end], sigdigits=4)) px²"
    end
    return Dict("obs_counts" => collect(obs_counts), "spacing" => spacing, "cost" => cost,
                "bias" => cost .- geom, "geom" => geom, "perimeter" => per)
end

"""
    run_metric_bias_analysis(; runs, stride, outdir) -> String

Run the contour-metric bias analysis and emit one figure per step.

The argument runs in three stages, each answering the question the next one depends on. All of
it is simulated (`sim_data`) ground truth — the mesh-convergence tree for the first two
stages, the constant-viscosity fits for the third.

1. **What does interpolation cost?** `interpolation_error_chamfer.pdf` — mesh_16 decimated to
   a coarse vertex count, splined back up, and scored against mesh_16 itself. Both sides come
   from the same solution, so no discretization error enters and the residual is the
   interpolation error alone. This is what licenses using an upsampled contour as a ruler.
2. **What does that ruler measure?** `chamfer_upsampling_mesh6_vs_mesh16.pdf` — mesh_6 against
   mesh_16, Chamfer, as the simulated side is upsampled. Chamfer without upsampling is
   dominated by the vertex-matching artefact; the curve shows it falling to the geometric
   difference once the sampling is matched.
3. **Do the costs converge to the same place?** `contour_metrics_iter.pdf` per experiment,
   from [`iteration_contour_metrics`](@ref): closest point, Chamfer without upsampling and
   Chamfer with upsampling against Gauss--Newton iteration.

# Returns
- `String`: the directory the first two figures were written to.
"""
function run_metric_bias_analysis(; runs=["1", "2", "3", "4", "5", "6"],
        stride::Int=MESH_STRIDE,
        outdir::String=plotpath(resolve_data_path(MESH_BIAS_OUT)))
    set_file(outdir)
    verts = [18, 26, 50, 98, 196]
    dens  = [501, 1000, 2000, 4000, 8000]

    # 1 — interpolation error, measured against the mesh it came from.
    # Spline only: the upsampling this licenses goes from the ne=6 vertex count up to ne=16's,
    # and `fit_curve` builds the observed contour the same way, so a linear arm would measure
    # a construction nothing in the pipeline uses.
    is = interpolation_error(n_verts=verts, method=:spline, stride=stride)

    p1 = _fig(margins=:all)
    Plots.plot!(p1, verts, [mean(is["ch"][n]) for n in verts], marker=:circle,
                color=def_red, label=false)
    Plots.plot!(p1, xscale=:log10, yscale=:log10)
    _label!(p1, L"\mathrm{Contour\;vertices}", L"\mathrm{Chamfer}\;[\mathrm{px}^2]")
    Plots.savefig(p1, joinpath(outdir, "interpolation_error_chamfer.pdf"))

    # 2 — Chamfer between the two meshes as the coarse side is upsampled.
    tests = read_mesh_contours(6;  kind="contour_data", stride=stride)
    refs  = read_mesh_contours(16; kind="contour_data", stride=stride)
    nf = min(length(tests), length(refs)); tests, refs = tests[1:nf], refs[1:nf]
    bias = mesh_contour_bias(kind="contour_data", densities=dens, stride=stride)
    ch_geom = mean(bias["ch_geom"])

    p2 = _fig(margins=:all, legend_column=2)
    Plots.plot!(p2, [], label=false)
    Plots.plot!(p2, dens, [mean(bias["ch_up"][N]) for N in dens], marker=:circle,
                color=def_red, label="Chamfer")
    Plots.hline!(p2, [ch_geom], color=:black, linestyle=:dash, label="Geometric difference")
    Plots.plot!(p2, xscale=:log10, yscale=:log10)
    _label!(p2, L"\mathrm{Simulated\;contour\;points}", L"\mathrm{Chamfer}\;[\mathrm{px}^2]")
    Plots.savefig(p2, joinpath(outdir, "chamfer_upsampling_mesh6_vs_mesh16.pdf"))

    # Per-frame Chamfer over the record, at mesh_16's own point count so both arms are scored
    # against the same reference sampling.
    ref_n  = _npts(refs[1])
    test_n = _npts(tests[1])

    # (a) interpolation error: mesh_16 decimated to mesh_6's count, splined back to mesh_16's.
    # Same solution on both sides, so this is what the interpolation costs and nothing else.
    round_trip = [upsample_contour(_decimate(c, test_n), ref_n) for c in refs]
    ch_interp = [chamfer_sq_distance_kdtree(a, b) for (a, b) in zip(round_trip, refs)]

    # (b) the measurement: mesh_6 splined up to mesh_16's count, against mesh_16.
    up6 = [upsample_contour(c, ref_n) for c in tests]
    ch_mesh6 = [chamfer_sq_distance_kdtree(a, b) for (a, b) in zip(up6, refs)]

    t = collect(0:length(ch_interp)-1) .* 0.1
    ylim = (0, 1.1 * maximum(vcat(ch_interp, ch_mesh6)))

    pa = _fig(margins=:all)
    Plots.plot!(pa, t, ch_interp, color=def_blue, label=false)
    _label!(pa, L"\mathrm{Time\;[s]}", L"\mathrm{Chamfer}\;[\mathrm{px}^2]";
            xlims=(0, t[end]), ylims=ylim)
    Plots.savefig(pa, joinpath(outdir, "chamfer_interpolation_error_time.pdf"))

    pb = _fig(margins=:all)
    Plots.plot!(pb, t, ch_mesh6, color=def_red, label=false)
    _label!(pb, L"\mathrm{Time\;[s]}", L"\mathrm{Chamfer}\;[\mathrm{px}^2]";
            xlims=(0, t[end]), ylims=ylim)
    Plots.savefig(pb, joinpath(outdir, "chamfer_mesh6_upsampled_time.pdf"))

    pc = _fig(margins=:all, legend_column=2)
    Plots.plot!(pc, [], label=false)
    Plots.plot!(pc, t, ch_interp, color=def_blue, label="Interpolation error")
    Plots.plot!(pc, t, ch_mesh6,  color=def_red,  label="mesh 6 upsampled")
    _label!(pc, L"\mathrm{Time\;[s]}", L"\mathrm{Chamfer}\;[\mathrm{px}^2]";
            xlims=(0, t[end]), ylims=ylim)
    Plots.savefig(pc, joinpath(outdir, "chamfer_interpolation_vs_mesh6_time.pdf"))

    # 3 — convergence of the three costs over the optimizer's iterations, per experiment.
    root = resolve_data_path(MESH_OPT_ROOT)
    for run in runs
        exp_path = joinpath(root, run, "Hex2_6", "simtime_5.0", "noise_0.0", "dt_0.1",
                            "gn", "view_1")
        isdir(joinpath(exp_path, "data", "opt_data")) || continue
        ep = read_json(datapath(exp_path, "experiment_parameters"))
        try
            iteration_contour_metrics(String(ep["data_type"]),
                                      resolve_data_path(String(ep["filepath_gt"])), exp_path)
            # `iteration_contour_metrics` writes into the experiment's own plots/ directory.
            # Copy rather than move, so the figure stays with the fit that produced it and a
            # re-run of either path stays consistent, while the analysis folder holds the set.
            src = joinpath(exp_path, "plots", "contour_metrics_iter.pdf")
            isfile(src) && cp(src, joinpath(outdir, "contour_metrics_iter_run$(run).pdf"),
                              force=true)
        catch err
            @warn "convergence metrics failed for run $run" exception=err
        end
    end

    @info "Wrote metric-bias figures to $outdir"
    return outdir
end

"""
    plot_beta_bounding(; run, λs, outdir) -> String

The figure for the claim that L2 regularization bounds β.

Per-window β against λ, with the unregularized fit as a reference line. Without a penalty β
runs far above its true value and scatters across windows; any λ above a threshold bounds it,
and the bound tightens monotonically. λ = 0 cannot sit on a log axis and is off the top of the
range anyway, so it is drawn as a horizontal line rather than a point.

Whiskers span the per-window minimum and maximum, so the figure shows the scatter collapsing
as well as the median coming down — the two things the penalty is there to do.

# Returns
- `String`: the directory the figure was written to.
"""
function plot_beta_bounding(; run::AbstractString="6",
                            λs=[1.0e-5, 1.16e-5, 1.35e-5, 3.5e-5, 1.0e-4, 3.5e-4],
                            outdir::String=plotpath(resolve_data_path(LAM_ROOT)))
    set_file(outdir)
    base = read_lambda_arm(run, 0.0)
    base === nothing && error("unregularized baseline missing for run $run")

    xs, med, lo, hi = Float64[], Float64[], Float64[], Float64[]
    for λ in λs
        a = read_lambda_arm(run, λ)
        a === nothing && continue
        push!(xs, float(λ)); push!(med, median(a.wβ))
        push!(lo, minimum(a.wβ)); push!(hi, maximum(a.wβ))
    end

    plt = _lam_fig(LAM_FIG_THIRD, legend_column=3)
    Plots.plot!(plt, [], label=false)
    Plots.hline!(plt, [median(base.wβ)], color=def_blue, linestyle=:dash,
                 label=L"\lambda_\beta = 0")
    Plots.hline!(plt, [base.β_gt], color=def_green, label=L"\beta_{\mathrm{gt}}")
    Plots.plot!(plt, xs, med, yerror=(med .- lo, hi .- med), marker=:circle,
                color=def_blue, label=L"\mathrm{Estimated}\;\beta")
    Plots.plot!(plt, xscale=:log10, yscale=:log10)
    _label!(plt, L"\lambda_\beta", L"\beta\;[\mathrm{MPa\,s\,m^{-1}}]")
    Plots.savefig(plt, joinpath(outdir, "beta_bounding_run$(run).pdf"))
    @info "Wrote β bounding figure to $outdir"
    return outdir
end
