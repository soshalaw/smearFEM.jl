# Validation metrics for the paper's experiments.
#
# These quantify how well a fit reproduces the observations; they are not part of the
# estimation method itself (that lives in `smearFEM`: `set_time_window`, `estimate_window`,
# `predict_window`). Kept separate from `pipeline.jl` so the pipeline reads as the paper's
# workflow rather than as workflow plus scoring.
#
# Included by `pipeline.jl`; `figures.jl` picks these up transitively through it. Both use
# every function here.

# Point count of one contour frame, whichever way round it is stored.
_npts(p) = max(size(p, 1), size(p, 2))

"""
    iteration_contour_metrics(data_type, filepath_gt, exp_path; view_folder, sim_view, write_results)
        -> Dict

Score every optimizer iterate against the observed contours.

Reads the per-iteration simulated contours written by `optimize` under
`<exp_path>/data/opt_data/iter_<n>/<sim_view>/2D_border_points` (produced when
`fit_model` is called with `store_border_pts=true`) and compares each against the
ground-truth contours, giving Chamfer and closest-point distances per frame, per iteration.

The Chamfer is reported twice: `chamfer` against the simulated contour as the optimizer
produced it, and `chamfer_up` against that contour resampled onto the observation's point
count. The gap between them is the sampling artefact — the two differ in nothing else — which is worth seeing rather than
silently correcting. Hausdorff distances are deliberately not reported: they are extremal,
so on this data they measure the interpolation choice as much as the fit.

Every metric is reported **squared**, in px², matching the optimizer's `ChamferCost` and
`ClosestPointCost` — which are squared for the same reason, Gauss-Newton needing `uᵀu`. So
the four share units outright and may be read against each other directly, unlike the
rooted forms `compare_pt_clouds` returns for three of them.

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
- `upsample_method::Symbol`: how the simulated contour is resampled for the upsampled
  Chamfer arm — `:spline` (default, the construction `fit_curve` uses for the observed
  contour) or `:linear`. See [`upsample_contour`](@ref); on real contours the two differ by
  up to ~1.2 px, so the choice is not cosmetic.
- `write_results::Bool`: also write the per-iteration series to
  `<exp_path>/data/opt_data/metrics/` (default: `true`).
- `plot_results::Bool`: also plot the time-averaged (frame-averaged) Chamfer (both arms)
  and closest-point distances against iteration, together on
  `<exp_path>/plots/contour_metrics_iter.pdf` (default: `true`). Each curve is
  normalized by its own value at the first iterate, so all three start at 1 and the axis
  reads as the fraction of the initial cost remaining; this affects the plot only, and the
  returned `Dict` and the written CSVs keep the raw px² distances. Squaring steepens every
  curve equally — a ratio of squares is the square of the ratio — so the ranking between
  metrics is unchanged and only the decay rate reads twice as fast on the log axis.
  Uses a log y-axis when every value is positive, falling back to linear otherwise.

# Returns
- `Dict` with `"iters"` (iteration indices) and, per iteration, `"chamfer_mean"`,
  `"chamfer_max"`, `"chamfer_up_mean"`, `"chamfer_up_max"` and `"closest_pt_mean"`, plus the
  per-frame series `"chamfer_frames"` and `"chamfer_up_frames"` (one vector per iteration).
  All in px². The `_up` entries are the upsampled arm; `closest_pt` is computed on the
  upsampled contour.
"""
function iteration_contour_metrics(data_type::String, filepath_gt::String, exp_path::String;
                                   view_folder::String="view_1", sim_view::String="view_1",
                                   upsample_method::Symbol=:spline,
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

    chamfer_mean    = Float64[]; chamfer_max    = Float64[]
    chamfer_up_mean = Float64[]; chamfer_up_max = Float64[]
    closest_pt_mean = Float64[]
    chamfer_frames    = Vector{Vector{Float64}}()
    chamfer_up_frames = Vector{Vector{Float64}}()

    for n in iters
        sim_pts, _, _ = read_csv(datapath(exp_path,"opt_data","iter_$n",sim_view,"2D_border_points"))
        nf = min(length(sim_pts), length(ObsDataList))
        nf > 0 || throw(ArgumentError("No overlapping frames for iteration $n"))
        obs = ObsDataList[1:nf]
        raw = sim_pts[1:nf]
        # The simulated contour carries ~160 mesh nodes against ~1300 observed points, and a
        # nearest-neighbour distance to those nodes overstates the distance to the contour
        # itself by an amount set purely by that sampling. `chamfer_up` resamples the coarse
        # side onto the observation's count and `chamfer` does not, so the gap between them is
        # that artefact alone. Metrics only — the optimizer's cost is untouched.
        up = [upsample_contour(s, _npts(o); method=upsample_method) for (s, o) in zip(raw, obs)]

        c    = [chamfer_sq_distance_kdtree(s, o) for (s, o) in zip(raw, obs)]
        c_up = [chamfer_sq_distance_kdtree(s, o) for (s, o) in zip(up, obs)]
        # `closest_point_distance_kdtree` is an RMSE; square it so every series is px² and
        # comparable with the Chamfer and with the optimizer's cost. RMSE² == mean(d²), so
        # this is the squared metric itself, not an approximation of it.
        cp  = [closest_point_distance_kdtree(s, o)^2 for (s, o) in zip(up, obs)]

        push!(chamfer_frames, c);             push!(chamfer_up_frames, c_up)
        push!(chamfer_mean, mean(c));         push!(chamfer_max, maximum(c))
        push!(chamfer_up_mean, mean(c_up));   push!(chamfer_up_max, maximum(c_up))
        push!(closest_pt_mean, mean(cp))
        @info "iter $n ($nf frames) [px²]: chamfer(mean)=$(round(mean(c), sigdigits=4)), " *
              "chamfer+upsampling(mean)=$(round(mean(c_up), sigdigits=4)), " *
              "closest_pt(mean)=$(round(mean(cp), sigdigits=4))"
    end

    metrics = Dict{String,Any}(
        "iters"              => iters,
        "chamfer_mean"      => chamfer_mean,
        "chamfer_max"       => chamfer_max,
        "chamfer_up_mean"   => chamfer_up_mean,
        "chamfer_up_max"    => chamfer_up_max,
        "closest_pt_mean"   => closest_pt_mean,
        "chamfer_frames"    => chamfer_frames,
        "chamfer_up_frames" => chamfer_up_frames,
    )

    if write_results
        for k in ("iters","chamfer_mean","chamfer_max","chamfer_up_mean","chamfer_up_max",
                  "closest_pt_mean")
            write_csv(datapath(exp_path,"opt_data","metrics",k), metrics[k])
        end
        # per-frame series: one file per iteration
        for (i, n) in enumerate(iters)
            write_csv(datapath(exp_path,"opt_data","metrics","chamfer_frames","iter_$n"), chamfer_frames[i])
            write_csv(datapath(exp_path,"opt_data","metrics","chamfer_up_frames","iter_$n"), chamfer_up_frames[i])
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
        chamfer_rel    = normalize_to_initial(chamfer_mean, "Chamfer")
        chamfer_up_rel = normalize_to_initial(chamfer_up_mean, "Chamfer + upsampling")
        closest_pt_rel = normalize_to_initial(closest_pt_mean, "closest-point")

        # Metrics decay over several decades as the fit converges, so prefer a log axis;
        # an exact zero (or a single iterate) would make that invalid.
        use_log = all(>(0), chamfer_rel) && all(>(0), chamfer_up_rel) && all(>(0), closest_pt_rel)
        yscale = use_log ? :log10 : :identity

        # The three metrics can sit almost on top of each other (a near-uniform contour
        # offset makes mean and max nearest-neighbour distance nearly equal), so vary the
        # line style as well as the colour to keep them separable.
        # No `²` on the individual labels: it is a property of all of them, so it belongs on
        # the axis, and on the Chamfer it would be a lie — the squared Chamfer is a mean of
        # squared distances, not the square of the Chamfer distance (see
        # `chamfer_sq_distance_kdtree`). The other two are genuine squares of their metrics.
        metric_plt = default_plot()
        Plots.plot!(metric_plt, iters, chamfer_rel, label=L"\mathrm{Chamfer}", marker=1,
                    linestyle=:solid, yscale=yscale, xminorgrid=:false,
                    legend=:outerbottom, legend_column=2)
        # The same metric with the simulated contour resampled onto the observation's point
        # count: the gap between these two curves is the sampling artefact, nothing else
        # about the fit differs between them.
        Plots.plot!(metric_plt, iters, chamfer_up_rel, label=L"\mathrm{Chamfer + upsampling}",
                    marker=1, linestyle=:dash, yscale=yscale, legend=:outerbottom, legend_column=2)
        Plots.plot!(metric_plt, iters, closest_pt_rel, label=L"\mathrm{Closest\;point}", marker=1,
                    linestyle=:dot, yscale=yscale, legend=:outerbottom, legend_column=2)
        _label!(metric_plt, L"\mathrm{Iterations}", L"d^{\imath}/d^{0}")
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
