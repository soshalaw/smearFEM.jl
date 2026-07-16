using HDF5
using Statistics
using ArgCheck
using Measures
using Plots
using LaTeXStrings
using smearFEM

global def_orange = RGB(245/255,118/255,0)
global def_blue = RGB(5/255,79/255,185/255)
global def_red = RGB(196/255,70/255,1/255)
global def_green = RGB(2/255,147/255,86/255)

const PLOT_CONFIG = Dict(
    :font_size => 11,
    :plot_height => 320,
    :plot_width => 480,
    :left_margin => 1pt,
    :right_margin => 1pt,
    :top_margin => -3pt
)

"""
    load_aruco_calibration(filepath)

Load the ArUco calibration data written by the validation script.

# Arguments:
- `filepath::String`: path to aruco_calibration.hdf5.

# Returns:
- `data::NamedTuple`: stage-marker separation time series and its statistics,
  together with the plate-marker scale validation against the reference molds.
"""
function load_aruco_calibration(filepath::String)

    if !isfile(filepath)
        throw(AssertionError("Calibration file not found: $filepath"))
    end

    return h5open(filepath, "r") do f
        gs = f["stage_separation"]
        gp = f["plate_separation"]

        (frame_index = read(gs["frame_index"]),
         separation  = read(gs["separation"]),
         mean_sep    = read_attribute(gs, "mean_mm"),
         std_sep     = read_attribute(gs, "std_mm"),
         fps         = read_attribute(f, "fps"),

         mold_a      = read(gp["mold_a"]),
         mold_b      = read(gp["mold_b"]),
         meas_diff   = read(gp["measured_difference_mm"]),
         true_diff   = read_attribute(gp, "true_height_difference_mm"),
         mean_diff   = read_attribute(gp, "measured_difference_mean_mm"),
         std_diff    = read_attribute(gp, "measured_difference_std_mm"),
         residual    = read_attribute(gp, "residual_mm"),
         scale       = read_attribute(gp, "scale_factor"))
    end
end

"""
    plot_aruco_calibration(filepath; savepath=nothing)

Plot the repeatability of the ArUco pose estimation.

The two stage markers are rigidly mounted, so their true separation is constant
throughout an experiment. The scatter in the recovered separation is therefore
attributable entirely to pose-estimation noise, and the absence of drift over
the recording indicates that no systematic degradation occurs as the specimen
deforms.

# Arguments:
- `filepath::String`: path to aruco_calibration.hdf5.
- `savepath::Union{Nothing,String}`: if given, the figure is written to this path.

# Returns:
- `plt::Plots.Plot`: the resulting figure.
- `data::NamedTuple`: the loaded calibration data, see [`load_aruco_calibration`](@ref).
"""
function plot_aruco_calibration(filepath::String)

    data = load_aruco_calibration(filepath)

    @argcheck !isempty(data.separation)

    t = data.frame_index ./ data.fps
    μ = data.mean_sep
    σ = data.std_sep

    μ_h = data.mean_diff
    σ_h = data.std_diff
    gt_h = data.true_diff

    savepath = resolve_data_path("experiments/physical_data/aruco_calibration/plots")

    set_file(dirname(savepath))
    plt = set_plot(PLOT_CONFIG[:font_size];
                   sz=(PLOT_CONFIG[:plot_width], PLOT_CONFIG[:plot_height]),
                   left_margin=PLOT_CONFIG[:left_margin],
                   right_margin=PLOT_CONFIG[:right_margin],
                   top_margin=PLOT_CONFIG[:top_margin],
                   legend_column=3)

    Plots.plot!(plt, t, fill(μ, length(t));
                ribbon    = σ,
                fillalpha = 0.20,
                lw        = 0,
                color     = def_blue,
                label     = L"\pm 1\sigma")
    Plots.plot!(plt, t, data.separation;
                lw    = 1.0,
                color = def_blue,
                label = L"d_{\mathrm{markers}}")
    Plots.hline!(plt, [μ];
                 ls    = :dash,
                 lc    = def_red,
                 lw    = 1.2,
                 label = L"\mathrm{mean}")
    Plots.ylims!(plt, 107, 109.1)
    Plots.Plots.xlims!(plt, 0.0, 10.0)
    Plots.xlabel!(plt, L"\mathrm{Time}\,[\mathrm{s}]")
    Plots.ylabel!(plt, L"d_{\mathrm{markers}}\;[\mathrm{mm}]")
    Plots.savefig(plt, joinpath(savepath,"stage_precision.pdf"))

    plt_h = set_plot(PLOT_CONFIG[:font_size];
                     sz=(PLOT_CONFIG[:plot_width], PLOT_CONFIG[:plot_height]),
                     left_margin=PLOT_CONFIG[:left_margin],
                     right_margin=PLOT_CONFIG[:right_margin],
                     top_margin=PLOT_CONFIG[:top_margin],
                     legend_column=4)
    Plots.plot!(plt_h, t, fill(μ_h, length(t));
                ribbon    = σ_h,
                fillalpha = 0.20,
                lw        = 0,
                color     = def_blue,
                label     = L"\pm 1\sigma")
    Plots.plot!(plt_h, t, data.meas_diff;
                lw    = 1.0,
                color = def_blue,
                label = L"h_{\mathrm{plate}}")
    Plots.hline!(plt_h, [μ_h];
                 ls    = :dash,
                 lc    = def_red,
                 lw    = 1.2,
                 label = L"\mathrm{mean}")
    Plots.plot!(plt_h, t, fill(gt_h, length(t));
                lw    = 1.0,
                color = def_red,
                label = L"\Delta h_{\mathrm{true}}")
    Plots.ylims!(plt_h, 9.0, 11.1)
    Plots.Plots.xlims!(plt_h, 0.0, 10.0)
    Plots.xlabel!(plt_h, L"\mathrm{Time}\,[\mathrm{s}]")
    Plots.ylabel!(plt_h, L"h_{\mathrm{plate}}\;[\mathrm{mm}]")
    Plots.savefig(plt_h, joinpath(savepath,"plate_accuracy.pdf"))

    drift = mean(data.separation[end-19:end]) - mean(data.separation[1:20])

    @info """ArUco pose estimation:
      stage separation : $(round(μ, digits=3)) ± $(round(σ, digits=3)) mm  (n=$(length(data.separation)))
      drift            : $(round(drift, digits=3)) mm
      plate scale      : $(round(data.scale, digits=4))  (residual $(round(data.residual, digits=3)) mm vs $(data.true_diff) mm)"""

    return plt, data
end

plot_aruco_calibration(resolve_data_path("ground_truth/physical_data/aruco_calibration/aruco_calibration.hdf5"))