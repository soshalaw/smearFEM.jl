using DelimitedFiles
using Plots
using Plots.PlotMeasures
using LaTeXStrings
using smearFEM


global y_lims = (10^-0.32, 10^0.32)


const PLOT_CONFIG = Dict(
    :font_size => 11,
    :plot_height => 320,
    :plot_width => 480,
    :left_margin => 1pt,
    :right_margin => 1pt,
    :top_margin => -3pt
)

# Result directories gained an `<opt_method>` level between `dt_*` and `view_*` on
# 2026-08-26. Resolve either layout so this script reads pre- and post-migration trees.
function _view_path(dt_path::String, opt_method::String, view::String="view_1")
    with_method = joinpath(dt_path, opt_method, view)
    isdir(with_method) && return with_method
    return joinpath(dt_path, view)   # pre-2026-08-26 layout
end

function plot_inv_mesh_convergence(filepath_res::String, filepath_gt::String; opt_method::String="gn")
    # Load the data from the CSV file
    _mesh_dir = readdir(filepath_res)
    mesh_dir = filter(x -> startswith(x, "Hex2_"), _mesh_dir)
    conv_η_list = Vector{Float64}(undef, length(mesh_dir))
    conv_β_list = Vector{Float64}(undef, length(mesh_dir))
    conv_cost_list = Vector{Float64}(undef, length(mesh_dir))
    mesh_size = Vector{Float64}(undef, length(mesh_dir))
    
    gt_data = read_json(joinpath(filepath_gt, "data", "sim_params.jld2"))
    gt_eta = gt_data["η"][1]
    gt_beta = gt_data["β"][1]
    r = gt_data["r"]
    h = gt_data["h"]
    V = π*r^2*h
    exp_h = (V / 192)^(1/3)

    for (i,mesh_folder) in enumerate(mesh_dir)
        if !startswith(mesh_folder, "Hex2_")
            continue
        end

        path = joinpath(_view_path(joinpath(filepath_res, mesh_folder, "dt_0.1"), opt_method), "data")
        exp_data = read_json(joinpath(path, "experiment_parameters.jld2"))
        cost_iter = readdlm(joinpath(path, "cost_iter.csv"), ',', Float64)
        β_iter = readdlm(joinpath(path, "β.csv"), ',', Float64)
        η_iter = readdlm(joinpath(path, "η.csv"), ',', Float64)

        ne = exp_data["num_ne"]
        effective_h = (V/ne)^(1/3)

        conv_cost_list[i] = abs(cost_iter[1]-cost_iter[end])/cost_iter[1]
        conv_β_list[i] = β_iter[end]/gt_beta
        conv_η_list[i] = η_iter[end]/gt_eta
        mesh_size[i] = effective_h
    end

    plot_path = joinpath(filepath_res, "post_analysis")
    set_file(plot_path)

    sorted_indices = sortperm(mesh_size)

    plt_cost = set_plot_from_config(PLOT_CONFIG)
    Plots.plot!(plt_cost, mesh_size[sorted_indices], conv_cost_list[sorted_indices], label=false, marker=:circle, markersize=4, xscale=:log10, yscale=:log10)
    Plots.hline!(plt_cost, [1.0], label=false, color=:black, linestyle=:dash)
    Plots.vline!(plt_cost, [exp_h], label=L"\mathrm{h}_{\mathrm{exp}}", line=:dash, color=:red)
    Plots.ylims!(plt_cost, y_lims)
    Plots.xlims!(plt_cost, 1,100)
    Plots.xlabel!(plt_cost, "Effective element size (h)")
    Plots.ylabel!(plt_cost, L"\Delta d_{\mathrm{rel}}")
    Plots.savefig(plt_cost, joinpath(plot_path, "mesh_convergence_cost.pdf"))

    plt_beta = set_plot_from_config(PLOT_CONFIG)
    Plots.plot!(plt_beta, mesh_size[sorted_indices], conv_β_list[sorted_indices], label=false,  marker=:circle, markersize=4, xscale=:log10, yscale=:log10)
    Plots.hline!(plt_beta, [1.0], label=false, color=:black, linestyle=:dash)
    Plots.vline!(plt_beta, [exp_h], label=L"\mathrm{h}_{\mathrm{exp}}", line=:dash, color=:red)
    Plots.ylims!(plt_beta, y_lims)
    Plots.xlims!(plt_beta, 1,100)
    Plots.xlabel!(plt_beta, "Effective element size (h)")
    Plots.ylabel!(plt_beta, L"\beta_{\mathrm{exp}}/\beta_{\mathrm{gt}}")
    Plots.savefig(plt_beta, joinpath(plot_path, "mesh_convergence_beta.pdf"))

    plt_eta = set_plot_from_config(PLOT_CONFIG)
    Plots.plot!(plt_eta, mesh_size[sorted_indices], conv_η_list[sorted_indices], label=false,  marker=:circle, markersize=4, xscale=:log10, yscale=:log10)
    Plots.hline!(plt_eta, [1.0], label=false, color=:black, linestyle=:dash)
    Plots.vline!(plt_eta, [exp_h], label=L"\mathrm{h}_{\mathrm{exp}}", line=:dash, color=:red)
    Plots.ylims!(plt_eta, y_lims)
    Plots.xlims!(plt_eta, 1,100)
    Plots.xlabel!(plt_eta, "Effective element size (h)")
    Plots.ylabel!(plt_eta, L"\eta_{\mathrm{exp}}/\eta_{\mathrm{gt}}")
    Plots.savefig(plt_eta, joinpath(plot_path, "mesh_convergence_eta.pdf"))

end

function plot_inv_time_convergence(filepath_res::String, filepath_gt::String; opt_method::String="gn")
     # Load the data from the CSV file
    _dt_dir = readdir(joinpath(filepath_res, "Hex2_6"))
    dt_dir = filter(x -> startswith(x, "dt_"), _dt_dir)
    conv_η_list = Vector{Float64}(undef, length(dt_dir))
    conv_β_list = Vector{Float64}(undef, length(dt_dir))
    conv_cost_list = Vector{Float64}(undef, length(dt_dir))
    dt_size = Vector{Float64}(undef, length(dt_dir))

    gt_data = read_json(joinpath(filepath_gt, "data", "sim_params.jld2"))
    gt_eta = gt_data["η"][1]
    gt_beta = gt_data["β"][1]
    r = gt_data["r"]
    h = gt_data["h"]
    V = π*r^2*h
    
    for (i,dt_folder) in enumerate(dt_dir)
        if !startswith(dt_folder, "dt_")
            continue
        end

        path = joinpath(_view_path(joinpath(filepath_res, "Hex2_6", dt_folder), opt_method), "data")
        exp_data = read_json(joinpath(path, "experiment_parameters.jld2"))
        cost_iter = readdlm(joinpath(path, "cost_iter.csv"), ',', Float64)
        β_iter = readdlm(joinpath(path, "β.csv"), ',', Float64)
        η_iter = readdlm(joinpath(path, "η.csv"), ',', Float64)

        dt = exp_data["dt"]

        conv_cost_list[i] = abs(cost_iter[1]-cost_iter[end])/cost_iter[1]
        conv_β_list[i] = β_iter[end]/gt_beta
        conv_η_list[i] = η_iter[end]/gt_eta
        dt_size[i] = dt
    end

    plot_path = joinpath(filepath_res, "post_analysis")
    set_file(plot_path)

    sorted_indices = sortperm(dt_size)

    plt_cost = set_plot_from_config(PLOT_CONFIG)
    Plots.plot!(plt_cost, dt_size[sorted_indices], conv_cost_list[sorted_indices], label=false,  marker=:circle, markersize=4, xscale=:log10, yscale=:log10)
    Plots.hline!(plt_cost, [1.0], label=false, color=:black, linestyle=:dash)
    Plots.vline!(plt_cost, [0.1], label=L"\Delta t_{\mathrm{exp}}", line=:dash, color=:red)
    Plots.ylims!(plt_cost, y_lims)
    Plots.xlims!(plt_cost, 0.08, 1.2)
    Plots.xlabel!(plt_cost, latexstring("Time step size \$(\\Delta t)\$"))
    Plots.ylabel!(plt_cost, L"\Delta d_{\mathrm{rel}}")
    Plots.savefig(plt_cost, joinpath(plot_path, "time_convergence_cost.pdf"))

    plt_beta = set_plot_from_config(PLOT_CONFIG)
    Plots.plot!(plt_beta, dt_size[sorted_indices], conv_β_list[sorted_indices], label=false,  marker=:circle, markersize=4, xscale=:log10, yscale=:log10)
    Plots.hline!(plt_beta, [1.0], label=false, color=:black, linestyle=:dash)
    Plots.vline!(plt_beta, [0.1], label=L"\Delta t_{\mathrm{exp}}", line=:dash, color=:red)
    Plots.ylims!(plt_beta, y_lims)
    Plots.xlims!(plt_beta, 0.08, 1.2)
    Plots.xlabel!(plt_beta, latexstring("Time step size \$(\\Delta t)\$"))
    Plots.ylabel!(plt_beta, L"\beta_{\mathrm{est}}/\beta_{\mathrm{gt}}")
    Plots.savefig(plt_beta, joinpath(plot_path, "time_convergence_beta.pdf"))

    plt_eta = set_plot_from_config(PLOT_CONFIG)
    Plots.plot!(plt_eta, dt_size[sorted_indices], conv_η_list[sorted_indices], label=false,  marker=:circle, markersize=4, xscale=:log10, yscale=:log10)
    Plots.hline!(plt_eta, [1.0], label=false, color=:black, linestyle=:dash)
    Plots.vline!(plt_eta, [0.1], label=L"\Delta t_{\mathrm{exp}}", line=:dash, color=:red)
    Plots.ylims!(plt_eta, y_lims)
    Plots.xlims!(plt_eta, 0.08, 1.2)
    Plots.xlabel!(plt_eta, latexstring("Time step size \$(\\Delta t)\$"))
    Plots.ylabel!(plt_eta, L"\eta_{\mathrm{est}}/\eta_{\mathrm{gt}}")
    Plots.savefig(plt_eta, joinpath(plot_path, "time_convergence_eta.pdf"))
end


"""
    main()

Plot the inverse mesh- and time-convergence results. Invoked automatically
when this file is run as a script (`julia inverse_convergence.jl`), but not
when it is `include`d.
"""
function main()
    filepath_res = resolve_data_path("experiments/sim_data/convergence_analysis/stokes_convergence/experiment_mesh_conv/cylinder/4")
    # Scenario A.IV (η = 100 kPa s, β = 100 MPa s m⁻¹); only sim_params (η, β, r, h) is read from here.
    filepath_gt = resolve_data_path("ground_truth/sim_data/Stokes/force/constant/Hex2_16/cylinder/4/")

    plot_inv_mesh_convergence(filepath_res, filepath_gt)
    plot_inv_time_convergence(filepath_res, filepath_gt)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
