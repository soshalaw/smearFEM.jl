using DelimitedFiles
using Plots
using Plots.PlotMeasures
using LaTeXStrings
using smearFEM

global def_orange = RGB(245/255,118/255,0)
global def_blue = RGB(5/255,79/255,185/255)
global def_red = RGB(196/255,70/255,1/255)
global def_green = RGB(2/255,147/255,86/255)

global y_lims = (0.7, 1.21)

const PLOT_CONFIG = Dict(
    :font_size => 11,
    :plot_height => 360,
    :plot_width => 330,
    :left_margin => 1pt,
    :right_margin => 1pt,
    :top_margin => -3pt
)

function plot_inv_mesh_convergence(filepath_res::String, filepath_gt::String)
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

    for (i,mesh_folder) in enumerate(mesh_dir)
        if !startswith(mesh_folder, "Hex2_")
            continue
        end

        path = joinpath(filepath_res, mesh_folder, "dt_0.1", "view_1", "data")
        exp_data = read_json(joinpath(path, "experiment_parameters.jld2"))
        cost_iter = readdlm(joinpath(path,"view_1","cost_iter.csv"), ',', Float64)
        β_iter = readdlm(joinpath(path, "view_1", "β.csv"), ',', Float64)
        η_iter = readdlm(joinpath(path, "view_1", "η.csv"), ',', Float64)

        ne = exp_data["ne_exp"]

        conv_cost_list[i] = abs(cost_iter[1]-cost_iter[end])/cost_iter[1]
        conv_β_list[i] = β_iter[end]/gt_beta
        conv_η_list[i] = η_iter[end]/gt_eta
        mesh_size[i] = ne
    end

    plot_path = joinpath(filepath_res, "post_analysis")
    set_file(plot_path)

    sorted_indices = sortperm(mesh_size)

    plt_cost = set_plot(PLOT_CONFIG[:font_size], 
                    sz=(PLOT_CONFIG[:plot_width], PLOT_CONFIG[:plot_height]), 
                    left_margin=PLOT_CONFIG[:left_margin],
                    right_margin=PLOT_CONFIG[:right_margin],
                    top_margin=PLOT_CONFIG[:top_margin],
                    legend_column=2)
    Plots.plot!(plt_cost, mesh_size[sorted_indices], conv_cost_list[sorted_indices], label=false, marker=:circle, markersize=4)
    Plots.hline!(plt_cost, [1.0], label=false, color=:black, linestyle=:dash)
    Plots.vline!(plt_cost, [6.0], label=L"\mathrm{h}_{\mathrm{exp}}", line=:dash, color=:red)
    Plots.ylims!(plt_cost, y_lims)
    Plots.xlabel!(plt_cost, "Effective element size (h)")
    Plots.ylabel!(plt_cost, "Relative cost reduction")
    Plots.savefig(plt_cost, joinpath(plot_path, "mesh_convergence_cost.pdf"))

    plt_beta = set_plot(PLOT_CONFIG[:font_size], 
                    sz=(PLOT_CONFIG[:plot_width], PLOT_CONFIG[:plot_height]), 
                    left_margin=PLOT_CONFIG[:left_margin],
                    right_margin=PLOT_CONFIG[:right_margin],
                    top_margin=PLOT_CONFIG[:top_margin],
                    legend_column=2)
    Plots.plot!(plt_beta, mesh_size[sorted_indices], conv_β_list[sorted_indices], label=false,  marker=:circle, markersize=4)
    Plots.hline!(plt_beta, [1.0], label=false, color=:black, linestyle=:dash)
    Plots.vline!(plt_beta, [6.0], label=L"\mathrm{h}_{\mathrm{exp}}", line=:dash, color=:red)
    Plots.ylims!(plt_beta, y_lims)
    Plots.xlabel!(plt_beta, "Effective element size (h)")
    Plots.ylabel!(plt_beta, L"\beta_{\mathrm{exp}}/\beta_{\mathrm{gt}}")
    Plots.savefig(plt_beta, joinpath(plot_path, "mesh_convergence_beta.pdf"))

    plt_eta = set_plot(PLOT_CONFIG[:font_size], 
                    sz=(PLOT_CONFIG[:plot_width], PLOT_CONFIG[:plot_height]), 
                    left_margin=PLOT_CONFIG[:left_margin],
                    right_margin=PLOT_CONFIG[:right_margin],
                    top_margin=PLOT_CONFIG[:top_margin],
                    legend_column=2)
    Plots.plot!(plt_eta, mesh_size[sorted_indices], conv_η_list[sorted_indices], label=false,  marker=:circle, markersize=4)
    Plots.hline!(plt_eta, [1.0], label=false, color=:black, linestyle=:dash)
    Plots.vline!(plt_eta, [6.0], label=L"\mathrm{h}_{\mathrm{exp}}", line=:dash, color=:red)
    Plots.ylims!(plt_eta, y_lims)
    Plots.xlabel!(plt_eta, "Effective element size (h)")
    Plots.ylabel!(plt_eta, L"\eta_{\mathrm{exp}}/\eta_{\mathrm{gt}}")
    Plots.savefig(plt_eta, joinpath(plot_path, "mesh_convergence_eta.pdf"))

end

function plot_inv_time_convergence(filepath_res::String, filepath_gt::String)
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

    for (i,dt_folder) in enumerate(dt_dir)
        if !startswith(dt_folder, "dt_")
            continue
        end

        path = joinpath(filepath_res, "Hex2_6", dt_folder, "view_1", "data")
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

    plt_cost = set_plot(PLOT_CONFIG[:font_size], 
                    sz=(PLOT_CONFIG[:plot_width], PLOT_CONFIG[:plot_height]), 
                    left_margin=PLOT_CONFIG[:left_margin],
                    right_margin=PLOT_CONFIG[:right_margin],
                    top_margin=PLOT_CONFIG[:top_margin],
                    legend_column=2)
    Plots.plot!(plt_cost, dt_size[sorted_indices], conv_cost_list[sorted_indices], label=false,  marker=:circle, markersize=4)
    Plots.hline!(plt_cost, [1.0], label=false, color=:black, linestyle=:dash)
    Plots.Plots.vline!(plt_cost, [0.1], label=L"\Delta t_{\mathrm{exp}}", line=:dash, color=:red)
    Plots.ylims!(plt_cost, y_lims)
    Plots.xlims!(plt_cost, (0, 1.1))
    Plots.xticks!(plt_cost, [0, 0.25, 0.5, 0.75, 1.0], ["0", "0.25", "0.5", "0.75", "1.0"])
    Plots.xlabel!(plt_cost, latexstring("Time step size \$(\\Delta t)\$"))
    Plots.ylabel!(plt_cost, "Relative cost reduction")
    Plots.savefig(plt_cost, joinpath(plot_path, "time_convergence_cost.pdf"))

    plt_beta = set_plot(PLOT_CONFIG[:font_size], 
                    sz=(PLOT_CONFIG[:plot_width], PLOT_CONFIG[:plot_height]), 
                    left_margin=PLOT_CONFIG[:left_margin],
                    right_margin=PLOT_CONFIG[:right_margin],
                    top_margin=PLOT_CONFIG[:top_margin],
                    legend_column=2)
    Plots.plot!(plt_beta, dt_size[sorted_indices], conv_β_list[sorted_indices], label=false,  marker=:circle, markersize=4)
    Plots.hline!(plt_beta, [1.0], label=false, color=:black, linestyle=:dash)
    Plots.vline!(plt_beta, [0.1], label=L"\Delta t_{\mathrm{exp}}", line=:dash, color=:red)
    Plots.ylims!(plt_beta, y_lims)
    Plots.xlims!(plt_beta, (0, 1.1))
    Plots.xticks!(plt_beta, [0, 0.25, 0.5, 0.75, 1.0], ["0", "0.25", "0.5", "0.75", "1.0"])
    Plots.xlabel!(plt_beta, latexstring("Time step size \$(\\Delta t)\$"))
    Plots.ylabel!(plt_beta, L"\beta_{\mathrm{est}}/\beta_{\mathrm{gt}}")
    Plots.savefig(plt_beta, joinpath(plot_path, "time_convergence_beta.pdf"))

    plt_eta = set_plot(PLOT_CONFIG[:font_size], 
                    sz=(PLOT_CONFIG[:plot_width], PLOT_CONFIG[:plot_height]), 
                    left_margin=PLOT_CONFIG[:left_margin],
                    right_margin=PLOT_CONFIG[:right_margin],
                    top_margin=PLOT_CONFIG[:top_margin],
                    legend_column=2)
    Plots.plot!(plt_eta, dt_size[sorted_indices], conv_η_list[sorted_indices], label=false,  marker=:circle, markersize=4)
    Plots.hline!(plt_eta, [1.0], label=false, color=:black, linestyle=:dash)
    Plots.vline!(plt_eta, [0.1], label=L"\Delta t_{\mathrm{exp}}", line=:dash, color=:red)
    Plots.ylims!(plt_eta, y_lims)
    Plots.xlims!(plt_eta, (0, 1.1))
    Plots.xticks!(plt_eta, [0, 0.25, 0.5, 0.75, 1.0], ["0", "0.25", "0.5", "0.75", "1.0"])
    Plots.xlabel!(plt_eta, latexstring("Time step size \$(\\Delta t)\$"))
    Plots.ylabel!(plt_eta, L"\eta_{\mathrm{est}}/\eta_{\mathrm{gt}}")
    Plots.savefig(plt_eta, joinpath(plot_path, "time_convergence_eta.pdf"))
end

function setup_plot_config(; left_margin=nothing, right_margin=nothing, top_margin=nothing)
    """Configure plot with standard settings."""
    figsize = (PLOT_CONFIG[:plot_width], PLOT_CONFIG[:plot_height])
    lm = left_margin !== nothing ? left_margin : PLOT_CONFIG[:left_margin]
    rm = right_margin !== nothing ? right_margin : PLOT_CONFIG[:right_margin]
    tm = top_margin !== nothing ? top_margin : PLOT_CONFIG[:top_margin]
    
    return set_plot(PLOT_CONFIG[:font_size], 
                          sz=figsize, 
                          left_margin=lm,
                          right_margin=rm,
                          top_margin=tm,
                          legend_column=2)
end

function collect_experiment_groups(filepath)
    groups = ExpGroup[]
    for elem_size_folder in readdir(filepath)
        elem_size_folder == "post_analysis" && continue

        for sim_time_folder in readdir(joinpath(filepath, elem_size_folder))
            startswith(sim_time_folder, "simtime_") || continue
            noise_base = joinpath(filepath, elem_size_folder, sim_time_folder)

            for noise_folder in readdir(noise_base)
                startswith(noise_folder, "noise_") || continue
                step_path = joinpath(noise_base, noise_folder)

                leaves = ExpLeaf[]
                for step_folder in readdir(step_path)
                    startswith(step_folder, "dt_") || continue
                    dt_path = joinpath(step_path, step_folder)

                    for view_folder in readdir(dt_path)
                        startswith(view_folder, "view_") || continue
                        push!(leaves, ExpLeaf(step_folder, view_folder, joinpath(dt_path, view_folder)))
                    end
                end
                push!(groups, ExpGroup(elem_size_folder, sim_time_folder, noise_folder, step_path, leaves))
            end
        end
    end
    return groups
end

function main()
    # Define the paths to the result and ground truth directories
    filepath_res = resolve_data_path("experiments/sim_data/convergence_analysis/stokes_convergence/experiment_mesh_conv/cylinder/3")
    filepath_gt = joinpath(resolve_data_path("ground_truth/sim_data/Stokes/force/constant/Hex_2/convergence_analysis/experiment_mesh_convergence_analysis/"))

    # Call the function to plot convergence
    plot_inv_mesh_convergence(filepath_res, filepath_gt)
    plot_inv_time_convergence(filepath_res, filepath_gt)
    
end

main()