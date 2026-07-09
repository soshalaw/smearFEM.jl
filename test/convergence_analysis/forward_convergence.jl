# Convergence analysis for 3D Stokes FEM solutions
# Analyzes both mesh convergence and time integration convergence

using smearFEM
using Plots
using LinearAlgebra
using DelimitedFiles
using CSV
using DataFrames
using Plots.PlotMeasures
using LaTeXStrings
using Dates
using Statistics
using Printf
using CurveFit
using Dierckx

# Plot configuration constants
const PLOT_CONFIG = Dict(
    :font_size => 11,
    :plot_height => 320,
    :plot_width => 480,
    :left_margin => 1pt,
    :right_margin => 1pt,
    :top_margin => -3pt
)

function get_camera_matrix()
    """Returns camera matrix for rendering."""
    return [[8*2048/7.07, 0.0, 2048/2] [0.0, 8*1536/5.3, 1536/2] [0.0, 0.0, 1.0]]'
end


function get_r(filepath::String)   
    if !isdir(filepath)
        throw(SystemError("Trying to read from $filepath, the directory does not exist."))
    end

    csv_files = readdir(filepath, join=true)        # get the list of the csv files in the directory
    border_r_list = AbstractArray[]
    border_z_list = AbstractArray[]
    for file in csv_files
        if !endswith(file, ".csv")  # check if the file is a CSV file
            continue
        end
        data = readdlm(file, ',', Float64, '\n', header=false)  # read the observation data
        data_sorted = data[sortperm(data[:, 3]), :]  # sort the data based on the z coordinate (3rd column)
        r = sqrt.(data_sorted[:, 1].^2 + data_sorted[:, 2].^2)  
        
        z_prev = data_sorted[1, 3]
        idx_list = [1]
        z_list = [z_prev]
        i = 1
        for zp in data_sorted[:, 3]
            if round(Int, z_prev) < round(Int, zp)
                push!(z_list, zp)
                push!(idx_list, i)
                z_prev = zp
            end
            i = i + 1
        end
        border_r = r[idx_list]
        push!(border_r_list, border_r)
        push!(border_z_list, z_list)
    end

    return border_r_list, border_z_list
end

function get_r_curves(filepath::String)
    csv_list = readdir(filepath, join=true)
    r_curves = []
    for file in csv_list
        name = basename(file)
        if !startswith(name, "free_curve") || !endswith(name, ".csv")
            continue
        end
        data = readdlm(file, ',', '\n')
        z = data[2:end, 1] .- minimum(data[2:end, 1])  # Normalize z to start from 0
        r = data[2:end, 2]
        idx = sortperm(z)
        spl = Spline1D(z[idx], r[idx]; k=1, bc="extrapolate")  # Cubic spline interpolation, clamp to boundary value outside range
        push!(r_curves, (z=z[idx], r=r[idx], spline=spl, file=name))
    end
    return r_curves
end

# Plotting utilities
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

function fit_convergence_rate(x_vals::Vector, y_vals::Vector)
    """Fit convergence rate using linear regression in log-log space.
    
    Returns: (intercept, convergence_rate)
    """
    log_x = log.(x_vals)
    log_y = log.(y_vals)
    A = [ones(length(log_x)) log_x]
    coeffs = A \ log_y
    return coeffs[1], coeffs[2]
end

"""
    plot_convergence_generic(x_vals, y_vals, x_label, x_label_latex=""; int_tol=0.2, shift=5.0)

Plot convergence data with a fitted power-law reference line.

# Arguments
- `x_vals::Vector`: x-axis values (h for mesh, Δt for time).
- `y_vals::Vector`: error values.
- `x_label::String`: label for the x-axis.
- `x_label_latex::String`: LaTeX symbol for the x-axis quantity, used to
  label the reference-line order (e.g. `"O(h^2)"`) (default: `""`).

# Keyword Arguments
- `int_tol::Float64`: if the fitted rate is within `int_tol` of the nearest
  integer, draw a single reference line at that integer; otherwise draw both
  the floor and ceil bounding lines (default: `0.2`).
- `shift::Float64`: multiplicative offset separating the reference line from
  the data at its anchor point (x_min for the single/floor line, x_max for
  the ceil line), so the reference line runs offset from the data rather
  than through it (default: `5.0`).

# Returns
- `plt`: the convergence plot.
"""
function plot_convergence_generic(x_vals::Vector, y_vals::Vector, x_label::String, x_label_latex::String="",;
                                  int_tol::Float64=0.2, shift::Float64=5.0)
    # Filter valid data
    valid_idx = findall(x -> !isnan(x) && !isinf(x) && x > 0, y_vals)
    y_vals_clean = y_vals[valid_idx]
    x_vals_clean = x_vals[valid_idx]

    if length(y_vals_clean) < 2
        @warn "Not enough valid data points for plotting"
        return
    end

    # Fit convergence rate
    intercept, convergence_rate = fit_convergence_rate(x_vals_clean, y_vals_clean)
    rate = convergence_rate

    x_min = minimum(x_vals_clean)
    x_max = maximum(x_vals_clean)
    fit_at(x0) = exp(intercept + rate * log(x0))

    ref_y_max = maximum(y_vals_clean)
    ref_y_min = minimum(y_vals_clean)

    ref_line(m, xs, C) = C .* xs .^ m
    get_x(y, m, C) = (y / C)^(1/m)
    # x-range that makes ref_line(m, ., C) span exactly [ref_y_min, ref_y_max] for this m, C.
    # For m == 0 the line is constant in y, so it can't be inverted to span a y-range —
    # just use the data's own x-range instead.
    x_range_for(m, C) = m == 0 ? exp.(range(log(x_min), log(x_max), length=100)) :
                                  exp.(range(log(get_x(ref_y_min, m, C)), log(get_x(ref_y_max, m, C)), length=100))

    conv_rate_str = @sprintf("%.3f", rate)
    plt = setup_plot_config()

    Plots.plot!(plt, x_vals_clean, y_vals_clean,
                label="rRMSE  ", mode="markers",
                xlabel=x_label, ylabel="Relative RMSE",
                marker=:circle, markersize=4, markerstrokewidth=1.5, color="#FF7F0E",
                yscale=:log10, xscale=:log10)
    # if rate == 0
    #     @warn "Convergence rate is 0 — data may not be converging; skipping reference line"
    # else
    #     nearest_int = round(Int, rate)
    #     slopes = abs(rate - nearest_int) <= int_tol ? [nearest_int] : [floor(Int, rate), ceil(Int, rate)]
    #     label_for(m) = if m == 0
    #         latexstring("O(1)")
    #     elseif m == 1
    #         latexstring("O(\$$(x_label_latex)\$)")
    #     else
    #         latexstring("O(\$$(x_label_latex)^{$m}\$)")
    #     end
    #     if length(slopes) == 1
    #         m = slopes[1]
    #         C = (fit_at(x_max) * 4) / x_max^m
    #         x_ref = x_range_for(m, C)
    #         Plots.plot!(plt, x_ref, ref_line(m, x_ref, C),
    #                     label=label_for(m),
    #                     line=2, linewidth=2.5, color=:teal,
    #                     yscale=:log10, xscale=:log10, linestyle=:dash)
    #     else
    #         # Floor-order line offset below the data (anchored near x_min), ceil-order line
    #         # offset above (anchored near x_max) — each spans exactly [ref_y_min, ref_y_max]
    #         # for its own slope, so the x-range differs between the two lines
    #         m_left, m_right = slopes[1], slopes[2]
    #         C_left = (fit_at(x_min) / 1.5) / x_min^m_left
    #         C_right = (fit_at(x_max) * 10) / x_max^m_right
    #         x_ref_left = x_range_for(m_left, C_left)
    #         x_ref_right = x_range_for(m_right, C_right)
    #         Plots.plot!(plt, x_ref_left, ref_line(m_left, x_ref_left, C_left),
    #                     label=label_for(m_left),
    #                     line=2, linewidth=2.5, color=:teal,
    #                     yscale=:log10, xscale=:log10, linestyle=:dash)
    #         Plots.plot!(plt, x_ref_right, ref_line(m_right, x_ref_right, C_right),
    #                     label=label_for(m_right),
    #                     line=2, linewidth=2.5, color=:teal,
    #                     yscale=:log10, xscale=:log10, linestyle=:dot)
    #     end
    # end
    @info "Convergence rate: O(h^$(conv_rate_str))"

    # Save plot
    return plt
end


function mesh_convergence_analysis(;radius::Float64=25.0, height::Float64=40.0, elem_sizes::Vector=[10, 8, 6, 4], 
                                  template_mesh_geo_path::String=joinpath(@__DIR__, "mesh.geo"))

    height_list = AbstractArray[]
    height_error_list = Float64[]
    rad_error_list = Float64[]
    effective_element_size_list = Float64[]
    μ_list = AbstractArray[]
    time_list = []

    β::Float64 = 100.0 # penalty parameter for the ground truth
    η::Float64 = 100.0

    viscosity_type::String = "constant" # "constant" or "bulk_viscosity"
    element_shape_x::Symbol = :Hex
    basis_order_x::Int = 2
    element_shape_u::Symbol = :Hex
    basis_order_u::Int = 2
    element_shape_p::Symbol = :Hex
    basis_order_p::Int = 1
    control::String = "force" # "force" or "velocity"

    geometry::Symbol = :cylinder # geometry type for the mesh generation

    gt_path = resolve_data_path(joinpath("ground_truth", "sim_data", "Stokes", "force", "constant", "Hex_2", "convergence_analysis", "mesh_convergence_analysis", "mesh_convergence_tfem_data"))

    F_ext::Float64 = 9518.61 # force applied to the cylinder in N
    sim_time::Float64 = 5.0 # simulation time in seconds
    step_size = 2.5e-3 # time step size in seconds
    steps = round(Int, sim_time/step_size)  

    obj_pose = [150, 0.0, height/2]
    camera_matrix = get_camera_matrix()
    nz_list = Union{Int,Float64}[2, 4, 6, 8, 10, 12, 14, 16] # number of elements for each mesh size
    # nz_list = Union{Int,Float64}[20, 14, 10, 9, 6, 3, 2] # number of elements for each mesh size
    volume = π*radius^2*height # approximate volume of the cylinder divided by number of elements for the coarsest mesh
    
    h_ref = readdlm(joinpath(gt_path, "data", "h.csv"), ',', Float64, '\n', header=false)[end] # reference height from the finest mesh solution
    r_ref = get_r_curves(joinpath(gt_path, "data"))[end]

    for (iter_index,nz) in enumerate(nz_list)
        filepath = resolve_data_path(joinpath("ground_truth", "sim_data", "Stokes", control, viscosity_type, "$(element_shape_x)_$(basis_order_x)", "convergence_analysis", "mesh_convergence_analysis", "mesh_sz_$nz"))
        # if nz == 4 || nz == 10 || nz == 12 
        #     continue
        # end
        # set_file(filepath)  # create the directory if it doesn't exist

        # Extract element size from directory name
        println("Running for element size = $nz")
        
        exp_params = Dict(
                "element_shape_u" => element_shape_u,
                "basis_order_u"   => basis_order_u,
                "element_shape_p" => element_shape_p,
                "basis_order_p"   => basis_order_p,
                "element_shape_x" => element_shape_x,
                "basis_order_x"   => basis_order_x,
                "ne_gt" => nz,
                "β_gt" => β,
                "η_gt" => η,
                "filepath_gt" => filepath,
                "control" => control,
                "viscosity_type" => viscosity_type,
                "obj_pose_gt" => obj_pose,
                "F_ext" => F_ext,
                "sim_time_gt" => sim_time,
                "steps_gt" => steps,
                "r" => radius,
                "h" => height,
                "camera_matrix" => camera_matrix,
                "animate" => false,
                "geometry" => geometry
            )
        
        try
            # write_gt_data(exp_params)
        catch e
            @error "Simulation failed for element size $nz" exception=(e, catch_backtrace())
            continue
        end
        
        # Read results
        try
            exp_params = read_json(joinpath(filepath, "data", "sim_params.jld2"))
            h_mesh = readdlm(joinpath(filepath, "data","h.csv"), ',', Float64)
            elapsed_time = readdlm(joinpath(filepath, "data","avg_time.csv"), ',', Float64)
            r, z = get_r(joinpath(filepath, "data", "sim_data", "surface_nodes"))

            in_range = (z[end] .>= minimum(r_ref.z)) .& (z[end] .<= maximum(r_ref.z))
            z_cmp = z[end][in_range]
            r_cmp = r[end][in_range]
            r_int = r_ref.spline(z_cmp)  # Evaluate reference curve at the z values of the current mesh

            ne = exp_params["ne"]
            effective_element_size = (volume / ne)^(1/3)

            δh = sqrt(mean(((h_mesh[end] - h_ref) / h_ref).^2))
            δr = sqrt(mean(((r_cmp - r_int) ./ r_int).^2))

            path =resolve_data_path("experiments/sim_data/convergence_analysis/stokes_convergence/mesh_convergence_analysis/plots")
            Plots.plot(z[end], r[end], label="Current mesh", xlabel="z", ylabel="r", title="Radius vs Height for element size $nz")
            Plots.plot!(z_cmp, r_int, label="Reference mesh", linestyle=:dash)
            Plots.savefig(joinpath(path, "radius_vs_height_mesh_sz_$nz.pdf"))

            println("Height error: ", δh)
            println("Radius error: ", δr)

            push!(effective_element_size_list, effective_element_size)
            push!(height_list, h_mesh)
            push!(height_error_list, δh)
            push!(rad_error_list, δr)
            push!(time_list, elapsed_time)
            iter_index += 1
        catch e
            @error "Simulation failed for element size $nz" exception=(e, catch_backtrace())
            continue
        end
    end
    write_csv(resolve_data_path("experiments/sim_data/convergence_analysis/stokes_convergence/mesh_convergence_analysis/effective_element_size"), effective_element_size_list)
    write_csv(resolve_data_path("experiments/sim_data/convergence_analysis/stokes_convergence/mesh_convergence_analysis/height_list"), height_list)
    write_csv(resolve_data_path("experiments/sim_data/convergence_analysis/stokes_convergence/mesh_convergence_analysis/height_error_list"), height_error_list)
    write_csv(resolve_data_path("experiments/sim_data/convergence_analysis/stokes_convergence/mesh_convergence_analysis/rad_error_list"), rad_error_list)
    write_csv(resolve_data_path("experiments/sim_data/convergence_analysis/stokes_convergence/mesh_convergence_analysis/elem_sizes"), elem_sizes)
    write_csv(resolve_data_path("experiments/sim_data/convergence_analysis/stokes_convergence/mesh_convergence_analysis/time_list"), time_list)
end

function time_intergration_convergence_analysis(; radius::Float64=25.0, height::Float64=40.0, dt_list::Vector=[2, 1, 0.5, 0.1, 0.05, 0.01, 0.005, 0.001], 
                                  template_mesh_geo_path::String=joinpath(@__DIR__, "mesh.geo"))    
    height_list = AbstractArray[]
    final_height_list = Float64[]
    height_error_list = Float64[]
    effective_element_size_list = Float64[]
    _dt_list = []

    β::Float64 = 100.0 # penalty parameter for the ground truth
    η::Float64 = 100.0

    viscosity_type::String = "constant" # "constant" or "bulk_viscosity"
    FunctionClass_x::String = "Q2" # Function space for the ground truth
    FunctionClass_p::String = "Q1"
    FunctionClass_u::String = "Q2"
    control::String = "force" # "force" or "velocity"

    F_ext::Float64 = 9518.61 # force applied to the cylinder in N
    sim_time::Float64 = 5.0 # simulation time in seconds

    obj_pose = get_object_pose(height)
    camera_matrix = get_camera_matrix()
    filepath = resolve_data_path("ground_truth/sim_data/Stokes/$control/$viscosity_type/$FunctionClass_x/convergence_analysis")
    gt_path = resolve_data_path("ground_truth","sim_data","Stokes","force","constant","Hex_2","convergence_analysis","mesh_convergence_analysis","mesh_convergence_tfem_data")
    
    volume = π*radius^2*height # approximate volume of the cylinder divided by number of elements for the coarsest mesh
    
    mesh_dir = dirname(filepath)
    set_file(mesh_dir)  # create the directory if it doesn't exist
    
    iter_index = 1
    
    mesh_filepath = joinpath(mesh_dir, "convergence_analysis", "mesh_convergence_analysis", "mesh_4.0")
    
    for (idx,step_size) in enumerate(Float64.(reverse(dt_list)))

        filepath = joinpath(mesh_dir, "convergence_analysis", "time_integration_convergence_analysis","step_$step_size", "simulation")
        println("Running for time step size = $step_size")
        steps = round(Int, sim_time/step_size)
        
        # Run simulation with the generated mesh
        exp_params = Dict("FunctionClass_x" => FunctionClass_x, "FunctionClass_u" => FunctionClass_u, "FunctionClass_p" => FunctionClass_p, 
                         "ne_gt" => 4, "β_gt" => β, "η_gt" => η, "filepath_gt" => filepath, 
                         "control" => control, "viscosity_type" => viscosity_type, "obj_pose_gt" => obj_pose, 
                         "F_ext" => F_ext, "sim_time_gt" => sim_time, "steps_gt" => steps, 
                         "r" => radius, "h" => height, "camera_matrix" => camera_matrix, "animate" => true, "mesh_path" => mesh_filepath)

        try
            write_gt_data(exp_params)
        catch e
            @error "Simulation failed for time step size $step_size" exception=e
            continue
        end
        
        # Read results
        try
            exp_params = read_json(joinpath(filepath, "data", "sim_params.jld2"))
            h_mesh = readdlm(joinpath(filepath, "data","h.csv"), ',', Float64)
            dt = exp_params["time_steps"]
            ne = exp_params["ne"]
            effective_element_size = (volume / ne)^(1/3)
            border_r = get_r(joinpath(filepath, "data", "sim_data", "surface_nodes"))

            # if idx == 1
            #     h_ref = h_mesh[end] # use the finest time step solution as reference for error calculation
            #     continue
            # end
            
            δh = abs(h_mesh[end] - h_ref) / h_ref
            _δh = (h_mesh[end] - h_ref) / h_ref
            
            println("Height error: ", _δh)
            push!(effective_element_size_list, effective_element_size)
            push!(height_list, h_mesh)
            push!(height_error_list, δh)
            push!(_dt_list, dt)
            push!(final_height_list, h_mesh[end])
            push!(border_r_list, border_r)
            iter_index += 1
        catch e
            @warn "Failed to read results for time step size $step_size: $e"
        end
    end
    println("Final height list: ", final_height_list)
    write_csv(resolve_data_path("experiments/sim_data/convergence_analysis/stokes_convergence/time_convergence_analysis/effective_element_size"), effective_element_size_list)
    write_csv(resolve_data_path("experiments/sim_data/convergence_analysis/stokes_convergence/time_convergence_analysis/height_list"), height_list)
    write_csv(resolve_data_path("experiments/sim_data/convergence_analysis/stokes_convergence/time_convergence_analysis/height_error_list"), height_error_list)
    write_csv(resolve_data_path("experiments/sim_data/convergence_analysis/stokes_convergence/time_convergence_analysis/t_steps"), _dt_list)
    write_csv(resolve_data_path("experiments/sim_data/convergence_analysis/stokes_convergence/time_convergence_analysis/final_height_list"), final_height_list)
end

function plot_convergence_mesh(file_path::String)
    """Plot mesh convergence analysis results."""
    try
        radius = 25.0
        height = 40.0
        height_error_list = readdlm(joinpath(file_path, "height_error_list.csv"), ',', Float64)
        height_list_ = readdlm(joinpath(file_path, "height_list.csv"), ',', Float64)
        elem_sizes = readdlm(joinpath(file_path, "effective_element_size.csv"), ',', Float64)
        time_list_ = readdlm(joinpath(file_path, "time_list.csv"), ',', Float64)
        rad_error_list_ = readdlm(joinpath(file_path, "rad_error_list.csv"), ',', Float64)

        list_end = size(height_error_list, 1)
        relative_error = vec(height_error_list[1:(list_end), end])
        rad_error_list = vec(rad_error_list_[1:(list_end), end])
        height_list = vec(height_list_[1:(list_end), end])
        time_list = vec(time_list_[1:(list_end), end])
        # println(time_list)
        elem_sizes_flat = vec(elem_sizes[1:(list_end), 1])
        
        selected_mesh_size =  (π*radius^2*height / 192)^(1/3) # effective element size for the experiment mesh.
        @info "Selected mesh size for experiment: $selected_mesh_size"
        plot_path = joinpath(file_path, "plots")

        plt1 = plot_convergence_generic(elem_sizes_flat, abs.(relative_error), "Effective element size (h)", "h")
        Plots.ylims!(plt1, 10^(-4.5), 10^(-2.8))
        Plots.xlims!(plt1, 1,100)
        Plots.vline!(plt1, [selected_mesh_size], label=L"h_{\mathrm{exp}}", line=:dash, color=:red, legend_column=4)

        # plt2 = plot_convergence_generic(elem_sizes_flat, time_list, "Effective element size (h)")
        plt2 = set_plot(PLOT_CONFIG[:font_size], 
                          sz=(PLOT_CONFIG[:plot_width], PLOT_CONFIG[:plot_height]), 
                          left_margin=PLOT_CONFIG[:left_margin],
                          right_margin=PLOT_CONFIG[:right_margin],
                          top_margin=PLOT_CONFIG[:top_margin],
                          legend_column=2)
        Plots.plot!(plt2, elem_sizes_flat, time_list, label=L"t_{\mathrm{step}}", xlabel="Effective element size (h)", ylabel="Time per step (ms)", marker=:circle, yscale=:log10, xscale=:log10)
        Plots.vline!(plt2, [selected_mesh_size], label=L"h_{\mathrm{exp}}", line=:dash, color=:red, legend_column=3)
        Plots.ylabel!(plt2, "Time per step (ms)")
        Plots.xlims!(plt2, 1,100)
        Plots.ylims!(plt2, 1e1, 10^5.05)

        plt3 = set_plot(PLOT_CONFIG[:font_size], 
                          sz=(PLOT_CONFIG[:plot_width], PLOT_CONFIG[:plot_height]), 
                          left_margin=PLOT_CONFIG[:left_margin],
                          right_margin=PLOT_CONFIG[:right_margin],
                          top_margin=PLOT_CONFIG[:top_margin],
                          legend_column=2)
        Plots.plot!(plt3, elem_sizes_flat, relative_error, label="RMSE", xlabel="Effective element size (h)", ylabel="RMSE", marker=:circle)
            
        plt4 = set_plot(PLOT_CONFIG[:font_size], 
                          sz=(PLOT_CONFIG[:plot_width], PLOT_CONFIG[:plot_height]), 
                          left_margin=PLOT_CONFIG[:left_margin],
                          right_margin=PLOT_CONFIG[:right_margin],
                          top_margin=PLOT_CONFIG[:top_margin],
                          legend_column=2)
        Plots.plot!(plt4, elem_sizes_flat, time_list, label="Time per step", xlabel="Effective element size (h)", ylabel="Time per step (s)", marker=:circle, yscale=:log10, xscale=:log10)

        plt5 = plot_convergence_generic(elem_sizes_flat, abs.(rad_error_list), "Effective element size (h)", "h")
        Plots.vline!(plt5, [selected_mesh_size], label=L"h_{\mathrm{exp}}", line=:dash, color=:red, legend_column=3)
        Plots.ylims!(plt5, 10^(-5), 10^(-2.05))
        Plots.xlims!(plt5, 1,100)

        plt6 = set_plot(PLOT_CONFIG[:font_size], 
                          sz=(PLOT_CONFIG[:plot_width], PLOT_CONFIG[:plot_height]), 
                          left_margin=PLOT_CONFIG[:left_margin],
                          right_margin=PLOT_CONFIG[:right_margin],
                          top_margin=PLOT_CONFIG[:top_margin],
                          legend_column=2)
        Plots.plot!(plt6, elem_sizes_flat, abs.(rad_error_list), label="RMSE", xlabel="Effective element size (h)", ylabel="RMSE", marker=:circle)


        Plots.savefig(plt1, joinpath(plot_path, "height_convergence.pdf"))
        Plots.savefig(plt2, joinpath(plot_path, "computational_cost_vs_elem_size.pdf"))
        Plots.savefig(plt3, joinpath(plot_path, "height_error_vs_elem_size.pdf"))
        Plots.savefig(plt4, joinpath(plot_path, "time_per_step_vs_elem_size.pdf"))
        Plots.savefig(plt5, joinpath(plot_path, "radius_convergence.pdf"))
        Plots.savefig(plt6, joinpath(plot_path, "radius_error_vs_elem_size.pdf"))

    catch e
        @warn "Failed to plot mesh convergence" exception=(e, catch_backtrace())
    end
end

function plot_convergence_time(file_path::String)
    """Plot convergence results for time integration analysis."""
    # try
        height_error_list = readdlm(joinpath(file_path, "height_error_list.csv"), ',', Float64)
        dt_list = readdlm(joinpath(file_path, "t_steps.csv"), ',', Float64)
        h_end_list = readdlm(joinpath(file_path, "final_height_list.csv"), ',', Float64)
        
        relative_error = vec(height_error_list[:])
        dt_sizes = vec(dt_list)
        
        selected_dt = 0.1 # time step size used in the experiment

        plot_path = joinpath(file_path, "plots")
        plt1 = plot_convergence_generic(dt_sizes, relative_error, "Time step size (Δt)", "\\Delta t")
        Plots.vline!(plt1, [selected_dt], label=L"\Delta t_{\mathrm{exp}}", line=:dash, color=:red, legend_column=3)
        Plots.xlims!(plt1, 10^(-3), 5)
        Plots.ylims!(plt1, 1e-6, 10^(-2.2))

        plt2 = set_plot(PLOT_CONFIG[:font_size], 
                          sz=(PLOT_CONFIG[:plot_width], PLOT_CONFIG[:plot_height]), 
                          left_margin=PLOT_CONFIG[:left_margin],
                          right_margin=PLOT_CONFIG[:right_margin],
                          top_margin=PLOT_CONFIG[:top_margin],
                          legend_column=2)
        Plots.plot(plt2, dt_list, h_end_list, label="Final height", xlabel="Time step size (Δt)", ylabel="Final height", yscale=:log10, xscale=:log10, marker=:circle)

        Plots.savefig(plt1, joinpath(plot_path, "time_convergence.pdf"))
        Plots.savefig(plt2, joinpath(plot_path, "final_height_vs_dt.pdf"))
        # catch e
    #     @warn "Failed to plot time integration convergence: $e"
    # end
end

# time_intergration_convergence_analysis()
plot_convergence_time(resolve_data_path("experiments/sim_data/convergence_analysis/stokes_convergence/time_convergence_analysis"))
# mesh_convergence_analysis()
plot_convergence_mesh(resolve_data_path("experiments/sim_data/convergence_analysis/stokes_convergence/mesh_convergence_analysis")) 