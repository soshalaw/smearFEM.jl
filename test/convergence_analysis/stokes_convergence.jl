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
    :font_size => 10,
    :plot_height => 360,
    :plot_width => 330,
    :left_margin => -6pt,
    :right_margin => 10pt,
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
    for file in csv_files
        if !endswith(file, ".csv")  # check if the file is a CSV file
            continue
        end
        data = readdlm(file, ',', Float64, '\n', header=false)  # read the observation data
        data_sorted = data[sortperm(data[:, 3]), :]  # sort the data based on the first column (time)
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
    end
    return border_r_list
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
        z = data[2:end, 1]
        r = data[2:end, 2]
        idx = sortperm(z)
        spl = Spline1D(z[idx], r[idx]; k=3, bc="extrapolate")
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

function plot_convergence_generic(x_vals::Vector, y_vals::Vector, x_label::String, x_label_latex::String="",; 
                                  reference_shift::Float64=0.4)
    """Generic function to plot convergence with fitted line.
    
    Arguments:
    - x_vals: x-axis values (h for mesh, Δt for time)
    - y_vals: error values
    - x_label: label for x-axis
    - reference_shift: factor to shift reference line
    """
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
    
    # Generate reference line
    x_min = minimum(x_vals_clean)
    x_max = maximum(x_vals_clean)
    x_line = range(x_min, x_max, length=100)
    y_ref_line = exp.(intercept .+ round(Int,convergence_rate) .* log.(x_line))
    x_line_offset = x_line .* reference_shift
    
    # Create plot
    conv_rate_str = @sprintf("%d", round(Int, convergence_rate))
    plt = setup_plot_config()
    
    Plots.plot!(plt, x_vals_clean, y_vals_clean, 
                label="Error", mode="markers",
                xlabel=x_label, ylabel="Absolute Relative Error",
                marker=:circle, markersize=4, markerstrokewidth=1.5, color="#FF7F0E",
                yscale=:log10, xscale=:log10)
    rate = round(Int, convergence_rate)
    if rate == 0
        @warn "Convergence rate is 0 — data may not be converging; skipping reference line"
    elseif rate == 1
        Plots.plot!(plt, x_line_offset, y_ref_line,
                    label=latexstring("O(\$$(x_label_latex)\$)"),
                    line=2, linewidth=2.5, color=:teal,
                    yscale=:log10, xscale=:log10, linestyle=:dot)
    else
        Plots.plot!(plt, x_line_offset, y_ref_line,
                    label=latexstring("O(\$$(x_label_latex)^{$conv_rate_str}\$)"),
                    line=2, linewidth=2.5, color=:teal,
                    yscale=:log10, xscale=:log10, linestyle=:dot)
    end
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
    step_size = 0.1 # time step size in seconds
    steps = round(Int, sim_time/step_size)  

    obj_pose = [150, 0.0, height/2]
    camera_matrix = get_camera_matrix()
    ne = Union{Int,Float64}[20, 19.75, 14, 10, 9, 6, 3, 2] # number of elements for each mesh size
    
    volume = π*radius^2*height # approximate volume of the cylinder divided by number of elements for the coarsest mesh
    
    h_ref = readdlm(joinpath(gt_path, "data", "h.csv"), ',', Float64, '\n', header=false)[end] # reference height from the finest mesh solution
    r_ref = get_r_curves(joinpath(gt_path, "data"))[end]

    for (iter_index,ne) in enumerate(ne)
        filepath = resolve_data_path(joinpath("ground_truth", "sim_data", "Stokes", control, viscosity_type, "$(element_shape_x)_$(basis_order_x)", "convergence_analysis", "mesh_convergence_analysis", "mesh_sz_$ne"))

        # set_file(filepath)  # create the directory if it doesn't exist

        # Extract element size from directory name
        println("Running for element size = $ne")
        
        exp_params = Dict(
                "element_shape_u" => element_shape_u,
                "basis_order_u"   => basis_order_u,
                "element_shape_p" => element_shape_p,
                "basis_order_p"   => basis_order_p,
                "element_shape_x" => element_shape_x,
                "basis_order_x"   => basis_order_x,
                "ne_gt" => ne,
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
            @error "Simulation failed for element size $ne" exception=(e, catch_backtrace())
            continue
        end
        
        # Read results
        try
            exp_params = read_json(joinpath(filepath, "data", "sim_params.jld2"))
            h_mesh = readdlm(joinpath(filepath, "data","h.csv"), ',', Float64)
            elapsed_time = readdlm(joinpath(filepath, "data","avg_time.csv"), ',', Float64)
            r, z = get_r(joinpath(filepath, "data", "sim_data", "surface_nodes"))

            r_int = r_ref.spline(z)  # Evaluate reference curve at the z values of the current mesh

            ne = exp_params["ne"]
            effective_element_size = (volume / ne)^(1/3)
            
            δh = (h_mesh[end] - h_ref) / h_ref
            δr = mean((r - r_int) ./ r_int)
            
            println("Height error: ", δh)
            println("Radius error: ", δr)

            push!(effective_element_size_list, effective_element_size)
            push!(height_list, h_mesh)
            push!(height_error_list, δh)
            push!(rad_error_list, δr)
            push!(time_list, elapsed_time)
            iter_index += 1
        catch e
            @error "Simulation failed for element size $ne" exception=(e, catch_backtrace())
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
            # write_gt_data(exp_params)
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
        Plots.ylims!(plt1, 1e-5, 10^(-2.7))
        Plots.xlims!(plt1, 1, 25)
        Plots.vline!(plt1, [selected_mesh_size], label=L"h_{\mathrm{exp}}", line=:dash, color=:red, legend_column=3)

        # plt2 = plot_convergence_generic(elem_sizes_flat, time_list, "Effective element size (h)", reference_shift=1.5)
        plt2 = set_plot(PLOT_CONFIG[:font_size], 
                          sz=(PLOT_CONFIG[:plot_width], PLOT_CONFIG[:plot_height]), 
                          left_margin=PLOT_CONFIG[:left_margin],
                          right_margin=PLOT_CONFIG[:right_margin],
                          top_margin=PLOT_CONFIG[:top_margin],
                          legend_column=2)
        Plots.plot!(plt2, elem_sizes_flat, time_list, label=L"t_{\mathrm{step}}", xlabel="Effective element size (h)", ylabel="Time per step (ms)", marker=:circle, yscale=:log10, xscale=:log10)
        Plots.vline!(plt2, [selected_mesh_size], label=L"h_{\mathrm{exp}}", line=:dash, color=:red, legend_column=3)
        Plots.ylabel!(plt2, "Time per step (ms)")
        Plots.xlims!(plt2, 1, 25)
        Plots.ylims!(plt2, 1e1, 10^6.05)

        plt3 = set_plot(PLOT_CONFIG[:font_size], 
                          sz=(PLOT_CONFIG[:plot_width], PLOT_CONFIG[:plot_height]), 
                          left_margin=PLOT_CONFIG[:left_margin],
                          right_margin=PLOT_CONFIG[:right_margin],
                          top_margin=PLOT_CONFIG[:top_margin],
                          legend_column=2)
        Plots.plot!(plt3, elem_sizes_flat, relative_error, label="Height error", xlabel="Effective element size (h)", ylabel="Height error", marker=:circle)
            
        plt4 = set_plot(PLOT_CONFIG[:font_size], 
                          sz=(PLOT_CONFIG[:plot_width], PLOT_CONFIG[:plot_height]), 
                          left_margin=PLOT_CONFIG[:left_margin],
                          right_margin=PLOT_CONFIG[:right_margin],
                          top_margin=PLOT_CONFIG[:top_margin],
                          legend_column=2)
        Plots.plot!(plt4, elem_sizes_flat, time_list, label="Time per step", xlabel="Effective element size (h)", ylabel="Time per step (s)", marker=:circle, yscale=:log10, xscale=:log10)

        plt5 = plot_convergence_generic(elem_sizes_flat, abs.(rad_error_list), "Effective element size (h)", "h")
        Plots.ylims!(plt5, 10^(-1.51), 10^(-1.29))
        Plots.xlims!(plt5, 1, 25)

        plt6 = set_plot(PLOT_CONFIG[:font_size], 
                          sz=(PLOT_CONFIG[:plot_width], PLOT_CONFIG[:plot_height]), 
                          left_margin=PLOT_CONFIG[:left_margin],
                          right_margin=PLOT_CONFIG[:right_margin],
                          top_margin=PLOT_CONFIG[:top_margin],
                          legend_column=2)
        Plots.plot!(plt6, elem_sizes_flat, rad_error_list, label="Radius error", xlabel="Effective element size (h)", ylabel="Radius error", marker=:circle)


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
# plot_convergence_time(resolve_data_path("experiments/sim_data/convergence_analysis/stokes_convergence/time_convergence_analysis"))
mesh_convergence_analysis()
plot_convergence_mesh(resolve_data_path("experiments/sim_data/convergence_analysis/stokes_convergence/mesh_convergence_analysis")) 