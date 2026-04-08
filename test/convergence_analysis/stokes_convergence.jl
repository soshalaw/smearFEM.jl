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
# ============================================================================
# Plot configuration constants
# ============================================================================
const PLOT_CONFIG = Dict(
    :font_size => 10,
    :plot_height => 360,
    :plot_width => 330,
    :left_margin => -6pt,
    :right_margin => 10pt,
    :top_margin => -3pt
)


# ============================================================================
# Simulation configuration helpers
# ============================================================================

function get_object_pose(height::Float64)
    """Returns object pose matrix."""
    obj_pose = zeros(Float64, 4, 4)
    obj_pose[1,1] = -1.0
    obj_pose[2,3] = -1.0
    obj_pose[3,2] = -1.0
    obj_pose[1:3,4] = [0.0, height/2, 150.0]
    return obj_pose
end

function get_camera_matrix()
    """Returns camera matrix for rendering."""
    return [[8*2048/7.07, 0.0, 2048/2] [0.0, 8*1536/5.3, 1536/2] [0.0, 0.0, 1.0]]'
end

# ============================================================================
# Plotting utilities
# ============================================================================

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

function plot_convergence_generic(x_vals::Vector, y_vals::Vector, x_label::String, 
                                  output_path::String, filename::String; 
                                  reference_shift::Float64=0.4)
    """Generic function to plot convergence with fitted line.
    
    Arguments:
    - x_vals: x-axis values (h for mesh, Δt for time)
    - y_vals: error values
    - x_label: label for x-axis
    - output_path: directory to save plot
    - filename: output filename
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
    y_ref_line = exp.(intercept .+ convergence_rate .* log.(x_line))
    x_line_offset = x_line .* reference_shift
    
    # Create plot
    conv_rate_str = @sprintf("%.2f", convergence_rate)
    plt = setup_plot_config()
    
    Plots.plot!(plt, x_vals_clean, y_vals_clean, 
                label="Error", mode="markers",
                xlabel=x_label, ylabel="Absolute Relative Error",
                marker=:circle, markersize=4, markerstrokewidth=1.5, color=:steelblue,
                yscale=:log10, xscale=:log10)
    Plots.plot!(plt, x_line_offset, y_ref_line, 
                label=latexstring("O(\$h^{$conv_rate_str}\$)"),
                line=2, linewidth=2.5, color=:darkorange, 
                yscale=:log10, xscale=:log10, linestyle=:dot)
    
    @info "Convergence rate: O(h^$(conv_rate_str))"
    
    # Save plot
    mkpath(output_path)
    Plots.savefig(plt, joinpath(output_path, filename))
end

# ============================================================================
# Main analysis functions
# ============================================================================

function generate_mesh_geo(radius::Float64, height::Float64, elem_size::Float64, output_path::String, template_path::String)
    """Generate mesh.geo file from template with specified parameters.
    
    Returns true if a new mesh.geo was generated, false if existing file has correct parameters.
    """
    # Check if file already exists with correct parameters
    if isfile(output_path)
        try
            content = read(output_path, String)
            has_radius = contains(content, "radius = $radius;")
            has_height = contains(content, "height = $height;")
            has_elem_size = contains(content, "elem_size_2d = $elem_size;")
            
            if has_radius && has_height && has_elem_size
                @info "mesh.geo already exists with correct parameters (radius=$radius, height=$height, elem_size=$elem_size), skipping generation"
                return false
            end
        catch
            # If error reading, regenerate
        end
    end
    
    # Read the template mesh.geo file
    template_content = read(template_path, String)
    
    # Replace parameters in the template
    geo_content = replace(template_content, 
        "radius = 25.0;" => "radius = $radius;",
        "height = 40.0;" => "height = $height;",
        "elem_size_2d = 10;" => "elem_size_2d = $elem_size;")
    
    # Write the modified content to the output path
    write(output_path, geo_content)
    @info "Generated mesh.geo at $output_path from template $template_path"
    return true
end

function run_gmsh(geo_path::String, msh_path::String, mesh_order::Int=2)
    """Run gmsh to generate mesh from .geo file"""
    # Check if gmsh is available
    gmsh_path = Sys.which("gmsh")
    if gmsh_path === nothing
        @warn "gmsh executable not found in PATH. Please install gmsh or add it to PATH."
        return false
    end
    
    cmd = `$gmsh_path $geo_path -3 -format msh -order $mesh_order -o $msh_path`
    @info "Running gmsh: $cmd"
    try
        run(cmd)
        @info "Mesh generated: $msh_path"
        return true
    catch e
        @error "Failed to run gmsh: $e"
        return false
    end
end

function mesh_convergence_analysis(; radius::Float64=25.0, height::Float64=40.0, elem_sizes::Vector=[10, 8, 6, 4], 
                                  template_mesh_geo_path::String="/home/soshala/SMEAR-PhD/smear-modules/smearFEM.jl/test/convergence_analysis/mesh.geo")

    height_list = AbstractArray[]
    height_error_list = Float64[]
    border_error_list = Float64[]
    effective_element_size_list = Float64[]
    μ_list = AbstractArray[]
    time_list = []

    β::Float64 = 100.0 # penalty parameter for the ground truth
    η::Float64 = 100.0

    viscosity_type::String = "constant" # "constant" or "bulk_viscosity"
    FunctionClass_x::String = "Q2" # Function space for the ground truth
    FunctionClass_p::String = "Q1"
    FunctionClass_u::String = "Q2"
    control::String = "force" # "force" or "velocity"

    F_ext::Float64 = 9518.61 # force applied to the cylinder in N
    sim_time::Float64 = 5.0 # simulation time in seconds
    step_size = 0.1 # time step size in seconds
    steps = round(Int, sim_time/step_size)  

    obj_pose = get_object_pose(height)
    camera_matrix = get_camera_matrix()
    filepath = string("/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/ground_truth/sim_data/Stokes/$control/$viscosity_type/$FunctionClass_x/convergence_analysis")

    volume = π*radius^2*height # approximate volume of the cylinder divided by number of elements for the coarsest mesh

    mesh_dir = dirname(filepath)
    set_file(mesh_dir)  # create the directory if it doesn't exist
    
    iter_index = 1
    h_ref = 37.514580952625970
    border_ref = []
    
    for elem_size in Float64.(elem_sizes)
        println("Running for element size = $elem_size")
        t_start = Dates.now()
        
        # Define mesh paths
        mesh_path = joinpath(mesh_dir, "convergence_analysis","mesh_convergence_analysis","mesh_$elem_size", "mesh.geo")
        mesh_msh_path_x = joinpath(dirname(mesh_path), "cylinder_x_$elem_size.msh")
        mesh_msh_path_p = joinpath(dirname(mesh_path), "cylinder_p_$elem_size.msh")
        
        # Create directory and check/generate mesh.geo (returns true if newly generated)
        set_file(dirname(mesh_path))
        mesh_generated = generate_mesh_geo(radius, height, elem_size, mesh_path, template_mesh_geo_path)
        
        # Only run gmsh if mesh.geo was newly generated or msh files don't exist
        if mesh_generated || !isfile(mesh_msh_path_x) || !isfile(mesh_msh_path_p)
            if !run_gmsh(mesh_path, mesh_msh_path_x, 2) || !run_gmsh(mesh_path, mesh_msh_path_p, 1) 
                @error "Mesh generation failed for element size $elem_size"
                continue
            end
        else
            @info "Mesh files already exist with correct parameters, skipping gmsh"
        end
        
        # Run simulation with the generated mesh
        mesh_filepath = dirname(mesh_path)
        exp_params = Dict("FunctionClass_x" => FunctionClass_x, "FunctionClass_u" => FunctionClass_u, "FunctionClass_p" => FunctionClass_p, 
                         "ne_gt" => elem_size, "β_gt" => β, "η_gt" => η, "filepath_gt" => joinpath(mesh_filepath, "simulation"), 
                         "control" => control, "viscosity_type" => viscosity_type, "obj_pose_gt" => obj_pose, 
                         "F_ext" => F_ext, "sim_time_gt" => sim_time, "steps_gt" => steps, 
                         "r" => radius, "h" => height, "camera_matrix" => camera_matrix, "animate" => true, "mesh_path" => mesh_filepath)

        try
            write_gt_data(exp_params)
        catch e
            @error "Simulation failed for element size $elem_size" exception=e
            continue
        end

        t_end = Dates.now()
        elapsed_time = t_end - t_start
        
        # Read results
        # try
            exp_params = read_json(joinpath(mesh_filepath, "simulation", "data", "sim_params.jld2"))
            h_mesh = readdlm(joinpath(mesh_filepath, "simulation", "data","h.csv"), ',', Float64)

            ne = exp_params["ne"]
            effective_element_size = (volume / ne)^(1/3)
            
            δh = abs(h_mesh[end] - h_ref) / h_ref
            
            println("Height error: ", δh)
            push!(effective_element_size_list, effective_element_size)
            push!(height_list, h_mesh)
            push!(height_error_list, δh)
            push!(time_list, elapsed_time.value/steps)
            iter_index += 1
        # catch e
        #     @warn "Failed to read results for element size $elem_size: $e"
        # end
    end
    write_csv("/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/experiments/sim_data/convergence_analysis/stokes_convergence/mesh_convergence_analysis/effective_element_size", effective_element_size_list)
    write_csv("/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/experiments/sim_data/convergence_analysis/stokes_convergence/mesh_convergence_analysis/height_list", height_list)
    write_csv("/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/experiments/sim_data/convergence_analysis/stokes_convergence/mesh_convergence_analysis/height_error_list", height_error_list)
    write_csv("/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/experiments/sim_data/convergence_analysis/stokes_convergence/mesh_convergence_analysis/elem_sizes", elem_sizes)
    write_csv("/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/experiments/sim_data/convergence_analysis/stokes_convergence/mesh_convergence_analysis/time_list", time_list)
end

function time_intergration_convergence_analysis(; radius::Float64=25.0, height::Float64=40.0, dt_list::Vector=[2, 1, 0.5, 0.1, 0.05, 0.01, 0.005, 0.001], 
                                  template_mesh_geo_path::String="/home/soshala/SMEAR-PhD/smear-modules/smearFEM.jl/test/convergence_analysis/mesh.geo")    
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
    filepath = string("/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/ground_truth/sim_data/Stokes/$control/$viscosity_type/$FunctionClass_x/convergence_analysis")
 
    volume = π*radius^2*height # approximate volume of the cylinder divided by number of elements for the coarsest mesh
    
    mesh_dir = dirname(filepath)
    set_file(mesh_dir)  # create the directory if it doesn't exist
    
    iter_index = 1
    h_ref = 37.514669827538519
    
    mesh_filepath = joinpath(mesh_dir, "convergence_analysis", "mesh_convergence_analysis", "mesh_4.0")
    
    for step_size in Float64.(dt_list)
        filepath = joinpath(mesh_dir, "convergence_analysis", "time_integration_convergence_analysis","step_$step_size", "simulation")
        println("Running for time step size = $step_size")
        steps = round(Int, sim_time/step_size)
        
        # Run simulation with the generated mesh
        exp_params = Dict("FunctionClass_x" => FunctionClass_x, "FunctionClass_u" => FunctionClass_u, "FunctionClass_p" => FunctionClass_p, 
                         "ne_gt" => 4, "β_gt" => β, "η_gt" => η, "filepath_gt" => filepath, 
                         "control" => control, "viscosity_type" => viscosity_type, "obj_pose_gt" => obj_pose, 
                         "F_ext" => F_ext, "sim_time_gt" => sim_time, "steps_gt" => steps, 
                         "r" => radius, "h" => height, "camera_matrix" => camera_matrix, "animate" => true, "mesh_path" => mesh_filepath)

        # try
            write_gt_data(exp_params)
        # catch e
        #     @error "Simulation failed for time step size $step_size" exception=e
        #     continue
        # end
        
        # Read results
        try
            exp_params = read_json(joinpath(filepath, "data", "sim_params.jld2"))
            h_mesh = readdlm(joinpath(filepath, "data","h.csv"), ',', Float64)

            dt = exp_params["time_steps"]

            ne = exp_params["ne"]
            effective_element_size = (volume / ne)^(1/3)

            δh = abs(h_mesh[end] - h_ref) / h_ref
            _δh = (h_mesh[end] - h_ref) / h_ref
            
            println("Height error: ", _δh)
            push!(effective_element_size_list, effective_element_size)
            push!(height_list, h_mesh)
            push!(height_error_list, δh)
            push!(_dt_list, dt)
            push!(final_height_list, h_mesh[end])
            iter_index += 1
        catch e
            @warn "Failed to read results for time step size $step_size: $e"
        end
    end
    println("Final height list: ", final_height_list)
    write_csv("/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/experiments/sim_data/convergence_analysis/stokes_convergence/time_convergence_analysis/effective_element_size", effective_element_size_list)
    write_csv("/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/experiments/sim_data/convergence_analysis/stokes_convergence/time_convergence_analysis/height_list", height_list)
    write_csv("/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/experiments/sim_data/convergence_analysis/stokes_convergence/time_convergence_analysis/height_error_list", height_error_list)
    write_csv("/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/experiments/sim_data/convergence_analysis/stokes_convergence/time_convergence_analysis/t_steps", dt_list)
    write_csv("/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/experiments/sim_data/convergence_analysis/stokes_convergence/time_convergence_analysis/final_height_list", final_height_list)
end

function plot_convergence_mesh(file_path::String)
    """Plot mesh convergence analysis results."""
    try
        height_error_list = readdlm(joinpath(file_path, "height_error_list.csv"), ',', Float64)
        elem_sizes = readdlm(joinpath(file_path, "effective_element_size.csv"), ',', Float64)
        
        relative_error = vec(height_error_list[:, end])
        elem_sizes_flat = vec(elem_sizes[:, 1])
        
        plot_path = joinpath(file_path, "plots")
        plot_convergence_generic(elem_sizes_flat, relative_error, "Effective element size (h)", 
                                plot_path, "height_convergence.pdf")
    catch e
        @warn "Failed to plot mesh convergence: $e"
    end
end

function plot_convergence_time(file_path::String)
    """Plot convergence results for time integration analysis."""
    # try
        height_error_list = readdlm(joinpath(file_path, "height_error_list.csv"), ',', Float64)
        dt_list = readdlm(joinpath(file_path, "t_steps.csv"), ',', Float64)
        h_end_list = readdlm(joinpath(file_path, "final_height_list.csv"), ',', Float64)
        
        Plots.plot(dt_list, h_end_list, label="Final height", xlabel="Time step size (Δt)", ylabel="Final height", yscale=:linear, xscale=:log10, marker=:circle)
        Plots.savefig(Plots.current(), joinpath(file_path, "plots", "final_height_vs_dt.pdf"))
        
        println("Height error list: ", height_error_list)
        relative_error = vec(height_error_list[:])
        dt_sizes = vec(dt_list)
        
        plot_path = joinpath(file_path, "plots")
        plot_convergence_generic(dt_sizes, relative_error, "Time step size (Δt)", 
                                plot_path, "time_integration_convergence.pdf")
    # catch e
    #     @warn "Failed to plot time integration convergence: $e"
    # end
end

time_intergration_convergence_analysis()
plot_convergence_time("/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/experiments/sim_data/convergence_analysis/stokes_convergence/time_convergence_analysis")
# mesh_convergence_analysis()
# plot_convergence_mesh("/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/experiments/sim_data/convergence_analysis/stokes_convergence/mesh_convergence_analysis") 