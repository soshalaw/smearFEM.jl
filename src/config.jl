"""
    config.jl

Centralized configuration and path resolution for smearFEM.jl.

Supports three-tier fallback:
1. Environment variables (SMEAR_DATA_DIR, SMEAR_MESH_DIR, SMEAR_SCRATCH_DIR)
2. Configuration file (config.toml in project root)
3. Relative paths from project root
"""

using TOML
import Logging

# Configuration Loading
"""
    load_config(config_path::String = "config.toml") -> Dict

Load configuration from TOML file. Returns empty dict if file not found.

# Arguments
- `config_path`: Path to config.toml (default: project root)

# Returns
- `Dict`: Configuration parameters, or empty dict if not found
"""
function load_config(config_path::String = "config.toml")
    # Try to find config in project root
    project_root = dirname(dirname(@__DIR__))
    full_path = joinpath(project_root, config_path)
    
    if isfile(full_path)
        try
            return TOML.parsefile(full_path)
        catch e
            @warn "Failed to parse config file at $full_path: $(e.msg)"
            return Dict()
        end
    end
    
    return Dict()
end

# Path Resolution Functions
"""
    get_data_dir() -> String

Resolve the data directory using hybrid fallback.

Priority:
1. Environment variable SMEAR_DATA_DIR
2. Config file: config.toml[data_dir]
3. Default: ~/SMEAR-Data

# Returns
- `String`: Absolute path to data directory
"""
function get_data_dir()::String
    # Try environment variable
    if haskey(ENV, "SMEAR_DATA_DIR")
        path = ENV["SMEAR_DATA_DIR"]
        isdir(path) && return abspath(path)
        @warn "SMEAR_DATA_DIR=$path does not exist; using fallback"
    end
    
    # Try config file
    config = load_config()
    if haskey(config, "data_dir")
        path = config["data_dir"]
        isdir(path) && return abspath(path)
        @warn "data_dir=$path from config does not exist; using fallback"
    end
    
    # Default fallback
    default = joinpath(homedir(), "SMEAR-Data")
    return default
end

"""
    get_mesh_dir() -> String

Resolve the mesh directory using hybrid fallback.

Priority:
1. Environment variable SMEAR_MESH_DIR
2. Config file: config.toml[mesh_dir]
3. Default: ~/SMEAR-Data/meshes

# Returns
- `String`: Absolute path to mesh directory
"""
function get_mesh_dir()::String
    # Try environment variable
    if haskey(ENV, "SMEAR_MESH_DIR")
        path = ENV["SMEAR_MESH_DIR"]
        isdir(path) && return abspath(path)
        @warn "SMEAR_MESH_DIR=$path does not exist; using fallback"
    end
    
    # Try config file
    config = load_config()
    if haskey(config, "mesh_dir")
        path = config["mesh_dir"]
        isdir(path) && return abspath(path)
        @warn "mesh_dir=$path from config does not exist; using fallback"
    end
    
    # Default fallback
    default = joinpath(get_data_dir(), "meshes")
    return default
end

"""
    get_scratch_dir() -> String

Resolve the scratch/temporary directory using hybrid fallback.

Priority:
1. Environment variable SMEAR_SCRATCH_DIR
2. Config file: config.toml[scratch_dir]
3. Default: /tmp/smear or ~/SMEAR-Scratch

# Returns
- `String`: Absolute path to scratch directory
"""
function get_scratch_dir()::String
    # Try environment variable
    if haskey(ENV, "SMEAR_SCRATCH_DIR")
        path = ENV["SMEAR_SCRATCH_DIR"]
        isdir(path) && return abspath(path)
        @warn "SMEAR_SCRATCH_DIR=$path does not exist; using fallback"
    end
    
    # Try config file
    config = load_config()
    if haskey(config, "scratch_dir")
        path = config["scratch_dir"]
        isdir(path) && return abspath(path)
        @warn "scratch_dir=$path from config does not exist; using fallback"
    end
    
    # Default fallback - prefer /tmp on Unix-like systems
    if Sys.isunix()
        default = "/tmp/smear"
    else
        default = joinpath(homedir(), "SMEAR-Scratch")
    end
    return default
end

"""
    resolve_data_path(relative_path::String) -> String

Resolve a data path relative to the data directory.

# Arguments
- `relative_path`: Path relative to data directory

# Returns
- `String`: Absolute path

# Example
```julia
gt_path = resolve_data_path("ground_truth/sim_data/Stokes/force/constant/Q2_16")
```
"""
function resolve_data_path(relative_path::String)::String
    return joinpath(get_data_dir(), relative_path)
end

"""
    resolve_mesh_path(mesh_name::String) -> String

Resolve a mesh file path relative to the mesh directory.

# Arguments
- `mesh_name`: Mesh filename or relative path

# Returns
- `String`: Absolute path to mesh file

# Example
```julia
mesh_path = resolve_mesh_path("channel.msh")
```
"""
function resolve_mesh_path(mesh_name::String)::String
    return joinpath(get_mesh_dir(), mesh_name)
end

"""
    create_output_dir(experiment_name::String; ensure_dir::Bool=true) -> String

Create and return a unique output directory for an experiment.

# Arguments
- `experiment_name`: Descriptive name for the experiment
- `ensure_dir`: Create directory if it doesn't exist (default: true)

# Returns
- `String`: Absolute path to output directory

# Example
```julia
out_dir = create_output_dir("stokes_convergence_2026-05-01")
```
"""
function create_output_dir(experiment_name::String; ensure_dir::Bool=true)::String
    scratch = get_scratch_dir()
    out_dir = joinpath(scratch, experiment_name)
    
    if ensure_dir && !isdir(out_dir)
        try
            mkpath(out_dir)
        catch e
            @warn "Failed to create output directory $out_dir: $(e.msg)"
        end
    end
    
    return out_dir
end

# ============================================================================
# Configuration Display
# ============================================================================

"""
    show_config()

Display current configuration settings (useful for debugging).
"""
function show_config()
    println("╔════════════════════════════════════════════════════════════╗")
    println("║           smearFEM.jl Configuration Summary              ║")
    println("╠════════════════════════════════════════════════════════════╣")
    println("║ Data Directory:      $(get_data_dir())")
    println("║ Mesh Directory:      $(get_mesh_dir())")
    println("║ Scratch Directory:   $(get_scratch_dir())")
    println("╚════════════════════════════════════════════════════════════╝")
end
