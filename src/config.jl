"""
    config.jl

Centralized configuration and path resolution for smearFEM.jl.

Supports three-tier fallback:
1. Environment variables (SMEAR_DATA_DIR, SMEAR_MESH_DIR, SMEAR_SCRATCH_DIR)
2. Configuration file (config.toml, searched from the project root upwards)
3. Relative paths from project root

The config file lives outside this repository (in the smear-modules directory
that holds both smearFEM.jl and smearPerception) so the Julia and Python
packages resolve the same paths. See src/config.py on the Python side.
"""

using TOML
import Logging

expand_env_vars(path::String) = replace(path, r"\$\{(\w+)\}" => m -> get(ENV, m[3:end-1], m))

# Configuration Loading
"""
    find_config(config_name::String = "config.toml") -> Union{String, Nothing}

Locate the config file by walking up from the project root. Honours
SMEAR_CONFIG_FILE if set.

# Returns
- `String`: Absolute path to the config file, or `nothing` if not found
"""
function find_config(config_name::String = "config.toml")
    if haskey(ENV, "SMEAR_CONFIG_FILE")
        path = ENV["SMEAR_CONFIG_FILE"]
        isfile(path) && return abspath(path)
        @warn "SMEAR_CONFIG_FILE=$path does not exist; searching for $config_name instead"
    end

    dir = dirname(@__DIR__)
    while true
        candidate = joinpath(dir, config_name)
        isfile(candidate) && return candidate
        parent = dirname(dir)
        parent == dir && return nothing  # reached the filesystem root
        dir = parent
    end
end

"""
    load_config(config_path::Union{String, Nothing} = nothing) -> Tuple{Dict, String}

Load configuration from TOML file. Returns an empty dict if no file is found.

# Arguments
- `config_path`: Explicit path to config.toml (default: search from project root upwards)

# Returns
- `Tuple{Dict, String}`: Configuration parameters and the directory holding the
  config file. Relative paths in the config are resolved against that directory.
"""
function load_config(config_path::Union{String, Nothing} = nothing)
    full_path = config_path === nothing ? find_config() : config_path

    if full_path !== nothing && isfile(full_path)
        try
            return TOML.parsefile(full_path), dirname(abspath(full_path))
        catch e
            @warn "Failed to parse config file at $full_path: $(e.msg)"
        end
    end

    return Dict(), dirname(@__DIR__)
end

"""
    config_path_value(config::Dict, base_dir::String, key::String) -> Union{String, Nothing}

Expand environment variables in `config[key]` and resolve it against `base_dir`
when it is relative. Returns `nothing` when the key is absent.
"""
function config_path_value(config::Dict, base_dir::String, key::String)
    haskey(config, key) || return nothing
    path = expand_env_vars(config[key])
    return isabspath(path) ? path : joinpath(base_dir, path)
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
    config, base_dir = load_config()
    path = config_path_value(config, base_dir, "data_dir")
    if path !== nothing
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
    config, base_dir = load_config()
    path = config_path_value(config, base_dir, "mesh_dir")
    if path !== nothing
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
    config, base_dir = load_config()
    path = config_path_value(config, base_dir, "scratch_dir")
    if path !== nothing
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
