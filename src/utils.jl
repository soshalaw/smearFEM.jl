using Logging
using ProgressMeter
import ProgressMeter: next!, finish!
using Dates
using CSV
using DataFrames

function mat_nan_inf_check(v::AbstractArray)
    row, col = size(v,1), size(v,2)
    for i in 1:row
        for j in 1:col
            if isnan(v[i,j])
                throw(ArgumentError("NaN value found in the matrix"))
            elseif isinf(v[i,j])
                throw(ArgumentError("Inf value found in the matrix"))
            end
        end
    end
end

"""
Configure a simple stderr-backed logger for the package.

This replaces the previous prog_* API. Call this from scripts/tests that
want logging routed to stderr.
"""
function setup_stderr_logger()
    Logging.global_logger(Logging.ConsoleLogger(stderr))
end

"""
Return a Progress-like object that is active when STDOUT is a TTY and a
no-op fallback otherwise. Use `progress_guard(total; desc=..., showspeed=...)`
in places that currently call `Progress(...)` so CI/non-interactive runs stay
clean.
"""
mutable struct NonTtyProgress
    count::Int
    total::Int
    interval_ns::Int    # reporting interval in nanoseconds
    last_report_ns::Int  # last report timestamp (ns)
    last_iteration_ns::Int  # last iteration timestamp (ns)
    desc::String
end

function next!(p::NonTtyProgress, n::Integer=1; kwargs...)
    p.count += n
    now_ns = time_ns()
    # Calculate time per iteration in seconds
    time_per_iter_ms = round((now_ns - p.last_iteration_ns) / 1e9, digits=3)
    p.last_iteration_ns = now_ns
    # report if enough time passed or on completion
    if (now_ns - p.last_report_ns) >= p.interval_ns || p.count >= p.total
        try
            @info "\e[32m$(isempty(p.desc) ? "Progress" : p.desc) | progress: $(round(p.count/p.total*100, digits=1))% of $(p.total) steps| time/iter: $(time_per_iter_ms)s/iter\e[39m"
        catch e
            @debug "Progress logging failed: $e"
        end
        p.last_report_ns = now_ns
    end
    return nothing
end

function finish!(p::NonTtyProgress; kwargs...)
    try
        @info "\e[32m$(isempty(p.desc) ? "Progress completed" : string(p.desc, " completed")) | progress: $(round(p.count/p.total*100, digits=1))%\e[39m"
    catch e
        @debug "Failed to log progress completion: $(e.msg)"
    end
    return nothing
end

function progress_guard(total; kwargs...)
    # Optional keyword to control reporting frequency in seconds
    report_interval = haskey(kwargs, :report_interval) ? float(kwargs[:report_interval]) : 5.0
        # Determine whether we have a TTY on a reasonable stdout object.
        function _get_default_stdout()
            # Try several likely places where a global stdout may be defined.
            if isdefined(Base, :STDOUT)
                return Base.STDOUT
            elseif isdefined(Base, :stdout)
                return Base.stdout
            elseif isdefined(Main, :STDOUT)
                return Main.STDOUT
            elseif isdefined(Main, :stdout)
                return Main.stdout
            else
                return nothing
            end
        end

        function _isatty_default()
            io = _get_default_stdout()
            if io === nothing
                return false
            end
            try
                return isatty(io)
            catch e
                @debug "isatty check failed: $e"
                return false
            end
        end

        if _isatty_default()
        return Progress(total; kwargs...)
    else
        desc = haskey(kwargs, :desc) ? string(kwargs[:desc]) : ""
        interval_ns = Int(max(1, floor(report_interval * 1e9)))
        now_ns = time_ns()
        return NonTtyProgress(0, Int(total), interval_ns, now_ns, now_ns, desc)
    end
end

# logging Functions

function write_time_log(start_time::Dates.DateTime, end_time::Dates.DateTime, params::Dict; dest_dir::String)
    log_filepath = joinpath(dest_dir, "time_log.txt")
    if !isdir(dest_dir)
        mkpath(dest_dir)
    end
    open(log_filepath, "a") do io
        println(io, "Start Time: ", Dates.format(start_time, "yyyy-mm-dd HH:MM:SS"))
        println(io, "End Time: ", Dates.format(end_time, "yyyy-mm-dd HH:MM:SS"))
        println(io, "Elapsed Time (seconds): ", canonicalize(end_time - start_time))
        println(io, "Parameters: ")
        for (k,v) in params
            println(io, "    ", k, " => ", v)
        end
        println(io, "----------------------------------------")
    end
end


function dataframe_2_vec(df::DataFrame)
    
    data_list::Vector{AbstractArray} = Vector{AbstractArray}()
    rows, cols = size(df)
    for i in 1:rows
        _data = Array(df[i,:])
        data = [x for x in _data if !ismissing(x)]
        push!(data_list, data)
    end
    return data_list
end

"""
    get_cMat(type; λ=nothing, μ=nothing, Young=nothing, ν=nothing)

Returns the constitutive matrix for a given material type.

# Arguments
- `type::String`: Type of constitutive matrix ("lame" or "standard").
- `λ::Float64`: Lame's first parameter.
- `μ::Float64`: Shear modulus.
- `Young::Float64`: Young's modulus.
- `ν::Float64`: Poisson's ratio.

# Returns
- `cMat::Matrix{Float64}`: Constitutive matrix.
"""
function get_cMat(c1::Float64, c2::Float64; type::String="lame")::Matrix{Float64}
    if type == "lame"
        cMat =  [[ 2*c2+c1  c1    c1    0  0  0]; 
                        [  c1   2*c2+c1  c1    0  0  0]; 
                        [  c1     c1   2*c2+c1 0  0  0]; 
                        [  0      0      0     c2 0  0]; 
                        [  0      0      0     0 c2  0]
                    [  0      0      0     0  0 c2]]  # constitutive matrix
        return cMat
    elseif type == "standard"
        cMat =  (c1/((1+c2)*(1-2*c2)))*[[1-c2 c2 c2   0      0      0]; 
                                                [c2 1-c2 c2   0      0      0]; 
                                                [c2 c2 1-c2   0      0      0]; 
                                                [0 0  0 (1-2*c2)/2 0      0]; 
                                                [0 0  0    0   (1-2*c2)/2 0]; 
                                                [0 0  0    0      0   (1-2*c2)/2]]  # constitutive matrix 
        return cMat
    else
        throw(ArgumentError("Type of cMat unknown it should be either 'lame' or 'standard'"))  
    end
end