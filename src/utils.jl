using Logging
using ProgressMeter
import ProgressMeter: next!, finish!

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
    desc::String
end

function next!(p::NonTtyProgress, n::Integer=1; kwargs...)
    p.count += n
    now_ns = time_ns()
    # report if enough time passed or on completion
    if (now_ns - p.last_report_ns) >= p.interval_ns || p.count >= p.total
        try
            @info (isempty(p.desc) ? "Progress" : p.desc) progress=p.count total=p.total
        catch
            # logging shouldn't break execution
        end
        p.last_report_ns = now_ns
    end
    return nothing
end

function finish!(p::NonTtyProgress; kwargs...)
    try
        @info (isempty(p.desc) ? "Progress completed" : string(p.desc, " completed")) progress=p.total
    catch
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
            catch
                return false
            end
        end

        if _isatty_default()
        return Progress(total; kwargs...)
    else
        desc = haskey(kwargs, :desc) ? string(kwargs[:desc]) : ""
        interval_ns = Int(max(1, floor(report_interval * 1e9)))
        return NonTtyProgress(0, Int(total), interval_ns, time_ns(), desc)
    end
end