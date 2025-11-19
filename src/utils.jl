using ProgressMeter
using Logging

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

# --- Simplified single-threaded progress helpers
# Remove channel/task-based manager; provide a small ProgressMeter-backed
# helper set for sequential runs. This keeps the `prog_*` API but avoids
# multithreading/channel plumbing.
function setup_stderr_logger()
    Logging.global_logger(Logging.ConsoleLogger(stderr))
end

const _PROG_TASKS = Dict{Any, Progress}()

function prog_init(id; total::Integer=1, desc::AbstractString="")
    try
        _PROG_TASKS[id] = Progress(total; desc=string(desc))
    catch
        # ignore if ProgressMeter not available or construction fails
    end
end

function prog_inc(id, n::Integer=1)
    try
        if haskey(_PROG_TASKS, id)
            next!(_PROG_TASKS[id], n)
        else
            # fallback: print a simple progress increment to stderr
            println(stderr, "[PROG_INC] ", string(id), ": +", string(n))
        end
    catch
    end
end

function prog_done(id)
    try
        if haskey(_PROG_TASKS, id)
            finish!(_PROG_TASKS[id])
            delete!(_PROG_TASKS, id)
        else
            println(stderr, "[PROG_DONE] ", string(id))
        end
    catch
    end
end

function log_via_channel(text; id=nothing)
    try
        if id === nothing
            println(stderr, string(text))
        else
            println(stderr, "(" * string(id) * ") " * string(text))
        end
    catch
    end
end