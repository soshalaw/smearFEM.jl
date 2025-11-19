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

# --- Option A helpers: progress manager + stderr logger
# Call `setup_stderr_logger()` at program startup to route Julia logs to stderr
function setup_stderr_logger()
    Logging.global_logger(Logging.ConsoleLogger(stderr))
end

# Channel that workers use to send progress/log messages to the single printer task
const PROG_CH = Channel{Any}(1024)

# ANSI color helpers (top-level so `const` is valid)
const _COLOR_CODES = Dict(
    :black => 30, :red => 31, :green => 32, :yellow => 33,
    :blue => 34, :magenta => 35, :cyan => 36, :white => 37
)
const _COLOR_PALETTE = [:green, :cyan, :magenta, :yellow, :blue, :red]

# Start the single-threaded progress manager. Returns the Task so caller can `wait` on it.
function start_progress_manager(ch::Channel{Any}=PROG_CH; bar_width::Integer=30, log_width::Integer=40, use_color::Bool=true)
    # Try ProgressMeter-backed manager when stdout is a TTY. If any error
    # occurs (or stdout isn't a TTY), fall back to the overlay/text manager.
    is_tty = begin
        try
            isatty(stdout)
        catch
            false
        end
    end

    if is_tty
        try
            # ProgressMeter-backed manager
            tasks = Dict{Any, Any}()
            return @async begin
                for msg in ch
                    try
                        if !(msg isa AbstractDict)
                            println(stderr, msg)
                            continue
                        end
                        typ = get(msg, :type, nothing)
                        if typ == :progress
                            id = msg[:id]
                            total = Int(get(msg, :total, 0))
                            total = total > 0 ? total : 1
                            inc = get(msg, :inc, 0)
                            if !(inc isa Integer)
                                try
                                    inc = parse(Int, string(inc))
                                catch
                                    inc = 0
                                end
                            end
                            desc = get(msg, :desc, "")
                            if !haskey(tasks, id)
                                tasks[id] = Progress(total; desc=string(desc))
                            end
                            if inc != 0
                                next!(tasks[id], inc)
                            end
                        elseif typ == :done
                            id = msg[:id]
                            if haskey(tasks, id)
                                try
                                    finish!(tasks[id])
                                catch
                                end
                                delete!(tasks, id)
                            end
                        elseif typ == :log
                            println(stderr, msg[:text])
                            flush(stderr)
                        else
                            println(stderr, "unknown PROG_CH message: ", msg)
                        end
                    catch e
                        println(stderr, "progress manager error: ", e)
                    end
                end

                # cleanup remaining progress bars
                for p in values(tasks)
                    try
                        finish!(p)
                    catch
                    end
                end
            end
        catch e
            # if ProgressMeter approach fails, fall through to overlay fallback
            println(stderr, "ProgressMeter manager failed, falling back: ", e)
        end
    end

    # Overlay/textual fallback (keeps previous overlay behavior), now supports
    # optional colors and truncation of last-log text.
    tasks = Dict{Any, Dict{Symbol, Any}}()
    order = Vector{Any}()

    function colorize(s::AbstractString, col::Symbol)
        if !use_color || !(col in keys(_COLOR_CODES))
            return s
        end
        code = _COLOR_CODES[col]
        return "\e[" * string(code) * "m" * s * "\e[0m"
    end

    function render_bar(curr, total; width=bar_width)
        if total <= 0
            return "[?]"
        end
        frac = clamp(curr/total, 0.0, 1.0)
        filled = Int(round(frac * width))
        bar = repeat('█', filled) * repeat(' ', max(0, width - filled))
        percent = round(frac*100; digits=1)
        return "|$bar| $percent% ($curr/$total)"
    end

    # helper to truncate log text
    function trunc(s::AbstractString, width::Integer)
        if length(s) <= width
            return s
        else
            return first(s, max(0, width-3)) * "..."
        end
    end

    # redraw the full overlay from the top of the terminal
    function redraw_overlay()
        if !is_tty
            return
        end
        # move to home and overwrite overlay lines
        print(stdout, "\e[H")
        for (i, id) in enumerate(order)
            info = tasks[id]
            last = get(info, :lastlog, "")
            desc = info[:desc]
            # pick a color based on worker index
            col = _COLOR_PALETTE[(i - 1) % length(_COLOR_PALETTE) + 1]
            desc_col = colorize(string(desc), col)
            line = string(desc_col, " ", render_bar(info[:curr], info[:total]; width=bar_width))
            if last != ""
                line *= "  | " * trunc(string(last), log_width)
            end
            print(stdout, "\r\e[K", line, "\n")
        end
        # ensure cursor is positioned just after the overlay
        print(stdout, "\e[", length(order) + 1, ";1H")
        flush(stdout)
    end

    return @async begin
        for msg in ch
            try
                if !(msg isa AbstractDict)
                    # treat raw messages as logs below the overlay
                    if is_tty
                        redraw_overlay()
                        println(stdout, string(msg))
                    else
                        println(stdout, string(msg))
                    end
                    continue
                end

                typ = get(msg, :type, nothing)
                if typ == :progress
                    id = msg[:id]
                    total = get(msg, :total, 0)
                    total = (total isa Integer) ? total : 0
                    inc = get(msg, :inc, 0)
                    if !(inc isa Integer)
                        try
                            inc = parse(Int, string(inc))
                        catch
                            inc = 0
                        end
                    end
                    desc = get(msg, :desc, "")
                    if !haskey(tasks, id)
                        push!(order, id)
                        tasks[id] = Dict(:curr => 0, :total => total, :desc => string(desc), :index => length(order))
                    end
                    if inc != 0
                        tasks[id][:curr] += inc
                    end
                    # update total/desc if changed
                    tasks[id][:total] = total > 0 ? total : tasks[id][:total]
                    if desc != ""
                        tasks[id][:desc] = string(desc)
                    end
                    redraw_overlay()

                elseif typ == :done
                    id = msg[:id]
                    if haskey(tasks, id)
                        tasks[id][:curr] = tasks[id][:total]
                        redraw_overlay()
                        # remove the finished worker after showing final state
                        delete!(tasks, id)
                        deleteat!(order, findfirst(==(id), order))
                        # reindex
                        for (i, k) in enumerate(order)
                            tasks[k][:index] = i
                        end
                        redraw_overlay()
                    end

                elseif typ == :log
                    text = get(msg, :text, "")
                    lid = get(msg, :id, nothing)
                    if lid !== nothing && haskey(tasks, lid)
                        # attach the last log to the worker's overlay line
                        tasks[lid][:lastlog] = string(text)
                        redraw_overlay()
                    else
                        # print general logs below overlay
                        if is_tty
                            redraw_overlay()
                            println(stdout, string(text))
                        else
                            println(stderr, string(text))
                        end
                    end

                else
                    if is_tty
                        redraw_overlay()
                        println(stdout, string(msg))
                    else
                        println(stderr, "unknown PROG_CH message: ", msg)
                    end
                end
            catch e
                println(stderr, "progress manager error: ", e)
            end
        end

        # cleanup
        if is_tty
            redraw_overlay()
        end
    end
    end

# Worker helpers -- send atomic messages to the manager
function prog_init(id; total::Integer=1, desc::AbstractString="")
    put!(PROG_CH, Dict(:type=>:progress, :id=>id, :total=>total, :inc=>0, :desc=>desc))
end

function prog_inc(id, n::Integer=1)
    put!(PROG_CH, Dict(:type=>:progress, :id=>id, :total=>0, :inc=>n))
end

function prog_done(id)
    put!(PROG_CH, Dict(:type=>:done, :id=>id))
end

function log_via_channel(text; id=nothing)
    if id === nothing
        put!(PROG_CH, Dict(:type=>:log, :text=>string(text)))
    else
        put!(PROG_CH, Dict(:type=>:log, :text=>string(text), :id=>id))
    end
end

# Usage notes (top-level):
#   setup_stderr_logger()
#   manager_task = start_progress_manager()
#   prog_init(:task1, total=100, desc="task1")
#   prog_inc(:task1)
#   prog_done(:task1)
#   close(PROG_CH)
#   wait(manager_task)