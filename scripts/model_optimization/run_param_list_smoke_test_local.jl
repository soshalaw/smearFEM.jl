using Base.Threads

println("Local smoke test: run_param_list threading check")

# dummy optimize to simulate work
function optimize(params)
    tid = Threads.threadid()
    println("[dummy optimize] thread=$tid params_id=$(get(params, "id", "?"))")
    sleep(0.12)
    return nothing
end

function run_param_list_local(params_list::Vector{<:AbstractDict}; max_workers::Int=8)
    nparams = length(params_list)
    if nparams == 0
        return
    end
    workers = max(1, min(Threads.nthreads(), max_workers))
    println("Starting run_param_list_local with $nparams params and $workers workers (Threads.nthreads()=$(Threads.nthreads()))")
    ch = Channel{Int}(nparams)
    @sync begin
        for i in 1:nparams
            put!(ch, i)
        end
        close(ch)
        tasks = Vector{Task}(undef, workers)
        for w in 1:workers
            tasks[w] = Threads.@spawn begin
                while true
                    idx = try
                        take!(ch)
                    catch
                        break
                    end
                    params = params_list[idx]
                    try
                        optimize(params)
                    catch err
                        @error "Dummy optimize failed for idx=$idx: $err"
                    end
                end
            end
        end
        for t in tasks
            wait(t)
        end
    end
    println("run_param_list_local completed")
end

param_list = [Dict("id"=>i) for i in 1:12]
run_param_list_local(param_list; max_workers=4)
println("Local smoke test finished")
