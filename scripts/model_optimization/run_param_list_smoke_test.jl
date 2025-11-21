println("Smoke test: run_param_list threading check")
# include the main script (local to this folder)
include("test_opt_stokes.jl")

# override optimize with a lightweight dummy to avoid heavy computations
optimize_old = optimize
optimize = function(params)
    tid = Threads.threadid()
    println("[dummy optimize] thread=$tid params_id=$(get(params, "id", "?"))")
    # small sleep to simulate work
    sleep(0.15)
    return nothing
end

# build small param list
param_list = [Dict("id"=>i) for i in 1:8]
println("Calling run_param_list with $(length(param_list)) params")
run_param_list(param_list; max_workers=4)
println("Smoke test finished")
