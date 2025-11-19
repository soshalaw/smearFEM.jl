# Minimal single-threaded test: the prog_* API was removed. Use logging or
# plain prints instead.
include("../src/utils.jl")

# route Julia logging to stderr
setup_stderr_logger()

id = :test_worker
@info "Starting test" id=id
for i in 1:20
    sleep(0.01 + 0.02*rand())
    @info "progress" id=id progress=i total=20
    if rand() < 0.15
        @info "logged at step $i" id=id
    end
end

@info "Completed test" id=id
println("thread_progress_test done")
