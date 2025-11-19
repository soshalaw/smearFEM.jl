# Minimal single-threaded test for the simplified progress helpers.
include("../src/utils.jl")

# route Julia logging to stderr (optional for this test)
setup_stderr_logger()

id = :test_worker
prog_init(id, total=20, desc="worker-test")
for i in 1:20
    sleep(0.01 + 0.02*rand())
    prog_inc(id, 1)
    if rand() < 0.15
        log_via_channel("logged at step $i"; id=id)
    end
end
prog_done(id)

println("thread_progress_test done")
