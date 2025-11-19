# Simple test for the Option A progress manager
# This test includes the utils file directly so it can be run in isolation.

include("../src/utils.jl")

# route Julia logging to stderr (optional for this test)
setup_stderr_logger()

manager = start_progress_manager()

nworkers = 4

Threads.@sync begin
    for w in 1:nworkers
        Threads.@spawn begin
            id = Symbol("w", w)
            prog_init(id, total=20, desc="worker $w")
            for i in 1:20
                sleep(0.005 + 0.02*rand())
                prog_inc(id, 1)
                if rand() < 0.15
                    # send the log together with the worker id so it appears in the overlay
                    log_via_channel("logged at step $i"; id=id)
                end
            end
            prog_done(id)
        end
    end
end

# shutdown manager
close(PROG_CH)
wait(manager)

println("thread_progress_test done")
