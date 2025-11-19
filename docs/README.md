# smearFEM: Thread-safe progress & logging helpers

This repository includes a small Option-A implementation to keep multithreaded
console output clutter-free: progress updates are serialized and logging is
routed to `stderr` so progress output on `stdout` isn't corrupted.

Quick usage
-----------

1. The package wires up the stderr logger automatically when the module is loaded.
   If you prefer to control this manually, comment out the auto-call in
   `src/smearFEM.jl` after the `include("utils.jl")` line.

2. Start the progress manager (optional — the test harness starts it for you):

   ```julia
   using smearFEM
   manager = start_progress_manager()   # returns a Task
   setup_stderr_logger()                # optional if you disabled auto-setup
   ```

3. From worker threads:

   ```julia
   prog_init(:task1, total=100, desc="Task 1")
   prog_inc(:task1)           # increment by 1
   prog_inc(:task1, 5)        # increment by 5
   prog_done(:task1)
   log_via_channel("worker message")
   ```

4. Shutdown gracefully:

   ```julia
   close(PROG_CH)      # signal manager to finish
   wait(manager)       # wait for manager Task to consume all messages
   ```

Notes & rationale
-----------------
- Progress updates are sent via `PROG_CH` so a single, dedicated manager task
  can serialize stdout writes. This prevents partial/interleaved lines.
- Logging is routed to `stderr` by default using `setup_stderr_logger()` to keep
  informational messages separate from progress output on `stdout`.
- The current manager prints simple atomic textual progress lines. If you want
  fancier terminal progress bars using ProgressMeter.jl, we can upgrade the
  manager to redraw bars from the manager task (this requires testing on your
  target terminal/TTY).

Example test harness
--------------------
See `test/thread_progress_test.jl` for an example that spawns multiple threads
and demonstrates non-interleaved progress and log messages.
