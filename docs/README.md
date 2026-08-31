# smearFEM: logging helpers

The small `prog_*` progress API was removed from this repository to simplify
the codebase and avoid hidden multithreading plumbing. Instead, use the
standard Julia `Logging` facilities or simple writes to `stderr`/`stdout`.

Quick usage
-----------

1. Configure logging to write to `stderr` (recommended for scripts/tests):

```julia
include("src/utils.jl")
setup_stderr_logger()
```

2. Use the `Logging` macros in your code or scripts:

```julia
@info "Starting optimization" params=some_params
@warn "Converged slowly" iter=100
println(stderr, "Important non-progress message")
```

Notes
-----
- The previous design serialized progress updates through a central manager
   task. That has been intentionally removed. If you need serialized output for
   a future multithreaded implementation, consider a simple channel/consumer
   printer or a ReentrantLock-protected print helper.
- `setup_stderr_logger()` is provided to make logging from library code appear
   on `stderr` by default; call it from your scripts or tests when desired.

Example test harness
--------------------
See `test/thread_progress_test.jl` for a tiny single-threaded example using
`Logging` to emit progress-like messages.
