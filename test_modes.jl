using smearFEM
include("scripts/stokes/write_gt_stokes.jl")

# This file will be used to test both modes
# Mode 1 (default - without calibration):
# main(use_parallel=true, calibrate=false, memory_per_experiment_mb=512.0)

# Mode 2 (with calibration):
# main(use_parallel=true, calibrate=true, memory_per_experiment_mb=512.0)

println("✓ Script loaded successfully and split functions are accessible")
println("\nAvailable modes:")
println("  1. WITHOUT calibration (default):  main(use_parallel=true, calibrate=false)")
println("  2. WITH calibration (optional):   main(use_parallel=true, calibrate=true)")
println("  3. Sequential (safest):            main(use_parallel=false, calibrate=false)")
