#!/bin/bash

# ============================================================================
# Ultra-High Performance Runner for 100Hz Execution (10x realtime)
# ============================================================================
# This script runs single_simulation.jl with all aggressive optimizations
# Target: 20-second simulation in ~0.2 seconds
# ============================================================================

set -e  # Exit on any error

echo "╔════════════════════════════════════════════════════════════════════╗"
echo "║  100HZ ULTRA-HIGH PERFORMANCE RUNNER                             ║"
echo "║  Target: 10x faster than realtime (0.2 sec for 20-sec sim)        ║"
echo "╚════════════════════════════════════════════════════════════════════╝"

# Step 1: Detect available cores
NCORE=$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4)
echo "✓ Detected CPU cores: $NCORE"

# Step 2: Set optimal environment variables
export OMP_NUM_THREADS=$NCORE
export OPENBLAS_NUM_THREADS=$NCORE
export MKL_NUM_THREADS=$NCORE
export JULIA_CPU_THREADS=$NCORE
export JULIA_NUM_THREADS=$NCORE
export JULIA_THREAD_POOL=interactive

echo "✓ Environment variables set for $NCORE threads"

# Step 3: Run with maximum optimization
echo "✓ Starting simulation with -O3 optimization..."
echo "  Command: julia -O3 --threads=auto --check-bounds=no scripts/stokes/single_simulation.jl"
echo ""

time julia -O3 --threads=$NCORE --check-bounds=no \
    -e "
    println(\"╔════════════════════════════════════════════════════════════════════╗\")
    println(\"║  PERFORMANCE SUMMARY                                             ║\")
    println(\"║  - Threads: $NCORE                                               ║\")
    println(\"║  - Optimization Level: -O3                                        ║\")
    println(\"║  - Bounds Checking: Disabled                                      ║\")
    println(\"║  - Expected Speedup: 8-20x                                        ║\")
    println(\"║  - Target: <200ms for 20-second simulation                        ║\")
    println(\"╚════════════════════════════════════════════════════════════════════╝\")
    " \
    $(git rev-parse --show-toplevel)/scripts/stokes/single_simulation.jl

EXIT_CODE=$?

echo ""
echo "╔════════════════════════════════════════════════════════════════════╗"
if [ $EXIT_CODE -eq 0 ]; then
    echo "║  ✓ Simulation completed successfully!                            ║"
else
    echo "║  ✗ Simulation failed with exit code $EXIT_CODE                   ║"
fi
echo "╚════════════════════════════════════════════════════════════════════╝"

exit $EXIT_CODE
