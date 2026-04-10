# Benchmarking & Profiling Suite

Comprehensive performance analysis tools for the Stokes solver.

## Quick Start

Run all benchmarks:
```bash
cd test/benchmarking
julia -O3 --threads=auto --check-bounds=no --project run_all_benchmarks.jl
```

## Available Tools

### 1. run_all_benchmarks.jl (Recommended)
**Purpose:** Master benchmarking orchestrator that runs all profilers in sequence

**Run:**
```bash
julia -O3 --threads=auto --project run_all_benchmarks.jl
```

**Output:**
- Comprehensive summary of all benchmarks
- Execution timing for each profiler
- Success/failure status
- Total performance metrics

---

### 2. detailed_profiler.jl
**Purpose:** Component-level performance analysis with detailed benchmark measurements

**Measures:**
- Assembly function performance
- Individual operation timing
- Memory allocation patterns
- Relative speedups

**Run:**
```bash
julia -O3 --project detailed_profiler.jl
```

---

### 3. profile_100hz.jl
**Purpose:** Validates that solver achieves 100Hz execution target (10x realtime)

**Validates:**
- End-to-end execution time
- Target: ~0.2s for 20-second simulation
- Performance scaling with threads
- Configuration verification

**Run:**
```bash
julia -O3 --threads=auto --project profile_100hz.jl
```

---

### 4. iterative_solver_profiler.jl
**Purpose:** Benchmarks iterative solver (GMRES + ILU) performance

**Measures:**
- GMRES iteration count per solve
- Solve time per iteration
- Comparison: Direct LU vs iterative
- Overall speedup

**Run:**
```bash
julia -O3 --project iterative_solver_profiler.jl
```

---

### 5. profileview_profiler.jl
**Purpose:** Generates flamegraph visualization of execution hotspots

**Features:**
- Visual call stack profiling
- Identifies bottlenecks
- Shows time distribution
- Requires ProfileView.jl

**Run:**
```bash
julia --project profileview_profiler.jl
```

**Note:** Produces interactive flamegraph for bottleneck analysis

---

## Recommended Workflow

1. **Full Suite (Recommended)**
   ```bash
   julia -O3 --threads=auto --project run_all_benchmarks.jl
   ```
   Get complete performance picture in one run

2. **Quick Validation**
   ```bash
   julia -O3 --project profile_100hz.jl
   ```
   Verify 100Hz target is met

3. **Deep Component Analysis**
   ```bash
   julia -O3 --project detailed_profiler.jl
   ```
   Identify individual component bottlenecks

4. **Solver-Specific Tuning**
   ```bash
   julia -O3 --project iterative_solver_profiler.jl
   ```
   Compare different solver strategies

5. **Visualization & Hotspot Analysis**
   ```bash
   julia --project profileview_profiler.jl
   ```
   See where time is actually spent

---

## Performance Targets

| Metric | Target | Status |
|--------|--------|--------|
| Simulation Speed | 100Hz (1 ms per ms of simulation) | Pass |
| Assembly Time | < 50% of total | Pass |
| Solver Time | < 40% of total | Pass |
| Memory Efficiency | No allocation in hot loops | Pass |

---

## Optimization Flags

For maximum performance, run with:
```bash
julia -O3 --threads=auto --check-bounds=no \
       --project run_all_benchmarks.jl
```

**Flags:**
- `-O3`: Aggressive compiler optimizations
- `--threads=auto`: Use all available cores
- `--check-bounds=no`: Disable array bounds checking

---

## Interpreting Results

### Detailed Profiler Output
- Times shown in seconds (s) or milliseconds (ms)
- Large allocation values indicate GC pressure
- Compare "min" time vs "median" to detect GC spikes

### 100Hz Validator Output
- Target: ~0.2s for 20-second simulation (100 timesteps)
- Scaling: Should see ~Nthreads speedup
- If below target: Architecture-dependent limits

### Solver Profiler Output
- GMRES iterations: < 50 typical, > 200 problematic
- Speedup factors: 1.0 = baseline (direct LU)

### Flamegraph
- Wider bars = more time spent
- Look for unexpected wide bars = optimization opportunity

---

## Troubleshooting

**"Module not found"**
```bash
julia --project  # Run from workspace root
```

**Very slow benchmarks**
- Ensure `-O3` flag is used
- Check `--check-bounds=no` is set
- May need `JULIA_NUM_THREADS=auto`

**Memory issues**
- Reduce problem size in profiler scripts
- Run one profiler at a time
- Check available RAM

---

## Directory Structure

```
test/benchmarking/
├── run_all_benchmarks.jl          Master orchestrator
├── detailed_profiler.jl            Component analysis
├── profile_100hz.jl                100Hz validation
├── iterative_solver_profiler.jl    Solver tuning
├── profileview_profiler.jl         Flamegraph visualization
└── README.md                       This file
```

---

## Performance References

Previous optimization work:
- `performance/OPTIMIZATION_SUMMARY.md` - Quick reference
- `performance/PERFORMANCE_100HZ_GUIDE.md` - Detailed guide
- `performance/run_100hz.sh` - Execution wrapper

---

**Last Updated:** April 10, 2026
