# smearFEM.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://soshalaw.github.io/smearFEM.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://soshalaw.github.io/smearFEM.jl/main/)
[![Build Status](https://github.com/soshalaw/smearFEM.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/soshalaw/smearFEM.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/soshalaw/smearFEM.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/soshalaw/smearFEM.jl)
[![Aqua](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

A model parameter estimation framework developed for the SMEAR (Soft Matter Exploring Autonomous Robot) project.

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/soshalaw/smearFEM.jl.git")
```

Or for local development:

```julia
using Pkg
Pkg.develop(path="/path/to/smearFEM.jl")
```

**Requirements:** Julia 1.6.7 or later. All dependencies are listed in `Project.toml` and installed automatically.

## Configuration

smearFEM resolves data paths via a three-tier system (highest priority first):

1. **Environment variable**: `export SMEAR_DATA_DIR=/path/to/data`
2. **Config file**: create `config.toml` in the package root (gitignored):
   ```toml
   data_dir    = "/path/to/SMEAR-DataFiles/Data"
   scratch_dir = "/path/to/tmp"
   ```
3. **Default**: `~/SMEAR-Data`

```julia
using smearFEM
resolve_data_path("ground_truth/sim_data/...")  # → full path via the resolver
get_scratch_dir()                                # → output/temp directory
```

## Reproducing Experiments

### Convergence Analysis

```bash
julia -O3 --threads=auto --project test/convergence_analysis/stokes_convergence.jl
```

**Output**: Convergence plots and rate analysis. **Duration**: 5–15 minutes.

### Parameter Optimization

```bash
julia -O3 --threads=auto --project scripts/model_optimization/test_opt_stokes.jl
```

**Output**: Optimization convergence, optimized parameters, model accuracy metrics.  
**Duration**: 30 minutes to 2 hours.

Key parameters (set at top of script): `β_gt`, `η_gt`, `control`, optimization tolerance.

### Single Forward Simulation

```bash
julia -O3 --threads=auto --project scripts/stokes/single_simulation.jl
```

### Performance Tips

```bash
julia -O3 --threads=auto --project script.jl   # multi-threaded BLAS + Gmsh
```

## Testing

```julia
using Pkg
Pkg.test("smearFEM")
```

Or individual test files:

```bash
julia --project -e "include(\"test/stokes_test.jl\")"
```

## License

MIT License — see [LICENSE](LICENSE) for details.

## Support

- **GitHub Issues**: [Report a bug or request a feature](https://github.com/soshalaw/smearFEM.jl/issues)
- **Documentation**: [API Reference](https://soshalaw.github.io/smearFEM.jl/dev/)
