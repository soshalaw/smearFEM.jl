# smearFEM.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://soshalaw.github.io/smearFEM.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://soshalaw.github.io/smearFEM.jl/dev/)
[![Build Status](https://github.com/soshalaw/smearFEM.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/soshalaw/smearFEM.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/soshalaw/smearFEM.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/soshalaw/smearFEM.jl)
[![Aqua](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

A finite element method (FEM) package for solving fluid flow and solid mechanics problems for the SMEAR (Soft Matter Exploring Autonomous Robot) project.

## Features

- **Stokes Flow**: Model viscous incompressible fluid flow in 2D and 3D
- **Shape Optimization**: Inverse design and parameter optimization for soft robotics applications
- **Mesh Generation**: Integration with Gmsh for robust mesh creation with thread-safe parallel support
- **Isogeometric Analysis (IGA)**: NURBS-based analysis and extraction capabilities
- **Parallel Computing**: Multi-threaded performance optimization for large-scale simulations
- **VTK Output**: Visualization-ready results export for post-processing
- **Linear Elasticity**: Support for 2D/3D linear elastic solid mechanics problems

## Installation

This package is distributed as supplementary material with the publication. To use it:

**1. Clone or download the repository:**
```bash
git clone https://github.com/soshalaw/smearFEM.jl.git
cd smearFEM.jl
```

**2. Install Julia dependencies:**
```julia
using Pkg
Pkg.instantiate()
```

**3. Load the package in your Julia session:**
```julia
include("src/smearFEM.jl")
using .smearFEM
```

Or run examples directly:
```julia
include("src/examples/squeeze_stokes.jl")
```

**Requirements:**
- Julia 1.6.7 or later
- All dependencies listed in `Project.toml` (automatically installed via `Pkg.instantiate()`)

## Quick Start

### Stokes Flow Simulation

```julia
using smearFEM

# Define fluid domain and boundary conditions
geometry = meshgrid_square(nx=20, ny=20)
model = Stokes(geometry, viscosity=1.0)

# Solve Stokes problem
velocity, pressure = simulate_stokes(model, boundary_conditions)

# Export results to VTK
write_vtk("stokes_results.vtu", solution)
```

### Shape Optimization for Soft Robotics

```julia
using smearFEM

# Define optimization parameters
initial_design = create_initial_shape()
objective = setup_stokes_objective()

# Run inverse design optimization
optimal_design = fit_model(initial_design, objective, constraints)
```

## Documentation

- **[Stable Documentation](https://soshalaw.github.io/smearFEM.jl/stable/)** — Released version docs
- **[Development Documentation](https://soshalaw.github.io/smearFEM.jl/dev/)** — Latest development version

Full API reference and usage examples are available in the [documentation](https://soshalaw.github.io/smearFEM.jl/dev/).

## Examples

Complete working examples are provided in the `src/examples/` directory:

- **`squeeze_stokes.jl`** — 2D Stokes flow between deformable surfaces (primary SMEAR application)
- **`fluid_flow_stokes.jl`** — 3D incompressible flow simulation
- **`squeeze_linear_elasticity.jl`** — 2D linear elasticity under compression
- **`shear.jl`** — Shear deformation example
- **`run_example.jl`** — Interactive example runner

Run a Stokes example:

```julia
include("src/examples/squeeze_stokes.jl")
```

## Testing

Run the full test suite:

```julia
using Pkg
Pkg.test("smearFEM")
```

Or run specific tests:

```julia
julia --project -e "include(\"test/linear_elasticity_test.jl\")"
julia --project -e "include(\"test/stokes_test.jl\")"
```

## Reproducing Experiments

This section describes how to reproduce the main computational experiments from the publication.

### 1. Convergence Analysis

Analyze mesh and time-step convergence of Stokes FEM solutions:

```bash
cd test/convergence_analysis/
julia -O3 --threads=auto --project stokes_convergence.jl
```

**Output**: Convergence plots and convergence rate analysis  
**Duration**: 5-15 minutes depending on mesh refinement levels

**Key Files**:
- `test/convergence_analysis/stokes_convergence.jl` — Main convergence study script
- `test/convergence_analysis/mesh.geo` — Test mesh definition

### 2. Shape Optimization (Model Inverse Design)

Run parameter identification and inverse design for Stokes flow:

```bash
cd scripts/model_optimization/
julia -O3 --threads=auto --project test_opt_stokes.jl
```

**Output**: Optimization convergence plots, optimized parameters, model accuracy metrics  
**Duration**: 30 minutes to 2 hours (depending on problem size and iterations)

**Key Parameters** (set in script):
- `β_gt`: Stiffness parameter (ground truth value)
- `η_gt`: Viscosity (ground truth value)
- `control`: "force" or "velocity" boundary conditions
- Optimization method: damped Newton with line search

### 3. Single Stokes Simulation

Run a single forward simulation of Stokes flow in squeeze geometry:

```bash
cd scripts/stokes/
julia -O3 --threads=auto --project single_simulation.jl
```

**Output**: Velocity and pressure field VTK files  
**Duration**: 1-5 minutes

### Parallel Execution

For optimal performance, use multi-threading:

```bash
julia -O3 --threads=auto --project script.jl
```

For distributed parallel optimization:

```bash
julia -O3 --threads=auto -p 4 --project test_opt_stokes.jl
```

where `-p 4` adds 4 worker processes.

### Performance Tips

- Use optimization level `-O3` for production runs
- Use `--threads=auto` to utilize all available CPU cores
- For large problems, enable distributed computing with `-p`
- Results scale approximately linearly with thread count

### Customizing Experiments

Key parameters are defined at the top of each script:

- **Mesh size**: `ne` (number of elements per dimension)
- **Simulation time**: `sim_time` (seconds)
- **Time step size**: `t_steps` (seconds per step)
- **Optimization tolerance**: Relative cost change threshold (~0.01%)
- **Viscosity and stiffness**: `η_gt`, `β_gt` parameters

## Testing

Run the full test suite:

```julia
using Pkg
Pkg.test("smearFEM")
```

Or run specific tests:

```julia
julia --project -e "include(\"test/linear_elasticity_test.jl\")"
julia --project -e "include(\"test/stokes_test.jl\")"
```

## Requirements

- **Julia**: 1.6.7 or later
- **Dependencies**: See `Project.toml` for the complete list
  - LinearAlgebra, SparseArrays (standard library)
  - Gmsh for mesh generation
  - IterativeSolvers for large-scale linear systems
  - WriteVTK for result export

## Performance

For optimal performance, use:

```bash
julia -O3 --threads=auto --project script.jl
```

Performance guidelines and benchmarking tools are available in the `performance/` directory.

## Contributing

We welcome contributions! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for:
- Development setup
- Code style guidelines
- Testing requirements
- Pull request process

## License

This project is licensed under the MIT License - see [LICENSE](LICENSE) for details.

## Citation

If you use smearFEM.jl in your research, please cite:

```bibtex
@software{smearFEM,
  title={smearFEM.jl: A Julia Package for SMEAR-Based Finite Element Analysis},
  author={Weerathunge, Soshala},
  url={https://github.com/soshalaw/smearFEM.jl},
  year={2024}
}
```

## Support

For issues, questions, or feedback:
- **GitHub Issues**: [Report a bug or request a feature](https://github.com/soshalaw/smearFEM.jl/issues)
- **Documentation**: [API Reference](https://soshalaw.github.io/smearFEM.jl/dev/)

## Project Status

This package is actively maintained. The `main` branch contains the stable release, while development occurs on feature branches.
