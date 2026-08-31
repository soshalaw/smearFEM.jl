# smearFEM.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://soshalaw.github.io/smearFEM.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://soshalaw.github.io/smearFEM.jl/main/)
[![Build Status](https://github.com/soshalaw/smearFEM.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/soshalaw/smearFEM.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/soshalaw/smearFEM.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/soshalaw/smearFEM.jl)
[![Aqua](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

A model parameter estimation framework developed for the SMEAR (Soft Matter Exploring Autonomous Robot) project.

smearFEM simulates the **squeeze flow** of a soft object compressed between two plates using a
mixed finite-element Stokes solver, renders the deforming silhouette through a camera model, and
then **inverts** that pipeline: given an observed border contour over time, it recovers the
material parameters `θ = [η, β]` (shear viscosity and a boundary slip/friction penalty).

The forward solve carries analytic sensitivities `∂q/∂η` and `∂q/∂β` alongside the solution, so
the inverse problem is solved with true Gauss–Newton / Levenberg–Marquardt steps rather than
finite differences.

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/soshalaw/smearFEM.jl.git")
```

Or for local development:

```julia
using Pkg
Pkg.develop(path="/path/to/smearFEM.jl")
Pkg.instantiate()
```

**Requirements:** Julia 1.6.7 or later. All dependencies are listed in `Project.toml` and installed
automatically. Gmsh ships as an artifact via `Gmsh.jl` — no system install needed.

## Package layout

| Path | Contents |
|------|----------|
| `src/fem/` | Element basis functions and quadrature (`fem.jl`), mesh types and generators (`meshes/`), model and scenario structs, geometry wrapper types, camera projection / border extraction (`PostProcess.jl`) |
| `src/solver/` | Stokes assembly, boundary conditions, time-stepping (`stokes_solver.jl`), problem construction (`stokes_setup.jl`), high-level drivers (`run_example.jl`) |
| `src/optimization/` | Gauss–Newton and Levenberg–Marquardt parameter fitting with Armijo / backtracking line search |
| `src/analysis/` | Point-cloud metrics (Hausdorff, Chamfer) and experiment-tree traversal |
| `src/io/` | CSV / JSON / HDF5 / VTK readers and writers, plotting and animation helpers |
| `experiments/` | Reproduction scripts — convergence studies, optimization batches, visualizations |
| `mesh_files/` | Gmsh `.geo` templates plus the generated-mesh cache (see [`mesh_files/USAGE.md`](mesh_files/USAGE.md)) |
| `test/` | Unit and integration tests |

## Core concepts

**Geometry** — the object being squeezed is described by a wrapper type, which selects the
meshing and `def_problem` method by dispatch:

```julia
Cylinder(r, h)                     # 3D
Cuboid(lx, ly, lz, edge_radius)    # 3D; edge_radius=nothing gives sharp vertical edges
Disk(r)                            # 2D
Square(lx, ly)                     # 2D
Segment(l)                         # 1D
```

**Meshing** — every geometry can be meshed structurally (pure Julia) or unstructurally (Gmsh,
`mesh_type=:unstructured`). Element shapes are `:Hex` / `:Tet` in 3D, `:Quad` / `:Tri` in 2D, and
`:Line` in 1D, at `basis_order` 1 or 2. Meshes are cached on disk under
`mesh_files/<geometry>/<Shape>/`; delete the cache after editing a `.geo` template.

**Mixed formulation** — a `Stokes` model holds three meshes: `mesh_u` (velocity), `mesh_p`
(pressure), and `mesh_x` (geometry / node tracking). The usual Taylor–Hood choice is quadratic
velocity over linear pressure (`:Hex` order 2 for `u`, `:Hex` order 1 for `p`).

**Scenario** — a `SqueezeFlow` bundles the model with the slip parameter `β`, Dirichlet boundary
operators, and the control history. Control is either `"force"` (a prescribed load `cParam`, with
top-plate displacement solved as an unknown) or `"velocity"` (a prescribed top-plate velocity).
Viscosity is `"constant"` or `"bulk_viscosity"` (time-varying, optionally from the power-law model).

**Conditions** — a `Conditions` struct carries the camera matrix, object pose, viewing angles,
output path, and the animate / render / VTK / contour output flags.

## Quick start

```julia
using smearFEM

camera_matrix = get_camera_matrix()
obj_pose      = [150.0, 0.0, 20.0]

sim_time, t_steps = 5.0, 2.5e-3
F = -9518.61 * ones(round(Int, sim_time / t_steps))   # applied load per time step

model, scene = def_problem(
    Cylinder(25.0, 40.0),         # geometry (mm)
    6, 100.0,                     # ne (target elements along the height), η
    :Hex, 2, 3,                   # velocity: shape, order, dofs/node
    :Hex, 1, 1,                   # pressure: shape, order, dofs/node
    :Hex, 2,                      # geometry mesh: shape, order
    100.0, F, "force", "constant",
    sim_time, t_steps)

conditions = Conditions(camera_matrix=camera_matrix, obj_pose=obj_pose,
                        ANIMATE=true, WRITEVTK=true, filepath="/tmp/smear/demo")

h, gradList, borderPts2D, displacement, surface_pts_3D,
    pos2D, pos3D, velocity, pressure, _ = simulate(model, scene, conditions)
```

To run the forward model and write the full result tree (fields, VTK, plots, `sim_params.json`) in
one call, use `write_sim_data(model, scene, camera_matrix, obj_pose, z_angles, filepath)`, or drive
it from a parameter dictionary with `write_gt_data(exp_params)`.

### Inverse problem

Given observed 2D border points per frame, recover `[η, β]`:

```julia
θ_init = [50.0, 10.0]

stats = fit_model(model, scene, conditions, obsBorderPts, θ_init;
                  method=:gn,                    # :gn (Gauss–Newton) or :lm (Levenberg–Marquardt)
                  line_search_method=:armijo)    # :armijo or :backtrack (:gn only)

stats["η"], stats["β"]        # fitted parameters
stats["cost_list"]            # cost per iteration
```

`fit_model` returns a `Dict` with keys `"η"`, `"β"`, `"ηList"`, `"βList"`, `"cost_list"`, and
`"iterList"`. Parameters are clamped to physical bounds each iteration — η to `[1e-3, 1e5]` Pa·s
and β to `[1e-3, 1e8]` L/mm. Pass `outliers` to skip frame indices.

## Configuration

smearFEM resolves data paths via a three-tier system (highest priority first):

1. **Environment variables**: `SMEAR_DATA_DIR`, `SMEAR_MESH_DIR`, `SMEAR_SCRATCH_DIR`
2. **Config file**: `config.toml`, searched for from the package root upwards. It lives in the
   `smear-modules` directory alongside this repo, so smearFEM.jl and smearPerception share one
   copy — start from `smear-modules/config.toml.example`:
   ```toml
   data_dir    = "../SMEAR-DataFiles/Data"
   scratch_dir = "/tmp/smear"
   ```
   Relative paths are resolved against the directory holding `config.toml`, and `${VAR}`
   references are expanded from the environment. Set `SMEAR_CONFIG_FILE` to point elsewhere.
3. **Defaults**: `~/SMEAR-Data` for data, `<data_dir>/meshes` for meshes, `/tmp/smear` for scratch
   (`~/SMEAR-Scratch` on Windows).

A tier is skipped with a warning if the directory it names does not exist.

```julia
using smearFEM
show_config()                                    # print all resolved directories
resolve_data_path("ground_truth/sim_data/...")   # → absolute path under data_dir
resolve_mesh_path("channel.msh")                 # → absolute path under mesh_dir
get_scratch_dir()                                # → scratch directory
create_output_dir("my_experiment")               # → created directory under scratch
```

## Reproducing experiments

Scripts in `experiments/` are run directly; each has a `main()` that fires only when the file is
executed as a script, so they can also be `include`d from the REPL. Some pull in packages beyond
`Project.toml` (e.g. `CurveFit`, `Dierckx` for the convergence plots) — add those to your
environment first.

### Convergence analysis

```bash
julia -O3 --threads=auto --project experiments/convergence_analysis/forward_convergence.jl
julia -O3 --threads=auto --project experiments/convergence_analysis/inverse_convergence.jl
```

**Output**: mesh- and time-convergence plots with fitted rates, under
`<data_dir>/experiments/sim_data/convergence_analysis/`. **Duration**: 5–15 minutes for the
forward study; `inverse_convergence.jl` only plots results already on disk.

### Parameter optimization

```bash
julia -O3 --threads=auto --project experiments/model_optimization/test_opt_stokes.jl
```

Runs the simulated, synthetic, and physical optimization batches, then plots results. Key
parameters (`β_gt`, `η_gt`, `control`, tolerances) are set at the top of the script.
**Duration**: 30 minutes to 2 hours.

### Ground-truth generation and single runs

```bash
julia -O3 --threads=auto --project experiments/stokes/write_gt_stokes.jl     # batch, parallel
julia -O3 --threads=auto --project experiments/stokes/single_simulation.jl   # one forward solve
julia -O3 --threads=auto --project experiments/visualizations/visualization_example.jl
julia -O3 --threads=auto --project experiments/verify_meshes.jl              # mesh sanity checks
```

`experiments/ParallelExecution.jl` is a shared, memory-aware work-stealing scheduler used by the
batch scripts; `run_parallel_tasks(param_list, task; max_workers=-1, memory_per_task_mb=512.0)`
sizes the worker pool against available RAM.

### Performance

```bash
julia -O3 --threads=auto --project script.jl
```

BLAS is set to `Threads.nthreads()` at module load, so `--threads=auto` benefits both the
assembly loops and the dense linear algebra.

## Testing

```julia
using Pkg
Pkg.test("smearFEM")
```

Or individual suites:

```bash
julia --project -e 'include("test/stokes_test.jl")'
julia --project -e 'include("test/fem_test.jl")'
julia --project -e 'include("test/mesh_templates_test.jl")'
julia --project -e 'include("test/stokes_optimization_test.jl")'
```

`test/qa.jl` runs [Aqua.jl](https://github.com/JuliaTesting/Aqua.jl) quality checks.

## Documentation

- [API reference](https://soshalaw.github.io/smearFEM.jl/dev/)
- [`mesh_files/USAGE.md`](mesh_files/USAGE.md) — mesh template usage, element-count calibration, caching
- [`CONTRIBUTING.md`](CONTRIBUTING.md) — development setup, code style, PR process
- [`docs/coding_guide.md`](docs/coding_guide.md) — internal conventions

Build the docs locally with `julia --project=docs docs/make.jl`.

## License

MIT License — see [LICENSE](LICENSE) for details.

## Support

- **GitHub Issues**: [Report a bug or request a feature](https://github.com/soshalaw/smearFEM.jl/issues)
- **Documentation**: [API Reference](https://soshalaw.github.io/smearFEM.jl/dev/)
