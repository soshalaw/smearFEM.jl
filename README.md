# smearFEM.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://soshalaw.github.io/smearFEM.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://soshalaw.github.io/smearFEM.jl/main/)
[![Build Status](https://github.com/soshalaw/smearFEM.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/soshalaw/smearFEM.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/soshalaw/smearFEM.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/soshalaw/smearFEM.jl)
[![Aqua](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

A real-time linear elasticity solver developed for the SMEAR (Soft Matter Exploring Autonomous
Robot) project.

This branch (`features/linear_elasticity`) holds the linear elastic model on its own, without the
Stokes squeeze-flow model and the parameter-optimization code of the main line. It provides:

- **Squeeze** (`simulate`): an object compressed between two rigid plates, one small-strain load
  increment per step on the current geometry (updated Lagrangian), under force or displacement
  control, with Navier slip on the plates and optional sensitivities of the node positions to
  `(λ, μ, β)`.
- **Volume update** (`volume_update!`): the mesh moved so that chosen surface nodes reach
  prescribed normal displacements, the interior following a compressible linear elastic solve —
  the volume update of the paper-2 pipeline.

Both keep one fixed-sparsity system per mesh and reuse its factorization, so a load step costs a
few milliseconds (see [Performance](#performance)).

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

**Requirements:** Julia 1.11 or later. All dependencies are listed in `Project.toml` and installed
automatically. Gmsh ships as an artifact via `Gmsh.jl` — no system install needed.

## Package layout

| Path | Contents |
|------|----------|
| `src/smearFEM.jl` | Module entry: the single export list and include order |
| `src/config.jl`, `src/utils.jl` | Data/scratch path resolution (`config.toml`), shared helpers (camera matrix, logging) |
| `src/fem/` | Element basis functions and quadrature (`fem.jl`), mesh types and generators (`meshes/`), the `LinearElasticity` model, camera projection / border extraction (`PostProcess.jl`) |
| `src/solver/` | Fixed-pattern stiffness assembly and refined solves (`elastic_system.jl`), the squeeze (`elastic_squeeze.jl`), the volume update (`volume_update.jl`), the sparse direct solver backend (`sparse_lu.jl`) |
| `src/analysis/` | Point-cloud metrics (Hausdorff, Chamfer) |
| `src/io/` | CSV / JSON / HDF5 / VTK readers and writers, plotting and animation helpers |
| `experiments/` | Scripts — the solver benchmark (`elasticity/benchmark.jl`) |
| `mesh_files/` | Gmsh `.geo` templates plus the generated-mesh cache (see [`mesh_files/USAGE.md`](mesh_files/USAGE.md)) |
| `test/` | Unit and integration tests |

## Core concepts

**Meshing** — a cylinder or box is meshed structurally (pure Julia, `:Hex`) or unstructurally
(Gmsh, `mesh_type=:unstructured`, `:Hex` or `:Tet`), at `basis_order` 1 or 2. Gmsh meshes are
cached in the repository under `mesh_files/<geometry>/<Shape>/` (not under the configured data
directory); delete the cache after editing a `.geo` template.

**Model** — `LinearElasticity(mesh=..., λ=..., μ=...)` is an isotropic linear elastic material on
one isoparametric 3D mesh. Lengths are in mm and forces in kg·mm/s² (1 N = 1e3), so stresses and
the Lamé parameters are in kPa; `lame_parameters(E, ν)` converts from Young's modulus and
Poisson's ratio. The mesh's `NodeList` is the current geometry: `simulate` starts from the
initial state (`reset_model!`) and leaves the mesh at the final one; `volume_update!` moves it
from wherever it is (`update_model!` makes the current geometry the new initial state).

**Squeeze** — the bottom plate holds `u_z = 0`, the top plate moves all its nodes rigidly in z, and
both plates carry Navier slip `β` on the tangential components. Each entry of `cParam` is one
load increment: of the plate force (`control="force"`) or of the plate displacement
(`control="displacement"`). `β = 0` leaves the horizontal rigid-body modes free; use a small
positive value for frictionless plates.

**Volume update** — `volume_update!(mdl, nodes, normals, d)` prescribes `u·n = d` at each listed
node (tangential motion free), `u_z = 0` on the table (the bottom face) and leaves the rest of the
surface traction-free. Observed and rigid-contact nodes take the same form, and a node may carry
several constraints. With displacement constraints only, the result depends on Poisson's ratio
alone. On an axisymmetric surface with radial normals nothing pins the rotation about z: pass a
small `β_table` (table slip, `[not decided]`), or the factorization warns that the stiffness is
singular and the result carries an arbitrary rotation.

## Quick start

```julia
using smearFEM

mesh = meshgrid_cylinder(25.0, 40.0; mesh_type=:unstructured, element_shape=:Hex, basis_order=2,
                         elem_size=6.0, mesh_path=joinpath(pkgdir(smearFEM), "mesh_files"))
λ, μ = lame_parameters(10.0, 0.35)                  # kPa
mdl = LinearElasticity(mesh=mesh, λ=λ, μ=μ)

# squeeze: 100 increments of -0.05 mm, with ∂positions/∂(λ, μ, β) per step
r = simulate(mdl, fill(-0.05, 100); control="displacement", β=1.0, sensitivities=true)
r.plate          # plate force increment per step
r.pos3D          # node positions, initial state first
r.grad3D         # (3, nNodes, 3) sensitivities per state

# volume update: push the lateral surface 0.5 mm outward
reset_model!(mdl)
X = mesh.NodeList
side = mesh.side_nodes
normals = vcat(X[1:2, side] ./ sqrt.(sum(abs2, X[1:2, side], dims=1)), zeros(1, length(side)))
u = volume_update!(mdl, side, normals, fill(0.5, length(side)); β_table=1e-3)
```

Sensitivities beyond the first step hold the geometry fixed within each step: exact at the first
step, and against finite differences of the full run about 0.5 % off after 0.6 % compression and
3 % after 12.5 % (25 steps, test cylinder). Write VTK output with `write_scene` / `write_vtk`.

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

`get_mesh_dir()` / `SMEAR_MESH_DIR` is resolved for external tools that share `config.toml`; the
Gmsh mesh cache itself always lives in the repository's `mesh_files/`.

A tier is skipped with a warning if the directory it names does not exist.

```julia
using smearFEM
show_config()                                    # print all resolved directories
resolve_data_path("ground_truth/sim_data/...")   # → absolute path under data_dir
resolve_mesh_path("channel.msh")                 # → absolute path under mesh_dir
get_scratch_dir()                                # → scratch directory
create_output_dir("my_experiment")               # → created directory under scratch
```

## Performance

```bash
julia -O3 --threads=auto --project experiments/elasticity/benchmark.jl [Hex|Tet] [order] [elem_size]
```

The sparsity pattern, the element-to-matrix gather lists and the factorization's symbolic
analysis are built once per mesh (0.15–0.5 s, cached across calls; the first call of a session
also compiles, ~1.5 s). A squeeze step then refills the stiffness at the current geometry in
parallel (~2 ms) and solves by iterative refinement against a factorization reused across steps.
The first step of every call factorizes afresh, so results never depend on earlier calls. The
factorization uses MKL Pardiso where MKL is available (x86-64 Linux and Windows), refactorizing
(~9 ms) once the mesh has deformed enough, typically every ~10 steps, and UMFPACK elsewhere —
every Mac and ARM machine — where a factorization takes ~150 ms.

Measured on an 8-thread laptop (other applications running), on the paper-1 cylinder (Gmsh,
element size 6 mm), 100 displacement-controlled steps of 0.05 mm (12.5 % compression); force
control (10.9 %) gives the same numbers within 0.5 ms. Medians and maxima are over steps 2–100:

| Mesh, backend | Squeeze step, median (max) | With sensitivities, median (max) | First step | Volume update |
|---------------|----------------------------|----------------------------------|------------|---------------|
| Hex Q2, 1,885 nodes, Pardiso | 6.9 ms (17 ms) | 15.0 ms (28 ms) | 17–27 ms | 20 ms |
| Tet P2, 2,555 nodes, Pardiso | 6.5 ms (17 ms) | 14.4 ms (27 ms) | 18–27 ms | 20 ms |
| Hex Q2, UMFPACK | 7.7 ms (15 ms) | 22 ms (40 ms) | 130–160 ms | 130 ms |
| Tet P2, UMFPACK | 7.8 ms (20 ms) | 24 ms (42 ms) | 155–175 ms | 160 ms |

With Pardiso and without sensitivities every step after the first stays within 20 ms (≥ 50 Hz);
with sensitivities ~7 % of steps exceed it (the steps that refactorize). With UMFPACK the steps
are as fast, but sensitivities push about half of them over 20 ms, and the first step and the
volume update each cost a full factorization.

## Testing

```julia
using Pkg
Pkg.test("smearFEM")
```

The suite takes about a minute. Individual files:

```bash
julia --project -e 'include("test/linear_elasticity_test.jl")'
julia --project -e 'include("test/fem_test.jl")'
julia --project -e 'include("test/mesh_templates_test.jl")'
```

`linear_elasticity_test.jl` checks the element matrices against `B C Bᵀ`, the squeeze against
the analytic uniaxial solution on Q2 hex and P1/P2 tet meshes (the force-controlled case also
guards the simplex quadrature weights), the sensitivities against finite differences, the
volume update against an exact linear field, and that a singular stiffness is reported. `test/qa.jl` runs
[Aqua.jl](https://github.com/JuliaTesting/Aqua.jl) quality checks; it needs the test environment,
so run it through `Pkg.test()`.

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
