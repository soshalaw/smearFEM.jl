# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

smearFEM.jl is a Julia FEM library for the SMEAR robotics project. It simulates **squeeze flow**
of a soft object compressed between two plates with a mixed-element Stokes solver, renders the
deforming silhouette through a camera model, and inverts that pipeline to recover material
parameters `θ = [η, β]` (shear viscosity, boundary slip/friction) from observed border contours,
via Gauss-Newton / Levenberg-Marquardt using analytic sensitivities.

`README.md` has the full quick-start, core-concept walkthrough, and config system. Read it before
making non-trivial changes — this file covers commands and architecture that isn't obvious from a
single file.

This repo is pure Julia — every convention below is Julia-specific. The sibling Python package,
`../smearPerception`, has its own `CLAUDE.md` with different (sometimes opposite, e.g. `print` vs
`@debug`/`@info`) style conventions; don't carry this file's rules over there. The two packages
share one `config.toml` (`../config.toml`) and mirror each other's path-resolution scheme
(`src/config.jl` here, `src/config.py` there) — check both sides if you touch config resolution.

## Two documentation layers (user, 2026-10-14)

- **Global, conceptual view:** `../claude_docs/` (i.e. `smear-modules/claude_docs/`, spelled
  `claude_docs`), shared by all modules and by the claude.ai Project — ONE file, `handoff.md`
  (merged 2026-10-14): status at a glance, reference sections §0–§8 (working rules, project
  context, Paper 2 perception concept and state), session log §9 (newest entry at the top).
  Concepts, decisions, outcomes and their implications for the papers only —
  no commit hashes, CI runs or code-level detail. Its writing rules: no fabrication, label
  claims sourced/inferred/speculative, US English, SI units.
- **Local, implementation view:** this `claud_docs/` (CLAUDE.md, memory.md, plans).
- At the end of each session update both: the details here, the concepts there.

## Branch `features/linear_elasticity` (2026-10-15) — read this first when on it

This `claud_docs/` is untracked and shared by every branch in the working copy. On
`features/linear_elasticity` the **Stokes model and all optimization code are deleted** (user
decision): no `Stokes`, `SqueezeFlow`, `def_problem`, `fit_model`, cost functions, replicate
statistics, geometry wrapper types, or their tests/experiments. Most sections below describe the
Stokes line (`main`/`features/realtime_solver`) and do NOT apply there. What does:
meshes, `fem.jl`, PostProcess, io, config, `sparse_lu.jl`, the Gmsh gotchas, units, simplex weights.

On that branch: `LinearElasticity(; mesh, λ, μ)` (kPa), `lame_parameters(E, ν)`;
`simulate(mdl, cParam; control="force"|"displacement", β, sensitivities=false)` (squeeze between
rigid plates, per-step load increments, updated Lagrangian) and `volume_update!(mdl, nodes,
normals, d; β_table=0)` (paper-2 volume update). Core: `src/solver/elastic_system.jl`
(`ElasticSystem`, fixed pattern + reused factorization). `BasisFunctionCache(mesh)` takes ONE mesh
there. Plate/table nodes come from `IEN_top`/`IEN_bottom` — the meshes' `top_nodes`/`bottom_nodes`
omit the nodes on the side. Tests: `test/linear_elasticity_test.jl` (also the Tet quadrature-weight
guard); `Pkg.test()` 225 checks, ~1 min. Benchmark: `experiments/elasticity/benchmark.jl`
(Hex Q2 paper-1 mesh, Pardiso: 7 ms/step, 15 ms with sensitivities; volume update ~20 ms. UMFPACK,
i.e. every Mac: same per step without sensitivities, but first step and volume update ~150 ms).
Neither sparse solver refuses a singular K — `_singular` (sparse_lu.jl) makes both warn. Full
history and measurements: `memory.md`, the 2026-10-15 entries.
Committed 2026-10-15 as `cc86cbf`, on top of `d19ff68` (the ArUco-script removal it shares with
`features/realtime_solver`, so it is that branch plus one commit). Pushed 2026-10-15 (tracks
origin). CI only via `gh workflow run CI.yml --ref features/linear_elasticity`; first run
38772998784 green (225/225 on Julia 1.11, ubuntu/macOS/windows; docs + doctests). `experiments/` holds only
`elasticity/benchmark.jl`: `aruco_accuracy.jl` (user: not needed on any branch; deferred #10),
`plot_style.jl` and `ParallelExecution.jl` (no users left here; user said remove) are deleted.

## Perception link (2026-10-23) — `src/io/perception_link.jl`

Two-process link to smearPerception (Python). User decisions 2026-10-23: offline and live both; live
over a Unix-domain socket (a shared-memory ring was built first and replaced the same day after the
review measured the socket faster and simpler); one record per FE step, keyed by the global step
index n (not the paper's within-window j; windows, fixed or sliding, are maps from n on this side).
- Offline: `read_perception_steps(h5)` → `(; grid, steps, t0_ns, dt_ns)` from `reconstruct_bag.py
  --save-field` output. Files written before the record change have no `n`: clear error, rerun.
- Live: `open_perception_link(path="/tmp/smear_perception.sock"; timeout=Inf)` + `next_step!(link)`;
  perception publishes with `fe_link.StepPublisher` / `reconstruct_bag.py --link [path]`.
  `PerceptionLink(io)` reads any IO carrying the same bytes (the test feeds a recorded stream).
- Record (`PerceptionStep`): `n`, `t_ns` (= t0 + n·Δt), `n_views`, `phi`, `count`, `points`, `normals`,
  `contact` (Int8: 1 = |F| > perception's threshold, 0 = below, −1 = unknown: no threshold or no
  F/T sample), `wrench` (F/T sensor frame, raw, every step), `probe_pose` (4×4 probe → world). User
  decisions: contact by force threshold (value not decided), wrench sent in and out of contact,
  pose every step. Wrench/pose are NaN where the recording has none (no OptiTrack yet).
- Units: the wire is SI, C order; this side converts once to FE units — lengths mm, forces
  kg·mm/s² (×1e3), torques kg·mm²/s² (×1e6). Fields are permuted back on read (user, 2026-10-28):
  φ and count are (nx, ny, nz), voxel (i,j,k) 0-based at `phi[i+1, j+1, k+1]`, same as Python;
  `next_step!` reads the wire bytes into `link.phi`/`link.count` and `permutedims!` into the step
  (no per-step allocation). Points arrive (3, P); `probe_pose` is transposed back.
  Message layout: docstring of smearPerception `src/smear_perception/fe_link.py` (the spec; keep the
  two files in step).
- `next_step!` returns the link's OWN reused buffer (deepcopy to keep); ~2–3 KB allocated per step.
  Reason: allocating ~5 MB per read caused full GCs of ~240 ms every ~8 steps (measured).
- Ends: `nothing` after perception's END message (complete run); `EOFError` when perception died.
  Any start order: `open_perception_link` retries until perception listens, perception's first
  step waits for the FE; perception blocks while the FE is behind, so no step is lost.
- `sample_field(phi, grid, pts_mm)`: trilinear by Interpolations.jl (`interpolate!` = no copy,
  `extrapolate(..., NaN)`), same rule as perception's SciPy `voxels.sample_field` (FE samples φ_s at
  its nodes; the field crosses every step anyway and it avoids a round trip). Interpolations names
  are imported explicitly: Distributions also exports `scale`. User guideline (2026-10-28): use
  off-the-shelf libraries wherever one fits.
- Crossings can carry NaN normals (undefined gradient; 9 per step on Carreau/5 steps 0–1).
- Test: `test/perception_link_test.jl` on a fixture written by perception's own writer, publisher and
  sampler (`test/data/perception_link/make_fixture.py` → steps.h5, stream.bin, expected.h5; grid
  5×7×9 so an axis swap fails; wrench/pose values make a missing conversion or transpose fail).
  Socket testset Unix only. `Pkg.test()` 275/275.

## Commands

```bash
# Install / update deps
julia --project -e "using Pkg; Pkg.instantiate()"

# Run an experiment / driver script (always with these flags — BLAS threading depends on --threads)
julia -O3 --threads=auto --project experiments/stokes/single_simulation.jl

# Run one test file directly (see "Pkg.test() is broken" below for why not `Pkg.test`)
julia --project -e 'include("test/stokes_test.jl")'
julia --project -e 'include("test/fem_test.jl")'
julia --project -e 'include("test/cost_functions_test.jl")'
julia --project -e 'include("test/stokes_optimization_test.jl")'
julia --project -e 'include("test/mesh_templates_test.jl")'
julia --project -e 'include("test/tikhonov_test.jl")'

# Build docs locally
julia --project=docs docs/make.jl
```

### `Pkg.test()` — FIXED 2026-09-20, use it

`julia --project -e 'using Pkg; Pkg.test()'` is now the canonical way to run the suite:
**97,987 tests, ~3.5 min** (2026-10-14, real-time solver branch; was 97,912 / ~22 min before
`test/squeeze_system_test.jl` and the faster `simulate`). Two faults had to be fixed together:

1. `test/Project.toml` listed only `smearFEM`, and when that file exists Pkg uses it *instead
   of* the root `[targets]` block. The root block was already complete, so the fix was to
   **delete `test/Project.toml` and `test/Manifest.toml`**, not to populate them. Do not
   reintroduce them.
2. `[targets]` was itself missing `Random`, which `cost_functions_test.jl` imports. Added to
   `[extras]` and `[targets]`.

Running a single file directly still works (`julia --project -e 'include("test/stokes_test.jl")'`)
for everything except `qa.jl`, which needs Aqua and therefore the test environment.

Do **not** work around it with `push!(LOAD_PATH, "@stdlib")` — it has poisoned the precompile
cache before, making a subsequent plain `julia --project=.` fail with a spurious "Package smearFEM
does not have Printf in its dependencies" error (Printf *is* a real dep). Recover with
`Pkg.precompile()`.

## Architecture

### Module structure

`src/smearFEM.jl` is a flat include list (no submodule nesting) — every file listed there shares
one namespace and one export list at the top of that file. When adding a new file, add both the
`include(...)` and the relevant `export` line there.

### The three-mesh Stokes model

A `Stokes <: AbstractModel` holds three separate meshes for a mixed (Taylor-Hood-style)
discretization, all built from the *same* geometry but with different element order:

- `mesh_u` — velocity, typically quadratic (`basis_order=2`)
- `mesh_p` — pressure, typically linear (`basis_order=1`)
- `mesh_x` — geometry / node tracking, matches the velocity mesh topology but `nDof=1`

Every mesh struct carries `volume_element_shape::Symbol` (`:Hex`/`:Tet` in 3D, `:Quad`/`:Tri` in 2D,
`:Line` in 1D) and `basis_order::Int` (1 or 2) as two independent fields — never a combined
"Q2"-style string. `BasisFunctionCache` precomputes basis-function values at all Gauss points
once per simulation and is passed into assembly — never recomputed inside a time-stepping or
optimization loop; build it once, before any control dispatch (see `simulate`).

(An earlier version of this file claimed `MeshgridCylinder` also stores precomputed quadrature
arrays `C_vol`/`C_top`/`C_btm`/`W`. **It has no such fields** — corrected 2026-09-17.)

### Geometry dispatch

The object being squeezed is a wrapper type (`Cylinder`, `Cuboid`, `Disk`, `Square`, `Segment` —
`src/fem/geometries.jl`) that dispatches `def_problem` and `set_model` (`src/solver/stokes_setup.jl`):
`def_problem` wires up physics (boundary conditions, DOF numbering) for that geometry, `set_model`
builds the three meshes and assembles them into a `Stokes` struct. Adding a new geometry means
adding a wrapper type plus a `def_problem`/`set_model` method pair, not branching inside existing
methods.

All five are `mutable struct` with a default for **every** field, and are constructed
**positionally** — `Cylinder(r, h)`, `Cuboid(lx, ly, lz, edge_radius)`. This is a deliberate,
documented exception to the guideline's "keyword-only inner constructors": Julia's strict
`function T(; r, h)` form rejects positional calls, and all 12 construction sites in the repo are
positional. Don't "fix" it to keyword-only — it breaks every caller.

`SqueezeFlow` (`src/fem/scenarios.jl`) got the same treatment, with defaults on 8 of its 9
parameters. `model` is the exception: it is typed `AbstractModel`, so no safe default value
exists. That field is also **never read anywhere** — a removal candidate that would dissolve the
exception.

### Meshing

Every geometry meshes either structurally (pure Julia, `meshgrid_*` functions in
`src/fem/meshes/`) or unstructurally via Gmsh (`mesh_type=:unstructured`, `src/io/gmsh_utils.jl`).
Gmsh meshes are generated from `.geo` templates in `mesh_files/templates/` (version-controlled)
into `mesh_files/<geometry>/<Shape>/` (git-ignored cache) — **delete the cache after editing a
`.geo` template**, it is not invalidated automatically. All Gmsh API calls must go through
`GMSH_LOCK` (the C library isn't thread-safe); never call `gmsh.*` outside that lock. Gmsh and
smearFEM use different node orderings for higher-order elements — `get_mesh_data` remaps them
(e.g. an explicit 27-entry permutation for Q2 hex); a new element type needs its own remap entry
in the `body_map` block there.

### Assembly conventions

Sparse system matrices are built in COO format (`push!` onto parallel `E`/`J`/`V` vectors during
the element loop, then a single `sparse(E, J, V)` call) rather than incremental sparse indexing.
Struct fields are pulled into local variables with `@unpack` before entering a hot loop, suffixed
`_cached` (e.g. `IEN_cached`) to mark them as snapshot locals rather than live struct fields.

### Optimization

`src/optimization/smearOptimize.jl` implements Gauss-Newton and Levenberg-Marquardt fits
(`fit_model`, dispatched by `method=:gn`/`:lm`, plus a Tikhonov-regularized `:gn_tikhonov`
variant), using the analytic `∂q/∂η`, `∂q/∂β` sensitivities the forward solve produces rather than
finite differences. Costs live in `src/optimization/cost_functions.jl` and are residual-based
(`ClosestPointCost`, `ChamferCost`, `SignedDistanceCost`) — a scalar metric can't be used directly
as an optimization cost; adding a new one means writing a `_residuals` method, not a scalar metric.

Internal helpers there are `_`-prefixed as of 2026-09-17: `_armijo_line_search`,
`_backtrack_line_search`, `_val_check`. Two exported symbols were also renamed repo-wide —
`readData` → `read_data` and `get_cMat` → `get_c_mat` — along with `plotting.jl`'s
`_animate_2d`/`_animate_3d`/`_normalize_field`/`_truncate_colormap` and `io.jl`'s
`_extract_p_from_u_nodes`/`_read_sparse_mat`. `_normalize_field` was called `normalize`, which
shadowed `LinearAlgebra.normalize`.

### Replicate statistics — three approaches (A/B/C)

`src/analysis/statistical_analysis.jl` reports cross-replicate agreement for a
normalized-per-run quantity (e.g. `h_est/h_m`) three ways, all returned by
`normalized_replicate_stats` and all plotted. They differ only in which axis (runs or time) a
replicate's `R × T` matrix is collapsed along first:

- **A** — `replicate_stats(per_run_mean)` (the `.mean` field): collapse each replicate to a
  scalar (its own time-average) *first*, then analyze across replicates. Discards within-run
  temporal spread entirely. `bias_vs_reference` deliberately still uses this — that test needs
  one independent sample per replicate, not pooled/correlated points.
- **B** — `pointwise_pooled_stats` (the `.pointwise_pooled` field): collapse across replicates
  *first*, separately at each time step (this is the plotted ribbon's `pointwise_mean`/
  `pointwise_sd`), then pool those `T` pointwise SDs over time (`sd = √(mean_t(sd_t²))`).
  Anything identical across every replicate at a given instant — e.g. a time-varying artifact
  shared by the fitting method itself, not genuine specimen noise — cancels out of `sd_t`
  exactly, before time is ever touched.
- **C** — `pooled_stats(series)` (the `.pooled` field): collapse each replicate over time
  *first* (its own mean **and** its own temporal SD), then pool those `R` per-replicate
  statistics across replicates via the standard pooled-variance formula. `sd_i` (each
  replicate's own temporal SD) cannot distinguish a pattern shared by every replicate from
  noise unique to that one — a shared pattern inflates `C`'s SD along with genuine noise, which
  `B` does not.

**B and C only agree when there is no time-pattern shared across replicates.** On this
project's data they diverge sharply for `..._pred` quantities specifically (2–4.5× on physical
`h_pred/h_m`, `rel. error pred [%]`, contour cost pred) because of the window-restart sawtooth,
which is common to (nearly) every replicate — B correctly treats that as *not* between-replicate
disagreement; C misattributes it to noise. `..._est` quantities, with little shared time
structure, keep A/B/C close together. `approach_stats(st, :A|:B|:C)` picks the right field;
never branch on the field names directly. Full derivation (generative model, why C is biased
and B isn't, the "two-stage formula that's secretly flattening" pitfall) in
`~/approaches_abc.md` (not checked into the repo).

Plotting (`experiments/model_optimization/test_opt_stokes.jl`): `_stat_fig(...; approach=:A|:B|:C)`
draws, per group (Estimation/Prediction), the raw per-run curves, an optional pointwise ribbon
(`_band!`, `ribbon=true` by default, approach-independent), and a flat reference band from
`approach_stats(st, approach)` (`_scalar_band!`, `bands=:sd`/`:ci`/`:both`) — a direct visual
check that the reported scalar matches what the curves underneath actually show. `ribbon=false`
drops the ribbon so a bad scalar (e.g. a silent `NaN`) isn't hidden behind an unrelated
pointwise band. `_stat_fig_approaches(groups, tgrid, dir, basename; kwargs...)` is the call-site
entry point every `post_analysis_*` function uses: for each of A/B/C it writes a `<dir>/
approach_<A|B|C>/<basename>_ci.pdf` (ribbon + nested SD/CI bands) and `..._scalar_only.pdf`
(runs + SD band only) — one self-contained folder per approach. `_write_approach_reports(stats,
dir, filename)` is the CSV counterpart, writing `<dir>/approach_<A|B|C>/<filename>` (mean/sd/
ci_lo/ci_hi/... — one row per quantity, via `write_approach_stats`) alongside the plots. The
older single combined CSV (`replicate_report`/`write_replicate_stats`, all of A/B/C/rmse/maxdev
as rows of one file in `plot_path_global` directly) is kept too, unchanged, for backward
compatibility.

**`plot_results`'s driver sweep does not cover every real result tree.** It only pairs
`data_type == "synthetic"` with `viscosity_type == "bulk_viscosity"` — `("synthetic",
"constant")` is never in its loop, even though that exact key is a real, configured entry in
`PLOT_PRESETS` and a genuine 6-replicate result tree exists at
`experiments/syn_data/optimization/Stokes/force/constant/Hex2_16/cylinder`. Don't infer "which
trees exist" from `plot_results`'s hardcoded config — check the filesystem. That tree also
predates multi-view support: its contour data lives at the flat `sim_data/2D_border_points/`
(60 points/row, one `x,y` pair per row) rather than the current `sim_data/view_1/
2D_border_points_est/` (50 points, two rows: all-x then all-y) — a real algorithm change, not
just a path/format rename, so don't attempt to convert one to the other. `post_analysis_const`
now tolerates a missing contour directory (`have_borders` flag): height/rel-error stats come
only from `est_h.csv`/`gt_h.csv` and never touch border-point data, so a tree in this old format
still contributes full height stats — only the `contour_plt*` diagnostic plots are skipped.

**Constant-viscosity `est_h` is mostly prediction, not estimation.** `optimize`'s
`gt_viscosity_type == "constant"` branch fits only against the first `sim_time_exp` seconds of
observations, then re-simulates with the converged, fixed `(η, β)` out to the full
`sim_time_gt` — e.g. a 5s fit produces a 301-row `est_h.csv` (≈30s), of which only the first 51
rows are in-sample. `post_analysis_const` now splits this properly at `n_est = round(Int,
sim_time_exp/t_steps) + 1` into `norm_est_runs`/`norm_pred_runs` (mirroring
`post_analysis_bulk`'s Estimation/Prediction groups, `vlines=[sim_time_exp]` marking the
boundary) — it used to push the whole unsplit series into `norm_est_runs`, so any
`h_est/h_gt` number computed before 2026-09-16 for a constant-viscosity tree is a mix of
in-sample agreement and out-of-sample drift, not a clean estimation-quality number. `replot`'s
own single-experiment `h_plt`/`error_plt` still don't make this split (still one "Estimation"
label) — only the cross-replicate stat-plot path was fixed.

### Config / path resolution

**`get_mesh_dir()` / `SMEAR_MESH_DIR` is resolved but never consumed by mesh code** — every Gmsh
generator and `def_problem` default builds its cache path from `@__DIR__` → repo `mesh_files/`.
Decided 2026-10-12 to *document* this rather than wire it in (wiring it would orphan existing
caches). README and coding_guide.md say so; don't "fix" one side without the other.

Never hardcode absolute paths. Use `get_data_dir()`, `get_mesh_dir()`, `get_scratch_dir()`,
`resolve_data_path(...)`, `resolve_mesh_path(...)` (`src/config.jl`), which resolve through env
vars (`SMEAR_DATA_DIR`/`SMEAR_MESH_DIR`/`SMEAR_SCRATCH_DIR`) → `config.toml` (lives one level up in
`smear-modules/`, shared with the sibling `smearPerception` package) → hardcoded defaults.

### Experiment scripts

`experiments/` scripts are run directly (each has a `main()` gated on
`abspath(PROGRAM_FILE) == @__FILE__`-style script detection so they're also `include`-able from
the REPL) and may depend on packages beyond `Project.toml` (e.g. `CurveFit`, `Dierckx`). Batch
scripts use `experiments/ParallelExecution.jl`'s memory-aware scheduler
(`run_parallel_tasks(param_list, task; max_workers=-1, memory_per_task_mb=512.0)`) rather than
plain `@threads`/`pmap`.

## Known gotchas

- **Two mesh generators, chosen by the machine.** `_run_gmsh` (`src/io/gmsh_utils.jl`) uses the
  `gmsh` CLI when it is on PATH — the user's machine has `/usr/bin/gmsh` 4.12.1 — and otherwise
  the Julia API (gmsh_jll, 4.15.x on CI) in ONE persistent session. `gmsh.clear()` keeps global
  options, so until `7c93012` (2026-10-14) `Mesh.RecombineAll` from the hex templates leaked
  into later tet/tri/line meshes: 13 mesh tests errored on every CI OS (likely why CI was red
  since May); locally invisible because the CLI path starts a fresh process per mesh. Test the
  API path with `env PATH=/nonexistent <julia> ...`. A fresh mesh can differ between machines
  (Gmsh 4.12 vs 4.15). After the fix CI is green on Ubuntu, Windows and macOS (2026-10-14).
- **Force units are kg·mm/s², not newtons**, because the solver works in mm (1 N = 1e3 in these
  units). Several comments in the codebase (`verify_meshes.jl`, `single_simulation.jl`,
  `forward_convergence.jl`, the `optimize_real` docstring) say "N" incorrectly — the numeric
  values are correct, trust those over the unit label in comments.
- **`nNodes^ndim` is always stale, never a convention.** `nNodes` meant nodes-*per-edge* on the
  old structured grids, so total DOFs were `nDof * nNodes^ndim`. Since unstructured meshing it is
  the *total* node count, and the correct sizing is `nDof * nNodes` — as in `simulate`
  (`stokes_solver.jl`). Seven surviving `^ndim` sites in `stokes_single_step_force` allocated a
  1885³ ≈ 6.7e9 matrix and died with `OutOfMemoryError`; fixed 2026-09-17. If you meet the
  exponent anywhere else, it is a bug.
- **Simplex quadrature weights must match the master element `∇λ` defines.** `_basis_tri`/
  `_basis_tet` take `∂λ/∂ξ` from Eq. (3.39) of the 4774 FEM notes, whose master element is the
  **equilateral/regular unit-edge** simplex — so weights must sum to its measure
  (`TRI_MASTER_AREA = √3/4`, `TET_MASTER_VOLUME = 1/(6√2)`), **not** the right-simplex 1/2 and
  1/6 of most published tables. They did not until 2026-09-19, over-integrating every Tet
  element by 6√2 and every Tri by 2/√3. `∇λ` is correct — do not "fix" it to unit-simplex
  gradients. Tensor-product elements (Hex/Quad/Line) have no such ambiguity and were never wrong.
- **RAM, not core count, limits parallelism.** Wide parameter sweeps should default to serial
  execution unless you've confirmed available memory against `run_parallel_tasks`'s
  `memory_per_task_mb` — a full `:gn` fit is CPU- and memory-heavy enough that a handful of
  parallel workers can exhaust a typical dev machine.
- **The skip-list lives in ONE place: `_window_subdirs(exp_path; keep_single_window=false)`**
  (`test_opt_stokes.jl`). It used to be copy-pasted into all three call sites, where editing two
  of three silently reintroduced the bug below. **`replot` passes `keep_single_window=true`** —
  it plots that condition and has the machinery for it (`pred_is_windowed`, and its
  `window_dir == "single_window"` branches); the two `post_analysis_*` functions must never see
  it. `replot` also guards each window with `_skip_incomplete` so an empty one is skipped, not
  thrown on.
- **FIXED (2026-09-16), keep the guard**: `window_dirs = readdir(leaf.exp_path)` in
  `post_analysis_bulk`/`post_analysis_real`/`replot`'s windowed branches must skip
  `single_window` alongside `Results`/`post_analysis_window`/`post_analysis_noise`.
  `single_window` is a *different experimental condition* per replicate (one global static fit
  evaluated against the same window grid, for comparison) — not another window of the
  `multi_window` re-fit — but it shares `data_ranges.csv`/`t_windows.csv` with `multi_window`,
  so nothing before this filter distinguished them. Left unfiltered it doubles `n` (replicates ×
  2 conditions) and its differently-shaped `pred_h.csv` breaks `_align_windowed`, producing
  `h_pred/h_gt`/`MOSD_pred/MOSD_gt` `NaN` — this was previously misdiagnosed as
  `_align_windowed`'s own (correct, intentional) `NaN`-gap behavior. Physical trees never
  exercised the bug only because they never had a `single_window` folder. If a future result
  tree adds another sibling directory under `view_1/` for a different comparison condition, it
  needs the same treatment.

## Real-time solver (2026-10-13/14), branch `features/realtime_solver` — committed

Commits (pushed 2026-10-14): `07594e9` fix (contour sensitivities of rotated views), `409421a`
perf (the solver), `1e2938e` refactor (projection out of `simulate`), `7c93012` fix (Gmsh
options in the Julia-API fallback). Then `ed36680` chore and `d19ff68` (ArUco-script removal,
2026-10-15), both pushed. CI runs on this branch only via `gh workflow run CI.yml
--ref features/realtime_solver` (push triggers only main/dev/publish_ready and PRs). CI is green on
`7c93012`. **Not merged into `main` — user, 2026-10-14: leave the FE as it is for now.** The
branch sits on `features/tet_mesh` (147 unmerged commits) and merges cleanly into `main`; merge the
two together, since `tet_mesh` alone still has the Gmsh option leak.

Target: **everything ≥ 50 Hz** — one dt = 0.1 s step of `simulate` in ≤ 20 ms, force AND
velocity control, sensitivities included, no loss of accuracy. Status (2026-10-14): median
9–12 ms/step, worst steady step ~18–22 ms, step 1 ~22 ms (per call), from ~586 ms; outputs
within ~2.5e-9 of the pre-branch code. Full suite 97,987/97,987 in ~3.5 min. **Living plan,
every measurement, gotchas: `claud_docs/realtime_solver_plan.md` (Phases A–F).** Scripts and
pre-branch reference outputs: `claud_docs/realtime/`. User decisions 2026-10-14: keep MKL
Pardiso with the UMFPACK fallback; the camera projection moved out of `simulate` (Phase F).

**API since Phase F:** `simulate` is 3D-only and returns 7 outputs — `output, displacement,
surface_pts_3D, pos3D, velocity, pressure, gradList_3d` (was 10). The camera lives in
`project_contours(mdl, pos3D, conditions; grad3D=nothing)` (stokes_solver.jl, after
`simulate`), which returns `borderPts2DList, gradList, pos2D` (gradList empty without grad3D),
threads over time steps, and owns `WRITECONTOUR` and the 2D GIFs of `ANIMATE` (simulate keeps
the 3D GIFs). Every fit goes through `_simulate_contours` (smearOptimize.jl). Results are
bit-identical to the old in-loop projection, except that the first frame's contour
sensitivities are now 2×nb×2 zeros (were 3×nb×2; every consumer reads rows 1:2).

How `simulate` solves now:
- `src/solver/squeeze_system.jl` (`SqueezeSystem`, cached per model in a WeakKeyDict): the
  force-control saddle matrix, unknowns `[q_f; p; μ]` (μ = top-plate velocity), on a **fixed
  sparsity pattern** refilled each step by a threaded gather of element matrices (kernel sums
  each node pair's 3×3 block over Gauss points in registers). Velocity control = the same
  matrix with μ's row/column set to the identity (`_prescribe_plate_velocity!`).
- `src/solver/sparse_lu.jl`: the factorization backend — **MKL Pardiso** (`Pardiso` dependency;
  ~10 ms factorization, ~1–2 ms solve) when `Pardiso.mkl_is_available()`, else UMFPACK
  (~150 ms / 3.6 ms; two solves at once on a `copy(F)` workspace). Pardiso gets K's lower
  triangle with every diagonal, Int32 indices, `iparm(34)` = thread count (bit-reproducible).
- Solves the unit force/velocity response `z` and its β-sensitivity `w` together
  (`_lagged_solve!`): cubic extrapolation in time (history kept at unit viscosity), iterative
  refinement against the reused factorization to STATE_TOL 1e-9 / SENS_TOL 1e-8. ∂/∂η comes
  from ∂/∂β by homogeneity (η∂x/∂η + β∂x/∂β = −[u;0] force, [0;p;0] velocity). Sensitivities
  remain the frozen-geometry approximation the old code used (exact at step 1 only).
- Refactorization: step 1 of every call (results never depend on earlier calls); Pardiso again
  at the start of the step after one shows staleness (`_stale`, ρ > 0.05); UMFPACK only after
  > 5 rounds. `simulate(...; sensitivities=false)` for prediction only.

Invariants / gotchas:
- **K must stay symmetric**: residuals are rhs − Kᵀx (column-parallel).
- **2D contour sensitivities were wrong for every view but 0° before 2026-10-14** (positions were
  rotated by the view angle, `dxdθ` was not) — fixed in `back_project`; the pre-branch
  reference `gradList` is wrong for views 30°/60°.
- **Velocity control was broken at HEAD** (`zero_matrix_p_q`); fixed by the unification;
  `output` is v·dt per step under velocity control. `μu_btm`/`μu_side` ≠ 0 throw.
- Pardiso.jl uses Int32 indices unless Julia's BLAS is MKL — anything else is converted (2 MB)
  per call. MKL symmetric matrices need every diagonal entry stored.
- MKL_jll 2025.2 ships for **x86-64 Linux (glibc) and Windows only** — every Mac and every ARM
  machine runs the UMFPACK fallback. Artifacts on disk: MKL 539 MB + IntelOpenMP 187 MB +
  oneTBB 2.5 MB (~730 MB; an earlier note said ~300). License: Intel Simplified Software
  License (free, redistributable, binary-only — not open source; the repo is MIT).
  `using Pardiso` adds ~0.55 s to load time.
- GC: a collection with the packages loaded costs ~60–80 ms (young) / ~320 ms (full); the loop
  allocates ~0.5 MB/step, mostly the returned history, so pauses are rare but possible.
- The projection never fed back into the solve: taking it out of the loop left every 3D
  output bit-identical and made a step ~0.5 ms cheaper (10.96 → 10.46 ms, 0.55 → 0.31 MB/step,
  quiet A/B 2026-10-14). Projected afterwards, threaded over steps, it costs ~0.06 ms/step.
- `contour_cost` reads only the FIRST view of each frame (`_as_2xN`/`_as_dudθ` take x[1]); the
  paper pipeline fits one view per call (`viewing_angles=[z_angle]`, ground truth at
  0/5/15/30/45°) — so the old `back_project` bug hit every per-view fit at 5–45°.
- `BLAS.set_num_threads(Threads.nthreads())` at the top of `src/smearFEM.jl` only runs at
  precompile; README's "set at module load" is wrong. Not changed (global policy).
- Benchmark on a quiet machine only (the user's browsers moved results by 30–100 %).

## Style (from `CONTRIBUTING.md` / `docs/coding_guide.md`)

- 4-space indentation; docstrings (standard Julia `# Arguments`/`# Returns` format) required on
  every exported function, not required on `_`-prefixed internal helpers.
- `@debug`/`@info`/`@warn`/`@error` only — no `println` in library code (`ENV["JULIA_DEBUG"] =
  "smearFEM"` to see debug output). The module installs a stderr logger at load time so log
  output doesn't collide with stdout data streams.
- Comment only the *why* (non-obvious invariants, library workarounds, constraints not visible
  from types) — never restate what the next line does or reference the PR/issue that introduced it.
- Types: `PascalCase`, mutable, keyword-only inner constructors with safe defaults for every
  field. Functions: `snake_case`. Constants: `UPPER_SNAKE_CASE`. Prefer keyword args with defaults
  over positional args.

Three conventions settled during the 2026-09-16 publish-readiness pass, because the source docs
were silent or self-contradictory on them:

- **`println` exception.** The "no `println`" rule targets *debug* output. A function whose
  entire purpose is formatted console output keeps its `println`s — `show_config`
  (`src/config.jl`) and `print_replicate_table` (`src/analysis/statistical_analysis.jl`) are the
  two in `src/`. `@info` prefixes each line with a log level and module, which destroys box and
  table layout. **A third category**: `write_time_log` (`src/utils.jl`) uses `println(io, ...)`
  against an open file handle — that is file formatting, not console logging, and a naive
  repo-wide grep-and-replace would corrupt `time_log.txt`. Every other stray `println` in `src/`
  was converted to `@debug` on 2026-09-17.
- **File size.** >600 lines invites a split review, >1000 lines is a strong flag. Neither
  `CONTRIBUTING.md` nor `docs/coding_guide.md` codifies a threshold; the only precedent was a
  1521-line `Meshes.jl` judged "too large" and split.
- **Docstring template.** Follow `docs/coding_guide.md`'s form (`# Arguments` / `# Returns` /
  optional `# Example`). `CONTRIBUTING.md` shows a different one with a separate
  `# Keyword Arguments` section — the two disagree, and coding_guide.md wins.

## Publish-readiness refactor — Phases 0-3 + 5a-5d complete (Phase 4 skipped)

A staged refactor of `src/` + `experiments/` for style conformance *and* structural cleanup.
**The plan is at `/home/soshala/.claude/plans/peppy-hopping-elephant.md`**; `claud_docs/memory.md`
carries the per-phase detail and the reasoning. Status summary only here.

**Done:** Phase 0 (commit hygiene), 1 (dead code, 2,605 deletions), 2a (`println` → `@debug`),
2b (renames + `_` prefixes), 2c (32 docstrings), 2d (22 misapplied comments + `smearFem` typo),
3a (`geometries.jl`), 3b (`SqueezeFlow` defaults), 3c (`stokes_single_step_force` repair),
5a-5d (`experiments/`: gates, docstrings, `println` conversions, two typo renames).
Plus two out-of-band changes: the `SqueezeFlow` boundary-key fix, and contour-cost replicate
statistics for synthetic bulk-viscosity runs (scenario B).

**Every exported symbol in the package now has a docstring** — 39 undocumented at the start of
2c, 0 now.

**Phase 4 (file splits) was consciously SKIPPED**, not forgotten. `stokes_solver.jl` is 1045
lines and `PostProcess.jl` 1103, both over the >1000 strong flag — but `stokes_solver.jl`'s
overage is *entirely the docstrings added in 2c*, and the plan itself marks the `PostProcess.jl`
split "optional, do last, droppable". Splitting a numerically-critical solver to satisfy an
advisory line count was judged poor value.

**Phase 5 resumed 2026-09-18: 5a-5d done, 5e-5h not started.** Remaining: 5e (extract the
duplicated `def_*` colours and `PLOT_CONFIG` to `experiments/plot_style.jl` — the one Phase 5
item that touches `test_opt_stokes.jl`), 5f (`write_gt_stokes.jl`), 5g (move `verify_meshes.jl`
to `test/`), 5h (`test/` consistency). **Phases 6-9 were never started.**

**`experiments/` has no test suite.** The plan's verification — "re-run each script and compare
to a pre-change run" — is infeasible for most files: the convergence sweeps are hours, the
batch scripts are multi-hour at 8 GB/task. Verify by the cheapest thing that exercises the real code path: re-run plot
functions against existing datasets, smoke-test library functions directly, and check the
script gate holds on `include`.

**Never `include` an `experiments/` file without checking it has a script gate first.**
`write_gt_stokes.jl` and `verify_meshes.jl` still lack one, and including `write_gt_stokes.jl`
immediately starts a full parallel ground-truth batch.

### Verification protocol used throughout — keep it

Every phase: package smoke-load, then all test files, compared against this baseline
(updated 2026-09-19 by Phase 5h) — **stokes 3770, stokes_optimization 93868, mesh_templates 144,
tikhonov 41, cost_functions 34, fem 9+4, verify_meshes 33**, plus (2026-10-13)
**squeeze_system 75** (`test/squeeze_system_test.jl`, ~30 s; 2026-10-14). All exit 0.

Two entries changed in 5h and the old values are *not* regressions to chase:
- `stokes_optimization` previously reported "clean, no testset" — its 93,868 assertions were
  running uncounted because the file had no enclosing `@testset`. That is 95% of the suite,
  which was effectively unmonitored.
- `tikhonov` previously exited 1 after all 41 passed (the line-101 `round`-on-`Vector` bug). It
  is fixed and now exits 0. **The suite has no expected-failure file any more** — a red file in
  the baseline masks genuine new failures, so keep it that way.

`verify_meshes_test.jl` is **slow (~8.5 min; 3m40s since the 2026-10-13 solver rewrite)**: six squeeze-flow solves across a β sweep. It is
the only guard against the simplex quadrature weights drifting out of step with the master
element `_basis_tet`/`_basis_tri` assume; see the 2026-09-19 entries in `memory.md`. Do not drop
it, or trim its β cases, to save time — β selects which half of the quadrature is tested. Anything touching
the solver additionally got a before/after numerical diff, not just a passing suite.

### Deferred — found during the refactor, deliberately NOT acted on

**Nine still open** (was thirteen). Closed 2026-09-20: #1 (dead `SqueezeFlow.model` field —
removing it also dissolved the "defaults for every field" exception), #7 (5 of 7 dead functions;
see the meshing caveat below), #8 (181 restating comments), #9 (`μu` knobs are now keyword
arguments defaulting to 0, not dead constants), #10 (`aruco_accuracy.jl` runs again; since
deleted, see the item), #11
(`plot_results` dispatch re-enabled). Full context in `claud_docs/memory.md`:

1. **`SqueezeFlow.model` is dead** — never read anywhere. Removing it also dissolves the
   "defaults for every field" exception, since it is the only field blocking it.
2. **`q_d` round-trips through dense.** Declared `Dict{Symbol, Matrix{Float64}}`, so the
   constructor densifies each `SparseMatrixCSC`; the solver then re-annotates them as
   `SparseMatrixCSC{Bool,Int64}` (`stokes_solver.jl:712-714`). Two pointless conversions per solve.
3. `get_c_mat`'s docstring signature is stale — shows `(type; λ, μ, Young, ν)`, actual is
   `(c1, c2; type="lame")`.
4. `smearFEM.jl`'s export comment places `plot_covariance!` in `analysis_plots.jl`; it is in
   `plotting.jl`.
5. `simulate`'s `output` is empty under `"velocity"` control — the `push!` is commented out.
6. `simulate_single_tstep_stokes` builds its mesh path from `dirname(dirname(@__DIR__))` instead
   of `get_mesh_dir()`, bypassing the `SMEAR_MESH_DIR` → `config.toml` resolver.
7. Three helpers with **zero call sites**, not exported: `_normalize_field`,
   `_truncate_colormap` (`plotting.jl`), `_read_sparse_mat` (`io.jl`).
8. ~30 comments that restate their own line but are *accurate* (`# solve the system of equations`
   on real `K\r` solves, and the `# assemble the system of equations` / `# assemble the solution`
   / `# apply the neumann boundary conditions` families). The style guide condemns them; removing
   them is a bigger, more opinionated sweep than 2d authorised.
9. The vestigial `μu_btm`/`μu_side` knobs — correctly labelled since the boundary-key fix, but
   still hardcoded `0` at every site.
10. `aruco_accuracy.jl` could not run (fixed 2026-09-20). **User, 2026-10-15: the script is not
    needed on any branch.** Deleted from all three branches that had it (pushed 2026-10-15):
    `features/tet_mesh` `fdf2989`, `features/realtime_solver` `d19ff68`, and
    `features/linear_elasticity` (sits on `d19ff68`). No other branch tip, local or remote, has it.
    The tet_mesh and realtime_solver commits delete the same file separately; they merge cleanly.
11. **`plot_results` calls no `post_analysis_*` function** — the dispatch at the end of its
    `opt_method` loop is commented out, so the driver regenerates only `replot`'s per-experiment
    figures and no cross-replicate statistics at all. Re-enabling changes what every run writes.
12. `const_vel` (`single_simulation.jl`) is a stale name — both methods set `control = "force"`,
    neither is velocity-controlled.
13. Six "Incomplete run, skipping" warnings during `post_analysis_bulk`. They cost no replicates
    (`n` is 6) but were never traced.

### `src/fem/meshes/` is off-limits for dead-code removal

Unused functions there are **held for 1D/2D/3D dimensional support**, not deletable.
`_meshgrid_ring` was deleted in the 2026-09-20 dead-code pass and had to be restored on the
user's instruction; `_second_derivative_pairs` (`fem.jl`, branches on `ndim == 1,2,3`) was
restored under the same rule. A zero-call-site scan tells you what is *unused*, not what is
*deletable* — ask before removing anything under `meshes/` or any `ndim`-dispatching helper.

### One pre-existing bug, still unfixed on purpose

(`test/tikhonov_test.jl`'s line-101 `round`-on-`Vector` crash was **fixed in Phase 5h**, together
with registering the file in `runtests.jl` — it had never been included, which is exactly why the
crash went unnoticed.)

- **`_get_F_ext` in `experiments/stokes/write_gt_stokes.jl`** tests `β_gt <= 5e3` before
  `β_gt <= 2e3`, making the 2e3 branch unreachable, so β in (1e3, 2e3] gets ×1.01 instead of
  ×0.995. Reordering changes the force that generated existing ground-truth data — check whether
  any synthetic tree used β in that range first.
