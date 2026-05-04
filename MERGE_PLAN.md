# Merge Plan: `features/gmsh_intergration` → `main`

**Date**: 2026-05-03  
**Branch**: `features/gmsh_intergration`  
**Goal**: Land all publication-readiness work on `main` with a clean, CI-passing state.

---

## Pre-conditions

- All 90,925 tests currently pass on this branch.
- ~40 files are modified but uncommitted — all changes are from Phases 1–4 of the publication readiness plan.
- Three IGA test files have been deleted from disk.

---

## Step 1 — Remove debug `println` calls from library source

**Files to edit:**

### `src/fem/PostProcess.jl` — inside `get_lagrange_proj`

Replace:
```julia
println("n_cp: ", n_cp)
println("n_lagrange: ", n_lagrange)
```
With:
```julia
@debug "get_lagrange_proj: n_cp=$n_cp, n_lagrange=$n_lagrange"
```

### `src/io/gmsh_utils.jl` — inside `get_mesh_data`

Replace:
```julia
println("NodeList size: ", size(NodeList))
println("IEN_volume size: ", size(IEN_volume))
println("Top nodes: ", length(nodes_top))
println("Bottom nodes: ", length(nodes_bottom))
println("Lateral nodes (including boundaries): ", length(nodes_lateral))
```
With:
```julia
@debug "get_mesh_data: NodeList size=$(size(NodeList))"
@debug "get_mesh_data: IEN_volume size=$(size(IEN_volume))"
@debug "get_mesh_data: Top nodes=$(length(nodes_top)), Bottom nodes=$(length(nodes_bottom)), Lateral nodes=$(length(nodes_lateral))"
```

**Note:** No `using Logging` import needed in either file — both are included inside the `smearFEM` module, which already imports `Logging` via `src/utils.jl`. The `@debug` macro is available in module scope.

**Scope note — out of scope for this step:**
- `src/io/io.jl:90` — `println(size(IEN))` inside `write_vtk` (debug artifact, not blocking)
- `src/config.jl:239-245` — `println(io, ...)` in `show_config()` (intentional user-facing output, leave as-is)

---

## Step 2 — Fix `Project.toml`

### 2a. Fix authors field (missing closing `>`)

Current line 4:
```toml
authors = ["Soshala Weerathunge <soshalaweerathunge@gmail.com"]
```
Change to:
```toml
authors = ["Soshala Weerathunge <soshalaweerathunge@gmail.com>"]
```

### 2b. Add missing `[compat]` entries

The following packages are in `[deps]` but absent from `[compat]`. Add these entries in alphabetical order within the existing `[compat]` block:

```toml
ArgCheck = "2"
ConvexHulls2d = "0.1"
DataInterpolations = "8"
DelimitedFiles = "1"
Distributions = "0.25"
ElasticArrays = "1"
HDF5 = "0.17"
LinearAlgebra = "1"
Plots = "1"
ProgressMeter = "1"
SparseArrays = "1"
TOML = "1"
WriteVTK = "1"
```

**Notes on pre-1.0 packages:** `ConvexHulls2d`, `Distributions`, and `HDF5` are pre-1.0. Their compat entries pin to the installed minor version (e.g. `"0.1"`, `"0.25"`, `"0.17"`) because pre-1.0 semver treats minor bumps as potentially breaking.

**Verify TOML is valid after editing:**
```
julia --project -e 'import TOML; TOML.parsefile("Project.toml"); println("TOML valid")'
```

---

## Step 3 — Fix CI to trigger on `main`

**File:** `.github/workflows/CI.yml`

Current:
```yaml
on:
  push:
    branches:
      - dev
    tags: ['*']
```

Change to:
```yaml
on:
  push:
    branches:
      - main
      - dev
    tags: ['*']
```

This is the only required change. The `docs` job has no branch filter and already runs on all pushes.

**Optional — `docs/make.jl`:** The current `devbranch="dev"` means only `dev` branch pushes deploy as the "dev docs". If you want `main` to serve as the dev docs source after merge, change `devbranch="dev"` to `devbranch="main"` in `docs/make.jl`. This is not required for the tests to pass.

---

## Step 4 — Run tests, commit, push, open PR

### 4a. Verify tests still pass after Steps 1–3

```
julia --project -e 'import Pkg; Pkg.test("smearFEM")'
```
Expected: **90,925 tests, 0 failures**.

### 4b. Stage all files

Stage modified source, config, docs, scripts, tests, and new community files. Do **not** stage `Manifest.toml`, `config.toml`, or any data/output files.

**Modified source files:**
```
src/fem/PostProcess.jl
src/fem/models.jl
src/io/gmsh_utils.jl
src/io/plotting.jl
src/optimization/smearOptimize.jl
src/examples/fluid_flow_stokes.jl
src/examples/run_example.jl
src/examples/squeeze_linear_elasticity.jl
src/examples/squeeze_stokes.jl
src/smearFEM.jl
src/utils.jl
```

**Modified config and CI:**
```
Project.toml
.gitignore
README.md
.github/workflows/CI.yml
```

**Modified scripts:**
```
scripts/cost_tests/cost_robustness_test.jl
scripts/cost_tests/cost_test.jl
scripts/cost_tests/cost_test_single_run.jl
scripts/cost_tests/cost_test_single_run_stokes.jl
scripts/cost_tests/model_accuracy_test.jl
scripts/elasticity/elasticity.jl
scripts/elasticity/single_simulation_lin_elasticity.jl
scripts/model_optimization/slice_and_plot.jl
scripts/model_optimization/test_opt_lin_elas_with_gt.jl
scripts/model_optimization/test_opt_stokes.jl
scripts/stokes/single_simulation.jl
scripts/stokes/stokes_model.jl
scripts/visualization/visualization_example.jl
```

**Modified test files:**
```
test/benchmarking/detailed_profiler.jl
test/benchmarking/iterative_solver_profiler.jl
test/benchmarking/profile_100hz.jl
test/benchmarking/profileview_profiler.jl
test/benchmarking/run_all_benchmarks.jl
test/convergance_test.jl
test/convergence_analysis/stokes_convergence.jl
test/intergration_tests/mesh_comparison_test.jl
test/intergration_tests/optimization_tuning_test.jl
test/linear_elasticity_optimization_test.jl
test/mesh_comparison/compare_init.jl
test/runtests.jl
test/stokes_optimization_test.jl
test/stokes_pressure_test.jl
test/vtk_test.jl
```

**Deleted files** (stage with `git rm` or `git add`):
```
test/IGA_tests/extraction_test.jl
test/IGA_tests/linear_elasticity_test.jl
test/IGA_tests/volume_test.jl
```

**New untracked files:**
```
CHANGELOG.md
CODE_OF_CONDUCT.md
CONTRIBUTING.md
memory/MEMORY.md
memory/project_simulate_return_todo.md
```

### 4c. Commit

```bash
git commit -m "$(cat <<'EOF'
Publication readiness: Phases 1-4 — bug fixes, path resolver, docstrings, CI

Phase 1 (critical bug fixes):
- Fix test file reference, function name typo (noramlize → normalize)
- Fix backwards isnothing guards (|| → &&) in 4 locations
- Add throw() to AssertionError patterns in squeeze_stokes/run_example/stokes_model
- Replace 9 empty catch blocks with catch e + @debug

Phase 2 (code quality):
- Add config.jl three-tier path resolver (SMEAR_DATA_DIR → config.toml → default)
- Replace all hardcoded /home/soshala/... paths with resolve_data_path/get_scratch_dir
- Replace 23 debug println calls in smearOptimize.jl with @debug

Phase 3 (documentation):
- Add docstrings to all public functions in PostProcess.jl, gmsh_utils.jl,
  squeeze_linear_elasticity.jl, smearOptimize.jl, squeeze_stokes.jl
- Add CONTRIBUTING.md, CHANGELOG.md, CODE_OF_CONDUCT.md
- Rewrite README.md with correct API examples and configuration docs
- Document Gmsh thread-safety in get_mesh_data and module-level GMSH_LOCK comment
- Replace remaining println debug calls with @debug in get_lagrange_proj and
  get_mesh_data

Phase 4 (testing & API):
- Fix LinearElasticity constructor new() field ordering (C, C_top, C_btm, W)
- Fix gaussian_quadrature positional args (19 occurrences)
- Add set_model/def_problem/simulate(::LinearElasticity) consistent with Stokes API
- Fix reset_model!/update_model! for LinearElasticity (model.NodeList_init)
- Remove three obsolete IGA test files

Project.toml:
- Fix authors field (missing closing >)
- Add missing [compat] entries for 13 packages

CI:
- Trigger test job on both main and dev branch pushes

Tests: 90,925 passing, 0 failures

Co-Authored-By: Claude Sonnet 4.6 <noreply@anthropic.com>
EOF
)"
```

### 4d. Push to origin

```bash
git push origin features/gmsh_intergration
```

### 4e. Open the PR

```bash
gh pr create \
  --base main \
  --head features/gmsh_intergration \
  --title "Publication readiness: Gmsh integration, path resolver, full docstrings, CI fixes" \
  --body "$(cat <<'EOF'
## Summary

This PR lands all publication-readiness work from the \`features/gmsh_intergration\` branch.

**Critical bug fixes (Phase 1):**
- Fixed backwards \`isnothing\` guards and missing \`throw()\` calls (silent error swallowing)
- Renamed \`noramlize\` → \`normalize\`, fixed test file reference
- Replaced 9 empty \`catch\` blocks with \`catch e\` + \`@debug\` logging

**Code quality (Phase 2):**
- Three-tier path resolver (\`SMEAR_DATA_DIR\` env → \`config.toml\` → \`~/SMEAR-Data\`)
- All hardcoded absolute paths replaced with \`resolve_data_path\`/\`get_scratch_dir\`
- 23 debug \`println\` → \`@debug\` in \`smearOptimize.jl\`; remaining \`println\`s in \`get_lagrange_proj\` and \`get_mesh_data\` likewise converted

**Documentation (Phase 3):**
- Docstrings on all public functions in PostProcess.jl, gmsh_utils.jl, squeeze_linear_elasticity.jl, smearOptimize.jl
- Gmsh thread-safety documented: \`GMSH_LOCK\` (ReentrantLock) serialises all Gmsh API calls
- CONTRIBUTING.md, CHANGELOG.md, CODE_OF_CONDUCT.md added
- README.md rewritten with correct working examples

**Testing & API (Phase 4):**
- Consistent \`def_problem\` / \`set_model\` / \`simulate(mdl, scene, conditions)\` API for both Stokes and LinearElasticity (multiple dispatch)
- Fixed \`LinearElasticity\` constructor field ordering, \`gaussian_quadrature\` positional args (19 occurrences), \`reset_model!\`/\`update_model!\`

**Project.toml / CI:**
- Authors field closing \`>\` fixed
- 13 missing \`[compat]\` entries added
- CI test job now triggers on both \`main\` and \`dev\` branch pushes

## Test plan

- [ ] CI passes on this PR (Ubuntu × macOS × Windows, Julia 1.11)
- [ ] Docs job builds without error
- [ ] After merge, verify CI triggers on first push to \`main\`
- [ ] Confirm \`Pkg.test("smearFEM")\` locally: **90,925 tests, 0 failures**

🤖 Generated with [Claude Code](https://claude.com/claude-code)
EOF
)"
```

### 4f. Merge (manual)

Review and merge on GitHub once all CI checks are green. Do not auto-merge.

---

## After Merge (next steps outside this plan)

Once on `main`:
1. **Zenodo DOI** — link the GitHub repo at zenodo.org, trigger a release for v1.0.0 to get a citable DOI
2. **GitHub release** — tag `v1.0.0` (TagBot may do this automatically once the registry is involved)
3. **Verify docs** — confirm the Documenter.jl badge links resolve correctly
4. **Julia General Registry** (optional) — submit via JuliaRegistrator to enable `Pkg.add("smearFEM")`
