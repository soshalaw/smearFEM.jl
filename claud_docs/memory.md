# memory.md

Running log of what was built, decisions made and why, and non-obvious findings/gotchas hit
along the way. See `CLAUDE.md` for the current-state architecture summary; this file is the
history/reasoning behind it.

## 2026-09-14 — 2026-09-15: replicate-statistics review, verification plots, pooled stats

Started from a review of `src/analysis/statistical_analysis.jl`'s handling of physical vs.
synthetic post-analysis (`post_analysis_real`/`post_analysis_bulk`/`post_analysis_const` in
`experiments/model_optimization/test_opt_stokes.jl`), which surfaced an uncommitted fix already
in progress: `bias_vs_reference` and `normalized_replicate_stats` now degrade to `NaN`/zero at
`n=1` replicate instead of crashing on `TDist(0)`.

**Discussion → build arc:**
1. Walked through exactly what `replicate_stats(per_run_mean)` computes (mean/SD/CI of each
   replicate's own time-average) and confirmed it discards within-run temporal spread by
   construction — a replicate that oscillates and one that's flat report the same
   `per_run_mean`.
2. Built `_scalar_band!` + wired it into `_stat_fig` (test_opt_stokes.jl): overlays the reported
   scalar as a flat band on top of the existing runs + pointwise ribbon, so the collapsed number
   can be checked against the curves it came from. Iterated on it live against the physical
   (5 replicates) and synthetic/simulated (6 replicates each: constant-viscosity `sim_data`,
   bulk-viscosity `syn_data/Stokes`, bulk-viscosity `syn_data/Carreau`) result trees.
3. Added nested SD/CI bands (`bands=:sd`/`:ci`/`:both` on `_scalar_band!`), then a `ribbon=false`
   mode on `_stat_fig` producing a `*_scalar_only.pdf` companion next to every `*_ci.pdf` — runs
   + scalar band only, no pointwise ribbon, so a bad scalar can't hide behind an unrelated
   time-varying band.
4. **That companion plot caught a real bug**: `h_pred/h_gt` and `rel. error pred [%]` are `NaN`
   across the board for both bulk-viscosity trees (Carreau and Stokes), `n` inflated to
   `2×n_runs` in the CSV. Root cause: `_align_windowed` intentionally leaves ungapped time steps
   as `NaN`; `mean()`/`std()` on a series containing any `NaN` returns `NaN` for the whole
   cross-replicate summary. The pointwise ribbon never showed this — it's computed fresh
   per-window from raw columns, unaffected. Not fixed yet (see `CLAUDE.md` gotchas); physical
   data doesn't hit it, reason not chased down.
5. User pushed further: the per-run collapse itself was flagged as lossy ("I see the time
   averaging resulting in a lot of info being lost"). Added `pooled_stats` to
   `statistical_analysis.jl` — grand mean/SD/CI over every pooled `(time, replicate)` value
   directly, no per-run collapse. Confirmed concretely on physical data: `h_pred/h_m` pooled SD
   (0.0252) came out ~2.75× the old per-run-collapsed SD (0.0092), because the sawtooth dips
   right after each window restart were being smoothed away by time-averaging each run first.
   `h_est/h_m`, a steadier series, barely moved (0.0090 vs 0.0089). `_scalar_band!` now plots
   `st.pooled` instead of `st.mean`; `.mean`/`per_run_mean` kept untouched since
   `bias_vs_reference` genuinely needs one independent sample per replicate, not pooled
   autocorrelated points.
6. Wrote up the pooled-vs-per-run derivation as a standalone LaTeX note
   (`pooled_stats.tex`, compiled and verified — not checked into the repo, lives in the
   session's scratch dir) for potential reuse in the paper.

**Non-obvious findings:**
- CI half-width (`t_{n-1,0.975}·sd/√n`) is *wider* than ±1 SD at small `n` (e.g. `n=5`:
  `t₄/√5 ≈ 1.24`), not narrower — small-sample `t` correction outgrows the `√n` benefit. Only
  once `n ≳ 7` does CI become the tighter of the two. `_scalar_band!`'s nested-band ordering
  handles both regimes by comparing actual half-widths rather than assuming.
- Pooling recovers `Var_total = Var_within + Var_between` (law of total variance) in a balanced
  design, but the *CI* on the pooled mean still must divide by `√n_runs`, not `√n_obs` — time
  steps within a run are autocorrelated, not independent draws, so using the full pooled count
  would understate the CI by roughly `√n_steps`.

**Process note:** this session did a lot of implementation before checking whether `memory.md`
existed — it didn't. Retroactively backfilled from conversation history on 2026-09-15. Going
forward: update this file and `CLAUDE.md` after each meaningful change, not just at session end.
Set up a global Stop-hook (`~/.claude/hooks/memory-reminder.sh`, registered in
`~/.claude/settings.json`) that nudges via `systemMessage` when a project's source files are
newer than its `memory.md`/`CLAUDE.md` — non-blocking, applies to every project, not just this
one. Caught two real bugs while building it: (1) `stat -f`/`stat -c` portability — GNU `stat -f`
doesn't error on the wrong flag, it silently switches to filesystem-info mode and "succeeds"
with garbage, so a BSD-first fallback chain never reaches the correct GNU form; fixed by trying
GNU (`-c %Y`) first. (2) `CLAUDE.md` is itself gitignored in this repo (confirmed intentional by
the user) — `git status` has no tracked baseline to diff an ignored+untracked file against, so a
git-diff-based detector would never see a CLAUDE.md-only edit; the hook uses file mtime instead,
which works regardless of git-tracking state.

## 2026-09-16: `pooled_stats` corrected to a two-stage (per-run-then-pooled) computation

User feedback: the original `pooled_stats` (`std(vcat(series...))`, flattening every raw
`(time, replicate)` value into one sample) still quietly treated autocorrelated time steps
within a run as independent contributions to that sum-of-squares — the CI already accounted for
this (dividing by `√n_runs` not `√n_obs`), but the `sd` point estimate itself didn't.

Fixed by computing each replicate's own mean *and* SD first (`per_run_mean[i]`,
`per_run_sd[i]` — purely descriptive of one trajectory, no independence assumption needed for
that), then combining across the `R` replicates with the standard pooled-variance formula:
within-run variance = `Σ(Tᵢ-1)·per_run_sd[i]² / Σ(Tᵢ-1)` (properly weighted by each replicate's
own degrees of freedom), between-run variance = `var(per_run_mean)` (the same quantity
`replicate_stats(per_run_mean).sd²` already reports), `sd = √(within + between)`. Verified
against hand-computed toy data (two 3-point series, known within/between variance) — exact
match. `_scalar_band!`/`_stat_fig` needed no changes since `st.pooled`'s field names are
unchanged, only the internal computation.

Practical impact on this repo's data was small (physical `h_pred/h_m` pooled SD: 0.02517 → 0.02553,
`h_est/h_m`: 0.00902 → 0.00987) — the flattened approach was already a decent numerical
approximation at ~300 steps/run, but the two-stage version is the methodologically correct
decomposition rather than an approximation to it, and matters more for shorter series. All
three post-analysis trees (physical, constant-viscosity, both bulk-viscosity) re-ran clean with
no errors after the change.

## 2026-09-16: third approach (B) added; full A/B/C plot + CSV output, per-approach folders

Follow-on from the pooled-stats correction above. User pointed out `pooled_stats` (approach C)
still has a real conceptual gap: it collapses each replicate over *time* first, so it cannot
tell a pattern shared by every replicate (e.g. the window-restart sawtooth — a property of the
fitting method, identical across replicates) from genuine idiosyncratic noise. Added approach
B, `pointwise_pooled_stats`: collapse across *replicates* first, at each fixed time step (this
is exactly the ribbon's `pointwise_mean`/`pointwise_sd`), then pool those `T` estimates over
time. Anything identical across replicates at a given instant cancels out of the per-timestep
SD before time is ever touched, so B is unbiased for genuine cross-replicate disagreement
regardless of shared time-structure, where C is biased upward by exactly the variance of that
shared structure. Verified the derivation with a formal generative model
(`v_{i,t} = μ + f(t) + b_i + e_{i,t}`) — full writeup in `~/approaches_abc.md` (LaTeX source
also produced, not kept).

**Concrete result (physical data)**: C is 2.2–4.5× larger than B specifically on `..._pred`
quantities (`h_pred/h_m`, `rel. error pred [%]`, contour cost pred) — matches the sawtooth being
shared across replicates. `..._est` quantities keep A/B/C within ~10% of each other, consistent
with little shared time-structure there. This means the previous C-only headline number for
`h_pred/h_m` reproducibility (0.0255) was overstating irreproducibility by roughly 2× relative
to B (0.0119) — B is the number to use when the claim is about genuine run-to-run
reproducibility, not about single-run wander.

**Also cross-checked and corrected a colleague's draft formula** that looked like a genuine
"per-timestep-first" two-stage derivation but was algebraically identical (via the standard
one-way ANOVA sum-of-squares identity, Total SS = Within SS + Between SS) to flattening every
raw value and computing `std` with `RT-1` "degrees of freedom" — silently reintroducing the
autocorrelation problem the whole A/B/C exercise was meant to fix, if that `RT-1` were ever used
for an actual CI rather than just a descriptive point estimate.

**Implementation**: `approach_stats(st, :A|:B|:C)` is the one dispatch point for picking a field
off a `normalized_replicate_stats` result — nothing else should branch on field names directly.
`_scalar_band!` now takes the already-selected stat NamedTuple directly (was `st.pooled`
hardcoded). `_stat_fig` gained an `approach` kwarg. New call-site helpers
`_stat_fig_approaches`/`_write_approach_reports` write one self-contained `approach_<A|B|C>/`
folder per approach (plots + a CSV with mean/sd/ci_lo/ci_hi per quantity) under every existing
`plot_path_global`, for all four trees (physical, sim-constant, synthetic-bulk-Stokes,
synthetic-bulk-Carreau) — replacing the manually-duplicated `_ci`/`_scalar_only` call pairs at
every quantity across all three `post_analysis_*` functions. The older single combined CSV
(all approaches as rows of one file) is kept unchanged alongside the new per-approach ones.
Verified: toy-data exact match for B same as was done for C earlier; all four trees re-ran with
zero errors; spot-checked CSVs and rendered plots for physical data match the hand-computed
comparison table from the design discussion exactly.

## 2026-09-16: found and fixed a 5th tree silently skipped — `syn_data/.../constant`

User asked why `experiments/syn_data/optimization/Stokes/force/constant/Hex2_16/cylinder` (6
replicates, real data) wasn't among the four trees processed above. Root cause: `plot_results`'s
driver only ever pairs `data_type == "synthetic"` with `viscosity_type == "bulk_viscosity"` —
`("synthetic", "constant")` is never in its sweep, even though that exact key already exists in
`PLOT_PRESETS` (confirms it's a real, intended configuration, not something that was never meant
to exist). My `run_post_analysis_synthetic.jl` driver script inherited the same blind spot by
copying `plot_results`'s path logic instead of checking the filesystem directly.

Running `post_analysis_const` against it directly hit `@warn "No simulated contours, skipping"`
for all 6 replicates — `sim_data/<view>/2D_border_points_est/` doesn't exist there. Investigated
by diffing against the working `sim_data/.../constant` tree before touching anything:
- **Not a simple rename.** The old tree's contour data lives at the flat, unqualified
  `sim_data/2D_border_points/` (no `view_1/` nesting, no `_est` suffix) — a genuinely older,
  pre-multi-view directory layout, not just a differently-named copy of the same thing.
- **Not even a reshape.** Old format: one `x,y` pair per row, 60 rows. New format
  (`view_1/2D_border_points_est/`): two rows (all-x, all-y), 50 values each. The point *count*
  differs — the border-extraction algorithm itself changed, not just the serialization. A
  from-scratch conversion script would have to guess how to resample 60→50 points, i.e.
  fabricate data, not convert it. Declined to write that; flagged it and asked the user first
  rather than guessing.
- **The actual fix, once the user pushed on "what do we actually need":** the height/rel-error
  stats that feed A/B/C come *only* from `est_h.csv`/`gt_h.csv` (`push!(norm_est_runs, ...)`
  at test_opt_stokes.jl, `post_analysis_const`) — both present and correctly formatted in the
  old tree. `_get_borders`'s output (needing the missing/incompatible contour directory) is used
  *only* for the `contour_plt`/`contour_plt_zoom*` diagnostic plots, never for height stats.
  The early `continue`-on-missing-borders gate was strictly more conservative than necessary.
  Fixed by making border availability a `have_borders` flag rather than a hard skip: when false,
  the six border variables are set to `nothing` and every `contour_plt*` line that reads them is
  individually guarded (`have_borders && Plots.plot!(...)`); everything else (height stats, η/β
  plots, cost plots) proceeds unchanged. No data touched, no resampling guessed at.

Re-ran clean, zero errors. The tree now has full `approach_A/B/C` output — and a genuinely
useful result: `h_est/h_gt`'s C-SD (0.00775) is ~3.8× A's (0.00204), because every one of the 6
runs shows the same steady downward drift over the fit window (visible directly in the raw
curves) — another real instance of the shared-time-pattern effect the whole A/B/C framework
exists to separate out, on top of the sawtooth case found in the bulk-viscosity trees earlier.

## 2026-09-16: `post_analysis_const` was mislabeling ~83% of its data as "Estimation"

Investigated at user's request: does `replot`/`post_analysis_const` correctly separate
in-sample fit agreement from out-of-sample forward prediction for constant-viscosity trees, the
way `post_analysis_bulk`/`post_analysis_real` already do for windowed fits? Answer: no — traced
`optimize`'s `gt_viscosity_type == "constant"` branch (lines ~1226-1289) and confirmed
`fit_model` only ever sees observations from the first `sim_time_exp` seconds (`ObsDataList`
sliced to `_range`, line 1229-1232), but the final `est_h` is produced by re-simulating with the
*converged, fixed* `(η, β)` over the *full* ground-truth duration `sim_time_gt` (line 1236,
1280). Confirmed on disk: `est_h.csv` has 301 rows (≈30s at dt=0.1) while the fit only used the
first 51 (5s) — **~83% of every constant-viscosity `est_h` series is forward prediction, not
estimation**, but `post_analysis_const` pushed the whole thing into `norm_est_runs` unsplit.
(The plotting code already half-knew this — `h_norm_plot_5`/`rel_height_error_glob_plot_5` drew
a vline and "Prediction" arrow annotation at a hardcoded `x=5` — but the underlying replicate
stats never acted on it.)

**Fix**: added `norm_pred_runs`/`rel_err_pred_runs` (plain vectors, for
`normalized_replicate_stats`) and `norm_pred_segs`/`rel_err_pred_segs` (the `(range, values)`
segment format `_stat_fig` already supports for windowed prediction, reused here for a single
segment) to `post_analysis_const`. Split point computed properly (`n_est = round(Int, sim_time
/ t_steps) + 1`) rather than the pre-existing hardcoded `51` littered through the plotting code.
`height_stats`/`_stat_fig_approaches` calls extended to carry both groups, exactly mirroring
`post_analysis_bulk`'s Estimation/Prediction pattern, with `vlines=[sim_time_exp]` marking the
boundary. Scope: only the `noise_level == 0.0` branch pushes into these arrays at all (checked);
the noise-repeatability branch is a separate feature, untouched.

**Numbers moved substantially** — this was not cosmetic. Syn-constant `h_est/h_gt`: mean/SD
went from 0.9855±0.0078 (C, mixed) to a *tight* 0.9970±0.0018 once genuinely restricted to the
in-sample region; the previously-reported spread was almost entirely prediction-region drift.
`h_pred/h_gt` (now separated): 0.9832±0.0063. Sim-constant same pattern, smaller magnitude
(no rendering noise there). Re-ran both affected trees (sim-constant, syn-constant — the only
two using `post_analysis_const`) clean, zero errors; visually confirmed both plots now show a
narrow flat orange Estimation band up to the vline at 5s and a widening blue Prediction band
after it, matching the underlying per-run curves exactly.

Did not also fix `replot`'s single-experiment `h_plt`/`error_plt` (still one unsplit
"Estimation" label) — out of scope for this pass, which targeted the cross-replicate stat
plots specifically; worth the same split later if those per-experiment diagnostics matter.

## 2026-09-16: CLOSED — the long-standing `h_pred/h_gt` NaN bug in bulk-viscosity trees

Root cause found and fixed, resolving the bug flagged repeatedly since 2026-09-15
([[project_replicate_stats_verification]]). It was never about `_align_windowed`'s NaN gap
(that was a red herring — a real, separate, correct behavior). The actual cause:
`post_analysis_bulk`/`post_analysis_real`/`replot`'s windowed branches all do
`window_dirs = readdir(leaf.exp_path); for window_dir in window_dirs`, expecting only genuine
per-window result directories. At some point a `single_window` sibling directory appeared next
to `multi_window` under each replicate's `view_1/` — **a different experimental condition (one
global static fit evaluated against the same window grid, for comparison), not another window
of the multi-window re-fit.** It shares `data_ranges.csv`/`t_windows.csv` with `multi_window`
(same window grid, different fit) but has a differently-shaped `pred_h.csv` (301 rows vs.
`multi_window`'s 306). The three occurrences of the `window_dir == "Results" || ...` skip-list
(lines ~2848, ~4069, ~4466) had no entry for `single_window`, so it silently got processed
alongside `multi_window`: doubled `n` (6 replicates × 2 conditions = 12, matching the exact
symptom reported) and its shape mismatch broke `_align_windowed`, producing the `h_pred/h_gt`/
`MOSD_pred/MOSD_gt` `NaN`. Physical data never showed the bug simply because its result trees
never had a `single_window` folder to begin with — the bug was latent there too, just
unexercised.

**Fix**: added `single_window` to all three skip-lists. Re-ran the synthetic bulk-viscosity
trees (Stokes and Carreau) clean, zero errors. Stokes bulk-viscosity now reports exactly `n=6`
with no `NaN` anywhere — and its numbers **exactly reproduce the existing paper table's Group B
entry** (`h_est/h_gt` 0.9960±0.0008, `h_pred/h_gt` 0.9949±0.0008, `MOSD_est/MOSD_gt`
1.0029±0.0005, `MOSD_pred/MOSD_gt` 1.0036±0.0005 — all match to 4 decimal places), confirming
the underlying data was correct all along and only the post-processing code had regressed
(`single_window` was presumably added by a later ablation-study run, after the table was
originally generated). Carreau now reports `n=5`, not 6 — also correct, not a new bug: replicate
5's `multi_window` folder is genuinely incomplete (`_skip_incomplete` correctly excludes it; a
data-completeness gap, not a code issue).

## 2026-09-16: Publish-readiness refactor — STARTED, Phase 1 done (uncommitted)

Full plan lives at `/home/soshala/.claude/plans/peppy-hopping-elephant.md` — read it first when
picking this up. It covers `src/` + `experiments/` for both style conformance and structural
cleanup, in 9 phases. What follows is status, not a restatement of the plan.

**Scope decisions already made (do not re-litigate):**
- `experiments/model_optimization/test_opt_stokes.jl` (5,247 lines): **dedup + style only, no
  file split.** It had just been verified against the paper's numbers; splitting it was judged
  not worth the risk. The 6-way split seams are documented in the plan if that ever changes.
- Latent bugs: delete the dead one, fix the live one.
- `experiments/stokes/stokes_model.jl`: legacy, delete.
- Commits: one per batch, diff reviewed before committing. User runs the git commands.

**Policies adopted** (now also in `CLAUDE.md`'s Style section): file-size flag at >600 lines /
strong flag >1000; `println` allowed in functions whose whole purpose is formatted console
output (`show_config`, `print_replicate_table`); docstring template follows
`docs/coding_guide.md` (Arguments/Returns/optional Example), not `CONTRIBUTING.md`'s variant.

**Phase 0 — done, committed** (`364118e`, `2677cc6`, `a6d7a0a`, `05bc5dd`). The previous
session's uncommitted work was split by concern so refactor diffs stay uncontaminated. Scratch
edits to `single_simulation.jl` (`t_steps` 2.5e-3→5.0, `ne` 16→6, swapped `const_vel` call)
were reverted, not committed.

**Baseline test counts** (captured before Phase 1, logs in that session's scratchpad): stokes
3770, mesh_templates 144, tikhonov 41, cost_functions 34, fem 9+4, stokes_optimization clean
(no testset, bare `@test`). Phase 1 re-run reproduced all of these exactly.

**Phase 1 — done, committed as `866677d`.** 2,605 deletions, zero insertions:
`stokes_solver.jl.bak` (1394), `src/fem/mesh_structured.jl` (662, byte-identical duplicate of
the `meshes/` copy), `assemble_system_A_dense`+`assemble_system_B_dense`+2 dead `println`s
(337), `stokes_model.jl` (87), `init_cylinder` (88), the dead 4-arg `plot_covariance!` (20),
`write_gt_stokes.jl`'s commented calibration table (17). Package loads clean; all 5 test files
pass identically to baseline.

Side effect worth knowing: `stokes_solver.jl` dropped 1347→1010 and `smearOptimize.jl` 936→848,
so **the planned Phase 4.1 file split is no longer warranted** — both are under the threshold.

**Plan corrections found by grepping rather than trusting the survey:**
- `apply_boundary_conditions_dense` is **not** dead — called at `run_example.jl:58`. Kept. Only
  the two `assemble_system_*_dense` functions were dead.
- The dead 4-arg `plot_covariance!` was confirmed genuinely unreachable: both call sites in
  `test_opt_stokes.jl` pass 3 positional args, which dispatch to the keyword method.

**Two pre-existing bugs found, neither fixed, both needing a decision:**
1. `test/tikhonov_test.jl` **fails at HEAD.** All 41 tests pass, then line 101's trailing
   `println` calls `round(v["η"], sigdigits=6)` on a `Vector{Float64}` and throws `MethodError`.
   Needs `round.(...)`. It goes unnoticed because `tikhonov_test.jl` is the one test file
   missing from `runtests.jl`'s include list — fix both together in Phase 5h.
2. `experiments/stokes/write_gt_stokes.jl`'s `_get_F_ext` ladder tests `β_gt <= 5e3` *before*
   `β_gt <= 2e3`, so the 2e3 branch is unreachable and β in (1e3, 2e3] silently gets the ×1.01
   multiplier instead of the intended ×0.995. Reordering changes the force used to generate
   ground-truth data, so it must not be fixed silently — check whether any existing synthetic
   tree used β in that range first.

**Next step**: Phase 2a (three stray `println` → `@info`/`@debug`), then 2b (the `readData` /
`get_cMat` renames and the internal-helper `_` prefixes). Open questions still unanswered are
listed at the end of the plan file.

## 2026-09-16: Claude docs moved to `claud_docs/`, and commit-history conventions

All Claude-facing documentation now lives in `claud_docs/` (gitignored): `CLAUDE.md`,
`memory.md` (this file) and `project_geometry_dispatch_refactor.md`. The repo root keeps a
one-line `CLAUDE.md` containing `@claud_docs/CLAUDE.md` — **that stub is load-bearing**, because
Claude Code only auto-loads `CLAUDE.md` from the project root; the `@` import pulls the real
file in. Both the stub and the folder are in `.gitignore`. Genuine project documentation
(`README.md`, `CONTRIBUTING.md`, `docs/`, `mesh_files/*.md`) deliberately stayed put — the
folder is only for material written for Claude rather than for the project's readers. The
global `~/.claude/CLAUDE.md` and the Stop hook `~/.claude/hooks/memory-reminder.sh` were both
updated to match; the hook checks `claud_docs/` paths first and falls back to the old root names
for projects not yet migrated.

Two git conventions settled at the same time, both now in the plan file:
- **Refactor commits get a bare conventional-commit subject line, no body.** The diffs are
  mechanical and self-explanatory. Behaviour-changing commits still explain why.
- **Never stage anything under `claud_docs/`.** Keep updating those files after every meaningful
  change as the global guideline requires — just never commit them, and never offer to.

Housekeeping done while applying this: the six unpushed commits had their message bodies
stripped via `git filter-branch --msg-filter 'head -n 1'`, and `memory.md` was purged from all
of them with an `--index-filter`, since it had been wrongly committed earlier in the session.
`memory/project_geometry_dispatch_refactor.md` was **not** purged — it entered history back in
the already-pushed `73c6f36`, so removing it would mean rewriting published commits; it is
simply deleted going forward. Commit hashes therefore changed twice; any hash written down
before 2026-09-16 is stale.

## 2026-09-17: Phase 2a — stray `println` → `@debug`

Three sites converted, all to `@debug` rather than `@info`: `run_example.jl:41` (mesh-load path,
fires inside a per-timestep helper that can be called in a loop), `io/io.jl:174` (node-count
summary sitting directly under an `@info`/`@warn` sanity-check branch that already reports the
meaningful result), `io/analysis_plots.jl:110` (dumps a whole `h_Vector` — a leftover debug
dump). In `analysis_plots.jl` the `@debug` was moved *above* the `# define a reference height`
comment, which describes the `h_ref = ...` line, not the print.

**`src/utils.jl`'s six `println`s are not in scope and must stay** — they are `println(io, ...)`
writing into an open `time_log.txt` handle inside `write_time_log`, i.e. file formatting, not
console logging. The "no `println`" rule is about stdout. This is a third exemption category
alongside `show_config`/`print_replicate_table`, and a naive repo-wide grep-and-replace would
break the time log.

Commented-out `println` lines remain in `stokes_solver.jl` (4), `run_example.jl` (4) and
`PostProcess.jl` (1) — dead debug scaffolding, deliberately left for the Phase 2d comment sweep
rather than widening 2a.

Verification: package loads clean; all six test files reproduce the baseline counts exactly
(stokes 3770, fem 9+4, cost_functions 34, mesh_templates 144, tikhonov 41-then-crash,
stokes_optimization clean). `tikhonov`'s `exit=1` is the known unfixed line-101 bug, unchanged.

Noted in passing, **not** changed (out of scope): `simulate_single_tstep_stokes` builds
`filepath_mesh` from `dirname(dirname(@__DIR__))` instead of `get_mesh_dir()`, bypassing the
`SMEAR_MESH_DIR` → `config.toml` resolver. Source-relative, so not an absolute-path violation,
but inconsistent with the rest of the codebase.

## 2026-09-17: Phase 2b — naming

Two exported renames (`readData` → `read_data`, `get_cMat` → `get_c_mat`) and nine internal
helpers given `_` prefixes. **Every call site was inside `src/`** — the pre-rename grep across
`src/`, `experiments/`, `test/` and `docs/` found nothing external, so no experiment script or
result-generating code was touched. Post-rename grep for the old names returns empty, and the
new-name hit counts match the old-name counts file-for-file.

| old | new | note |
|---|---|---|
| `readData` | `read_data` | exported; export list updated |
| `get_cMat` | `get_c_mat` | exported; export list updated |
| `animate2D` / `animate3D` | `_animate_2d` / `_animate_3d` | see judgement call below |
| `normalize` | `_normalize_field` | also stops shadowing `LinearAlgebra.normalize` |
| `truncate_colormap` | `_truncate_colormap` | |
| `extract_p_from_u_nodes` | `_extract_p_from_u_nodes` | |
| `read_sparse_mat` | `_read_sparse_mat` | |
| `armijo_line_search` | `_armijo_line_search` | |
| `backtrack_line_search` | `_backtrack_line_search` | |
| `val_check` | `_val_check` | |

**Judgement call:** the plan said only "add `_` prefixes", which for `animate2D` would give
`_animate2D` — still not `snake_case`. Since 2b *is* the naming phase and the guideline mandates
`snake_case` for functions, they became `_animate_2d`/`_animate_3d` instead. Zero call sites
outside `plotting.jl`, so this is trivially revertible if the literal plan wording is preferred.

**Plan correction:** `get_cMat` had only one call site (`stokes_solver.jl:42`), not the two the
plan lists — Phase 1 deleted the `:168` one along with the dense assembly functions. Don't trust
the plan's line numbers; re-grep.

**Three renamed helpers have zero call sites** and are not exported, i.e. unreachable dead code
that Phase 1's sweep missed: `_normalize_field` (`plotting.jl:400`), `_truncate_colormap`
(`plotting.jl:427`), `_read_sparse_mat` (`io.jl:404`). Renamed rather than deleted — deletion is
a decision, still open.

Verification: package loads clean, `read_data`/`get_c_mat` both resolve in the module, all six
test files reproduce baseline counts exactly (stokes 3770, fem 9+4, cost_functions 34,
mesh_templates 144, tikhonov 41-then-the-known-crash, stokes_optimization clean).

## 2026-09-17: Phase 2c — docstrings for undocumented exported symbols

32 docstrings across 10 files. Undocumented exported bindings went **39 -> 7**, the remaining 7
being exactly `geometries.jl`'s (`AbstractGeometry`, `Cylinder`, `Cuboid`, `Disk`, `Square`,
`Segment`, `ndim`), deferred to Phase 3a so that file is touched once.

**How to enumerate them** (the plan's per-file counts were close but not right — it missed
`run_example.jl` entirely, plus `Model` and `AbstractMeshgrid`):

```julia
for n in sort(names(smearFEM))
    d = string(Docs.doc(Docs.Binding(smearFEM, n)))
    occursin("No documentation found", d) && println(n)
end
```

**`@doc(getfield(smearFEM, n))` does NOT work** — `@doc` is a macro operating on the
*expression*, so it documents `getfield` itself and reports zero undocumented symbols. Use
`Docs.doc(Docs.Binding(mod, name))`.

**`write_gt_data` was not actually undocumented** — it had a docstring that Julia never attached:
indented onto a stray continuation line, naming `write_data` instead of `write_gt_data`, and
separated from the `function` line by a blank line. Repaired in place rather than adding a
second one. Worth grepping for the same shape elsewhere.

**Docstring style choices made here:**
- For the 5 `set_model` / 2 `def_problem` methods, the `Cylinder` method carries the full
  `# Arguments` block and the others say "see the `Cylinder` method for the shared arguments"
  plus what differs. Repeating all ~12 arguments five times would have added ~125 lines of
  duplication to a 167-line file.
- The 5 `Meshgrid*` structs share one paragraph naming the common `AbstractMeshgrid` fields, then
  list only the geometry-specific ones.
- Documented actual behaviour, not intent, where they differ: `plot_data`'s `n` kwarg is
  immediately overwritten by `n = 10` internally, and `arrow0!`'s `color_u`/`color_v`/`color` are
  unused (the colour-mapped plot calls are commented out).

**FINDING — `SqueezeFlow`'s boundary-condition keys are swapped (NOT fixed, deliberately).**
`set_boundary_cond` returns `(q_upper, q_side, q_lower, C_uc)`; `def_problem` passes them on as
`[q_tp, q_side, q_btm]`; the `SqueezeFlow` constructor maps position 2 -> `:bottom` and position
3 -> `:border`. So `q_d[:bottom]` holds the **side** nodes and `q_d[:border]` the **bottom** ones,
and the solver's `q_d .= (μu_btm*q_d_cached_btm + μu_side*q_d_cached_brdr)` applies the
bottom-plate multiplier to side nodes and vice versa. **No current number changes**, because both
`μu_btm` and `μu_side` are hardcoded `0` (`stokes_solver.jl:703-704`, `run_example.jl:180-181`),
making the whole expression zero either way. It only bites if the bottom plate or side wall is
ever driven. The `SqueezeFlow` docstring documents the real mapping, not the intended one.

**Also noticed, not acted on** (deferred to the end-of-refactor sweep by the user's instruction):
- `get_c_mat`'s pre-existing docstring signature line is stale — it shows
  `get_c_mat(type; λ=..., μ=..., Young=..., ν=...)` but the function is `(c1, c2; type="lame")`.
- `smearFEM.jl`'s export comment puts `plot_covariance!` in `analysis_plots.jl`; it is in
  `plotting.jl`.
- `CLAUDE.md` claims `MeshgridCylinder` precomputes and stores `C_vol`/`C_top`/`C_btm`/`W`. It has
  no such fields — the claim is stale.
- `simulate`'s `output` vector is only populated under `"force"` control; the `push!` is commented
  out in the `"velocity"` branch, so it returns empty there.

Verification: package loads clean; all six test files reproduce baseline counts exactly.

## 2026-09-17: SqueezeFlow boundary-key mismatch — FIXED (relabel only)

**What was wrong.** `set_boundary_cond` (`stokes_solver.jl`) returns `q_upper, q_side, q_lower,
C_uc`, but `def_problem` passed them to `SqueezeFlow` as `[q_tp, q_side, q_btm]` while the
constructor maps position 1/2/3 to `:top`/`:bottom`/`:border`. So `q_d[:bottom]` held the *side*
vector and `q_d[:border]` the *bottom-plate* one.

**The sharper detail:** `q_side` is allocated in `set_boundary_cond` and **never written to** —
it is always all-zero. Only `q_upper` (nodes at `z == z1Bound`) and `q_lower` (nodes at `z == 0`)
are ever marked. So `μu_btm * q_d[:bottom]` multiplied a zero vector and could never do anything,
while `μu_side * q_d[:border]` was the term that would have driven the bottom plate.

**Why no result was ever affected.** Both multipliers are hardcoded `0` (`stokes_solver.jl:738-739`,
`run_example.jl:202-203`), and every call site of `simulate_single_tstep_stokes` — the only
function taking them as arguments — passes `0` for both. `q_d` is therefore identically zero in
the force branch and exactly `cParam[iter]*q_top` in the velocity branch. The bottom plate is
pinned correctly regardless, because its z-DOFs enter `rCol` and are eliminated from the free set
by `C_uc = C[:, setdiff(1:size(C,2), rCol)]` — the constraint comes from the elimination, not
from these terms.

Note the direct path was always wired correctly: `run_example.jl`'s
`q_d = (μu_btm*q_btm + μu_tp*q_tp + μu_side*q_side)` uses the raw tuple with matching names. Only
the `SqueezeFlow` dict route shifted.

**Fix chosen: relabel at the call sites, not in the constructor.** The constructor's contract
(`:top`, `:bottom`, `:border` = top plate, bottom plate, lateral surface) is coherent as written;
the two callers simply passed the middle two in the wrong order. Both `def_problem` methods now
pass `[q_tp, q_btm, q_side]`. Two lines, one file. The wider options (deleting the dead side slot,
or removing the vestigial `μu_btm`/`μu_side` knobs) were considered and deliberately not taken.

**Numerical verification** (scratch script `bc_check.jl`, not kept): 20-step force-controlled
`Cylinder(10, 20)`, `Hex2/Hex1`, `ne=6`, `η=β=100`. Ran with the new order and the old, wrote
`output`, final `pos3D`, `velocity`, `pressure` and `gradList_3d[:,:,1]`, and diffed — all five
**byte-identical**. Full test suite reproduces baseline.

## 2026-09-17: Phase 2d — misapplied comments and the `smearFem` typo

22 comments removed, not the plan's estimated ~8. Two families, both matching the plan's
"misapplied to unrelated calls" criterion:

- **6 x `assemble_system_B(...)  # assemble the stiffness matrix`** — flatly wrong. That
  function's own docstring says it returns a "sparse coupling matrix relating velocity and
  pressure degrees of freedom". The 6 sitting on `assemble_system_A` are accurate but pure
  restatement, so the whole named family went.
- **10 x `drdη/drdβ = [...]  # solve the system of equations`** — those lines assemble
  right-hand-side derivative vectors; they solve nothing. Careful when grepping for these: a
  pattern like `drd[ηβ].*# solve` also matches the 10 **legitimate** `dsoldη = K\(drdη + ...)`
  lines, which genuinely do solve and were correctly left alone.

Also fixed `smearFEM.jl:94`'s `end # module smearFem` -> `smearFEM`.

**Deliberately left:** roughly 30 comments that are redundant but *accurate* — `# solve the
system of equations` on the real `K\r` solves, plus the `# assemble the system of equations`,
`# assemble the solution` and `# apply the neumann boundary conditions` families. The style
guide's "never restate what the next line does" would condemn all of them, but that is a much
larger and more opinionated sweep than Phase 2d authorises. Added to the end-of-refactor list.

**Verification:** a diff check confirmed comments-only — for all 23 changed line pairs the text
before the `#` is byte-identical. Package loads; all six test files reproduce baseline.

**Process note / mistake to avoid repeating:** Phase 2d was started before the Phase 2c batch had
been committed, and both touch `run_example.jl`, `stokes_solver.jl` and `smearFEM.jl`. The two
phases are now intermixed in those files and **cannot be separated by `git add <path>`**. Commit
each batch before starting the next one when the file sets overlap.

## 2026-09-17: Phase 3a — geometries.jl, and the keyword-only question RESOLVED

`src/fem/geometries.jl` rewritten (30 -> 108 lines): the five geometry wrappers are now
`mutable struct` with a default for every field, and all 7 exported symbols
(`AbstractGeometry`, the 5 types, `ndim`) have docstrings. **Undocumented exported bindings are
now 0**, down from 39 at the start of Phase 2c.

**RESOLVED: defaulted-positional, NOT keyword-only.** `docs/coding_guide.md` mandates
"keyword-only inner constructors with safe defaults for every field", but Julia's strict
`function T(; r, h)` form **rejects positional calls**, and the call-site grep found **12
construction sites, every one positional** (`Cylinder(r, h)`, `Cuboid(lx, ly, lz, edge_radius)`
in `run_example.jl`, `single_simulation.jl`, `test_opt_stokes.jl`, `tikhonov_test.jl`,
`stokes_optimization_test.jl`) and **zero** keyword-style. Keyword-only would have broken all
twelve for no gain. So `geometries.jl` satisfies "mutable" and "defaults on every field" but
deliberately **not** "keyword-only" — the one documented exception to that rule.

`Disk`, `Square` and `Segment` have **no construction sites anywhere** and could have been
keyword-only for free, but splitting the convention across five sibling types is worse than one
documented exception.

**Silent semantic change, checked first:** these were immutable `struct`s, so `==` compared
field-wise; `mutable struct` switches `==` to identity. Verified inert — no `==` comparison, no
`Dict`/`Set` keying and no `hash` of a geometry anywhere in `src/`, `experiments/` or `test/`.

Parametric typing is preserved: `Cylinder(10, 20.0)` still infers `Cylinder{Int64, Float64}`.
Verified by construction that both caller forms work, all five types construct bare from
defaults, and fields mutate. Full suite reproduces baseline.

## 2026-09-17: Phase 3b — SqueezeFlow defaults

`SqueezeFlow`'s inner constructor now defaults **8 of its 9 parameters**. Defaults are empty or
zero-sized (`Float64[]`, `spzeros(Float64, 0, 0)`, three empty sparse blocks for `q`) except
`control="force"` and `viscosity_type="constant"`, which take the most common configuration.
A bare `SqueezeFlow(mdl)` still cannot run a simulation — `sim_time=0.0, t_steps=0.0` makes the
time range degenerate — so no default can silently produce a meaningless run.

**`model` deliberately has no default.** It is typed `AbstractModel`, an abstract type, so there
is no value to put there. Widening the field to `Union{AbstractModel,Nothing}` or deleting it are
both more than "add defaults". This is the second documented exception to the guideline's
"defaults for every field", after 3a's keyword-only one.

**Blast radius was narrower than the plan predicted, not wider:** only **2 construction sites**,
both in `stokes_setup.jl`. Everything else referencing `SqueezeFlow` is a type annotation, a
docstring, or a `scene.β = [...]` mutation in `smearOptimize.jl` (fine — the struct is mutable).

**FINDINGS for the end-of-refactor list:**
- **The `model` field is dead.** Never read anywhere: both consumers' `@unpack` lists omit it
  (`run_example.jl:172`, `stokes_solver.jl:708`), there is no `scene.model` access in `src/`,
  `experiments/` or `test/`, and every function needing the model takes `mdl::Stokes` separately.
  Removing it is the honest fix and would make "defaults for every field" fully satisfiable — the
  exception documented above exists *only* because of a field nothing uses.
- **`q_d` round-trips through a dense representation.** The field is declared
  `Dict{Symbol, Matrix{Float64}}`, so the constructor silently converts each
  `SparseMatrixCSC{Float64,Int64}` to a dense `Matrix{Float64}`; the solver then annotates them
  straight back to `SparseMatrixCSC{Bool,Int64}` on read (`stokes_solver.jl:712-714`) — a
  re-sparsification plus a `Float64`->`Bool` conversion. Harmless in cost (these are n x 1
  columns) but two pointless conversions per solve; declaring the field as what it actually holds
  would remove both.

Verification: package loads; minimal (`SqueezeFlow(mdl)`), partial and full-positional
construction all work; full suite reproduces baseline.

## 2026-09-17: Phase 3c — stokes_single_step_force was broken FOUR ways, now ported

The plan framed 3c as "pass an explicit cache at the call sites". It was far deeper.
`stokes_single_step_force` has **zero callers** (exported, never invoked in `src/`,
`experiments/` or `test/`) and had never run: four independent defects, each masking the next.

| # | Defect | Fix |
|---|---|---|
| 1 | `@unpack ... SIDES = conditions` — no such field on `Conditions`, never used after unpacking | dropped from the unpack |
| 2 | 6 x 1-arg `assemble_system_A/B` / `apply_boundary_conditions` — the plan's target | `cache = BasisFunctionCache(mdl)` before the control dispatch, passed at all six |
| 3 | calls single-view `extract_borders` but indexes `[1]` at 14 sites as if per-view | all 3 sites -> `_get_2D_data(..., angles=rot_angle_cached)` |
| 4 | `nDof * nNodes^ndim` sizing -> 1885^3 ~ 6.7e9 -> `OutOfMemoryError` | exponent dropped at 7 sites |

**Defect 4 is the dateable one.** `nNodes^ndim` is only correct when `nNodes` means
nodes-*per-edge*, as on a structured grid. It has been the total node count since unstructured
meshing landed, and `simulate` has no exponent at its equivalent lines (`stokes_solver.jl:764-767`).
**If you meet `(nNodes_*)^ndim` anywhere else, it is stale — not a convention.**

Defect 3 got simpler than feared: `_get_2D_data` already computes `obs_border_pts =
vcat(pi', qi')` internally, which is exactly what the function was doing by hand, so the manual
`fit_curve` calls were deleted rather than rewritten.

**Correctness cross-validation** (scratch `3c_xval.jl`, not kept): since nothing calls this
function, "no longer errors" proves nothing. Ran it against `simulate` on an identical one-step
force-controlled `Cylinder(10,20)`, `Hex2/Hex1`, `ne=6`, `sim_time == t_steps == 2.5e-3`:
`output` agrees to **1.2e-15 relative** (round-off, different summation order) and the border
contours are **bit-identical** (`max|Δ| = 0.0`, both `(2, 86)`). The ported path reproduces the
maintained path's physics.

**Two plan premises were wrong**, both now corrected in the plan file:
- The 1-arg bug is **not** confined to the `"velocity"` branch — the primary `"force"` branch has
  it too (was `run_example.jl:255`).
- There is **no** "cache-construction pattern used elsewhere in the same function" to mirror. One
  had to be introduced, copying `simulate`'s placement (build once, before the control dispatch —
  consistent with `BasisFunctionCache` never being rebuilt inside a loop).

Full suite reproduces baseline.

## 2026-09-17: Refactor PAUSED after Phase 3 — status snapshot

Stopping point by request. Phases 0, 1, 2a-2d, 3a-3c are complete, plus the out-of-band
`SqueezeFlow` boundary-key fix. **Phases 5-9 never started.**

**Phase 4 (file splits) was consciously SKIPPED, not forgotten.** Re-measured on 2026-09-17:
`stokes_solver.jl` 1045 lines, `PostProcess.jl` 1103 — both over the >1000 strong flag. But
`stokes_solver.jl` had been brought *under* the threshold by Phase 1 (1010) and was pushed back
over **purely by the `simulate` docstring added in Phase 2c**. Splitting a numerically-critical
solver to satisfy an advisory line count, when the overage is good content, was judged poor
value; 4.2 the plan itself already marked "optional, do last, droppable"; 4.3 was always
"no split". Note the earlier claim in this log that "Phase 4.1 is no longer warranted" is
**stale** — it was true at 1010 lines and is not now.

**Headline outcome: every exported symbol in the package has a docstring.** 39 undocumented
when Phase 2c began, 0 now.

**The verification baseline, unchanged through every phase** — stokes 3770, mesh_templates 144,
tikhonov 41, cost_functions 34, fem 9+4, stokes_optimization clean. `tikhonov` exits 1 *after*
all 41 pass; that is the known line-101 bug, not a regression. Anything touching the solver also
got a before/after numerical diff rather than just a passing suite — that is what caught nothing
and proved two changes inert (the boundary-key relabel: byte-identical; the 3c port:
bit-identical contours vs `simulate`).

**Nine deferred items** are listed in `claud_docs/CLAUDE.md` under the refactor section, each
with its reasoning in the per-phase entries above. Two pre-existing bugs (`tikhonov_test.jl:101`,
`_get_F_ext`'s branch order) remain deliberately unfixed — both need a decision, the second
because it changes ground-truth-generating force.

**Process lesson worth keeping:** commit each phase before starting the next when the file sets
overlap. Phase 2d was started while 2c was uncommitted and both touched `run_example.jl`,
`stokes_solver.jl` and `smearFEM.jl`, briefly making them unsplittable by `git add <path>`.

## 2026-09-17: Contour-cost replicate statistics for scenario B (synthetic bulk-viscosity)

**The question:** why does the results section report est/pred for height and MOSD but not for
contour matching, on scenario B?

**Answer — two independent layers, both omissions, neither a data limitation.**

**Layer 1: `post_analysis_bulk` never computed the contour cost.** Those stats (`cpd_runs`,
`cpd_pred_runs`) existed *only* in `post_analysis_real`. The giveaway that it was an omission
rather than a decision: `post_analysis_bulk` already read `sim_border_pt_lst` and
`obs_border_pt_lst` and truncated both to the aligned window length — then **never referenced
them again**. It did exactly the loading the contour cost needs and dropped the result.

The data was never the blocker. Scenario B has strictly *more* contour data than the physical
tree that does report the metric:

| tree | `_est` | `_pred` | `_gt` | `cost_windows.csv` |
|---|---|---|---|---|
| syn bulk-viscosity (B) | 301 | 306 | **301** | present |
| physical | 301 | 308 | **0** | present |

The probable original rationale is that for synthetic runs MOSD (3D surface distance vs the GT
mesh) is a stronger accuracy measure, so the 2D contour distance read as redundant — the comment
at `post_analysis_bulk`'s head documents exactly that asymmetry in the other direction ("the
physical data has no ground-truth mesh ... `post_analysis_real` does not compute this"). Contour
cost was the physical path's *substitute* for MOSD. Defensible, but it left B without the metric
the paper reports for A.

**Layer 2 — bigger: `plot_results` does not call ANY `post_analysis_*` function.** The
`if viscosity_type == "constant" / elseif ... bulk_viscosity / elseif ... physical` dispatch at
the end of its `opt_method` loop is **commented out**. The driver currently regenerates only
`replot`'s per-experiment figures — no cross-replicate statistics for any scenario or quantity,
height and MOSD included. Left commented (re-enabling changes what every run writes); regenerate
by calling `post_analysis_bulk(gt, res, ["post_analysis_global"]; method="gn")` directly.

**What was added** (mirrors `post_analysis_real` exactly, reusing `_windowed_contour_cost`):
accumulators `cpd_runs`/`cpd_pred_runs`/`cpd_pred_segs`; the two reads it was missing
(`cost_init` from `window_data/cost_windows.csv`, and `2D_border_points_pred` guarded by `isdir`
like `pred_h`); and a third reporting block alongside height and MOSD, writing
`contour_cost_statistics.csv` plus the A/B/C folders and a `closest_point_distance_error` figure
with the same `ylims=(0,3)` clamp `replot` uses.

**Verified:** `n = 6` on both new quantities (no `single_window` leak), and existing numbers
unchanged — height still reproduces the paper's Group B entry exactly (`h_est/h_gt`
0.9960±0.0008, `h_pred/h_gt` 0.9949±0.0008). Purely additive.

**New numbers** (approach B, the one to quote): contour cost est/initial **0.386, 95% CI
[0.257, 0.515]**; contour cost pred/initial **0.473, 95% CI [0.325, 0.621]**. A/B/C diverge
exactly as predicted for `..._pred`: C's SD is 1.74x B's (0.245 vs 0.141), the window-restart
sawtooth again; on `..._est` the ratio is only 1.10x.

**Watch when writing up:** at `n = 6` the 95% CI half-width is *wider* than ±1 SD
(`t₅/√6 ≈ 1.05`), visible in every row. Only at `n ≳ 7` does CI become the tighter of the two.

**Still open:** six "Incomplete run, skipping" warnings appear during the run. They cost no
replicates (`n` is 6 everywhere) but were not traced — check before quoting these numbers.

## 2026-09-18: Refactor resumed — Phase 5a-5d done (experiments/)

Phase 4 stays skipped. Resumed at Phase 5. **There is no test suite for `experiments/`** — the
plan's verification is "re-run each script and compare to a pre-change run", which turned out to
be infeasible for most files (see each item).

**5a — `visualization_example.jl` + `single_simulation.jl`.** Removed the dead `ndim = 3`
(`initialize_mesh` has no such parameter); 4 docstrings; script gates on both. Nothing in the
repo `include`s either file, so no side effect was lost. Verified both: gate suppresses `main()`
on `include` and fires as a script; `single_simulation.jl` ran end-to-end (~10.6 min, exit 0).

**5b — `aruco_accuracy.jl`.** Fixed `Plots.Plots.xlims!` at 2 sites and added the gate. The
qualifier fix is **provably inert** — `Plots.Plots === Plots` and `Plots.Plots.xlims! ===
Plots.xlims!`, since a module binds its own name.

> **`aruco_accuracy.jl` cannot run to completion, and could not before this change either.**
> Two pre-existing path defects: (1) it reads
> `ground_truth/physical_data/aruco_calibration/aruco_calibration.hdf5`, but the file is at
> `ground_truth/aruco_calibration/...` — the `physical_data/` segment is spurious; (2) `savepath`
> is `experiments/physical_data/aruco_calibration/plots`, which does not exist, and unlike
> sibling scripts this one never calls `set_file(savepath)` before `Plots.savefig`. The real
> output dir in the data tree is `experiments/aruco_calibration/plots` — no `physical_data`. Both
> paths carry the same spurious segment, so the data layout was likely reorganised without
> updating this script. Left unfixed (behaviour-changing, outside 5b).

**5c — `forward_convergence.jl`.** 7 docstrings (file now has zero undocumented functions);
renamed `time_intergration_convergence_analysis` → `time_integration_...` (only 2 occurrences,
both in-file, no external callers); 6 `println` → 2 `@info` (per-sweep progress) + 4 `@debug`
(per-iteration error dumps), matching the file's own 9 existing logging calls. Verified by
re-running **both** plot functions against the existing convergence datasets — all 8 PDFs
regenerated. The two sweep functions were not run: `main()` is an 8-step time sweep plus a 4-size
mesh sweep, each step a full simulation, i.e. hours.

**5d — `inverse_convergence.jl` (2 docstrings), `ParallelExecution.jl` (4), `slice_and_plot.jl`
(3 + 12 `println` + `run_test` → `main`).** All three now have zero undocumented public
functions.
- `ParallelExecution.jl` really is **4**, as the plan said — an initial scan said 5 because
  `my_task` appears inside the *module docstring's usage example*, not as a real function.
  Smoke-tested live rather than by re-running consumers: 8/8 trivial tasks across 2 workers, with
  the spinner and per-worker timings rendering, which exercises all of `get_available_memory_mb`,
  `allocate_workers`, `print_progress_spinner`, `display_batch_info`, `run_parallel_tasks`;
  `format_task_info` tested separately including its path truncation. Re-running the real
  consumers is not feasible (multi-hour batches at 8 GB/task) and the changes are docstrings only.
- `slice_and_plot.jl`'s 11 diagnostics became `@debug`; the final "done" marker became `@info`,
  being genuinely user-facing when run as a script.
- **Plan correction:** 5d says "add script gate" for `slice_and_plot.jl` — it already had one.
  Only its call site needed updating to `main(plotflag)`.

**MISTAKE, worth not repeating:** I ran `include("experiments/stokes/write_gt_stokes.jl")` to
check it still loaded. That file **has no script gate**, so including it immediately started a
full parallel ground-truth batch. Killed it; it had written **zero files**, only three empty dirs
(`ground_truth/sim_data/Stokes/force/bulk_viscosity/Hex2_16.0/cube/4`), which were removed. The
pre-existing `Hex2_16` tree was untouched (35,586 files intact). This empirically confirms the
warning in Phase 5f — never `include` an `experiments/` file before checking it has a gate. Files
still lacking one: `write_gt_stokes.jl`, `verify_meshes.jl`.

## 2026-09-18: Phase 5e — shared plot style, and why PLOT_CONFIG was NOT extracted

`experiments/plot_style.jl` now holds the four `def_orange`/`def_blue`/`def_red`/`def_green`
colours; `aruco_accuracy.jl` and `test_opt_stokes.jl` include it instead of defining them.
`paper_1.jl` inherits them through its `paper_1 -> test_opt_stokes -> plot_style` chain.

**`PLOT_CONFIG` was deliberately left duplicated in all four files.** The plan called the four
dicts "near-identical" and told me to verify before extracting. They are **not** identical —
three distinct variants:

| file | font | h x w | margins l/r/t |
|---|---|---|---|
| `inverse_convergence.jl` | 11 | 320 x 480 | 1pt / 1pt / -3pt |
| `aruco_accuracy.jl` | 11 | 320 x 480 | 1pt / 1pt / -3pt |
| `forward_convergence.jl` | **12** | **360 x 330** | 1pt / 1pt / -3pt |
| `test_opt_stokes.jl` | **12** | 320 x 480 | 1pt / **5pt** / **1pt** |

A shared dict would silently resize figures and change font sizes in at least two scripts —
including the one that produces the paper's figures. `plot_style.jl` carries a comment saying so,
so a future pass doesn't "finish the job".

Note the colours are `global`, not `const`, so re-including `plot_style.jl` in one session (e.g.
`paper_1.jl` pulling it in transitively) reassigns without a redefinition warning.

Verified: colours compare equal to the pre-extraction literals from `test_opt_stokes.jl` and
`aruco_accuracy.jl`; aruco's own `PLOT_CONFIG` still reads font 11 / 320x480; `paper_1.jl`
resolves `def_blue` correctly through the chain.

## 2026-09-18: Phase 5f — write_gt_stokes.jl documented+gated, worker error handling consolidated

Done in two commits because the second half touches `test_opt_stokes.jl`.

**Part 1 (`4bfbffb`):** `write_gt_stokes.jl` went from **zero** docstrings to zero undocumented
functions (6 added), one `println` → `@info`, and gained the script gate.

**The gate is the substantive change.** `include`ing this file used to launch a full parallel
ground-truth batch; it is now inert (3.6 s, nothing written). Nothing in the repo includes it, so
no dependent broke — but it was a live footgun for interactive use, and it is the thing that bit
this session on 2026-09-17.

Two docstrings worth their ink: `_get_F_ext`'s records that the force is **kg·mm/s², not
newtons**, and `calibrate_multipliers`' records that it is the procedure that *produced* the
table hardcoded in `_get_F_ext` — not obvious from either function alone. The `_get_F_ext`
unreachable-branch bug was left alone **by explicit decision** ("no effect at the moment"); the
docstring says nothing the mis-ordered ladder would contradict.

**Part 2: `_handle_worker_error` + `write_error_log` moved into `ParallelExecution.jl`** and
exported, deleting the divergent copy in `test_opt_stokes.jl`. The shared version keeps
`write_gt_stokes.jl`'s log-persisting behaviour, as the plan wanted, with two adjustments:
- **Signatures loosened** to `(err, idx::Int, params)` and `write_error_log(err, bt::Vector;
  params=Dict(), ...)`. The `write_gt` original annotated `err::Exception` / `params::Dict`, but
  `test_opt_stokes.jl`'s call sites pass them untyped.
- Message generalised from `"write_gt_data failed for params index $idx"` to `"Task $idx failed"`.

**This is a deliberate behaviour change for `test_opt_stokes.jl`:** its workers previously only
`@error`'d; they now also persist a log file. Checked first that its param dicts carry
`"filepath_res"` at all three parallel call sites (lines ~4871, ~4985, ~5117), so logs land in
each run's own `Results/logs` rather than dumping into the cwd via the `"."` fallback.

**Subtle risk that had to be checked:** `_handle_worker_error` calls `catch_backtrace()`. That is
dynamically scoped, so moving it into a module function still called from inside a `catch` block
works — verified by a live synthetic failure: the log was written with both params and a
stacktrace. Also verified both consumers still load, the gate still holds, and
`test_opt_stokes.jl`'s colours and `PLOT_CONFIG` (font 12, right margin 5pt = 1.7639mm) are
unchanged.

## 2026-09-19: SIMPLEX QUADRATURE BUG — every Tet/Tri integral was wrong (FIXED)

Found via Phase 5g. `experiments/verify_meshes.jl` compares a Tet and a Hex solve of the same
squeeze-flow problem; it had not run since the geometry-dispatch refactor (stale `def_problem`
signature) and was never in `runtests.jl`. Once ported it failed 20/20: Tet and Hex disagreed by
**87%** (cumulative compression −2.13 mm vs −13.10 mm).

**Root cause — one bug, in the quadrature weights only.**

Numerical integration pulls back to a *master element*: `∫_Ω f dx ≈ Σ w_k f(ξ_k)|det J|` with
`J = coords · ∇λ`. So **`∇λ` defines the master element**, and the weights must satisfy
`Σ w_k = that master element's measure`. The two are halves of one convention.

`_basis_tri`/`_basis_tet` take `∇λ` from **Eq. (3.39) of the 4774 FEM lecture notes**
(`/media/soshala/Shared Vol/SMEAR-Notes/Literature/Lecture Notes/FEM lecture notes 4774.pdf`,
§3.3.5), whose master element is the **equilateral / regular unit-edge simplex**:

    λ1 = 1/3 − ξ1 − (1/√3)ξ2 ;  λ2 = 1/3 + ξ1 − (1/√3)ξ2 ;  λ3 = 1/3 + (2/√3)ξ2

Differentiating gives `[-1,-1/√3; 1,-1/√3; 0,2/√3]` — the code's matrix verbatim. **`∇λ` is
correct; do not "fix" it to unit-simplex gradients.** (An earlier diagnosis in this session
wrongly called it a bug — the lecture-notes cross-check overturned that.)

The **weights**, however, came from published tables built on the *right/unit* simplex, and from
two different conventions at that:

| rule | was | should be | error |
|---|---|---|---|
| tri n=1, n=3 | 0.5 (unit-simplex area) | √3/4 = 0.433013 | 2/√3 = 1.1547x |
| tri n=4 | 1.0 (fractional table) | √3/4 | 2.3094x |
| tet n=1,4,5 | 1.0 (fractional table) | 1/(6√2) = 0.117851 | 6√2 = 8.4853x |

Over-integrated stiffness ⇒ far too stiff ⇒ under-deformation. The signature was unmistakable:
Tet crept at a near-constant −0.11 mm/step while Hex physically decelerated −0.95 → −0.79 as the
gap closed.

**Fix:** `TRI_MASTER_AREA = √3/4` and `TET_MASTER_VOLUME = 1/(6√2)` in `src/fem/fem.jl`, with all
six rules written as normalised fractions times their measure. `∇λ` untouched.

**Results:** unit tet volume 1.41421 → 0.1666667 (exact 1/6); unit tri area 0.57735 → 0.5;
cylinder-mesh volume ratio 8.48528 → 1.00000. Tet vs Hex solve: 87% apart → **1.8e-6 apart**
(−13.101336 vs −13.101359). `verify_meshes_test` 0/20 → **20/20**.

**No published number is affected.** Hex/Quad/Line are tensor-product rules on [-1,1] and were
always exact (volume ratio 1.0000); `_surface_element_shape` maps Hex → Quad, so a Hex pipeline
never touches `_basis_tri` either. Every paper result tree is `Hex2_16`/`Hex2_6`. Tet converged
*onto* the Hex answer, which independently confirms Hex was right all along.

**Why it survived so long:** only simplices have the master-element ambiguity (tensor-product
elements don't); `mesh_templates_test.jl` exercises Tet *mesh generation* (node counts,
connectivity) but never Tet *integration*, so its 144 tests passed throughout; and the one test
that would have caught it lived in `experiments/`, was absent from `runtests.jl`, and could not
even run. **`test/verify_meshes_test.jl` is now registered in `runtests.jl`** (~6 min) as the
standing guard.

## 2026-09-19: verify_meshes_test — β sweep, and why the tolerance is relative

Extended `test/verify_meshes_test.jl` from one β to three, because **β selects which half of the
simplex quadrature is exercised**: the system is `A = η·A_bar + β·b`, where `A_bar` is the
*volume* integral (Tet rule) and `b` the *surface* integral (Tri rule for Tet meshes, Quad for
Hex). At the original β=1e-5 the surface term vanishes, so the test verified only the Tet volume
fix and left the Tri fix completely unchecked.

Cases are `((1e-5, 20.0), (100.0, 5.0), (1e4, 5.0))` — (β, sim_time). The β→0 case keeps 20 s as
the only one able to catch Tet and Hex *drifting apart* over time; the surface-dominated cases
run 5 s because a weight mismatch shows ~88% error on step one. 33 tests, ~8m30s.

**The added cases initially failed** (0.1-0.4% disagreement vs `atol=1e-3`). A convergence sweep
settled whether that was a second bug or discretization error:

| ne | β=1e-5 | β=100 | β=1e4 |
|---|---|---|---|
| 4 | 0.0023% | 0.7891% | 0.7994% |
| 6 | 0.0040% | 0.4091% | 0.3596% |
| 8 | 0.0001% | 0.3381% | 0.4449% |
| 10 | 0.0002% | 0.2319% | 0.3148% |

**Discretization error, not a bug.** It converges under refinement (~0.8% → ~0.3%), vanishes
entirely at β→0 (the volume rule is exact to 1e-4 %), and scales with how much the surface term
matters — i.e. Tri-vs-Quad faceting of the curved side wall. Non-monotonicity at β=1e4/ne=8 is
Gmsh retopologising; unstructured meshes do not nest, so convergence is noisy.

**Tolerance is now relative** (`VM_RTOL = 0.02`), not `atol=1e-3`. Step magnitudes vary with β
(≈0.95 mm at β→0, ≈0.65 mm at β=1e4) so a fixed absolute tolerance means a different strictness
per case. 2% clears the measured 0.79% worst case by ~2.5x and sits ~40x below the 87% signal of
a quadrature mismatch. Hex is the denominator.

**`ne` is elements-per-direction, NOT a characteristic length** — bigger is finer. At ne=11 the
Hex mesh has 1430 elements ≈ 11³. I initially read it backwards and ran a "refinement" sweep that
was actually coarsening. Production trees: ground truth `Hex2_16` (ne=16), fits `Hex2_6` (ne=6).

**Hex is more accurate than Tet per DOF**, not less — worth knowing when reading a disagreement.
A Q2 hex spans the triquadratic space (27 nodes); a P2 tet only the complete quadratic (10). At
ne=11 Tet needs 7896 elements against Hex's 1430 for the same ~12.5k velocity DOFs, and tets are
the more locking-prone choice in the near-incompressible Stokes limit. Tets' advantage is
meshing complex geometry, which a cylinder does not need. So when Tet and Hex disagree, expect
the Tet side to carry most of the error — consistent with Tet converging *onto* Hex once the
quadrature was fixed.

## 2026-09-19: Phase 5h — test/ consistency, and two things it exposed

Three changes, necessarily one unit: registering `tikhonov_test.jl` without fixing line 101 turns
the suite red, and fixing line 101 without registering it leaves the next such bug equally
invisible.

1. **`stokes_optimization_test.jl` wrapped in a `@testset`** — the only test file without one.
   Verified whitespace-and-wrapper only: exactly one `@testset` line and one `end` added, nothing
   removed (compared as a `Counter` of stripped lines).
2. **`tikhonov_test.jl:101` `round(...)` → `round.(...)`** — the bug deferred since the start of
   the refactor. Some entries store η/β as a 1-element `Vector`, others as a scalar, so the
   scalar method threw `MethodError` on the vector ones. **Now exits 0.**
3. **`tikhonov_test.jl` registered in `runtests.jl`.** All seven test files are now registered
   (audited, not assumed).

**Two discoveries worth keeping:**

- **`stokes_optimization_test.jl` runs 93,868 assertions**, not "a few bare `@test`s". Without an
  enclosing `@testset` they reported into the ambient context, so the baseline recorded the file
  as "clean" and uncounted — **95% of the suite was unmonitored**. It now reports a number that a
  regression would move.
- **The suite no longer has an expected-failure file.** The old "tikhonov exits 1, that is
  expected" caveat is deleted from `CLAUDE.md`. A permanently-red file is precisely what hides a
  genuine new failure.

**`runtests.jl` still cannot execute here**: `ArgumentError: Package Aqua not found in current
path` — the documented `Pkg.test()` breakage (`test/Project.toml` lists only `smearFEM` in
`[deps]`, shadowing the root `[targets]`). So registering `tikhonov_test.jl` fixes the bookkeeping
but is **latent** until that is resolved; verification remains per-file, as `CLAUDE.md` prescribes.
The proper fix — adding `Test`/`Aqua` to `test/Project.toml` — was deliberately **not** bundled
into 5h: it is a dependency-resolution change that can affect precompilation, and it deserves its
own commit. Do **not** use `push!(LOAD_PATH, "@stdlib")`; it has poisoned the precompile cache
before.

**New baseline (all seven exit 0):** stokes 3770, stokes_optimization 93868, mesh_templates 144,
tikhonov 41, cost_functions 34, fem 9+4, verify_meshes 33.

## 2026-09-19: Phase 6 — optimiser consolidation (partial, by measurement)

The plan gated this as "the highest-risk `src/` change — the code that produced the paper's
fits", with "skip and record why" as a legitimate outcome. A longest-common-run analysis
(`difflib.SequenceMatcher` over comment/blank-stripped, whitespace-normalised lines) split the
answer rather than confirming either premise:

| pair | identical | longest contiguous run |
|---|---|---|
| GN vs Tikhonov | **62%** | **33 lines** |
| GN vs LM | 38% | 9 |
| Tikhonov vs LM | 44% | 9 |

**LM is legitimately different — left alone.** Trust-region gain ratio `ρ` for step acceptance,
rejects and restores `θ_prev`, adaptive `λ`, no line search. Its stats dict also returns only 6
keys (**no `simBorderPtsList`**), so routing it through the shared builder would have silently
added a key it never returned. Only the 3-line accumulation loop was shared with it.

**GN and Tikhonov genuinely duplicated** three blocks: a 16-line setup, a 17-line iteration
header, and a 33-line convergence/termination/stats block.

**Extracted (4 helpers):** `_accumulate_derivatives!` (all three), `_warn_ill_conditioned`,
`_check_stop`, `_fit_stats`.

**Deliberately NOT extracted:** the step computation and step acceptance — the parts that
actually differ. `p = t∂2d \\ t∂d` (GN, with a damped `+reg·I` fallback) vs
`(t∂2d/len_d + ∇²R) \\ (t∂d/len_d + ∇R)` (Tikhonov — note it normalises the data term by `len_d`
and GN does not) vs LM's adaptive damping. Forcing these together would be a correctness
regression dressed as cleanup.

**One case needed care:** GN's ill-conditioning `if` does not merely warn — it computes the
damped step inside it, while Tikhonov's only warns. So `_warn_ill_conditioned` returns
`(condition_number, ill_conditioned)` and each function keeps its own branch; output ordering
and arithmetic are unchanged.

**Verification — bit-for-bit.** Baseline captured *before* any edit from a 7-configuration
harness (`:gn`, five `:gn_tikhonov` variants incl. λ_scale=0/1/1e4/vector and a prior θ_p, and
`:lm`) mirroring `tikhonov_test.jl`'s setup (ne=6, 5 steps), dumping η, β, iteration count and
every `cost_list` entry at `%.17g`. Post-refactor output is **identical, all 47 lines**. Full
suite also unchanged: stokes 3770, stokes_optimization 93868, mesh_templates 144, tikhonov 41,
cost_functions 34, fem 9+4, verify_meshes 33.

Net 848 → 892 lines: four documented helpers cost more than the duplication they remove. **The
gain is not line count** — it is that the convergence thresholds now live in one place instead of
two, the same single-point-of-truth argument as the `window_dirs` skip-list.

## 2026-09-19: Phase 7 — test_opt_stokes.jl dedup (7.1, 7.2, 7.5 done; 7.3 skipped; 7.4 open)

**7.1 — `_window_subdirs(exp_path)` extracted (the plan's highest-value item).** The identical
skip-list plus its 8-line explanatory comment was copy-pasted in `replot`, `post_analysis_bulk`
and `post_analysis_real`; editing two of three copies would silently reintroduce the
`single_window` NaN bug. Now one helper, pure extraction, filtering logic untouched.
**Verified on real data**: across 16 experiment directories in the synthetic-bulk and physical
trees, `_window_subdirs` returns byte-identical lists to the inline filter — 0 mismatches.

> **Finding, not acted on:** `replot` still has two branches at (now) ~2948 and ~3125 testing
> `window_dir == "single_window"` / `!= "single_window"`, **inside the very loop that filters
> `single_window` out**. They are unreachable — and were already unreachable before 7.1, since
> the old skip-list `continue`d at the top of the same loop. They clearly once plotted
> `single_window` with a special range, which suggests the 2026-09-16 skip may have been
> **over-applied to `replot`**, silently dropping those per-experiment plots. The skip is
> unambiguously right for the two `post_analysis_*` functions (it fixed the n=12 / NaN bug);
> whether `replot` should also skip it is a separate question worth answering.

**7.2 — docstring misattachment fixed.** `_npts(p) = ...` sat between
`iteration_contour_metrics`'s docstring and the function, so Julia bound that 3636-char
docstring to `_npts`. Moved `_npts` above it; verified the docstring now binds to
`iteration_contour_metrics`.

**7.3 — SKIPPED, premise wrong.** The plan calls these "26 leftover `printstyled` calls". They
are not leftovers. The file already uses **57 `@info`, 51 `@debug`, 41 `@warn` and zero
`println`**, and all three tiers appear within three consecutive lines — `printstyled` is a
deliberate always-visible tier, colour-coded by nesting depth (green=directory, blue=element
size, cyan=sim time, magenta=experiment) for watching multi-hour batches. 13 of the 26 are in
`optimize` itself. Converting them makes batch progress invisible by default and destroys the
hierarchy — exactly what `CLAUDE.md`'s `println` exception exists to prevent.

**7.5 — partial; the plan over-reached.** Removed the two commented-out alternate `PLOT_CONFIG`
blocks. **Only 3 of the 9 globals could be `const`-ified** — `end_obs_win`,
`eta_beta_plt_width`, `eta_beta_plt_height`. The other six (`fs`, `plt_height`, `plt_width`,
`plt_lft_margin`, `plt_right_margin`, `plt_top_margin`) are **reassigned by `set_plot_config`**
at runtime; `const` would hard-error. The plan said "const-ify the load-once globals (lines
27-65)" without distinguishing them.

**7.4 — NOT acted on, needs a decision** (the plan says "ask before acting"). This is deferred
item 11: `plot_results`'s entire per-viscosity `post_analysis_*` dispatch is commented out, and
`main()`'s docstring claims it runs "the full pipeline" while its batch calls are commented out
and only `plot_results()` runs.

## 2026-09-19: single_window restored for replot (only), with a completeness guard

Follow-up to the Phase 7.1 finding. `_window_subdirs` gained `keep_single_window::Bool=false`,
and **only `replot` passes `true`**.

**The asymmetry is deliberate and must stay.** In `post_analysis_bulk`/`post_analysis_real` the
exclusion is load-bearing — it is what fixed the n=12 / `h_pred/h_gt` NaN bug, and re-including
`single_window` there would reintroduce it. `replot` is different: it draws one experiment at a
time and already had the machinery — the `window_dir == "single_window"` branches plotting over
the first window's range, and `pred_is_windowed` detecting the differing `pred_h.csv` length.
Those branches were not vestigial, they were **orphaned by an over-broad filter** applied to all
three call sites on 2026-09-16.

**Empty-directory safety:** `replot`'s window loop now calls `_skip_incomplete(win_exp_path) &&
continue` — reusing the guard the two post-analysis loops already had rather than adding a new
one. Without it an empty `single_window` throws on `read_json(datapath(…,
"experiment_parameters"))`, the loop's first read. Verified both branches: empty dir → listed
but skipped; populated dir → processed.

**Where `single_window` actually lives:** 12 directories, **all under `Carreau`** (plus a
separate `syn_data/single_window_test` tree, unexamined). None in the Stokes bulk tree — which
is why the Stokes bulk post-analysis reported a clean n=6 earlier while Carreau reported n=5.

**End-to-end verification:** `replot` on Carreau replicate 3 exits 0, processes *both* windows at
501 time steps each, and writes PDFs into both `multi_window/plots` and `single_window/plots`.
Before this change the `single_window` set was never regenerated.

## 2026-09-19: `Pkg.test()` FIXED, and 7.4 resolved (batches re-enabled)

**`Pkg.test()` works.** It had been documented as broken for the whole refactor. The recorded
fix ("add `Test`/`Aqua` to `test/Project.toml`") was not quite right — there were **two** faults:

1. **Duplicate mechanisms.** The root `Project.toml` already had a complete `[extras]`/`[targets]`
   block. But a `test/Project.toml` containing only `smearFEM` also existed, and when that file
   is present Pkg uses it *instead of* `[targets]`. Populating it would mean maintaining the dep
   list twice; **removing** it (and its `Manifest.toml`) lets the already-correct root block work.
2. **`[targets]` was itself incomplete.** `cost_functions_test.jl` does `using ... Random ...`,
   and `Random` was in neither `[deps]` nor `[extras]`. Even after removing the shadow,
   resolution would still have failed — later and less obviously. Added to both.

Verified: `Pkg.test()` resolves a temp env with **Aqua v0.8.18**, precompiles it, and runs the
suite. **Do not reintroduce `test/Project.toml`**, and do not use
`push!(LOAD_PATH, "@stdlib")` — that has poisoned the precompile cache before.

**7.4 resolved — user chose "re-enable everything".** In `plot_results` the
`post_analysis_const/bulk/real` dispatch is uncommented, so cross-replicate statistics now
regenerate alongside the plots. In `main()` the three `optimize_*(false)` batches are
uncommented and `data_types` is actually passed to `plot_results` (it had been assigned and
ignored); the dead `avoid_dirs` local is gone. Docstring updated to state plainly that this is a
**multi-hour** run and that `plot_results()` alone regenerates figures without refitting.

> **The script gate is now load-bearing.** `paper_1.jl` does
> `include(".../test_opt_stokes.jl")`. With `optimize_*` live, an ungated include would launch
> hours of refitting. Verified: include takes 9.4 s, `main` is defined but unrun. **Never remove
> that gate.**

**Tet/Tri data sweep — nothing to discard.** Every result tree in `SMEAR-DataFiles` is Hex
(`Hex2_6` x13, `Hex2_16` x9, plus single `Hex2_8/4/2/10` and `Hex_2`); **zero Tet trees exist**.
The 21 paths matching "tet|tri" are false positives — `me-tri-cs` and `mesh_metric_bias`. So the
simplex quadrature bug corrupted no stored data at all: it was purely latent, consistent with
Tet having been unreachable since the geometry-dispatch refactor.

## 2026-09-20: `Pkg.test()` passes — 97,912 tests, 22m32s

First successful `Pkg.test()` in this repo. See the 2026-09-19 entry for the two faults. The
number is worth recording because ~96% of it is `stokes_optimization_test.jl`'s 93,868
assertions, which ran uncounted until Phase 5h gave that file a `@testset`.

Breakdown: stokes 3770 · stokes_optimization 93868 · mesh_templates 144 · tikhonov 41 ·
cost_functions 34 · fem 9+4 · verify_meshes 33 · qa (Aqua) the remainder.

## 2026-09-20: Q4 (`isatty`) resolved — it was hiding a live defect

The plan's open question asked whether the `const isatty = Base.isatty` alias in `smearFEM.jl`
was worth replacing by fixing "the caller that references `smearFEM.isatty` incorrectly".
Investigating it answered both halves and found a real bug:

1. **There is no such caller.** Grepped `smearFEM.jl` *and* the sibling `smearPerception` —
   nothing references `smearFEM.isatty`. The alias protected code that does not exist.
2. **`isatty` does not exist in Julia 1.12 at all** — not `Base.isatty`, not `Base.Libc.isatty`,
   not bare. So `isdefined(Base, :isatty)` was **false**, the guard's `else` branch always won,
   and `const isatty = (io->false)` meant `_isatty_default()` returned `false` unconditionally.
   **`progress_guard` therefore always returned `NonTtyProgress`, even in a live terminal** — the
   redrawing progress bar had been silently dead.

**Fix:** deleted the alias, and replaced `progress_guard`'s ~30 lines of probing
(`_get_default_stdout` tried `Base.STDOUT` and `Main.STDOUT`, both Julia-0.6 names that do not
exist; `_isatty_default` wrapped the call in try/catch; both were nested *inside* `progress_guard`
so they were redefined on every call) with the one check that is correct on current Julia:

```julia
if stdout isa Base.TTY
```

**Verified in both contexts**, which is the only way to test this: piped → `NonTtyProgress`
(unchanged); real terminal (forced with `script -qec`) → `ProgressMeter.Progress` (**restored**).
Zero `isatty` references remain in `src/`.

**Lesson worth keeping:** a `isdefined(Base, :x)` feature-guard with a silent fallback hides the
failure permanently once upstream drops `x`. This one degraded a UI feature, but the same shape
around a numerical fallback would be far worse. There are no other `isdefined(Base, ...)` guards
in `src/` — checked.

## 2026-09-20: deferred-list pass — Q2, #10, #7, #8, #9 closed

**Q2 (`figures.jl` script gate) — resolved by measurement, not judgement.** The plan guessed the
partial gate was deliberate. In fact **exactly 1 of 19 public functions refits**:
`run_lambda_sweep`, which calls `optimize`. Everything else in all four parts — including Part
2's own `lambda_sweep_table`/`plot_lambda_sweep`/`plot_lambda_contour_error` — only reads. So
"Part 1 only" had no rationale. The gate now runs 8 entry points across all four parts and
excludes `run_lambda_sweep` alone, with the rule stated inline.
- It had to **move to the end of the file**: at line 533 it could not reference
  `plot_beta_bounding` (defined at 1539). That constraint is probably why it only ever
  covered Part 1.
- Each call is wrapped so one missing dataset warns and the rest still regenerate.
- Verified end-to-end: all 8 parts ran, zero failures.

**#10 `aruco_accuracy.jl` — was broken three layers deep, now runs.** My deferred note was
wrong in one detail and incomplete in another:
1. Read path had a spurious `physical_data/` segment. (as recorded)
2. `set_file` **was** present — but written `set_file(dirname(savepath))`, creating the *parent*
   of `plots` rather than `plots`. My note said it was absent; it was misapplied.
3. **A third bug the first two hid:** once it could read the file, `BoundsError: 115-element
   Vector at index [1:162]`. The script builds one time axis from `stage_separation/frame_index`
   (162 samples) and used it for the plate figure too — but `plate_separation` is a separate,
   shorter acquisition (115 samples; `mold_b` is 146 again) with no `frame_index` of its own.
   It now derives its own axis as index÷fps, the same way the stage axis is built.
   **Caveat:** that assumes the 115 plate samples are consecutive frames at the recorded fps.
   The file does not state it. The reported scalars do not depend on it; only the x-axis does.
Runs clean now, first time since 2026-07-14: stage separation 107.663 ± 0.138 mm (n=162),
drift −0.109 mm, plate scale 1.018 (residual 0.18 mm vs 10.0 mm).

**#7 dead code — the list said 3, a fresh scan found 7.** Deleted 5: `_normalize_field`,
`_truncate_colormap` (`plotting.jl`), `_read_sparse_mat` (`io.jl`), `closest_point_contour`,
`fit_curve_2D` (`PostProcess.jl`). **Restored 2 on the user's instruction — "1D 2D 3D meshing
functions should not be touched":** `_meshgrid_ring` (`meshes/mesh_structured.jl`) and
`_second_derivative_pairs` (`fem.jl`, branches on `ndim == 1,2,3`). **`src/fem/meshes/` is
off-limits for dead-code removal** — unused there means held for dimensional support, not
deletable.

**#1 `SqueezeFlow.model` deleted** — and with it the documented guideline exception. It was the
only abstractly-typed field, so it forced "defaults for every field" to carry a caveat.
`SqueezeFlow()` now constructs bare with all 8 fields defaulted.

**#8 comments — 181 restating comments removed** from `run_example.jl` (91),
`stokes_solver.jl` (84), `io.jl` (5), `plotting.jl` (1). Verified **comment-only: zero code-text
changes** in every file. `src/fem/meshes/` excluded. **Deliberately kept** anything the code does
not say: the `# Jacobian matrix [dx/dxi ...]` layout notes, `# *2 because we have two surfaces`,
`# second derivatives vanish for linear basis`, `# slip boundary condition only affects the
tangential components`, the `(for ndof = 1) / ...` dual-meaning notes, `# x direction` on loop
indices, and the `# Meshes.jl` export-list file tags.

**#9 `μu_btm`/`μu_side` — now configurable, not deleted** (user's call). They are keyword
arguments on `simulate` and `stokes_single_step_force`, defaulting to `0.0`, documented on
`simulate`. Verified: defaults are **bit-identical** to the old hardcoded zeros, and setting
`μu_btm=0.05` genuinely changes the solution (step 1: −0.022295 → −0.021670, less net
compression from a driven bottom plate — the right direction). This also gives the
boundary-key fix teeth: `q_d[:bottom]`/`q_d[:border]` now matter because they can be non-zero.

**Deferred list: 13 → 9 open.** Closed: #1, #7 (partly — 5 of 7, 2 restored by instruction),
#8, #9, #10, #11. Still open: #2 (`q_d` dense round-trip), #3 (`get_c_mat` docstring), #4
(`plot_covariance!` export comment), #5 (`simulate` output empty under velocity control),
#6 (`simulate_single_tstep_stokes` bypasses `get_mesh_dir()`), #12 (`const_vel` stale name),
#13 (six untraced "Incomplete run" warnings).

## 2026-09-20: estimate/predict extraction — the paper's method moved into src/

User's architecture call, verified before acting: `estimate` and `predict` are **core library**
functions (every symbol they need — `fit_model`, `update_model!`, `reset_model!`, `simulate`,
`get_height`, and the three types — is already exported by `smearFEM`; neither touches
`datapath`, `exp_params` or anything paper-local), while orchestration and I/O stay in
`paper_1/`. `set_time_window` moved too, as Eq. 20 of the same formulation.

**New: `src/optimization/window_estimation.jl`** — the complete moving-window formulation of
Section 3:
- `set_time_window` — Eq. 20, `w_k = [t_k, t_{k+1}]`, moved verbatim (109 lines, was
  self-contained).
- `estimate_window` — Eq. 25, `θ_k* = argmin_θ d_k(θ)`, **plus step (iii.c)**: it overwrites
  `θ` in place so it seeds the next window (`θ̂_{k+1} ← θ_k*`) and calls `update_model!` so
  `M_{k+1,0} ← M_{k,N_k}`. Both mutations are documented — a caller needing the pre-fit state
  must copy first.
- `predict_window` — the k → k+1 prediction under a fixed preceding estimate; `reset_model!`
  first so it starts from the committed reference rather than wherever the last prediction left
  the mesh.
All three exported.

**Duplication removed.** The k→k+1 prediction previously existed **twice** — inline in
`optimize`'s multi-window loop and again in `run_window_predictions` (which `predict` calls from
disk). Both now call `predict_window`. That was the same single-point-of-truth hazard as the
`window_dirs` skip-list.

`pipeline.jl`: 5020 → 4911 lines.

**NOT done, deliberately:** splitting `optimize` (776 lines) into per-concern functions. That is
the riskier half and needs baselines for the `constant` and `physical` branches too — the
baseline below covers only `synthetic`/`bulk_viscosity`/`multi_window`/`:gn`.

**Verification approach — and why the first attempt was abandoned.** `optimize` at the
production `ne=6` ran **1h25m without finishing** (~108 full 50-step simulations; GN iterations
plus every line-search evaluation each cost one). Two baselines plus a post-change run would
have been 6+ hours, so that run was killed and the harness re-cut at **`ne=4`** — same branch,
same loop, ~34 min per run. Both runs therefore sit on the *same* `src/`, which also removes the
staleness problem the original before/after pair existed to solve.

Before-run: 2014 s, **4553 CSVs**. The after-run is compared file-by-file; see the next entry
for the result.

## 2026-09-20: extraction VERIFIED — 4553/4553 CSVs byte-identical

Before-run 2014 s, after-run 2002 s, **4553 CSVs each, all byte-identical**. The only differing
file is `logs/time_log.txt`, and only in its wall-clock fields and the `filepath_res` path that
the harness varies deliberately. `set_time_window`/`estimate_window`/`predict_window` are
therefore a pure extraction on the `synthetic`/`bulk_viscosity`/`multi_window`/`:gn` path.

**The verification earned its cost — the first after-run failed.** `predict_window` was typed
`θ_prev::AbstractVector{<:Real}` but both call sites pass a tuple, and `Tuple` is not an
`AbstractVector`; the run died 535 s in with a `MethodError`. Signature widened to
`Union{AbstractVector{<:Real},Tuple{Real,Real}}`.

Note the bug was in the *prediction* path — the one place two implementations were being merged
into one, i.e. exactly where a mistake was most likely.

**Process lesson:** run a ~1 min smoke test exercising the new functions with the *real*
argument types before committing 34 min to a full run. A cheap gate first would have caught this
immediately. The smoke test now lives at `scratchpad/smoke_ep.jl` in shape: build a tiny
`def_problem`, call `estimate_window` (checks it recovers η=β=100 and mutates `θ` in place per
step iii.c), then `predict_window` with both a tuple and a vector.

**Coverage caveat, unchanged:** this covers one path only. The `constant`-viscosity branch,
`single_window` mode and `physical` data are NOT exercised by this baseline. They share
`optimize`'s setup code, which the extraction did not touch — but if the deferred per-concern
split of `optimize` ever happens, it needs its own baselines for those branches.

## 2026-10-12: publish docs pass — README / CONTRIBUTING / coding_guide brought in line with code

Audit of README.md, CONTRIBUTING.md and docs/coding_guide.md against the actual source. Fixed:
- README: `def_problem` exists only for `Cylinder`/`Cuboid` (README implied all 5 geometries);
  clamp bounds are `[1e-12, 1e12]` for both η and β (README and `_val_check`'s own docstring
  said `[1e-3,1e5]`/`[1e-3,1e8]`); `experiments/verify_meshes.jl` → now a test; test list
  completed; layout table covers config/utils/cost_functions/window_estimation/statistics.
- **Julia floor raised `1.6.7` → `1.11`** in Project.toml + README. Stdlib compat pins
  (Dates/Logging 1.11.0, Statistics 1.11.1) already made <1.11 unresolvable; CI only tests 1.11.
- `docs/coding_guide.md` rewritten: it described a long-gone layout (`src/examples/`,
  `Meshgrid1D/2D/3D`, `MeshgridCube`, `meshgrid_cube`, `init_cylinder`, `C_vol` on
  MeshgridCylinder, Makie). Mesh field is `volume_element_shape` (+ `surface_element_shape`),
  not `element_shape` — the `AbstractMeshgrid` docstring had the same error, also fixed.
- `docs/src/index.md`: `pkg "add smearFEM"` was invalid syntax and the package is unregistered.
- **35 `printstyled` calls** in `smearOptimize.jl` (33) and `run_example.jl` (2) converted —
  the 2026-09-17 `println` sweep missed them. Mapping: `[WARNING]` → `@warn`; start /
  iteration / converged / complete → `@info`; per-step `[ACCEPT]`/`[REJECT]`/`[DAMPED]` and the
  Tikhonov prior/data-precision line → `@debug`. Plus 11 commented-out `println`/`display` lines
  deleted.
- Four internal notes moved `docs/` → `claud_docs/` and untracked (`git rm --cached`):
  meshes_refactor_plan.md, refactor_element_basis_split.md, data_migration_view_level.md,
  webserever.bash.txt.

`StaticArrays` (declared, zero uses in src/, experiments/, test/) removed via `Pkg.rm`; Aqua 7/7.

## 2026-10-12: GitLab mirror synced — gitlab.tue.nl/robotics-lab/soft-robotics/smear-fem

Pushed `6ab3553` to GitLab `main` (fast-forward from `466186f`). Local checkout of that repo is
`../smear-fem`. It is NOT a plain copy of this repo — keep these GitLab-only differences on every
sync: `.gitlab-ci.yml` + `docs/make.jl` (GitLab Pages, no `deploydocs`), `config.toml.example` at
repo root, README/CONTRIBUTING/docs-index links pointing at gitlab.tue.nl, no `.github/`.
Sync method that worked: `git archive HEAD` (tracked files only — a raw `cp` leaks `paper_1.tex`,
`CLAUDE.md`, mesh caches, because GitLab's `.gitignore` lacked those rules) into a fresh clone.
README docs/API links point at GitHub Pages (`soshalaw.github.io/smearFEM.jl/dev/`) on the user's instruction (`c67aaea`) — unsure gitlab.tue.nl can host Pages. Everything else links to GitLab.

## 2026-10-12: GitLab CI OOM — measured test memory (Julia 1.12 local, Pkg.test)

GitLab runner is a Kubernetes executor ("Error in container build ... OOMKilled", exit 137).
`JULIA_NUM_PRECOMPILE_TASKS: "2"` (a9c78d7) got precompilation through; the kill moved to the
test run. Suite split via `SMEAR_TEST_GROUP=fast|slow|all` (runtests.jl, default all).
Peak RSS: fast group + coverage **1.90 GB** (15 min; 7.4 min without coverage);
slow group (`verify_meshes_test.jl`) **~4.8-5.0 GB** (12-17 min). `--heap-size-hint=2G` changes
neither — it's live mesh/solver data, not uncollected garbage, so don't retry GC tuning.
Unused test extras BenchmarkTools/Images/ProfileView dropped (4578f67 / a9c78d7).

## 2026-10-13: real-time solver — IN PROGRESS on branch `features/realtime_solver` (uncommitted)

**Goal (user):** forward prediction at **50 Hz with dt = 0.1 s**, i.e. one simulated step in
**≤ 20 ms** wall time, and the fast path must live **inside the fitting code** (`simulate`, so
`fit_model` gets it too, sensitivities included). User asked to stay in Julia after discussion
(C++ port buys nothing: UMFPACK is already C; the slowness was allocation, fixable in Julia).

**Baseline (Hex2_6 cylinder, 1885 u-nodes, 287 p-nodes, 192 el, saddle 5653 DOF, 3 views):
~600 ms/step.** Breakdown: sparse LU ~165 ms, `assemble_system_A` 155 ms (139 MB alloc/call),
`_B` 35 ms, M via hvcat + `C'*A*C` products ~11%, `_get_2D_data(GRAD=true)` ~20 ms.

**Done so far (state when stopped): ~157 ms/step before the last two tweaks, not re-timed since.**
- New file `src/solver/stokes_force_system.jl` (included after stokes_solver.jl, nothing exported):
  `ForceSystem` struct, `_force_system` (fixed sparsity pattern + element→nzval position maps),
  `_assemble_force_system!` (threaded element kernels, scatter into `nz_visc`/`nz_slip`/`nz_div`,
  `K.nzval = η nz_visc + β nz_slip + nz_div`), `_lagged_solve!` (iterative refinement against a
  stale UMFPACK factorization, `lu!` refactor when > `refactor_after` corrections or no
  convergence), `_extrapolate_guess!`.
- Key fact that makes it work: `C_uc` and `q_top` are disjoint 0/1 selection columns (checked
  numerically), so M = [T'AT T'B; B'T 0] with T = [C_uc q_top] is assembled directly; unknown
  order kept as [q_f; p; μ]. dM/dη = pattern with `nz_visc`, dM/dβ = `nz_slip`; RHS is
  [0;…;F], θ-independent, so sensitivities are M ∂sol/∂θ = −(∂M/∂θ) sol (one 2-column solve).
- **Gotcha: `BasisFunctionCache` is `NTuple{7, Any}`** — first version was 1.9 s/10 GB alloc
  from dynamic dispatch. Fixed by function barriers with concretely typed kernel args
  (`_element_matrices!`, `_slip_faces!`). Assembly now **10 ms, 129 kB** (was 359 ms).
- Viscous element block written as K[(a,i),(b,j)] = w(δij ∇Na·∇Nb + ∂jNa ∂iNb) — equals
  2wBCB' with the "standard" C = diag(1,1,1,½,½,½).
- `simulate` force branch rewired (old M/dMdη/dMdβ build, the `@debug` block that used A/B/M,
  and the per-step `lu` removed). Throws `ArgumentError` if `μu_btm`/`μu_side` ≠ 0 (all call
  sites pass 0; fast path assumes q_d = 0). Velocity-control branch untouched.
- **Aliasing gotcha:** `sol` is now reused in place, so `p = copy(p_f)'` — a bare view would
  make every stored pressure in the history equal the last step.

**Verification done:** M vs old construction max rel diff 1.1e-15 (initial) / 1.9e-15 (deformed
mesh). Full `simulate` outputs (all 10 return values, 40 steps, constant + bulk_viscosity) vs
pre-change reference: ≤ 1.2e-11 rel. **That diff was with tol=1e-12 and no extrapolation.**

**NOT yet done / next steps:**
1. Re-run the before/after diff after the last change (tol 1e-10 + `_extrapolate_guess!`).
   Scripts in the session scratchpad are gone after reboot: reference = run old `simulate`
   (git stash) on Cylinder(25,40), ne=6, η=β=100, F=-9518.61, 40 steps dt=0.1, views 0/30/60,
   `serialize` the return tuple; diff each element.
2. Run the test suite (`Pkg.test()`, or fast group `SMEAR_TEST_GROUP=fast`) — never run on this
   branch yet. `stokes_optimization_test.jl` exercises fits through the new path.
3. Remaining per-step cost: ~3 ms per UMFPACK triangular solve (memory-bound, 3.7M L+U nnz;
   UMFPACK IRSTEP already 0), 0.5 ms per K*x. With extrapolation the main solve takes 2–3
   corrections, sensitivities 3–4 (initial residual ~1e-6 vs ~1e-3 without). Refactor (`lu!`,
   146 ms, reuses symbolic) is synchronous → 50 Hz jitter; consider a background refactor or
   Krylov (GMRES, would need e.g. Krylov.jl — ask before adding deps) or MKL Pardiso.
4. `back_project` GRAD path is O(n²): `R*dxdθ[:,:,i]` is computed inside the per-node loop
   (PostProcess.jl ~line 240) — hoist it out of the j-loop; pure perf, same numbers.
5. Prediction-only flag to skip sensitivity solves + history pushes for the real-time loop.
6. Then decide whether the β/η scaling lookup or a ROM is still needed: under force control with
   fixed β/η, velocity ∝ 1/η exactly (pressure unchanged), so shape-vs-height depends only on
   κ = β/η for constant params — not yet checked numerically.
- `steps_since_factor` counts `_lagged_solve!` calls (2 per step), not time steps.
- Two `@debug` lines in `_lagged_solve!` log residual/corrections (`JULIA_DEBUG=smearFEM`).

## 2026-10-13 (later): real-time solver — 586 → ~17 ms/step, velocity control fixed

Continuation of the entry above; the step-by-step record with every measurement is
`claud_docs/realtime_solver_plan.md` (B1–B18, C1–C3). Goal restated by the user: **everything
≥ 50 Hz** (≤ 20 ms per 0.1 s step), force and velocity control, sensitivities included.

What moved the number, in order of impact:
- **One sensitivity solve, not two**: M = ηĀ + βb + D with a θ-free RHS is homogeneous, so
  η∂x/∂η + β∂x/∂β = −[u;0] (force) or [0;p;0] (velocity). Checked to 1e-14.
- **Better initial guess**: solve for the *unit* force/velocity response and scale by the
  control (immune to noisy measured force); cubic extrapolation of a unit-viscosity history.
  Offline the cubic guess is 5e-11 (constant η); in situ ~1e-8, because stored solutions
  carry tolerance-level residuals that a cubic amplifies ×15.
- **Tolerances** 1e-9 (state) / 1e-8 (sensitivities) relative residual: one correction per
  solve instead of two/three; outputs still ≤ 2.6e-9 from the pre-branch code.
- **Memory traffic, not flops**: the threaded scatter/gather was bandwidth-bound (no gain from
  threads until Int32 lists and fewer stored arrays); symmetric half of the element kernel;
  residuals as column-parallel Kᵀx (4×); ∂M/∂β on its own 8k-entry pattern (100×).
- Smaller items: cached SqueezeSystem per model (setup 766 → ~170 ms/call), hull-boundary and
  `back_project` loops in scalars (2D projection 2.6 → 0.77 ms), nodal velocity by index
  gather instead of an n×1 *sparse* matrix (3.5 → 0.13 ms).

Dead ends: GMRES (≤ 1 iteration saved), UMFPACK METIS/symmetric/"best" orderings (≤ 10 %),
η-rescaling of the stale factorization for bulk viscosity (kept — exact and cheap — but only
2.21 → 2.0 corrections; the poor *guess* was the real problem), and a **background
refactorization** on a pinned thread, built and then removed: A/B showed it slower in 40-step
runs (contention for cores/bandwidth/OpenBLAS threads; GC stalls while a long UMFPACK ccall
holds no safepoint) and only break-even in 200-step runs. Details in the plan's B1 entry.

Found along the way:
- **Velocity control was broken at HEAD** — `zero_matrix_p_q` only ever assigned in the force
  branch; never ran since v1.0.0. Reference for the rewrite came from a patched HEAD export.
- Julia 1.12 `-t auto` = 16 default + 1 interactive thread, and the main task runs on the
  interactive one, so `Threads.@spawn` work never lands there. On 1.11 or `-t N,0` a long task
  picked up by thread 1 while the sticky main task waits would block it.
- A GC requested while another thread is inside a long ccall (UMFPACK numeric factorization)
  waits for that call to return — background C work and GC-heavy Julia code do not mix.
- Benchmarks on this laptop drift ±2 ms run to run and badly when the user's browsers are
  active (load average 5.7 seen) — always compare configurations back to back.
- Top-level `BLAS.set_num_threads(Threads.nthreads())` in `src/smearFEM.jl` only runs at
  precompile time; runtime BLAS threads stay at 8. Not changed (global policy).

## 2026-10-14: Phase E — Pardiso backend (in progress, uncommitted)

User asked to reach 20 ms everywhere, "open to anything" without losing accuracy. Done: paired
z/w solve (concurrent on UMFPACK), register-accumulated element kernel (assembly 6.1 → 2.8 ms),
MKL Pardiso backend with UMFPACK fallback (factor 9.3 vs 150 ms, solve 1.3 vs 3.6 ms),
staleness-triggered refactor. Median ~9–12 ms/step, accuracy unchanged (≤ 2.5e-9). Gotchas:
Pardiso.jl Int32 indices (convert-per-call → GC spikes), explicit diagonal for symmetric MKL,
iparm(34) needed for bit-reproducible multithreaded results. Remaining: full test run, GC
pauses from `_get_2D_data` allocations, docs. Details: plan, Phase E.

## 2026-10-14: Phase E finished — MKL Pardiso, 2D projection, a multi-view gradient bug

- Pardiso backend completed: refactor at the start of the step after staleness is seen (no
  wasted round), `iparm(34)` for bit-reproducible multithreaded results, Int32 lower triangle.
- 2D projection 1.15 → 0.31 ms/step and 0.93 → 0.25 MB/step: no spline resampling unless
  WRITECONTOUR, scalar `back_project`, hull from `eachcol`, views in parallel.
- **Bug (pre-existing): `back_project` never rotated `dxdθ` by the view angle**, so 2D contour
  sensitivities were wrong (~12 % at 45°) for every view except 0°; multi-view Gauss-Newton fits
  used them. Found while rewriting `back_project`, confirmed by finite differences, fixed, test
  added. The pre-branch reference `gradList` is wrong for views 30°/60°.
- Median 9–12 ms/step, worst steady ~18–22 ms, step 1 ~22 ms; full suite 97,982 green in ~3 min.
- GC with Plots & co. loaded: young ~60–80 ms, full ~320 ms — explicit GC before the loop is
  not a fix; fewer returned outputs would be.

## 2026-10-14: what dropping the 2D projection would do; MKL facts corrected

- User asked what happens if the 2D camera projection is dropped (not needed for the next
  project steps). A/B (scratch copy of `simulate` without the three `_get_2D_data` calls,
  force/constant, 40 steps, 3 views, quiet machine): 3D results bit-identical; 10.96 → 10.46
  ms/step with sensitivities, 8.11 → 7.59 without; allocations 0.55 → 0.31 MB/step (GC
  pauses roughly half as often). Breaks: `fit_model` and all line searches (`contour_cost` on
  `simBorderPts`/`gradList`), cost functions, `window_estimation.jl`, paper_1 pipeline/
  figures/metrics, WRITECONTOUR/ANIMATE output, and the tests carrying ~96 % of the suite's
  assertions (stokes_optimization, tikhonov, cost_functions, the rotated-view testset).
  The projection only needs `pos3D` and `gradList_3d`, which `simulate` already returns, so it
  can move out of `simulate` (or behind a keyword) without losing anything. Awaiting decision.
- Corrected earlier claims about MKL: MKL_jll 2025.2 covers only x86-64 Linux/Windows (no
  macOS), artifacts ~730 MB on disk (not ~300), Intel Simplified Software License,
  `using Pardiso` +0.55 s. README said "(x86-64)" — now "(x86-64 Linux and Windows)".

## 2026-10-14: camera projection moved out of `simulate` (Phase F); Pardiso kept

- User: keep MKL Pardiso with the UMFPACK fallback; move the projection out of `simulate`
  (the next project steps need no camera).
- `simulate` now returns 7 outputs (`output, displacement, surface_pts_3D, pos3D, velocity,
  pressure, gradList_3d`). New exported `project_contours(mdl, pos3D, conditions; grad3D)` →
  `borderPts2DList, gradList, pos2D`; owns WRITECONTOUR and the 2D GIFs; threaded over time
  steps (~2.3 ms for 41 states vs ~0.5 ms/step inside the loop before). Fits call
  `_simulate_contours`. Updated: write_sim_data, paper_1 pipeline (9 sites), tests, README,
  coding guide.
- Verified bit-identical against references saved from the old code (sims in 3 modes, gn/lm/
  Tikhonov fits, written files). Only intended difference: the first frame's contour
  sensitivities are 2×nb×2 zeros instead of 3×nb×2 (consumers read rows 1:2 only).
- Found: `contour_cost` uses only the first view of each frame; the paper pipeline fits one
  view per call at 0/5/15/30/45°, so the back_project gradient bug (fixed earlier today)
  affected every per-view fit at 5–45°, not "multi-view fits" as first reported.

## 2026-10-14: next phase keeps the solver as is

- User: next phase = same solver, a list of applied forces as in paper 1, fitting as a
  post-processing step; online use comes later. Step-by-step API and model reduction are
  parked (plan, "Status after Phase F").
- The back_project gradient bug does not touch paper 1's numbers as far as the data shows:
  every fit tree is view_1 except one synthetic cube run (view_2–5, July 2026) that no
  paper script reads.

## 2026-10-14: real-time work committed as three commits

`07594e9` fix (rotated-view contour sensitivities, minimal fix on the old back_project),
`409421a` perf (solver; pre-move state rebuilt from the transcript, bit-identical to saved
references), `1e2938e` refactor (projection out of simulate). Built and tested in a scratch git
worktree (copy the gitignored Manifest.toml and mesh caches in), then `git reset` (mixed) moved
the branch without touching the working tree. Not pushed. The user's own edits stay uncommitted.

## 2026-10-14: first CI run of the branch — 13 mesh errors, pre-existing, fixed (7c93012)

CI (workflow_dispatch; pushes to feature branches don't trigger it) failed on Ubuntu, Windows
and macOS with the same 13 errors, all Gmsh mesh generation (mesh_templates_test, tet cases
of verify_meshes_test); every solver test passed, incl. Pardiso on Windows and the UMFPACK
fallback on macOS (97,856 passed each). Cause: no gmsh CLI on CI → Julia-API fallback in one
persistent Gmsh session, and gmsh.clear() keeps global options, so Mesh.RecombineAll from hex
templates leaked into tet/tri/line meshes. Reproduced locally with PATH hiding /usr/bin/gmsh;
fix `gmsh.option.restoreDefaults()` before each mesh → full suite 97,987/97,987 without CLI
and with an empty mesh cache. Old CI logs (May) expired (HTTP 410), so "red since May for the
same reason" is likely, not confirmed.

## 2026-10-14: parent-folder handoff docs; two documentation layers

`smear-modules/claude_docs/` holds handoffs from the claude.ai Project: `handoff.md` (dated log,
newest first) and `HANDOFF.md` (a Paper 2 handoff for Claude Code; perception notes belong in
`../smearPerception/claud_docs/`, not here — user, 2026-10-14: smearFEM.jl is for modeling). User:
the parent docs are the global, conceptual view; each module's `claud_docs/` the local view;
update both each session. Added a conceptual 2026-10-14 FE entry + status line to handoff.md and
corrected HANDOFF.md §1/§5/§6 Q3/§7 (FE = smearFEM.jl, ~10 ms/step on the paper-1 hex mesh,
tets not timed, forward model 3D-only). Paper 2 uses tet elements and a 3D closest-point cost —
the next FE-side items are a tet benchmark and that residual.

Same day: the user asked to merge `handoff.md` and `HANDOFF.md` into one. Now a single
`smear-modules/claude_docs/handoff.md`: intro, status at a glance, reference §0–§8 (the former
HANDOFF.md, numbering kept so its internal §/2.5 references hold), session log §9 (the former
handoff.md entries, newest first). Every line carried over except the old headers and four
deliberate edits; pre-merge copies were in the session scratchpad only.

CI on `7c93012` (run 38547028134, workflow_dispatch): all green — Ubuntu 10 min, Windows 13 min,
macOS 31 min (x64 under Rosetta, UMFPACK fallback), docs 3 min. First green CI on record since
at least May 2026.

The user confirmed the green CI themselves. Branch layout for the merge (2026-10-14):
`features/realtime_solver` = `features/tet_mesh` + our 4 commits; `tet_mesh` = `publish_ready`
+ 147 commits; neither is merged and no PR exists. `main` has 5 commits the branch lacks (the
PR #21 merge, a Codecov-flags CI change, a docs/make.jl fix, a .gitignore line, a dependabot
merge); a dry-run `git merge-tree` of the branch into `main` is clean. Merge the two together:
`tet_mesh` alone still has the Gmsh option leak, so its CI would be red.

2026-10-14 (afternoon), user: do not merge `features/realtime_solver` into `main` yet; the FE solver
stays as it is for now.

## 2026-10-15: Branch `features/linear_elasticity` — real-time elasticity, Stokes/optimization stripped

User request: find the branch with the linear elasticity implementation, make it a separate
branch independent of the Stokes model and the optimization functions, optimized for real-time.
Decisions (user, via questions): **strip** Stokes + optimization (not just decouple); implement
**both** the legacy squeeze and the paper-2 volume update; sensitivities **optional** (off by
default).

**Where the old code was.** `src/examples/squeeze_linear_elasticity.jl` (+ `test/linear_elasticity_test.jl`,
`scripts/elasticity/`) exists only on the OLD lineage (dev, IGA*, gmsh_intergration, gpu_optimization,
model_optimization, …). Newest: `features/gmsh_intergration` @ `2b066a7` (2026-05-01);
`gpu_optimization`'s copy differs only by a kwarg name. The current lineage (`main` → `tet_mesh` →
`realtime_solver`) never had it — `git merge-base` with gmsh_intergration is `ca8a2ce`, where the
file did not exist. Old code used the pre-refactor API (`FunctionClass` strings, `ne^ndim` sizing,
`meshgrid_cylinder` 6-tuple, camera contours inside `simulate`), so it was **re-implemented on the
current FEM core**, not merged. The legacy analytic test (frictionless uniaxial: u_z = εz,
u_r = −νεr) was kept as the oracle.

**Branched from `features/realtime_solver` @ `ed36680`.** Removed: `src/solver/{stokes_solver,
stokes_setup,squeeze_system,run_example}.jl`, `src/optimization/`, `src/fem/{scenarios,geometries}.jl`,
`src/analysis/{experiment_tree,statistical_analysis}.jl`, `src/io/analysis_plots.jl`,
`write_stokes_scene`/`_extract_p_from_u_nodes` (io.jl), `plot_matches`/`plot_covariance!`
(plotting.jl), the Stokes/optimization tests and experiments (stokes/, convergence_analysis/,
paper_1/, model_optimization/, visualizations/). Kept: meshes, fem.jl, PostProcess.jl (camera
observation model — neither Stokes nor optimization), io, config, pointcloud_metrics,
`aruco_accuracy.jl` (deleted later that day, see the last 2026-10-15 entry). Deps dropped (Aqua stale-deps): Colors, IterativeSolvers, LinearSolve,
Measures, Printf. `pointcloud_metrics.jl` had been getting `KDTree` from the deleted
cost_functions.jl — now imports NearestNeighbors/Statistics itself.

**New code.**
- `LinearElasticity(; mesh, λ, μ)` (models.jl) replaces the stale struct; `lame_parameters(E, ν)`.
  Units kPa (mm, kg·mm/s²).
- `BasisFunctionCache(mesh)` now one mesh, concrete fields `N, ΔN, w, N_s, ΔN_s, w_s` (was a
  Stokes u/p/x `NTuple{7,Any}`).
- `src/solver/elastic_system.jl` — `ElasticSystem`: fixed pattern on a `dof` map (3 × nNodes,
  0 = eliminated, shared unknowns allowed), gather lists, threaded element kernel writing μA+λL
  and L (`A = 2ε:ε`, same node-pair form as the Stokes viscous kernel; `L = div div`), slip b on
  its own small pattern (`dKdβ`), `_refined_solve!` (single-rhs port of `_lagged_solve!`, same
  `_stale`/`_refactor_now` policy), `_symv!` column-parallel product (K must stay symmetric).
- `src/solver/elastic_squeeze.jl` — `simulate(mdl, cParam; control="force"|"displacement", β,
  sensitivities)`: unknowns = all components except bottom z (eliminated) and top z (→ one plate
  unknown, last). Displacement control = identity row/col (`_prescribe_plate!`), reaction from the
  saved plate column. cParam = per-step INCREMENT (legacy semantics). Unit-response solve + cubic
  extrapolation, refactorize at step 1 of every call. Sensitivities: solves for ∂/∂λ (L) and ∂/∂β
  (b); ∂/∂μ from homogeneity (Σθ∂z/∂θ = −z force, 0 displacement). Returns a NamedTuple
  `(plate, displacement, pos3D, grad3D)`; grad3D (3, nNodes, 3) = [λ, μ, β].
- `src/solver/volume_update.jl` — `volume_update!(mdl, nodes, normals, d; β_table=0)`: u·n = d
  (tangential free), table u_z = 0 auto-added, rest traction-free. Constrained nodes rotated into a
  frame spanning their normals (Gram–Schmidt, several constraints per node OK), prescribed rotated
  unknowns → identity rows. Pattern = full 3N unconstrained K (all node-pair blocks are full 3×3),
  so the Pardiso symbolic analysis is reused whichever nodes are observed. `β_table` is
  `[not decided]` regularization.

**Gotchas found.**
- **`top_nodes`/`bottom_nodes` omit the nodes shared with the side** (structured cuboid Q2 ne=4:
  49 of 81 top-face nodes). First smoke test was 3.5 % off the analytic solution because the
  plate rims were unconstrained. Plate/table nodes now come from `unique(vec(IEN_top/IEN_bottom))`.
- **Gmsh cylinder `IEN_sides` is unset** (4×1 garbage) — use `side_nodes`.
- **`reset_mesh!` aliases `NodeList = initial_state`** (no copy): solvers must replace NodeList
  (`X = X + u`), never mutate it. Tested (`initial_state == X` after runs).
- `mesh.nNodes` through the abstractly typed `mdl.mesh` made `volume_update!`'s loops dynamic:
  34 → 19 ms after `nN::Int = ...`. Generic `mul!` on 3×3 views was 18 % of a call → explicit loops.
- β = 0 leaves plate rigid modes free (x, y translation, z rotation): use small β. Pardiso then
  perturbs pivots silently, so `_refactor!` warns when `iparm(14) > 0`.

**Verification (all by running code).** Kernel = B C Bᵀ to 2e-16 (Hex1/2, distorted). Global K =
dense assembly to 4e-17. Uniaxial squeeze: u_z to 1e-11 (Hex2), plate force = EAε to 1e-10 (Hex2),
1e-8 (Tet2 — the 1e-9 residual tol × near-frictionless conditioning; test uses rtol 1e-6).
Sensitivities vs central FD, one step: 1e-7…1e-5 rel. Volume update vs exact linear field
(normal strains + rigid z-rotation, zero shear on faces): 1e-15 on Hex2, Tet1, Tet2. Mutation
test: `TET_MASTER_VOLUME = 1/6` makes the Tet force tests fail (guard works; replaces the deleted
`verify_meshes_test.jl`). `Pkg.test()`: **221/221, ~1 min**.

**Performance (Pardiso, 8 threads, paper-1 cylinder elem 6 mm, 12.5 % compression in 100 steps).**
Hex Q2 1885 nodes (n = 5366, nnz 862k): step median 7.2 ms / max 16 ms; with sensitivities
15.7 / 29 ms (17 % of steps > 20 ms — the refactorization steps). Tet P2 2555 nodes (n = 7128):
6.9 / 19 ms; sens 14.3 / 37 ms. Breakdown: assembly 2.2 ms (kernel 0.8, gather ~1.3), Pardiso
refactor ~9 ms every ~10 steps, solve ~1 ms, ~1 refinement round per solve. Combining μA+λL in
the kernel (instead of gathering A and L separately) and slip on its own pattern took assembly
7.3 → 2.2 ms. Volume update ~20 ms (half of it the fresh numeric factorization). Pattern build
~1.5 s once per mesh (cached in a WeakKeyDict).

**Open / not done.** Sensitivity-step spikes need a background refactorization (not built).
Nothing committed yet at the time of writing (user had not asked).

## 2026-10-15 (later): verification pass on `features/linear_elasticity` — what it found

User asked "did you check everything?". Ran the checks that had been skipped:

- **UMFPACK fallback** (every Mac/ARM machine, macOS CI): forced via
  `@eval smearFEM _sparse_lu(K::SparseMatrixCSC{Float64,Int}) = _Umfpack(nothing)` and ran the
  elasticity tests: all pass.
- **Gmsh Julia-API path** (what CI uses): `env PATH=/nonexistent <abs julia path>`, test file
  `include_string`-ed with `MESH_DIR` swapped for a scratch copy of `mesh_files/templates`, so the
  repo cache was untouched and every Tet box was generated by gmsh_jll: all pass.
- **Docs**: `docs/make.jl` + `doctest(smearFEM)` in a scratch env (Documenter + dev'd smearFEM, not
  the repo's docs/Project.toml, which pulls PyCall/Conda): both exit 0; only warnings are 54
  docstrings not in the manual (warnonly) and the skipped local deploy.
- **BUG (fixed): UMFPACK hides a singular stiffness.** The earlier comment "UMFPACK fails loudly"
  was false: `lu`/`lu!` keep rounding-level pivots, so a volume update with radial normals and
  `β_table = 0` silently returned an arbitrary 0.27 mm rigid rotation (Pardiso: 0.004 mm, and it
  warned). Now `_singular(f)`: Pardiso `iparm(14) > 0`; UMFPACK
  `F.info[LibSuiteSparse.UMFPACK_RCOND + 1] < 1e-12` (min|diag U|/max|diag U|: ~1e-15 singular,
  2.5e-9 for β = 1e-9 squeeze, ~1e-2 volume update with any table slip). `_refactor!` warns on
  either; new testset "singular stiffness is reported" (`@test_logs`, both backends). Relies on the
  `UmfpackLU.info` field (SparseArrays internals, present since 1.9).
- **BUG (fixed): `experiments/aruco_accuracy.jl` broke** — `using Measures` (for `pt`) after I
  dropped Measures from Project.toml. Now `using Plots, Plots.PlotMeasures`. Include + gate verified;
  `main()` not run (needs the calibration data, writes into the data tree). Moot: the script was
  deleted later that day (next entry).
- **Wrong number (fixed in README):** "pattern build ~1.5 s" was compilation. Compiled: 0.15–0.5 s.
- **Wrong scenario label (fixed):** the README timing table said "12.5 % over 100 steps" but its
  numbers came from a force run of unknown compression. Re-measured (browser active), 100 × 0.05 mm
  displacement steps = 12.5 %; force −25/step = 10.9 %, same within 0.5 ms:
  Pardiso Hex2 6.9 (max 17) ms, sens 15.0 (28) ms, 7/99 steps > 20 ms; Tet2 6.5 (17), sens 14.4
  (27); first step 17–27 ms; volume update 20 ms. UMFPACK: 7.7–7.8 ms (max 15–20; it never
  refactorizes after step 1, refinement keeps up), sens 22–24 ms (48–76/99 > 20 ms), first step
  130–175 ms, volume update 130–160 ms.
- **Multi-step sensitivity accuracy** (frozen geometry) vs FD of the full run: force 1/10/25 steps
  (0/0.2/0.6 % compression): ~1e-6 / 0.2 % / 0.5 %; displacement (0.5/5/12.5 %): ~1e-7 / 1.3 % /
  3.4 % (worst of λ, μ, β). Now stated in the `simulate` docstring and README.
- Julia 1.11 (the CI version) is NOT installed locally (only 1.12.6) — untested until CI runs.
- `Pkg.test()`: **225/225**, exit 0, no stray warnings (~1 min).

## 2026-10-15 (later): `aruco_accuracy.jl` removed — user: not needed on any branch

- Only three branch tips had it: `features/tet_mesh` (4578f67), `features/realtime_solver`
  (ed36680) and `features/linear_elasticity`. Checked with `git ls-tree` on every local and remote
  ref, after fetching the first two (their remote tips matched the local ones). Not on `main`,
  `dev`, `publish_ready` or any other branch.
- `features/linear_elasticity`: `git rm -f` (my `Plots.PlotMeasures` fix went with it); the README
  layout row no longer lists it. No other tracked file referred to it, so neither the package nor
  the tests change (Aqua stale-deps already passed: every dep is used by `src/`).
- The other two branches: the one-file removal commits, made from temporary worktrees so this
  checkout's uncommitted work stays untouched, were **denied by the auto-mode classifier ("Modify
  Shared Resources")**. Nothing was applied (worktree list and both tips checked); waiting for the
  user. Nothing else on those branches refers to the script. `plot_style.jl` stays there: there it
  is also included by `paper_1/pipeline.jl` (corrected: I first wrote `test_opt_stokes.jl`, which
  is not on those branches).
- Here, `experiments/plot_style.jl` (its only remaining user was the aruco script) and
  `experiments/ParallelExecution.jl` (its only users were the deleted batch scripts; the README
  still documents it) now have no users. Kept, flagged to the user.

## 2026-10-15 (later): committed on all three branches — user: "remove the unused scripts and commit on all branches"

- Checked first, by `git grep` on each branch: on `features/tet_mesh`/`features/realtime_solver`
  `plot_style.jl` is included by `paper_1/pipeline.jl`, and `ParallelExecution.jl` by
  `paper_1/pipeline.jl` and `stokes/write_gt_stokes.jl`. So on those two only the ArUco script goes.
- `features/tet_mesh` `fdf2989`, `features/realtime_solver` `d19ff68`: "chore: remove the ArUco
  accuracy experiment" (one file each), from temporary worktrees; auto mode allowed it this time.
  `git merge-tree` of the two: no conflicts (both delete the same file).
- `features/linear_elasticity`: the branch pointer moved `ed36680` → `d19ff68` with
  `git reset --soft` (ancestry checked first; index and working tree untouched), so the branch
  shares that commit and is `realtime_solver` + 1. Then `git rm` of `experiments/plot_style.jl` and
  `experiments/ParallelExecution.jl`, the README paragraph on ParallelExecution dropped, and
  everything committed as `cc86cbf` "feat: real-time linear elasticity, without the Stokes model
  and the optimization code" (46 files, +1231/−15644; the body lists the parts). Staged list
  checked: nothing from claud_docs/, CLAUDE.md, memory.md or Manifest.toml (all gitignored).
- Nothing pushed. tet_mesh and realtime_solver are 1 ahead of origin; linear_elasticity has no
  upstream yet.
- Clean checkout of `cc86cbf` (detached worktree in the scratchpad, Manifest.toml copied in,
  meshes regenerated by the gmsh CLI): `Pkg.test()` 225/225, 1m07s, exit 0. So the commit holds
  every file it needs.
- Then (user: "yes") pushed all three; `features/linear_elasticity` is new on origin, tracking
  it. No open PR has these branches as head (only dependabot PRs), so the pushes triggered no
  CI. Dispatched CI run 38772998784 on `features/linear_elasticity` (tests on Julia 1.11:
  ubuntu-22.04, macOS-latest, windows-latest; docs job on Julia 1, which does not deploy from a
  feature-branch dispatch).
- **CI run 38772998784: green.** 225/225 on every OS — test step ubuntu 1m14s, macOS 3m24s (x64
  Julia on the macOS runner; no MKL there, so UMFPACK), windows 1m25s; jobs 430/1017/685 s
  including setup. Documentation (build + doctests) 277 s. So Julia 1.11 and the UMFPACK path in
  CI are now verified, not just assumed.

## 2026-10-23: perception link (FE ↔ smearPerception), branch `features/linear_elasticity`, uncommitted

User: "we need a way to communicate between smearFEM and smearPerception"; chose (questions) to build
the offline and the live link now, shared memory for live (handoff §5 proposal); asked which side
should sample φ_s at FE nodes "computationally efficient" → FE side (field crosses each step anyway
for the cost; sampling on the perception side adds a per-window round trip; sampling itself is µs).

Built: `src/io/perception_link.jl` (+ `Mmap` stdlib dep, exports in smearFEM.jl), its Python twin
`smearPerception/src/smear_perception/fe_link.py` (holds the spec), `reconstruct_bag.py --shm`.
Design: one /dev/shm file per direction (plain mmap on both sides — NOT Python
multiprocessing.shared_memory, whose resource tracker in 3.12 also unlinks segments a process
only attached to — known CPython issue, `track=False` only from 3.13; from memory, not tested here),
64 KB header with self-describing JSON layout, ring of fixed-size slots, per-slot seqlock
(2k+1 writing, 2k+2 done), records-written counter; reader errors when overrun (never skips);
reader-attached flag so the writer waits while Julia loads (first test lost 16+ records to the
~5 s `using smearFEM`). Julia atomics: `unsafe_load(p, :acquire)` / `unsafe_store!(p, v, :release)`
(Julia ≥ 1.10). Python side relies on x86-64 TSO.

Verified (all by running code):
- Python writer → Julia reader and HDF5 file → Julia on a 5×7×9 grid: fields identical,
  `sample_field` = Python's to 1e-7 mm, NaN pattern identical (incl. last-layer and corner points).
- Julia writer → Python reader (FE → perception direction), overrun detection both ways.
- Live, real size (71×71×51, 10 k crossings, 100 k capacity, 16 slots), Python at 30 Hz, Julia
  concurrent: 600/600 in order; publish 0.63 ms; publish → available in Julia median 2.0 ms, p99
  3.9, max 4.6 ms (incl. 1 ms poll sleep). First read ~0.5 s (compilation); backlog drains.
- GC: first version allocated 5.1 MB/read → full GCs ~240 ms every ~8 steps (`@timed` gctime);
  now reads into preallocated buffers, 0 MB/read, no GC.
- End to end: Carreau/5 bag, steps 0–2, cap 7 mm, real driver with `--shm`: what Julia got over shm
  is bit-identical to what it reads back from the HDF5 file of the same run (NaN normals included).
- `Pkg.test()` 261/261 (~1 min); Python `tests/test_fe_link.py` and `test_imports.py` pass.
Not done: FE-side consumer (volume-update correspondences, closest-point cost), the window-start
mesh message FE → perception (channel exists, record not defined), live rclpy node.

## 2026-10-23 (later): launch order, stale segments, end of stream

User: "can I launch everything at the same time with the live method?" Tested before answering:
both together worked, but (B) perception > 10 s slower → FE `open` timed out (default 10 s) and
perception then waited forever in `wait_for_reader`; (C) a segment left by a crashed run → the FE
attached to it, the new Python writer truncated and rewrote THE SAME inode (O_TRUNC), the FE errored,
perception waited forever. Fixed on both sides: writer PID at header offset 48; readers wait only for
a segment whose writer is alive (`kill(pid, 0)`); the Python writer now unlinks before creating (as
Julia's already did); default open timeouts Inf; a reader whose writer is gone returns
`nothing`/None once every record is read (end of bag or crash), so
`while (st = next_step!(link)) !== nothing` ends. Re-tested: A together, B perception 15 s late,
C stale segment, D perception SIGKILLed after 3 steps (FE gets 0–2, then ends), E FE 30 s early —
all correct. Test added (stale PID ignored; `getpid` needs the process still running, so
`run(...; wait=false)`). `Pkg.test()` 262/262. Residual: a stale file whose PID was reused by an
unrelated process would look live (unlikely; delete /dev/shm/smear_* if a reader attaches to nothing).

## 2026-10-23 (review): is the shared-memory ring the right live transport? Measured: no

User asked for a critical re-check of the Julia↔Python link. Measured on the laptop, Carreau/5 grid
(71×71×51, 10 k crossings), 30 Hz, concurrent Julia reader, latency publish → step in Julia, steps 31+:
- ring as built (1 ms sleep poll): median 2.0–2.5 ms, p99 3.9–5.7, max 4.6–7.1 ms
- ring, reader spinning (`sched_yield`, burns a core; @eval'd `_read_into!`): median 1.1, p99 3.6, max 6.0 ms
- Unix-domain socket prototype (scratch: Python `sendall` of a 64-byte header + φ + count + n×3 points
  + normals; Julia `listen`/`accept`, `read!`/`unsafe_read` into preallocated buffers):
  1.53 MB/msg, median 0.44, p99 1.22–1.35, max 1.24–1.74 ms; ~1.4 KB allocated per message (one
  2.2 MB stream-buffer growth once); back-to-back 600 msgs: median 0.14 ms. Excludes the ×1000 mm
  conversion (~0.1 ms, inferred).
Why the ring loses: it copies twice per step (Python → shm, shm → Julia buffer; the seqlock check needs
the copy), the same count as a socket (→ kernel → user), but moves 3.7 MB (fixed capacity) vs 1.53 MB
and does the copies in series; the socket pipelines them. Zero-copy would need in-place consumption +
slot ownership (more complexity).
Socket robustness, verified: sender first or receiver first both work (client retries connect, listener
`rm`s a stale socket file); sender SIGKILLed after 54 msgs → Julia `EOFError` immediately.
Ring fragility found: first `next_step!` compiles 0.43–0.78 s; 16 slots at 30 Hz = 0.53 s → one run
overflowed (the earlier 600/600 runs passed by a small margin; at 10 Hz the margin is 1.6 s).
Checked and NOT a problem: a dtype change on the writer (e.g. φ Float64) fails loudly (MethodError in
`open_perception_link`), not silently.
Remaining review items (transport-independent unless noted): no step index j in the record or the HDF5;
offline reader loads every step at once and has a different API/type from the live one; the FE must drain
the link in its own task with window-sized buffering before online estimation (a GN pass takes seconds;
16 slots = 1.6 s at 10 Hz); ring-only: missing seqlock fences (formally racy; x86 TSO + opaque calls
make it work in practice; unsafe on ARM), no close flag (publisher closed in a live REPL → reader waits
forever), late joiner errors at once, pid 0 counts as alive and VERSION was not bumped for the PID field,
fixed max_points; `sample_field` allocates two small arrays per point.
Recommendation given to the owner: keep two processes and the HDF5 offline path; replace the ring by a
socket stream (UDS, same code over TCP for a second machine); add j; one Julia iterator for file and live.
Owner's decision pending; nothing changed in the code by this review.

## 2026-10-23 (latest): socket link and the per-step record (owner decisions); ring removed

Owner, after the review: socket stream (the recommended option). Then: "each index should consist of
surface observation, the fact if the probe is in contact or not, if in contact the force
measurement"; chose (questions): contact by force threshold (value not decided), the wrench at every
step (not only in contact), the probe pose at every step.
Built: `perception_link.jl` rewritten — `ShmChannel` and its API removed, `Mmap` dep → `Sockets`;
`PerceptionLink(io)`, `open_perception_link(path)` (retries connect), `next_step!(link)` (no timeout
argument any more); `PerceptionStep` gains `n`, `contact`, `wrench`, `probe_pose`;
`read_perception_steps` returns `(; grid, steps, t0_ns, dt_ns)`. Python twin (`fe_link.py`) rewritten
the same way; driver `reconstruct_bag.py --link [path] --contact-threshold F`.
Design choices (Claude, owner may veto): `contact` Int8 with −1 = unknown (no threshold set, or NaN
wrench), so missing data is never read as "no contact"; pose as perception's 4×4 T (no quaternion-
order ambiguity), transposed back on read; END message so a crash (EOF) differs from a complete run;
hello sent with the first step, so t0 = t_n − n·Δt is known; "STEP"/"END." tags catch a desync; FE
units at the boundary (forces ×1e3 → kg·mm/s², torques ×1e6, lengths ×1e3).
Verified by running code:
- `Pkg.test()` 275/275; link testset 50 checks: HDF5 fixture, the byte stream Python's publisher sent
  (`stream.bin`), socket (stale file → retries then timeout error; stream cut inside the last step →
  `EOFError`). `read!` on the contiguous point views takes Base's StridedArray path: 0 allocations.
- Carreau/5 steps 2–4, real driver with `--link` and `--save-field`, FE reading concurrently: live steps
  bit-identical to the HDF5 read-back (17/17), n = 2, 3, 4, t0 = t_2 − 2Δt, END received, socket file gone.
- Latency, real classes, 71×71×51 grid + 10 k crossings (1.53 MB per step), send → converted step in
  Julia, steps 31+: at 30 Hz median 0.86 ms, p99 1.59, max 2.05 ms; back to back median 0.33, p99 0.55,
  max 0.8 ms; 2.1–2.5 KB allocated per step. First step 3–5 s (FE load + compile; perception waits,
  nothing lost). Above the prototype's 0.44 ms at 30 Hz: float64 → float32 conversion in Python and the
  unit conversion here (inferred); idle wake-up between steps would explain the back-to-back figure being
  lower (inferred).
- Start order and crashes: FE 15 s before perception ✓; perception SIGKILLed → FE `EOFError` after 128
  steps (socket file left) ✓; with that stale file, FE first and perception 5 s later ✓ (file replaced,
  removed after the clean run); FE SIGKILLed → perception `BrokenPipeError`, socket file removed ✓.
- F/T and probe path: a TEST bag (scratch `make_sensor_bag.py`: 4 Carreau/5 frames plus made-up 100 Hz
  signals, WrenchStamped F = (0.1, 0, −3t) N, M = (0.01, 0.02, 0.5t) N·m, and /tf world → probe, 30° about
  z, (0.09 + 0.01t, 0, 0.05) m) → driver with `--contact-threshold 0.5 --link` → Julia: wrench, contact
  [0, 0, 1, 1] and pose exact in FE units (13/13); the HDF5 file holds the same. Bag deleted (385 MB).
Open: FE consumer of the record; FE → perception message (window-start mesh; the socket is two-way);
for a live camera node, send from a background thread with a bounded queue so a slow FE cannot stall
the camera loop (blocking is right for bag replay); one Julia iterator for file and live (the offline
reader still loads a whole file).
- Buffering between the processes, measured (user asked whether there is a buffer if the FE is late):
  effectively none. Kernel socket buffer = SO_SNDBUF 212 992 B (`net.core.wmem_default`; `wmem_max`
  4 MB) ≈ 1/7 of a Carreau/5 step; Julia reads the socket only while a task waits in `next_step!` and
  calls `stop_reading` when the wait ends (`base/stream.jl:422`, 1.12.6), so nothing accumulates on the
  FE side while it computes. FE paused 5 s after 3 steps (sleep, or a non-yielding busy loop):
  perception completed 0 steps during the pause, resumed ~1 ms after the FE did, 30/30 delivered. So
  an FE delay goes straight back to perception: fine for bag replay; live it would stall the camera
  loop. A Julia reader thread alone would not fix it (GC stops every thread; pauses ~60–320 ms measured
  in this package). Proposed: perception-side queue + sender thread, sized for the longest transient FE
  delay (~1.5 MB per step; 10 s at 10 Hz = 150 MB); overflow policy is the owner's call.
- Window memory for the estimation (discussed with the user, not built): each GN/LM iteration re-simulates
  the window, so every step of the window stays in memory until its estimation ends — per step the
  observed surface (points + normals, ~0.25 MB at Carreau/5) for the closest-point cost, wrench/contact/
  pose for the forward run; at the window start the field for the volume update. Dense φ_s must also be
  kept until the window-start FE shape exists (combined field, §2.5), per window if windows overlap
  (handoff §6 q17). Proposed API: a preallocated store of `PerceptionStep` slots addressed by n and
  `next_step!(link, slot)` reading straight into a slot (no deepcopy → no GC churn).

## 2026-10-27 — Integrity check after the file deletion/recovery

- `git fsck` clean; every tracked file present, no NUL bytes or invalid UTF-8 (text files incl. claud_docs).
- `test/data/perception_link/` had come back as `test/data/#9977220` (e2fsck's lost+found name for an
  orphaned directory); contents intact (valid HDF5 headers), renamed back.
- `Pkg.test()`: 225/225 pass (2.5 min incl. precompile). `test/perception_link_test.jl` is not in
  `runtests.jl`; run alone it passes 50/50 (also inside Claude Code's sandbox).

## 2026-10-28 — perception link: fields permuted on read; `sample_field` by Interpolations.jl (uncommitted)

- User: permute the fields to perception's axis order (`phi[i,j,k]` = Python's), and use
  off-the-shelf libraries. `next_step!` reads into `link.phi`/`link.count` (wire order) and
  `permutedims!` into the step; `read_perception_steps` permutes per step.
- `sample_field` now Interpolations.jl (`extrapolate(scale(interpolate!(phi, BSpline(Linear())),
  centres...), NaN)`): no copy of the field (656 B per call vs 78 KB with `interpolate`), 0 B per
  sample. Interpolations 0.16.2 was already in the Manifest (indirect); added to Project.toml
  [deps] + compat "0.16". Manifest NOT re-resolved: the only Julia installed is 1.13.3 and a
  resolve rewrote the whole 1.12.6 Manifest (stdlib bumps, format 2.1) — reverted; loading works
  without it. User to run `Pkg.resolve()` when they choose to move the Manifest to 1.13.
- **Bug found by the swap (both languages):** the hand-written sampler's "last layer" special case
  caught every point in [n−1, n) grid units, so points up to one voxel past the last voxel centre
  got the last layer's value instead of NaN. Interpolations.jl/SciPy return NaN there. Elsewhere
  identical (7e-16 on 20 000 random points, same NaN pattern). Low impact (grid edges are far
  from the specimen).
- Fixture: `expected.h5` recomputed with the SciPy sampler from the stored fields (as
  `make_fixture.py` does; the script itself could not run: Unix socket blocked in the sandbox,
  bypass declined). 110 samples changed, all past the last centre. Test: size check now (5, 7, 9);
  "defined points > 100" lowered to > 50 (now 97/123/151 of 302 by design).
- Run (sandbox, scratch depot stacked over ~/.julia): 45/46; the socket testset errors only
  because /tmp is read-only in the sandbox. Full `Pkg.test()` not run.
