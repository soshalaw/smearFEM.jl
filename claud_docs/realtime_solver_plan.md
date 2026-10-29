# Real-time solver plan — branch `features/realtime_solver`

Living plan: tick items off and add findings as they land, so a session that runs out of
context can resume from here. History and reasoning: the 2026-10-13 entries at the end of
`memory.md`.

## Goal (user, 2026-10-13)

**Everything runs at ≥ 50 Hz**: one dt = 0.1 s step of `simulate` in **≤ 20 ms wall time**,
for **both force and velocity control**, inside `simulate` itself so `fit_model` gets it too
(sensitivities included). Faster than 50 Hz is welcome; 50 Hz is the floor. Stay in Julia.

Acceptance, on the benchmark case below:
- mean ≤ 20 ms/step **with** sensitivities (the fitting path), and worst step ≤ 20 ms too —
  a synchronous refactorization stall breaks 50 Hz even if the mean is fine;
- `simulate` return values vs the pre-branch code (HEAD `4578f67`): state outputs within
  **1e-8 relative**, sensitivity outputs (`gradList`, `gradList_3d`) within **1e-6** — the
  sensitivities only feed the Gauss-Newton step and already ignore ∂K/∂x·∂x/∂θ;
- full test suite green (`Pkg.test()`; fast group `SMEAR_TEST_GROUP=fast` while iterating).

## Benchmark / verification harness (kept here, NOT in the scratchpad — it is wiped on reboot)

- Case: `Cylinder(25, 40)`, ne = 6, Hex2/Hex1/Hex2 (1885 u-nodes, 287 p-nodes, saddle
  5653 DOF), η = β = 100, F = −9518.61, dt = 0.1, 40 steps, views 0/30/60.
- `realtime/ref_run.jl <out.jls> [constant|bulk_viscosity] [force|velocity]` — runs
  `simulate` once to warm up, times the second run, serializes the 10-tuple.
  `julia -O3 --threads=auto --project=<repo> claud_docs/realtime/ref_run.jl out.jls constant force`
- `realtime/diff.jl a.jls b.jls` — per-return-value relative diff.
- References from the pre-branch code: `realtime/ref/old_{constant,bulk_viscosity}.jls` (force)
  and `old_velocity_{constant,bulk_viscosity}.jls` (velocity, from a patched export — see C1).
  To regenerate: `git archive 4578f67 | tar -x -C <dir>`, copy `Manifest.toml`, run
  `ref_run.jl` with `--project=<dir>`.
- `realtime/step_cost.jl [vtype] [control]` — steady per-step cost (slope over 2/10/40-step
  runs) and per-call overhead; `BLAS_THREADS=n` env to set BLAS threads.
- `realtime/instrument.py <repo> <out>` + `instrument_run.jl [vtype] [control] [nsteps]` —
  copy of `src/` with lap timers: per-section ms and MB per step, per-step wall times
  (step 1, steps 2–5, median/max of the rest), corrections/solve, refresh/swap/wait counts,
  GC total. Markers are matched by code text; it fails loudly if one is not found once.
- `realtime/umfpack_concurrency.jl [rounds]` — stress test of concurrent UMFPACK use on
  distinct UmfpackLU objects (no crash, bit-identical to serial); evidence for B6.
- **Timing on this laptop (i7-11800H, 8C/16T, 24 MB L3, powersave governor) drifts ±2 ms
  between runs and much more when the user's browsers are busy** — compare configurations
  back to back, repeat, and check `uptime`/`ps` before trusting a number.

## Status

### Phase A — validate the force-control rewrite
- [x] A1 before/after diff with tol 1e-10 + extrapolation: worst rel. diff 5.1e-10 (constant),
      4.3e-10 (bulk). Speed 586 → 105 ms/step (constant), 560 → 114 (bulk).
- [x] A2 fast group green on this branch: 97,879/97,879, 3m26s (slow group not yet run).
- [x] A3 breakdown (`realtime/profile_step.jl`, `profile_flat.jl`, `convergence_study.jl`):
      93 ms/step after B2. **70 % is the two lagged solves** (sens solve, 2 RHS, 1575 samples;
      main 865). 1 RHS triangular solve 3.5 ms, 2 RHS 6.2 ms, K*x 0.5 ms. Only the initial
      factorization happens in 40 steps; corrections/solve mean 3.74, max 5.
      `_squeeze_system` **527 ms + first `lu` 239 ms per `simulate` call** (≈ 19 ms/step
      amortized over 40 steps, and paid on every fit evaluation). Assembly 9.8 ms =
      kernels 3.6 (16 threads; 21 serial, scales to ~8 threads) + **serial scatter 5.5** + slip
      0.4. `_get_2D_data` GRAD 2.8 ms after B2.
      Convergence study (benchmark trajectory): initial relative residual for extrapolation
      order 0/1/2 = 7e-4 / 1.5e-6 / **6.9e-9**; one refinement step contracts by
      **ρ ≈ 1.8e-3 × s** (s = factorization age in steps); GMRES saves ≤ 1 iteration over
      plain refinement — not worth it. UMFPACK ordering: default already picks the symmetric
      strategy; METIS no better; "best" only −10 % fill — not a lever.
      ⇒ quadratic guess + factor age ≤ ~7 steps ⇒ 1 correction at tol 1e-10; at tol 1e-9 one
      correction holds up to s ≈ 75 (no refactor at all inside a typical horizon).

### Phase B — force control to ≤ 20 ms/step
Expected budget after B7–B10: assembly ~4 + main solve ~4.5 + sens solve ~5 + 2D 2.8 +
bookkeeping ~2 ≈ 18 ms.
- [x] B7 homogeneity ⇒ ONE sensitivity solve: M is ηĀ + βb + D with a θ-free RHS, so
      η ∂sol/∂η + β ∂sol/∂β = −[u; 0] (u = q_f and μ, pressure part 0). Checked to 1e-14
      (`realtime/check_homogeneity.jl`). Solve for ∂/∂β, derive ∂/∂η.
- [x] B8 quadratic extrapolation (`_extrapolate!`/`_push_history!`). Also: the loop now solves
      for the **unit-force response z** and scales by F_t (robust to noisy measured force,
      no F = 0 special case), and keeps the history **at unit viscosity**
      (`_scale_unknowns!`: velocity × η, ∂/∂β velocity × η², pressure × η) — only 3.56 → 3.28
      corrections in bulk mode; remove if it never pays off.
- [x] B9 `_cached_squeeze_system` (WeakKeyDict per model) + `_refactor!` at step 1 of every
      call (keeps each call's result independent of earlier calls). `_squeeze_system`
      527 → 348 ms with concrete locals, now once per model.
- [x] B10 threaded gather (`_gather_lists`/`_gather_entries!`, same summation order as the
      serial scatter). **No speed-up**: the pass is memory-bandwidth bound (~70 MB/pass at
      ~15 GB/s; 4.7 ms on 1 or 8 threads). Next lever: fewer bytes (Int32 lists, don't store
      nz_visc/nz_div, symmetric half of Ke_visc).
- [x] B11 tolerances: `STATE_TOL = 1e-9`, `SENS_TOL = 1e-8` ⇒ 1.06 corrections/solve
      (was 2.02). Acceptance revised: state outputs ≤ 1e-8, sensitivity outputs ≤ 1e-6 vs the
      old code (now 3.7e-9 / 2.7e-9).
- Instrumented breakdown (copy of `src/` with lap timers, 2026-10-13, before B11):
  assembly 11.2, main solve 9.3, sens solve 10.0, 2D data 4.9 (in situ), bookkeeping 1.6,
  step-1 refactor 3.8 amortized. `realtime/step_cost.jl` fits per-call vs per-step cost.
- [x] B12 `_correction`: F factorized at η_factor, K at η ⇒ rescale around `F \ res` by
      diag(√(η_factor/η) velocity, √(η/η_factor) pressure) (K(η,β) = D K(1,β/η) D). Bulk only
      2.21 → 2.0 corrections: the factorization was not the main problem —
- [x] B16 — the **guess** was: bulk r0 ~3e-7 even with exact history (η(t) varies, so z(t)
      is less smooth). Cubic extrapolation of the unit-viscosity history: bulk r0 1.7e-7 →
      1.3e-8 offline (unit scaling is worth 7× there, keep it), constant 7e-9 → 5e-11
      offline. In situ the stored solutions carry ~tol-level residuals that a cubic
      amplifies up to 15× (|4|+|6|+|4|+|1|), so constant mode still takes ~1 correction per
      solve; bulk drops 2.0 → 1.1. **Both modes now 20.0 ms/step steady** (worst diff vs
      old: constant 2.6e-9, bulk 4.5e-10).
- [x] B13 `get_convex_hull_boundary_points` rewritten with scalar arithmetic (same evaluation
      order): 0.42 → 0.021 ms/view; `_get_2D_data` 3 views 2.6 → 1.4 ms. Same border sets.
- [x] B14 bookkeeping: `q = C_uc*q_f + μ q_top` into an n×1 *sparse* matrix + hcat per
      component replaced by `_nodal_velocity(sol, sys)` (gather through `sys.node_dofs`; the
      products only select, so values are identical): 3.5 → 0.13 ms/step.
- [x] B15 element kernel computes only node pairs a ≤ b of the symmetric viscous block; the
      gather lists (`src_of`) read the mirror for the rest. **Steady 20.0 ms/step** (constant).
- Breakdown at 21.1 ms (before B15), `instrument.py` + `instrument_run.jl`: assembly 7.7,
  main solve 5.6, sens solve 5.8, 2D 2.2, bookkeeping 0.13; per call ~175–185 ms (step-1
  `_refactor!` ~150 ms + setup).
- [x] B17 `_residual!`: rhs − Kᵀx split by column over threads (K symmetric to 3.5e-18
      relative; 0.50 → 0.12 ms per product). **Invariant: K must stay symmetric.**
- [x] B18 `dKdβ` on its own compact pattern (8,452 nonzeros vs K's 1.0M: 0.5 → ~0.01 ms per
      product); `slip_pos` maps it into K, β·b is added after the gather.
- **Steady 16.9 / 16.1 (force), 17.5 / 17.9 (velocity) ms/step** (constant / bulk) on a quiet
  machine; all outputs within 2.6e-9 of the references. Under the user's browser load later
  the same code measured 16.5–22.8 (force) and 17.4–25.3 (velocity).
- Remaining ideas for margin: one 2-RHS correction for z and w together (6.0 vs 2×3.6 ms),
  `back_project` per-node allocations (~0.4 ms), skip `fit_curve` unless WRITECONTOUR,
  kernel BLAS calls on 3×27 matrices (`mul!(Jac, …)`, `mul!(G, dN, dN')`) per Gauss point.
Then, from the earlier list:
- [✗] B1 background refresh — **built, measured, removed** (2026-10-13). Design was: when the
      first correction contracts by less than REFRESH_RHO, snapshot K and `lu!` it into a
      second UmfpackLU on the last thread (pinned via `jl_set_task_tid`), swap it in a fixed
      number of solves later (deterministic). Back-to-back A/B, same load:
      40 steps, ρ-trigger 0.02 — force 21–23 vs 16.5–19 ms/step (on vs off), velocity 33–35
      vs 22–25; trigger 0.1 — force equal (never fires), velocity 33 vs 17 (fires twice:
      assembly 6.2 → 12.6 ms/step from contention, GC 15 → 113 ms, a 151 ms step);
      200 steps — medians equal, corrections 1.39 → 1.12 (force) and synchronous refactors
      5 → 2 (velocity). The factorization competes for cores, memory bandwidth and 8 OpenBLAS
      threads, and its ~84 MB of counted mallocs plus a long ccall (no safepoint) feed GC
      stalls. Not worth an internal runtime API and an unlocked UMFPACK call. If revisited:
      Julia ≥ 1.12 `@ccall gc_safe=true` around the numeric factorization, BLAS at 1 thread,
      and a trigger that only fires in long runs. Kept from it: `_correction` (exact η
      rescaling), `realtime/umfpack_concurrency.jl` (no crash, bit-identical — 1 refresh vs
      solves, and 4 workers each doing lu/lu!/solve). Also learned: Julia 1.12 `-t auto`
      gives 16 default + 1 interactive thread and the main task runs on the interactive one;
      on 1.11 / `-t N,0` a long task picked up by thread 1 would block the sticky main task.
- [x] B2 `_get_2D_data` GRAD path was O(n²): `R*dxdθ[:,:,i]` inside the per-node loop
      (`PostProcess.jl` `back_project`) — hoisted to `dxdθ_cam`; same arithmetic.
- [x] B3 allocations: B14 + `back_project` per-node 3×3 algebra in scalars (0.88 → 0.16 MB
      and 0.27 → ~0.1 ms per view; `_get_2D_data` 3 views 1.5 → 0.77 ms, 3.1 → 0.94 MB).
- [x] B4 ordering — no lever (A3).
- [x] B5 `simulate(...; sensitivities=false)`: skips the w solve, the 3D sensitivity update
      and the GRAD projection; `gradList`/`gradList_3d` keep only their initial entry.
- [~] B6 `SPARSE_LU_LOCK` kept around every UMFPACK call (its 2026-07-07 "exp" commit records no
      cause for the segfault it guards against). `umfpack_concurrency.jl` suggests distinct
      objects are safe concurrently, so parallel fits could drop it — user's call.

### Phase C — velocity control on the same machinery
Key fact (checked 2026-10-13): `set_boundary_cond` builds `C_uc`/`q_top` independent of
`control`, so with `μu_btm = μu_side = 0` velocity control is the **leading
(n−1)×(n−1) block of the force system** with μ (the plate velocity) prescribed:
`M_ff x = −M_fμ v`, and `M_ff ∂x/∂θ = −(∂M/∂θ · [x; v])[1:n−1]`. Plan: same `SqueezeSystem`,
replace row/column n by the identity after reading the coupling column (keeps the pattern
and the factorization machinery), RHS from that column.
**Velocity control is broken at HEAD** (found 2026-10-13): `UndefVarError: zero_matrix_p_q`.
It was only ever assigned inside the *force* branch (an `if` block is not a scope, but the
assignment never runs under velocity control); the dead-code commit 866677d removed even
that. So it has not run since at least the v1.0.0 squash, and no test covers it.
- [x] C1 references `realtime/ref/old_velocity_{constant,bulk_viscosity}.jls`, from the HEAD
      export with ONE patch: `zero_matrix_p_q = zeros(size(B,2), size(q_d,2))` at the top of
      the velocity branch (the force branch's definition). cParam = −1 mm/s. Old code
      ~584 ms/step. Old velocity outputs: `output` empty (push commented out), `gradList_3d`
      only the initial entry — compare the other 8 return values.
- [x] C2 one loop for both modes in `simulate` (old velocity branch and its dead setup
      deleted; `ForceSystem` renamed `SqueezeSystem`, file `src/solver/squeeze_system.jl`).
      `_prescribe_plate_velocity!` zeroes row/column n (keeps the diagonal = 1, pattern
      unchanged) and moves the column to the RHS for a unit plate velocity; sol = v_t z.
      Velocity-control homogeneity is different: η∂x/∂η + β∂x/∂β = [0; p; 0]. Unit-viscosity
      history scaling: velocity part η^(a-1), pressure η^a, w one power lower (a = 1 velocity,
      0 force). Diff vs C1: all 8 comparable outputs ≤ 1.04e-9 (gradList 6.3e-10 validates
      the identity). New behaviour: `output` = v·dt per step (was empty — deferred item #5),
      `gradList_3d` filled every step. `μu_btm`/`μu_side` ≠ 0 now throw in both modes.
      Steady: force 19.6 / 19.2, velocity 21.6 / 21.8 ms/step (constant / bulk).
- [x] C3 `test/squeeze_system_test.jl` (fast group, 45 tests, ~50 s): fixed-pattern K and
      dKdβ equal the old sparse-product construction; K symmetric; homogeneity in both
      modes; force run's plate velocity replayed under velocity control reproduces positions
      (1e-8) and pressure (1e-6); velocity sensitivities vs finite differences at step 1
      (exact there — later steps differ ~1 % at 25 % compression because the sensitivities
      omit ∂K/∂x·∂x/∂θ, as they always have); repeat calls bit-identical; prediction-only
      state matches; nonzero `μu_*` throws.

### Phase D — wrap-up
- [x] D1 docs: `simulate` docstring (control modes, `sensitivities`, μu restriction, solver
      note), README Performance paragraph, `claud_docs/CLAUDE.md` real-time section rewritten,
      `memory.md` entry "2026-10-13 (later)". README's "BLAS is set at module load" sentence
      left as is — it is false (top-level `BLAS.set_num_threads` only runs at precompile) and
      reported to the user rather than "fixed" one way or the other.
- [x] D2 final code (background refresh removed): fast group 97,924/97,924 in 2m57s, slow
      group 33/33 in 3m40s (was ~8.5 min) — **97,957 total, all green**. Final timings below
      are under browser load; the quiet-machine numbers (B18) predate only cleanups that remove
      work (back_project scalars, B5 restructure, refresh removal).
- [x] D3 not needed for 50 Hz: the full solve now meets it. For the record, homogeneity makes
      it exact that under force control with constant θ, (sη, sβ) traces the same shapes as
      (η, β) at 1/s the speed (continuous in time; forward Euler at fixed dt only
      approximately) — so shape-vs-height is a one-parameter family in κ = β/η, should a
      lookup/ROM ever be wanted for much finer meshes.
- Open for the user: commit; global BLAS thread policy; dropping SPARSE_LU_LOCK for parallel
  workers (B6); `sensitivities=false` in prediction-only calls of `experiments/paper_1/
  pipeline.jl`. Still above 20 ms: the per-call setup (step 1, ~170–200 ms: the step-1
  refactorization), steps 2–3 (~30–40 ms: low-order guess, 2–3 corrections), occasional GC
  pauses (incremental ~2–7 ms; a 50–90 ms one every few hundred steps), and in long runs a
  synchronous refactorization (~150–200 ms) whenever a solve needs > 5 corrections — about
  2 per 200 steps in force mode, 5 in the velocity case that compresses 50 %.

### Phase E — push below 20 ms everywhere (user, 2026-10-13: "open to anything as long as
the code works without losing accuracy") — STATUS AT SESSION END (2026-10-14 ~00:00)
Constraint kept: same tolerances (STATE_TOL 1e-9, SENS_TOL 1e-8), same discretization.
- [x] E1 `_lagged_solve!(sys, z, rhs_z, w, rhs_w!)` solves state + sensitivity together; w's RHS
      is rebuilt from the corrected z before the convergence check (same acceptance). UMFPACK
      backend runs both corrections at once (`copy(F)` = second workspace; 2 solves 6.1 → 3.5 ms).
- [x] E2 element kernel rewritten (`_element_matrices!`, `_inverse_jacobian`): gradients of all
      Gauss points first (gp fastest), node-pair 3×3 blocks summed in registers, no BLAS on
      3×27. Assembly 6.1 → 2.8 ms (3.3 in situ).
- [x] E3 per-chunk kernel scratch in `sys.scratch`; in-place `_nodal_velocity!` + `dNodeList .+=`;
      2D projection (2026-10-14): `_get_2D_data(...; FIT)` skips the spline resampling unless
      `conditions.WRITECONTOUR` (it only feeds the contour file), `back_project` computes each
      point with scalars (`_to_camera`, `_rotation_z`), the hull takes `eachcol` instead of one
      vector per point. 2D 1.15 → 0.31 ms/step, 0.93 → 0.25 MB/step; GC per 200 steps ~170 →
      ~6 ms (one 65 ms pause still seen in one run).
- [x] E4 views of `_get_2D_data` projected in parallel (`Threads.@spawn` per view).
- **BUG FOUND AND FIXED (2026-10-14, pre-existing, not from this branch):** `back_project`
  rotated node positions by the view angle but not their sensitivities `dxdθ`, so every view
  other than 0° got wrong 2D contour sensitivities (FD at step 1, 45° view: error 1.1e-3 on
  9.5e-3, ~12 %; now 1.4e-10). The 0° view was right (single-view tests never saw it).
  Multi-view fits used inconsistent Gauss-Newton gradients. Consequence for this harness:
  the reference `gradList` for views 2–3 (30°, 60°) is wrong; compare view 1 only (≤ 2e-9).
  New test: "contour sensitivities of a rotated view vs finite differences".
- [x] E5 **MKL Pardiso backend, in the repo, UNCOMMITTED, full suite NOT yet run.**
      `src/solver/sparse_lu.jl` (included before squeeze_system.jl): `_Pardiso` (default when
      `Pardiso.mkl_is_available()`) and `_Umfpack` (fallback, keeps SPARSE_LU_LOCK + concurrent
      pair solves). `Pardiso` added to Project.toml/Manifest (+ MKL_jll, IntelOpenMP_jll,
      oneTBB_jll). Measured on K (n = 5653): Pardiso numeric **9.3 ms** (8 threads) vs UMFPACK
      150; solve **1.27 ms** vs 3.6; analysis ~150 ms once per model (cached); residual 1e-14.
      Gotchas fixed: (1) Pardiso.jl uses Int32 indices here (MklInt = Int32 unless Julia's
      BLAS is MKL) → store the lower triangle as SparseMatrixCSC{Float64,Int32}, else 2 MB is
      converted per call (GC spikes); (2) MKL symmetric needs every diagonal present → added
      zeros (`from` = 0); (3) multithreaded Pardiso was NOT bit-reproducible run to run
      (~3e-14) → `set_iparm!(ps, 34, nprocs)` fixes it (12 runs identical); (4) iparm[8]=0
      still means "2 internal refinement steps if pivots were perturbed" (0 perturbed here).
      Refactor policy `_refactor_now`: Pardiso as soon as a round contracts by less than 20×
      (ρ > 0.05); UMFPACK after 5 rounds. Loop throws if > 2·max_rounds.
      Verified: squeeze_system tests 45/45; all four modes ≤ 2.5e-9 vs references.
      Prototype timings (40 steps, quiet): median 9.0–10.8 ms/step, steps 2–5 13–18 ms,
      step 1 ~30 ms (fresh factorization + 2 rounds); 200 steps: median 9.4 / 12.1, no
      refactor stalls; GC pauses ~100 ms still occur 1–2× per 200 steps.
      KMP_BLOCKTIME=0: inconclusive (noise) — not applied.
- [x] E6 refactor at the START of the next step once a round shows staleness (`sys.stale`,
      `_stale`), instead of mid-step after a wasted round (`_refactor_now` keeps an emergency
      path: Pardiso ρ > 0.5, UMFPACK > 5 rounds). Pardiso refactorization in situ ~10–12 ms;
      KMP_BLOCKTIME and MKL thread count make no difference; no iparm variant factorizes
      faster (single precision: same numeric time, only the solve is faster, residual 1e-5).
- Timings 2026-10-14, quiet machine, `instrument_run.jl`: median force 9.0–9.1 ms, velocity
  10.3–11.8 ms; worst steady step 18–21 ms (velocity refactor steps late in 200-step runs);
  step 1 ~33 ms instrumented, of which ~10 ms is first-use compilation of the `@debug` line in
  `_refactor!` (warm-up runs have no debug logger) — ~23 ms real; per-call preamble ~17 ms.
- [x] Full suite 2026-10-14: **97,982/97,982 in 2m58s** (fast + slow; squeeze_system now 70
      tests incl. rotated-view FD and a forced-UMFPACK run).
- Final timings 2026-10-14 (quiet, instrumented; step 1 includes ~10 ms debug-compile artifact):
  | mode (steps)            | median | worst steady | steps 2–3 |
  |-------------------------|--------|--------------|-----------|
  | force, constant (40)    | 9.0    | 18.4         | 15–18     |
  | velocity, constant (40) | 10.0   | 18.2         | 16–22     |
  | force, bulk (40)        | 11.9   | ~19 (+ a 62 ms GC in one run) | 18–20 |
  | velocity, bulk (40)     | 9.0    | ~19 (+ a 63 ms GC in one run) | 16–20 |
  | force, constant (200)   | 8.9    | 21.3         |           |
  | velocity, constant (200)| 11.8   | 22.1         |           |
  A 1-step `simulate` (preamble + step 1) takes 22.4 ms; the preamble itself is ~0.25 ms.
- Tried and rejected 2026-10-14: `GC.gc(false)` before the loop (it costs 60–80 ms itself —
  the young generation holds the previous calls' discarded outputs); KMP_BLOCKTIME; Pardiso
  iparm variants and single precision (no faster numeric factorization).
- OPEN / ideas for the user: (1) GC pauses — `simulate` returns ~0.35 MB/step of history the
  fit mostly discards; an option to skip unused outputs (the fit uses only output, gradList,
  borderPts2DList) would make collections rarer; (2) steps 2–3 cost 3–4 rounds because the
  extrapolation has little history; (3) step 1 (~22 ms) is the fresh factorization kept for
  call-to-call reproducibility; (4) `sensitivities=false` in prediction-only calls of
  `experiments/paper_1/pipeline.jl`; (5) Pardiso/MKL adds ~730 MB of artifacts on disk (corrected from ~300), x86-64 Linux/Windows only.
- (superseded) NEXT (in order): run full `Pkg.test()`; re-time with `claud_docs/realtime/instrument.py`
  (markers updated for the paired solve; `instrument_run.jl` prints per-step ms/rounds);
  cut `_get_2D_data` allocations (E3) and parallelize views (E4) to remove GC pauses; check
  step 1 (~30 ms); add a test that forces the UMFPACK backend (`_Umfpack(nothing, nothing)`)
  since CI only exercises Pardiso; update README Performance paragraph (still says
  UMFPACK-style numbers) and the CLAUDE.md real-time section; tell the user Pardiso/MKL is
  x86-64 only (fallback elsewhere) and ~300 MB of artifacts.

## Timing log (benchmark case, ms/step, warm second run)

| date       | change                                   | constant | bulk |
|------------|------------------------------------------|----------|------|
| baseline   | HEAD 4578f67                             | 586      | 560  |
| 2026-10-13 | fixed pattern + lagged LU, tol 1e-10, extrap | 105  | 114  |
| 2026-10-13 | B2 hoist                                  | 93       |      |
| 2026-10-13 | B7 one sens solve + B8 quadratic guess    | 52       | 64   |
| 2026-10-13 | B9 cached SqueezeSystem                     | 41       | 51   |
| 2026-10-13 | B11 tol 1e-9 / 1e-8 (steady 28.0 + 181/call) | 34    | 43   |
| 2026-10-13 | B13–B18 (steady 16.9 / 16.1, velocity 17.5 / 17.9) | 22 | 23 |
| 2026-10-13 | final, under browser load (load avg 5–8): steady force 16.5–22.8, velocity 17.4–25.3 (constant) | 22–28 | 23–28 |

(ms/step here = total `simulate` time / 40, so it includes the per-call refactor; the
`step_cost.jl` "per step" number is the steady-state slope.)

### Phase F — projection out of `simulate` (user, 2026-10-14: "moving the projection out of simulate sounds good"; Pardiso kept with the UMFPACK fallback)
Why: the next project steps need no camera; the projection costs ~0.5 ms/step and ~45 % of the
loop's allocations, and never feeds back into the solve (A/B: 3D results bit-identical).
Design: `simulate` returns the 3D results only — `output, displacement, surface_pts_3D, pos3D,
velocity, pressure, gradList_3d` (7, was 10). New exported `project_contours(mdl, pos3D,
conditions; grad3D=nothing)` (PostProcess.jl) returns `borderPts2DList, gradList, pos2D` (gradList
empty without grad3D), threads over time steps, owns `WRITECONTOUR` and the 2D half of `ANIMATE`
(simulate keeps the 3D GIFs). First frame: projected with grad3D[1] = 0 → 2×nb×2 zeros (was a
3-row zeros array; every consumer reads rows 1:2 only).
- [x] F0 references from the old code: scratchpad `proj_ref_save.jl` → `proj_ref/` (sims: 3 modes ×
      sensitivities on/off, 3 views, 10 steps; fits gn/lm/tikh from test/tikhonov_test.jl setup;
      `write_sim_data` files with WRITECONTOUR).
- [x] F1 `project_contours` + export + docstring — in stokes_solver.jl after `simulate`
      (PostProcess.jl is included before types.jl, so `Conditions` is not defined there).
- [x] F2 `simulate` 3D-only, docstring.
- [x] F3 smearOptimize.jl: `_simulate_contours` helper at the 8 sites.
- [x] F4 run_example.jl `write_sim_data` (projects without sensitivities; avg_time includes
      the projection as before); window_estimation.jl `predict_window` docstring.
- [x] F5 tests: squeeze_system (indices, rotated-view test, new testset "contours projected
      after the solve": 5 tests), stokes_optimization, tikhonov.
- [x] F6 experiments/paper_1/pipeline.jl: 7 simulate + 2 predict_window sites (none used
      gradList, so they project without sensitivities); parses.
- [x] F7 docs: README quick start, docs/coding_guide.md (file map now also lists
      squeeze_system.jl and sparse_lu.jl), claud_docs.
- [x] F8 verified: scratchpad `proj_ref_check.jl` → ALL IDENTICAL (3D outputs, contours,
      contour sensitivities for steps ≥ 1, initial ones zero, in force/velocity/bulk with and
      without sensitivities; gn/lm/tikh fits — every iterate, cost, H, λ, stored contours —
      compared with `isequal`; `write_sim_data` files incl. contour_data identical except
      avg_time.csv). Full suite 97,987/97,987 in 3m39s. `project_contours` over 41 states
      (3 views, with sensitivities) ~2.3 ms per call.
Finding while reading: `contour_cost` uses only the FIRST view of each frame (`_as_2xN`/`_as_dudθ`
take x[1]); the paper pipeline fits one view per call (`viewing_angles=[z_angle]`, ground truth
at 0/5/15/30/45°), so the back_project gradient bug hit every per-view fit at 5–45°.

### Status after Phase F — solver frozen for the next phase (user, 2026-10-14)
"For the next phase the solver would stay the same, at least for now": a prescribed list of
forces as in paper 1, fitting as a post-processing step. This will change later (online use),
not now. Parked until then, NOT to be started without the user:
- step-by-step interface (`step!(state, force)`, buffers allocated once): would remove the
  per-call step-1 refactorization (~22 ms), the slow steps 2–3 and the history-driven GC pauses;
- model-order reduction / lookup for many rollouts per control step (shape history depends on
  (η, β) only through β/η under constant force and parameters; η sets the speed).
In batch use only total time matters (~0.5 s per 50-step simulate on the benchmark mesh).
Paper impact of the back_project bug, checked on disk 2026-10-14: 470 result trees are view_1
(0°, unaffected); the only rotated-view fits are one synthetic cube run (bulk viscosity,
Hex2_16 gt, replicate 4, Hex2_6, simtime 5, dt 0.1, gn, view_2–view_5, July 2026), and no
paper_1 script reads anything but view_1 by default.

### 50 Hz status after Phase F (2026-10-14 10:16, scratchpad `step_times.jl`: lap timer per loop step, 5 calls × 40 steps, 2 × 200)
| mode                         | step 1 | steps 2–3 ≤ | steady median | p95  | steps > 20 ms |
|------------------------------|--------|-------------|---------------|------|---------------|
| force, constant              | 22.0   | 23.3        | 8.8           | 14.5 | 10/200        |
| velocity, constant           | 21.3   | 28.6        | 10.3          | 16.1 | 9/200         |
| force, bulk viscosity        | 19.4   | 33.4        | 11.9          | 18.9 | 16/200 (126, 70 ms GC) |
| velocity, bulk viscosity     | 20.3   | 26.2        | 9.3           | 16.1 | 9/200         |
| force, constant, no sens     | 18.3   | 16.1        | 6.7           | 10.9 | 0/200         |
| force, constant, 200 steps   | 23.8   | 23.0        | 8.9           | 15.1 | 8/400 (101 ms GC) |
| velocity, constant, 200 steps| 26.5   | 20.9        | 12.0          | 19.6 | 20/400 (78 ms GC) |
A run two minutes earlier on the same code was 1.5–2× slower across the board (medians
14.8–22.0) — machine state; benchmark only when the machine is idle and repeat.

### Where a step goes now + further options (2026-10-14 10:20, scratchpad `breakdown.jl`, load ~2)
Steady step without refactorization: total 14.9 ms = assembly 3.8 + lagged solve 10.7 + rest 0.4.
Rounds per steady step: 1 round in 160/180, 2 in 15, 0 in 5. A 0-round step costs 6.2 ms in
total, a 1-round step 14.8. Isolated building blocks: assembly 3.1, refactorization 12.5, one
solve 1.55, pair solve 2.57, residual 0.15 ms; 2 RHS in ONE Pardiso call 2.96 ms (slower than
2 calls, same result to 3.5e-16 — not worth it in isolation).
**Finding:** one refinement round costs ~9 ms inside the step but ~3 ms in isolation (pair solve
+ 2 residuals). Suspects: the assembly evicting the factor from L3 (24 MB) so the solves stream
it from DRAM, or Julia threads and MKL OpenMP threads (KMP_BLOCKTIME 200 ms spin) contending.
Next experiments (cheap): `--threads=8`; KMP_BLOCKTIME=0 in situ; 2-RHS solve in situ; pin
threads. Potential: ~5 ms of a 10–15 ms step.
Other options, not started: warm-start each fit iteration from the previous iterate's z, w
(+ sensitivity correction) → more 0-round steps, ~1e-9 call-to-call differences; several fits
in parallel processes with fewer threads each (batch throughput; RAM-limited); `--gcthreads=8,1`
for shorter GC pauses; step-by-step API (parked); MOR / β/η lookup (constant parameters only).

### Commit split IN PROGRESS (user, 2026-10-14: "do the commits as 3") — stopped at usage limit
- [x] C1 `07594e9 fix: rotate contour sensitivities with the camera view` on branch
      `realtime-split` (from HEAD 4578f67): minimal fix in HEAD's back_project
      (`R = obj_trans[1:3,1:3] * _rotate_round_z(Matrix{Float64}(I, 3, 3), rot_angle)`).
      FD check on HEAD's solver: 45° view 1.4e-10 (0° unchanged). Fast test group NOT run yet
      (no HEAD test uses a non-zero view, so the change is a no-op for the suite).
- [ ] C2 fast solver = the pre-Phase-F state: sparse_lu.jl, squeeze_system.jl, final
      PostProcess.jl, Project.toml (Pardiso), runtests.jl, smearFEM.jl includes only, README perf
      paragraph only, coding_guide file-map lines for squeeze_system.jl/sparse_lu.jl only,
      stokes_solver.jl and squeeze_system_test.jl with the Phase F edits reversed (10 outputs,
      70 tests; keep the corrected "keeps no state between calls" comment). Full suite.
- [ ] C3 `refactor: move the camera projection out of simulate` (one line): everything else;
      the tree must equal the main working tree minus the user's own edits.
- [ ] Then in the main tree: `git reset realtime-split` (mixed; working tree untouched) →
      `git status` must show only forward_convergence.jl, single_simulation.jl and the deleted
      `,` (the user's, not to be committed); `git branch -d realtime-split`.
Build in a worktree (copy the gitignored Manifest.toml and mesh_files caches into it).

### Commit split DONE (2026-10-14) — on `features/realtime_solver`, not pushed
- `07594e9 fix: rotate contour sensitivities with the camera view` — minimal fix in HEAD's
  back_project; FD check 1.4e-10 at 45°; fast group 97,879/97,879 (old solver, 8.5 min).
- `409421a perf: real-time simulate — fixed-pattern system, reused factorization, MKL Pardiso` —
  the pre-move state, rebuilt from the session transcript's file reads; reproduces the saved
  pre-move references bit-identically (sims, fits, written files); full suite 97,982/97,982.
- `1e2938e refactor: move the camera projection out of simulate` — tree equals the working tree
  except the user's own edits; full suite 97,987/97,987.
Still uncommitted (the user's, left alone): forward_convergence.jl, single_simulation.jl, the
deleted file `,`. The "Commit split IN PROGRESS" section above is superseded.
