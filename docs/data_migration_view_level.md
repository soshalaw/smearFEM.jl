# Data migration: `view_*` level in `experiment_mesh_conv`

Status: **closed — resolved by regeneration, not by migration.**
Written 2026-07-17. Re-verified against disk and rewritten 2026-08-26.

The migration described here was never run. Between 2026-07-17 and 2026-08-26 the
`experiment_mesh_conv` tree was regenerated from scratch: the directories holding
the defects no longer exist, and the ones that replaced them were written by the
corrected code. Nothing needs moving. This note is kept for the invariants in
[What caused it](#what-caused-it), which still constrain the writers.

Everything below concerns one tree:

```
$DATA/experiments/sim_data/convergence_analysis/stokes_convergence/experiment_mesh_conv/
```

where `$DATA` = `~/SMEAR-PhD/SMEAR-DataFiles/Data`. Its only reader is
`experiments/convergence_analysis/inverse_convergence.jl`.

## Current layout (verified 2026-08-26)

```
experiment_mesh_conv/cylinder/4/<elem>/dt_<Δt>/gn/view_1/data/
├── experiment_parameters.jld2
├── est_h.csv
├── η.csv  β.csv  gt_h.csv  cost_iter.csv  stats.jld2
└── sim_data/view_1/…        <- correct; written by write_2d_data, not by the splice
```

11 leaves, every one in that shape and every one holding `η.csv`:

| sweep | leaves |
|---|---|
| mesh (`dt_0.1`) | `Hex2_2`, `Hex2_4`, `Hex2_6`, `Hex2_8`, `Hex2_10`, `Hex2_16` |
| time step (`Hex2_6`) | `dt_0.1`, `dt_0.2`, `dt_0.5`, `dt_0.7`, `dt_0.8`, `dt_1.0` |

Plots land in `cylinder/4/post_analysis/`.

The `gn/` component is the optimizer that produced the run, added 2026-08-26 so
that `:gn`, `:lm` and `:gn_tikhonov` results cannot overwrite each other. It sits
between `dt_*` and `view_*`.

## What was wrong, and what became of it

All three defects were confined to `cylinder/3` and `cylinder/6`, which no longer
exist — the tree is `cylinder/4` only.

| # | Defect (2026-07-17) | State 2026-08-26 |
|---|---|---|
| 1 | Redundant `data/view_<i>/` level — 32 dirs | 0 remain |
| 2 | Inverted `view_*/dt_*` ordering — 32 `dt_` dirs | 0 remain |
| 3 | Stale flat July-9 copies shadowing July-11 data — 4 leaves | gone with `cylinder/3` |

The two collisions that step 2 flagged as needing a human decision
(`cylinder/3/hex2_12/dt_0.1/view_1`, where two runs six minutes apart disagreed,
and `cylinder/6/Hex2_4/dt_0.1/view_1`) are moot: neither path exists. **No
decision is outstanding.**

Verification commands, all returning 0 as of 2026-08-26:

```bash
cd "$DATA/experiments/sim_data/convergence_analysis/stokes_convergence/experiment_mesh_conv"
find . -type d -path '*/data/view_*'          | wc -l   # redundant level
find . -type d -name 'dt_*' -path '*/view_*/dt_*' | wc -l   # inverted ordering
for d in $(find . -type d -name data); do [ -e "$d/η.csv" ] || echo "PARTIAL: $d"; done
```

## What caused it

Worth keeping, because both invariants still bind the current writers.

**The `view_$i` splice is positional.** `optimize` builds each result directory as

```julia
exp_path = joinpath(dirname(filepath_res), "view_$i", basename(filepath_res))
```

so `view_$i` is inserted before whatever the *last* component of `filepath_res`
happens to be. Defect 2 was `optimize_sim`'s `:conv_exp_mesh` branch omitting the
trailing `window` component, which pushed the splice one level too high
(`<elem>/view_1/dt_0.1/` instead of `<elem>/dt_0.1/view_1/`). The same rule is why
the `<opt_method>` segment added on 2026-08-26 goes *before* `window` rather than
after it — appending it would put `view_*` above the method instead of below.

**`sim_data/view_<i>/` is not the splice.** It is written by `write_2d_data`,
which splices `view_$a` relative to its own `filepath` argument. It is correct and
must stay. Defect 1 was `optimize` additionally writing `datapath(exp_path,
"view_$i", …)` while `exp_path` already ended in `view_<i>`; the writer no longer
does this.

## Notes on readers

- `collect_experiment_groups` **cannot walk this tree** — it requires
  `simtime_*`/`noise_*` levels, which `experiment_mesh_conv` does not have.
  `inverse_convergence.jl` reads it by constructing paths directly, via the
  `_view_path` helper that accepts both the `gn/` layout and the pre-2026-08-26
  one.
- `inverse_convergence.jl`'s `main()` reads `cylinder/4`. The 2026-07-17 note said
  `cylinder/3`; that is stale.
- `optimize` no longer skips views after the first — the `if i != 1; continue`
  guard the original note referred to is gone, so a ground truth with several
  `z_angle_list` entries produces several `view_*` leaves again. `_single_leaf` in
  `test_opt_stokes.jl` errors on multi-leaf groups by design, but it is not
  reachable from this tree (see above). It *is* reachable from the optimization
  trees once a second optimizer is run over one experiment — pass
  `collect_experiment_groups(…; method=…)` there.

## Adjacent tree: constant optimization — also resolved

`$DATA/experiments/sim_data/optimization/Stokes/force/constant/` was recorded in
2026-07-17 as holding only `Hex2_3.25/` and `Q2_16/` in a pre-migration shape,
with `plot_results`' expected `Hex2_16/cylinder` missing entirely.

It now holds exactly `Hex2_16/cylinder/<1…6>` in the target hierarchy
(`<dir>/<elem>/simtime_*/noise_*/dt_*/gn/view_1`), with `post_analysis/`,
`post_analysis_time/` and `post_analysis_global/` present. `post_analysis_const`
can be verified against it.
