# Data migration: `view_*` level in `experiment_mesh_conv`

Status: **pending — no data has been moved or deleted yet.**
Written 2026-07-17, against the tree as it stood that day.

Two code fixes changed where `optimize` writes constant-viscosity results. The
code is now self-consistent, but the data on disk was written by three different
past versions and does not match it. This note records what is stale, what is
valid, and where the valid data has to move.

Everything below concerns one tree:

```
$DATA/experiments/sim_data/convergence_analysis/stokes_convergence/experiment_mesh_conv/
```

where `$DATA` = `~/SMEAR-PhD/SMEAR-DataFiles/Data`. Readers of this tree are
`experiments/convergence_analysis/inverse_convergence.jl` (which reads
`cylinder/3`, and `Hex2_6` for the dt sweep).

## Target layout

```
experiment_mesh_conv/<geom>/<dir>/<elem>/dt_<Δt>/view_<i>/data/
├── experiment_parameters.jld2
├── est_h.csv
├── η.csv  β.csv  gt_h.csv  cost_iter.csv  stats.jld2
└── sim_data/view_<i>/…        <- keep; written by write_2d_data, not the splice
```

Note `sim_data/view_<i>/` is **correct and must stay**. Only the `data/view_<i>/`
level directly holding `η.csv` etc. is redundant.

## What is wrong, and why

### 1. Redundant `data/view_<i>/` level — 32 directories

`optimize` used to write `datapath(exp_path, "view_$i", …)` while `exp_path`
already ended in `view_<i>`, producing `dt_0.1/view_1/data/view_1/η.csv`. The
writer no longer does this. The inner level must be flattened up one.

The inner view name always matches its parent (`view_2/data/view_2`), verified —
so flattening can never collide between different views.

### 2. Inverted `view_*/dt_*` ordering — 32 `dt_` directories

`optimize_sim`'s `:conv_exp_mesh` branch omitted the trailing `window` component
from `filepath_res`, so the `view_$i` splice landed above `dt_` instead of below:

```
<elem>/view_1/dt_0.1/…     <- wrong
<elem>/dt_0.1/view_1/…     <- right
```

Fixed in code. The 32 already-written directories are still inverted:

| elem dir | inverted `dt_` count |
|---|---|
| `cylinder/3/hex2_12` | 3 |
| `cylinder/6/Hex2_4` | 1 |
| `cylinder/6/Hex2_16` | 3 |
| `cylinder/6/Hex2_10`, `Hex2_12`, `Hex2_14`, `Hex2_6`, `Hex2_8` | 5 each |

`cylinder/6` is inverted in its entirety. `cylinder/3` is correct except `hex2_12`.

Inverted leaves are also doubly wrong: they contain `dt_0.1/data/view_1/` too.

### 3. Stale flat copies — 4 leaves

Four correct-shape leaves have **both** `data/η.csv` (2026-07-09) and
`data/view_1/η.csv` (2026-07-11). The July-9 flat files are leftovers from a
still older layout. Content was byte-identical in every file sampled, but the
July-11 copy is authoritative and must win.

Since the code now reads the flat path, these four currently resolve to the
**stale** copy. Flattening (step 1) overwrites them correctly.

## Census

Counts verified 2026-07-17.

| shape | `data/` dirs | holds `η` data | partial (no `η`) |
|---|---|---|---|
| correct `dt_*/view_*` | 21 | 18 (4 also have stale flat copies) | 3 |
| inverted `view_*/dt_*` | 32 | 2 | 30 |

- `dt_` dirs in correct position: 23 · inverted: 32
- `data/view_*` dirs to flatten: 32 (spans both shapes)

**Leaves holding real `η` data (18)** — these are the ones that must survive:

```
cylinder/3/Hex2_2/dt_0.1/view_1   view_2   view_3   view_4
cylinder/3/Hex2_4/dt_0.1/view_1
cylinder/3/Hex2_6/dt_0.1/view_1   dt_0.2   dt_0.5   dt_0.7   dt_0.8   dt_1.0  (view_1)
cylinder/3/Hex2_8/dt_0.1/view_1   dt_0.3/view_1
cylinder/3/Hex2_10/dt_0.1/view_1
cylinder/3/Hex2_16/dt_0.1/view_1
cylinder/3/hex2_12/dt_0.1/view_1
cylinder/3/hex2_14/dt_0.1/view_1
cylinder/6/Hex2_2/dt_0.1/view_1
```

## Migration

Run in order. **Dry-run first** — every step prints before it moves.

### Step 1 — flatten `data/view_<i>/` up into `data/`

```bash
cd "$DATA/experiments/sim_data/convergence_analysis/stokes_convergence/experiment_mesh_conv"
for v in $(find . -type d -path '*/data/view_*'); do
  d=$(dirname "$v")
  echo "mv $v/* -> $d/"
  # mv -f "$v"/* "$d"/ && rmdir "$v"
done
```

Uncomment the `mv` to execute. This overwrites the 4 stale July-9 flat copies
with their July-11 counterparts, which is the intent.

### Step 2 — un-invert `view_*/dt_*` → `dt_*/view_*`

**Resolve the two collisions below first.** 30 of the 32 move cleanly; 2 have a
target that already exists. A bare `mv` into an existing directory would nest the
source inside it (`dt_0.1/view_1/dt_0.1`) and silently corrupt both — hence the
`-d` guard.

```bash
for old in $(find . -type d -name 'dt_*' -path '*/view_*/dt_*'); do
  dt=$(basename "$old")                      # dt_0.1
  view=$(basename "$(dirname "$old")")       # view_1
  elem=$(dirname "$(dirname "$old")")        # …/hex2_12
  new="$elem/$dt/$view"
  if [ -d "$new" ]; then
    echo "SKIP (target exists, resolve by hand): $old -> $new"
    continue
  fi
  echo "mv $old -> $new"
  # mkdir -p "$(dirname "$new")" && mv "$old" "$new"
done
```

Then remove the emptied `view_*` shells:

```bash
find . -type d -name 'view_*' -empty -delete
```

#### Collision 1 — `cylinder/3/hex2_12/dt_0.1/view_1` — **needs a decision**

Both copies exist and **their contents differ** — two distinct runs six minutes
apart, not duplicates:

| copy | `η.csv` mtime |
|---|---|
| inverted `hex2_12/view_1/dt_0.1/…` | 2026-07-11 02:52 (newer) |
| correct `hex2_12/dt_0.1/view_1/…` | 2026-07-11 02:46 |

Picking by "newest wins" keeps the inverted copy and discards the one already in
the right place. **Do not automate this** — confirm which run is the one you want
before overwriting, then move it into `dt_0.1/view_1/` manually.

#### Collision 2 — `cylinder/6/Hex2_4/dt_0.1/view_1` — safe

Both hold only `experiment_parameters.jld2` and are byte-identical; the inverted
copy is older. Delete the inverted one:

```bash
# rm -rf cylinder/6/Hex2_4/view_1
```

### Step 3 — delete partial runs (destructive, do last)

33 `data/` dirs contain no `η.csv` under any shape — abandoned/partial runs
(30 of them inverted). Only delete after steps 1-2 confirm the 18 real leaves
survived. Identify them with:

```bash
for d in $(find . -type d -name data); do
  [ -e "$d/η.csv" ] || echo "PARTIAL: $d"
done
```

**Run this only after step 1.** Beforehand it reports 49, not 33 — the extra 16
are valid leaves whose `η.csv` still sits in `data/view_*/` and has not been
flattened yet. Deleting on the pre-migration list would destroy real data.

Note `cylinder/3/Hex2_2/dt_0.1/view_4` is real but incomplete (only
`experiment_parameters.jld2`); decide per-leaf rather than by rule.

### Step 4 — verify

```bash
find . -type d -path '*/data/view_*' | wc -l          # expect 0
find . -type d -name 'dt_*' -path '*/view_*/dt_*' | wc -l   # expect 0
find . -name 'η.csv' -path '*/dt_*/view_*/data/η.csv' | wc -l  # expect >= 18
```

Then run `experiments/convergence_analysis/inverse_convergence.jl`.

## Open items (decide before/while migrating)

- **Lowercase elem dirs.** `cylinder/3/hex2_12` and `cylinder/3/hex2_14` are
  lowercase; every other is `Hex2_*`. No uppercase twin exists, so renaming is
  collision-free — but confirm nothing indexes them by the lowercase name.
- **Multi-view leaves.** `cylinder/3/Hex2_2/dt_0.1` has `view_1`…`view_4` and
  `cylinder/3/Hex2_10/dt_0.1` has `view_1`,`view_2`. `optimize` now skips every
  view but the first (`if i != 1; continue` in `optimize`), so this data predates
  that guard and cannot be regenerated without lifting it. Keep or discard — but
  note `_single_leaf` in `test_opt_stokes.jl` errors on multi-leaf groups by
  design. It does not fire here, because `collect_experiment_groups` cannot walk
  this tree at all (no `simtime_*`/`noise_*` levels).

## Unrelated but adjacent: the old constant optimization tree

`$DATA/experiments/sim_data/optimization/Stokes/force/constant/` holds only
`Hex2_3.25/` and `Q2_16/` — pre-migration layout throughout (`Q2_6` element
naming, data at `noise_*/Results/<n>/data/`). `plot_results` expects
`Hex2_16.0/cylinder`, which **does not exist**.

`post_analysis_const` has been migrated to the new hierarchy but cannot be
verified against this data. It needs regenerating, not moving — the shape is too
far from the target to migrate mechanically. Out of scope for this note.
