# Initial FE mesh from φ_s — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** `mesh_from_field(phi, grid; table_z, elem_size, basis_order)` builds the paper-2 initial tet
mesh from perception's fused field, with the table faces as `IEN_bottom`.

**Architecture:** close the TSDF (seal below the table, fill NaN by connected components, cut at the
table) → marching cubes → cap/free triangles in an MSH 2.2 file → Gmsh `.geo` template through the
existing `_get_mesh_data`/`_run_gmsh` → `MeshgridObserved`. A Julia experiment plus a smearPerception
tool compare the result with the Carreau/5 GT mesh.

**Tech Stack:** Julia 1.13, Gmsh (CLI 4.12 locally, gmsh_jll API fallback), MarchingCubes.jl 0.1.11,
ImageMorphology.jl 0.4.7; Python (env `smear_ros`): Open3D 0.19, pytest.

**Spec:** `claud_docs/specs/2026-10-30-mesh-from-field-design.md` (revised 2026-10-31 after the prototype).

## Global Constraints

- Lengths in mm; world frame; `table_z` default `0.0`.
- Tets only (`:Tet`), `basis_order` 1 or 2; `Mesh.SecondOrderLinear = 1`.
- All Gmsh calls through `_get_mesh_data`/`_run_gmsh` (they hold `GMSH_LOCK`); no new `gmsh.*` calls.
- Style (`claud_docs/CLAUDE.md`): docstrings with `# Arguments`/`# Returns` on exported functions,
  `@info`/`@warn` not `println` in `src/`, comments only for the why, keyword args.
- Julia commands need the sandbox off (`~/.julia` is read-only inside it). Run Julia as
  `julia --project=. ...` from `smearFEM.jl/`.
- Commit locally on `features/volume_update`; do not push code (owner reviews first). Never stage
  `claud_docs/` together with code in one commit.

## Review Focus

1. Table on a voxel-centre plane (the rig: grid bottom = table = 0) → exact zeros; the mesh must still
   build with bottom nodes at z = 0. Test in Task 3.
2. Spurious blobs of material (perception leaves spurious crossings up to 19 mm out) → only the
   largest piece is meshed, with a warning. Test in Task 1.
3. Material in the outermost voxel layer (specimen larger than the grid) → warning that the mesh is
   cut flat there. Test in Task 1.
4. A lip ring detached from the wall by a one-voxel gap → a separate piece of material, so it is
   dropped as spurious with a warning (not meshed, not an error). Test in Task 1.
5. Material that does not reach the table (wrong `table_z`) → `ArgumentError` before Gmsh, not a
   Gmsh topology error. Test in Task 2.

---

### Task 1: Close the field

**Files:**
- Create: `src/fem/meshes/mesh_observed.jl`
- Modify: `src/io/perception_link.jl:254-259` (`sample_field` uses the new `_grid_centres`)
- Modify: `src/smearFEM.jl` (include after line 49 `include("io/perception_link.jl")`, with the comment
  `# after io/perception_link.jl: uses PerceptionGrid`)
- Modify: `test/runtests.jl` (include `mesh_observed_test.jl` after `linear_elasticity_test.jl`)
- Test: `test/mesh_observed_test.jl`
- Commit also: `Project.toml` (the owner's dep additions, uncommitted)

**Interfaces:**
- Produces: `_grid_centres(grid::PerceptionGrid) -> NTuple{3,StepRangeLen}` (voxel centres per axis, mm),
  in `perception_link.jl`.
- Produces: `_close_field(phi::AbstractArray{<:Real,3}, grid::PerceptionGrid, table_z::Real)
  -> (field::Array{Float64,3}, coords::NTuple{3,Vector{Float64}})`: padded by one voxel per side,
  no NaN, no exact zeros, negative inside.
- Produces (test helper, top of the test file): `tsdf_cylinder(; r=30.0, h=40.0, tz=4.685, Δ=1.0,
  μ=3.0, depth=7.0, nz=53, extra=nothing) -> (phi::Array{Float32,3}, grid::PerceptionGrid)` and
  `cylinder_sd(x, y, z; r=30.0, h=40.0, tz=4.685)`.

- [ ] **Step 1: Write the test helper and failing tests**

```julia
using Test, LinearAlgebra
using smearFEM
using smearFEM: _close_field, _grid_centres

cylinder_sd(x, y, z; r=30.0, h=40.0, tz=4.685) =
    let d = (hypot(x, y) - r, max(tz - z, z - (tz + h)))
        norm(max.(d, 0)) + min(max(d...), 0)
    end

# A cylinder standing on the table at z = tz as perception leaves it: NaN more than μ outside or
# `depth` inside (the distance cap) and below the table. `extra(x, y, z)` may return a value to
# overwrite the voxel with (NaN for a hole, a negative value for a blob), or `nothing`.
function tsdf_cylinder(; r=30.0, h=40.0, tz=4.685, Δ=1.0, μ=3.0, depth=7.0, nz=53, extra=nothing)
    n = round(Int, 36 / Δ)
    grid = PerceptionGrid([0.0, 0.0, 0.0], Δ, [n, n, 0], [2n + 1, 2n + 1, nz])
    xs, ys, zs = _grid_centres(grid)
    phi = Array{Float32}(undef, length(xs), length(ys), length(zs))
    for (i, x) in enumerate(xs), (j, y) in enumerate(ys), (k, z) in enumerate(zs)
        d = cylinder_sd(x, y, z; r, h, tz)
        v = (z < tz || d > μ || d < -depth) ? NaN : d
        e = isnothing(extra) ? nothing : extra(x, y, z)
        phi[i, j, k] = isnothing(e) ? v : e
    end
    return phi, grid
end
```

Testsets (`@testset "observed mesh" begin ... end` wrapping all):
- `"closing a TSDF-like cylinder"`: `F, (px, py, pz) = _close_field(tsdf_cylinder()..., 4.685)`;
  `@test size(F) == (75, 75, 55)`; `@test !any(isnan, F)`; `@test !any(iszero, F)`;
  `@test px[1] == -37.0 && pz[1] == -1.0`; for every voxel with `abs(cylinder_sd(x, y, z)) > 1.0`,
  `sign(F) == sign(cylinder_sd)` (`@test all(...)`).
- `"a hole through the side wall is not closed"`: hole = `extra=(x, y, z) -> (21 ≤ x ≤ 35 && abs(y) ≤ 3 &&
  20 ≤ z ≤ 30) ? NaN : nothing`; `err = try _close_field(phi, grid, 4.685); nothing catch e; e end`;
  `@test err isa ErrorException && occursin("not closed", err.msg)`.
- `"a detached lip ring is dropped"`: lip = `extra=(x, y, z) -> (32 ≤ hypot(x, y) ≤ 33 && 41 ≤ z ≤ 43) ?
  -0.5 : (30 < hypot(x, y) < 32 && 41 ≤ z ≤ 43) ? 0.5 : nothing` (a negative ring outside the wall, a
  positive gap between); `@test_logs (:warn, r"spurious") match_mode=:any`, and the voxel nearest
  (32.5, 0, 42) is positive afterwards.
- `"a spurious blob is dropped"`: blob = `extra=(x, y, z) -> norm((x, y, z - 49)) ≤ 1.5 ? -0.5 : nothing`;
  `@test_logs (:warn, r"spurious") match_mode=:any F, c = _close_field(...)`; then the voxel nearest
  (0, 0, 49) is positive.
- `"material at the grid side"`: `tsdf_cylinder(r=36.5)` (at r = 36 the outermost layer only touches
  φ = 0); `@test_logs (:warn, r"grid boundary")
  match_mode=:any _close_field(...)`.
- `"bad input"`: `@test_throws ArgumentError _close_field(phi, grid, -5.0)` (below the grid);
  `@test_throws ArgumentError _close_field(fill(NaN32, size(phi)), grid, 4.685)` (no material).

- [ ] **Step 2: Run, verify they fail**

Run: `julia --project=. -e 'include("test/mesh_observed_test.jl")'`
Expected: error, `_close_field` / `_grid_centres` not defined.

- [ ] **Step 3: Implement `_grid_centres` and `_close_field`**

`_grid_centres(grid) = ntuple(a -> grid.x0[a] .+ grid.delta .* ((0:grid.shape[a]-1) .- grid.offset[a]), 3)`;
`sample_field` uses it (`scale(..., _grid_centres(grid)...)`).

`mesh_observed.jl` header: `using MarchingCubes: MC, march` and
`using ImageMorphology: label_components`, then a short comment on what the file does (paper-2 initial
mesh from φ_s). `_close_field`, in this order (the order is the design):

```julia
Δ = grid.delta
xs, ys, zs = _grid_centres(grid)
zs[1] <= table_z < zs[end] || throw(ArgumentError(...))         # table inside the grid's z range
px, py, pz = (vcat(c[1] - Δ, collect(c), c[end] + Δ) for c in (xs, ys, zs))
F = fill(NaN, length(px), length(py), length(pz)); F[2:end-1, 2:end-1, 2:end-1] .= phi
plane = reshape(table_z .- pz, 1, 1, :)                          # signed distance to the table
F .= ifelse.(plane .>= 0, plane, F)                               # 1. seal at and below the table
# 2. NaN components touching the padded border are outside; the rest are enclosed: inside
# 3. F .= max.(F, plane)
# 4. keep the largest component of F .< 0 (others := +Δ, @warn "… spurious …" with their voxel count);
#    no component at all → ArgumentError("no material above the table")
# 5. leak check (`_ringed_slices`) → error("observed surface is not closed: … z = a … b mm …
#    unobserved patch, or table_z too low"); @warn "… grid boundary …" if any F < 0 in the unpadded
#    outermost layer of the sides (i, j) or the top (k)
# 6. F[F .== 0] .= 1e-6Δ;  @info voxels filled inside / outside
```

`_ringed_slices(F, pz, table_z) -> Vector{Float64}`: z of every slice above the table where
`label_components(F[:, :, k] .>= 0)` has a component with no voxel on the slice border.
Message values rounded to 0.1 mm.

- [ ] **Step 4: Run, verify they pass**

Run: `julia --project=. -e 'include("test/mesh_observed_test.jl")'` → all pass.
Also `julia --project=. -e 'include("test/perception_link_test.jl")'` → 50/50 (sample_field refactor).

- [ ] **Step 5: Commit**

```bash
git add Project.toml src/fem/meshes/mesh_observed.jl src/io/perception_link.jl src/smearFEM.jl test/mesh_observed_test.jl test/runtests.jl
git commit -m "[ADD] Close perception's fused field for meshing: seal at the table, fill unobserved voxels, leak check"
```

### Task 2: Surface by marching cubes

**Files:**
- Modify: `src/fem/meshes/mesh_observed.jl`
- Test: `test/mesh_observed_test.jl`

**Interfaces:**
- Consumes: `_close_field` (Task 1).
- Produces: `_field_surface(field::Array{Float64,3}, coords::NTuple{3,Vector{Float64}}, table_z::Real,
  delta::Real) -> (V::Matrix{Float64} (3, nv), T::Matrix{Int} (3, nt), on_cap::BitVector (nt))`.

- [ ] **Step 1: Failing tests** — testset `"surface of the closed field"`:
  `V, T, cap = _field_surface(_close_field(phi, grid, 4.685)..., 4.685, 1.0)`;
  - closed: every undirected edge `minmax(a, b)` of every triangle appears exactly twice;
  - `vol = abs(sum(det(V[:, T[:, t]]) for t in axes(T, 2))) / 6`; `@test vol ≈ π * 30^2 * 40 rtol=0.005`;
  - `@test all(abs.(V[3, unique(vec(T[:, cap]))] .- 4.685) .< 1e-3)`; `@test count(cap) > 0`;
  - no zero-area triangle: `all(norm(cross(...)) > 0)`.
- Testset `"floating material"`: a fully defined sphere field (`phi[i,j,k] = norm((x, y, z - 26)) - 10`,
  no NaN) with `table_z = 4.685` → `@test_throws ArgumentError _field_surface(...)` (no cap triangle).

- [ ] **Step 2: Run, verify they fail** (`_field_surface` not defined).

- [ ] **Step 3: Implement `_field_surface`**: `mc = MC(field, Int; x=coords[1], y=coords[2],
  z=coords[3]); march(mc, 0.0)`; `V = reduce(hcat, mc.vertices)`, `T = reduce(hcat, mc.triangles)`;
  `on_cap` = all three vertices within `1e-3 * delta` of `table_z`; no cap triangle →
  `ArgumentError("the material does not reach the table at z = …; check table_z")`.

- [ ] **Step 4: Run, verify pass.**

- [ ] **Step 5: Commit** — `[ADD] Surface of the closed field by marching cubes; the cap is the triangles on the table`

### Task 3: Gmsh volume, `MeshgridObserved`, `mesh_from_field`

**Files:**
- Create: `mesh_files/templates/observed_tet.geo`
- Modify: `src/fem/meshes/mesh_types.jl` (append `MeshgridObserved`)
- Modify: `src/fem/meshes/mesh_observed.jl` (`_write_surface_msh`, `mesh_from_field` ×2)
- Modify: `src/smearFEM.jl` (export `MeshgridObserved` on the line-17 export list; new export line
  `export mesh_from_field # fem/meshes/mesh_observed.jl`)
- Modify: `src/solver/volume_update.jl:120` docstring: "`IEN_bottom` (cylinder, cuboid or observed)"
- Test: `test/mesh_observed_test.jl`

**Interfaces:**
- Consumes: `_close_field`, `_field_surface`, `_get_mesh_data(path; params, template_path,
  mesh_order, mesh_dim=3, body_group, body_dim=3, face_groups, node_set_groups)`.
- Produces: `MeshgridObserved(; NodeList, IEN, IEN_bottom, IEN_free, ID, volume_element_shape=:Tet,
  surface_element_shape=_surface_element_shape(volume_element_shape), basis_order=2, nNodes=0, ne=0,
  bottom_nodes=Int[], free_nodes=Int[], effect_elem_sze=nothing)`; `initial_state = copy(NodeList)`.
- Produces: `mesh_from_field(phi::AbstractArray{<:Real,3}, grid::PerceptionGrid; table_z::Real=0.0,
  elem_size::Real, basis_order::Int=2, out_dir::String=mktempdir()) -> MeshgridObserved` and
  `mesh_from_field(step::PerceptionStep, grid::PerceptionGrid; kwargs...)`.

Template (exact):

```
// Paper-2 initial mesh: volume inside the closed surface from smearFEM's `mesh_from_field`
// (surface.msh: entity 1 = cap on the table, 2 = free surface, shared nodes).
elem_size = 3.0;
Merge "surface.msh";
CreateTopology;
CreateGeometry;
Surface Loop(1) = {1, 2};
Volume(1) = {1};
Physical Surface("Bottom") = {1};
Physical Surface("Free") = {2};
Physical Volume("Volume") = {1};
Mesh.MeshSizeMax = elem_size;
Mesh.MeshSizeMin = elem_size / 4;
// mid-edge nodes on straight edges: curved onto the reparametrized surface they fold P2 tets
Mesh.SecondOrderLinear = 1;
```

- [ ] **Step 1: Failing tests**
  - `"observed mesh of a cylinder (Tet$order)" for order in (1, 2)`:
    `m = mesh_from_field(tsdf_cylinder()..., table_z=4.685, elem_size=6.0, basis_order=order)`;
    `@test m isa MeshgridObserved && m.volume_element_shape == :Tet && size(m.IEN, 1) == (4, 10)[order]`;
    corner-tet signed volumes `det([X[:, b] - X[:, a] X[:, c] - X[:, a] X[:, d] - X[:, a]]) / 6` all `> 0`
    and their sum `≈ π * 30^2 * 40 rtol=0.02`;
    `@test all(abs.(m.NodeList[3, m.bottom_nodes] .- 4.685) .< 1e-9)`;
    vertex nodes of the free faces (`unique(vec(m.IEN_free[1:3, :]))`) satisfy
    `abs(cylinder_sd(X[:, a]...)) < 1.0`;
    cap area (sum over `IEN_bottom` corner triangles) between `π * 29^2` and `π * 30^2`.
  - `"table on a voxel plane"`: `tsdf_cylinder(tz=0.0)` meshed with `table_z=0.0` → bottom nodes at `z == 0`
    within 1e-9.
  - `"Gmsh API path"`: `withenv("PATH" => "/nonexistent") do mesh_from_field(...; basis_order=1) end`
    → `MeshgridObserved` with non-empty `IEN_bottom` and `IEN_free`.
  - `"volume update on the observed mesh"` (P2 mesh from above): side vertex nodes with
    `abs(hypot(x, y) - 30) < 1` and `tz + 5 < z < tz + 35`, radial unit normals, `d = fill(0.5, n)`;
    `mdl = LinearElasticity(mesh=m, λ=λ, μ=μ)` with `λ, μ = lame_parameters(1.0, 0.3)`;
    `u = volume_update!(mdl, nodes, normals, d; β_table=1e-3)`;
    `@test all(abs(dot(u[:, a], n) - 0.5) < 1e-8 ...)` and `@test maximum(abs.(u[3, m.bottom_nodes])) < 1e-10`.

- [ ] **Step 2: Run, verify they fail** (`mesh_from_field` not defined).

- [ ] **Step 3: Implement.** `_write_surface_msh(path, V, T, on_cap)`: MSH 2.2 ASCII (`$MeshFormat`
  `2.2 0 8`, `$Nodes` with 1-based ids, `$Elements` type 2 with 2 tags = entity, entity: 1 cap, 2 free).
  `mesh_from_field`: `elem_size > 0` and `basis_order ∈ (1, 2)` else `ArgumentError`; `Float64` copy of
  phi happens in `_close_field`; write `out_dir/surface.msh`; `rm(joinpath(out_dir, "observed.msh"),
  force=true)`; `_get_mesh_data(joinpath(out_dir, "observed.msh"); params=Dict("elem_size" =>
  Float64(elem_size)), template_path=<repo>/mesh_files/templates/observed_tet.geo (same
  `dirname(dirname(dirname(@__DIR__)))` rule as `mesh_gmsh.jl`), mesh_order=basis_order,
  body_group="Volume", face_groups=["Bottom", "Free"], node_set_groups=["Bottom", "Free"])`; `ID` as in
  `_get_gmsh_cylinder` (ndof = 3). Docstrings on both methods and the type (template of
  `MeshgridCylinder`'s docstring; say the cap is ~one voxel smaller in radius than the base).

- [ ] **Step 4: Run, verify pass**; then the whole suite:
  `julia --project=. -e 'using Pkg; Pkg.test()'` → previous 225 + new, all pass (Aqua included).

- [ ] **Step 5: Commit** — `[ADD] mesh_from_field: paper-2 initial tet mesh from the fused field (Gmsh), MeshgridObserved`

### Task 4: Ground-truth comparison tools

**Files:**
- Create: `smearFEM.jl/experiments/mesh_init/gt_comparison.jl`
- Create: `smearPerception/tools/analysis/compare_mesh_gt.py`
- Test: `smearPerception/tests/test_compare_mesh_gt.py`

**Interfaces:**
- Consumes: `read_perception_steps`, `mesh_from_field` (Task 3); smearPerception `sdf.load_mesh(path,
  scale=MM_TO_M)`, `sdf.signed_distance(scene, points_m)`.
- Produces: `gt_comparison.jl`: `export_boundary_obj(mesh::MeshgridObserved, path)` (corner nodes of
  `IEN_bottom` and `IEN_free`, re-indexed, OBJ `v`/`f` lines, mm) and
  `main(h5, out_dir; step=0, table_z=4.685, elem_size=3.0, basis_order=2)` gated by
  `if abspath(PROGRAM_FILE) == @__FILE__` (args: `<recon.h5> <out_dir> [elem_size]`), as
  `experiments/elasticity/benchmark.jl`; prints node/tet counts, tet volume, the OBJ path and the
  `compare_mesh_gt.py` command.
- Produces: `compare_mesh_gt.py`: `compare(observed_obj, sim_dir, step=0) -> dict` with keys
  `obs_to_gt`, `gt_to_obs` (signed distances, mm, numpy arrays), `volume_obs`, `volume_gt` (mm³);
  CLI `python tools/analysis/compare_mesh_gt.py <observed.obj> <sim data dir> [--step 0]` prints
  `observed -> GT [mm]: median ±x, p95 |.| y, max |.| z (n = …)`, the same for `GT -> observed`, and
  `volume [mm^3]: observed …, GT …, difference …%`.

- [ ] **Step 1: Failing pytest** (`tests/test_compare_mesh_gt.py`): a fake sim dir in `tmp_path` with
  `sim_params.json` `{"cyl_T": [90.0, 0.0, 4.685], "h0": 40.0}` and `Blender_obj/play_dough_000.obj` =
  `o3d.geometry.TriangleMesh.create_cylinder(radius=25, height=40, resolution=96, split=8)` shifted by
  +20 in z (mesh frame, bottom at 0). Tests: (a) observed = the same cylinder shifted by `cyl_T` →
  `max |obs_to_gt| < 1e-6`, `max |gt_to_obs| < 1e-6`, `volume_obs ≈ volume_gt` (rel 1e-9);
  (b) observed shifted 1 mm further up → `max |obs_to_gt| ≈ 1.0` (abs 1e-6), volumes equal.

- [ ] **Step 2: Run, verify it fails**: `/opt/anaconda3/envs/smear_ros/bin/python -m pytest
  tests/test_compare_mesh_gt.py -q` → import error.

- [ ] **Step 3: Implement `compare`**: GT vertices from `o3d.io.read_triangle_mesh`; distances via
  `load_mesh`/`signed_distance` in metres, converted to mm; volume = |Σ a·(b×c)|/6 after
  `merge_close_vertices(1e-9)` and `orient_triangles()` (our OBJ's two Gmsh entities may face opposite
  ways). Then `gt_comparison.jl`.

- [ ] **Step 4: Verify.** pytest passes. Julia end to end on the analytic field (no reconstruction on
  disk): in the REPL, `include("experiments/mesh_init/gt_comparison.jl")` (gate holds),
  `export_boundary_obj(mesh_from_field(<tsdf_cylinder r=25, tz=4.685+... >), …)`, fake sim dir with an
  r = 25 Open3D cylinder → `compare_mesh_gt.py` reports sub-voxel distances and a volume difference near
  −1 %. Report the numbers.

- [ ] **Step 5: Commit** (two repos, two commits; neither pushed):
  smearFEM.jl `[ADD] Experiment: observed mesh vs the Carreau/5 GT mesh`;
  smearPerception `[ADD] compare_mesh_gt.py: FE mesh vs GT mesh, both directions and volume`.

### Task 5: Docs

- [ ] smearFEM.jl `claud_docs/CLAUDE.md`: section "Initial mesh from φ_s (2026-10-31)": API, pipeline in
  one paragraph, the four prototype findings, test count; `claud_docs/memory.md` entry with the
  measurements. smearPerception `claud_docs/CLAUDE.md`/`memory.md`: the tool. `../claude_docs/handoff.md`
  §9 entry (concepts only) and fix the stale "plate hides the top face" line in §9 2026-10-21 GT notes
  only if it is a statement of current fact (leave dated entries as they are; add a correction line).
- [ ] Commit the notes in each repo (`[DOCS] ...`); do not push (the code commits precede them on the
  same branch; the owner reviews first).
