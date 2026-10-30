# Initial FE mesh from the fused field φ_s — design

Branch `features/volume_update`. Agreed with the owner 2026-10-30.

## Goal

Paper 2, window 1: build the FE tet mesh at t = 0 from the surface perception observes, on the FE
side (Julia), from the fused field φ_s the FE already receives (`PerceptionStep.phi` + `PerceptionGrid`,
offline from `read_perception_steps` or live from `next_step!`). Later windows do not remesh; they
move this mesh with `volume_update!`, which needs its table faces as `IEN_bottom`.

It runs once per experiment, offline: robustness matters, speed does not.

## Decisions

- **No Poisson reconstruction** (owner, 2026-10-30). The old plan (handoff 2026-10-22: Open3D Poisson
  from points + normals on the perception side) went through points only because the Python toolchain
  could not triangulate our field. In Julia, marching cubes extracts the zero set of φ_s directly.
- **The table closes the bottom.** Table height `table_z`, default 0 (rig: ChArUco origin). Keyword,
  because synthetic bags put the table at the plate top: Carreau/5 `cyl_T[3]` = 4.685 mm.
- **Unobserved patches other than the bottom are not filled.** A hole in the observed surface is an
  error. The static, unobstructed initial period is expected to show the side and the top; on
  synthetic Carreau/5 two cameras see the top at ~9° elevation (smearPerception memory 2026-10-20).
- **Input field:** reconstructed with the free-space veto on (`--free-space-veto`), so the top-rim lip
  (smearPerception memory 2026-10-30) does not become a bump on the mesh.
- **Tets only**, P1 or P2. **New mesh type** `MeshgridObserved`; plate squeeze (`simulate`) is not
  supported on it (no top plate in paper 2: a probe).
- **Gmsh through the existing template path** (`_get_mesh_data` → `_run_gmsh`), so the CLI/API choice,
  `GMSH_LOCK`, option reset and P2 node remap are reused unchanged.
- **New dependencies:** MarchingCubes.jl (surface extraction), ImageMorphology.jl
  (`label_components`). Both pure Julia; owner's guideline: off-the-shelf where one fits.

## Pipeline

`mesh_from_field(phi, grid; table_z=0.0, elem_size, basis_order=2, out_dir) -> MeshgridObserved`,
and `mesh_from_field(step::PerceptionStep, grid; kwargs...)`. Lengths in mm (the link's FE units).

1. **Seal below the table.** Pad the grid by one voxel on every side. Every voxel with
   z ≤ table_z gets φ = table_z − z (outside). This cuts the NaN region under the material off from
   the NaN region below the table.
2. **Fill NaN.** `label_components` on the NaN mask. A component touching the padded border is
   outside: φ = +Δ. Every other component is enclosed by the observed surface: inside, φ = −L with L
   the grid's diagonal (any value below every real one).
3. **Cut at the table.** φ ← max(φ, table_z − z) everywhere, so the bottom cap lies exactly at
   table_z (the inside fill never wins over the plane).
4. **Leak check** (next section).
5. **Surface.** MarchingCubes.jl on the padded field at iso 0, vertices mapped to world mm with the
   grid's `x0`, `delta`, `offset`. Written as `out_dir/surface.stl`.
6. **Volume.** Template `mesh_files/templates/observed_tet.geo`, after Gmsh tutorial t13:
   `Merge "surface.stl"` (fixed name, resolved next to the generated `.geo`), `ClassifySurfaces`
   (angle splits the flat cap from the free surface), `CreateGeometry`, `Surface Loop` of all surfaces,
   `Volume`, mesh size from `elem_size`. Physical groups: `"Bottom"` = `Surface In BoundingBox` of a
   slab of half-thickness Δ/2 around table_z; `"Free"` = all other surfaces; `"Volume"`. Parameters go
   through `_generate_mesh_geo`'s `key = value;` substitution (numbers only).
   `_get_mesh_data` generates only when the `.msh` is missing, and every observation differs, so the
   `.msh` in `out_dir` is deleted first: the mesh is always rebuilt.
7. **Mesh struct.** `_get_mesh_data(...; face_groups=["Bottom","Free"],
   node_set_groups=["Bottom","Free"], body_group="Volume")` → `MeshgridObserved`.

## `MeshgridObserved <: AbstractMeshgrid`

Shared fields as every mesh type (`NodeList`, `IEN`, `ID`, `volume_element_shape` = `:Tet`,
`surface_element_shape`, `basis_order`, `nNodes`, `ne`, `initial_state`, `effect_elem_sze`) plus
`IEN_bottom`, `IEN_free`, `bottom_nodes`, `free_nodes`. Keyword constructor with a default for every
field, like the other mesh types in `src/fem/meshes/mesh_types.jl`.

## Errors

- **Leak** (step 4): a hole in the observed surface lets the interior NaN region reach the border, so
  it is filled as outside and the mesh would be a hollow shell. Check: in every horizontal slice above
  the table, the inside set (φ < 0) must not enclose outside voxels (`label_components` on the slice's
  outside set; a component not touching the slice border is a ring's interior). Error message: the
  z-range of the ringed slices and the ring interior's centroid, i.e. near the hole. Assumes the
  specimen's horizontal cross-sections have no holes (true for a lump of material at t = 0; a ring
  shape would be refused).
- `table_z` outside the grid's z range; no zero crossing in the field: `ArgumentError`.
- Gmsh failures pass through with the path of `surface.stl` to inspect.
- `@info`: the share of voxels filled as inside and as outside.

## Tests — `test/mesh_observed_test.jl` (in `Pkg.test()`)

Analytic field: cylinder r = 30 mm, h = 40 mm on a table at table_z ≠ 0 (catches a hard-coded 0),
grid Δ = 1 mm, made TSDF-like: NaN beyond a band μ = 3 mm outside, NaN deeper than 7 mm inside, NaN
below the table.

- Mesh volume within 2 % of πr²h; free-surface nodes within one voxel of the true surface.
- Every `Bottom` node at z = table_z; `Bottom` faces cover the base (area within 2 % of πr²).
- P1 and P2 both build; `volume_update!` runs on the result.
- A hole cut through the band of the side wall: error naming a z-range that contains the hole.

## Ground-truth comparison (experiment, not in `Pkg.test()`: the data is not in the repo)

- Input: step 0 of `reconstruct_bag.py <synthetic Carreau/5 bag> out.h5 --gt <sim data dir>
  --save-field --free-space-veto --steps 0 1` (no saved reconstruction exists on disk now; run it first).
- `experiments/mesh_init/gt_comparison.jl`: `read_perception_steps` → `mesh_from_field(step,
  grid; table_z=4.685)` → writes the mesh's free surface as STL.
- `smearPerception/tools/analysis/compare_mesh_gt.py`: GT `Blender_obj/play_dough_000.obj` shifted by
  `cyl_T` (mm), scored with Open3D `RaycastingScene` as `--gt` scoring does, both ways: observed surface
  nodes → GT signed distance; GT vertices → distance to the observed surface. Reports median, p95, max,
  and the relative volume difference. GT data:
  `SMEAR-DataFiles/_keep_carreau5/ground_truth/sim_data/Carreau/5/data/`.

## Files

- New: `src/fem/meshes/mesh_observed.jl`, `mesh_files/templates/observed_tet.geo`,
  `test/mesh_observed_test.jl`, `experiments/mesh_init/gt_comparison.jl`,
  `smearPerception/tools/analysis/compare_mesh_gt.py`.
- Changed: `src/fem/meshes/mesh_types.jl` (`MeshgridObserved`), `src/smearFEM.jl` (include, export),
  `Project.toml` (deps + compat), `test/runtests.jl`.

## Not decided (defaults to pick during implementation, marked `[not decided]` in code)

- `elem_size` default; `ClassifySurfaces` angle (40° as in t13; the cap/side edge is 90°).
- Whether to smooth the marching-cubes surface before Gmsh (start without; the GT comparison decides).

## Out of scope

Filling unobserved patches; time-averaging the initial field (perception side, question 18); sending
the mesh back to perception; plate squeeze on `MeshgridObserved`.
