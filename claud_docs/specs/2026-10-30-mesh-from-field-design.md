# Initial FE mesh from the fused field φ_s — design

Branch `features/volume_update`. Agreed with the owner 2026-10-30, approved 2026-10-31. Revised
2026-10-31 after a prototype (scratch, not committed): see "Changed by the prototype".

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
  (`label_components`). Both pure Julia; owner's guideline: off-the-shelf where one fits. Added to
  `Project.toml` by the owner 2026-10-31.

## Pipeline

`mesh_from_field(phi, grid; table_z=0.0, elem_size, basis_order=2, out_dir=mktempdir()) ->
MeshgridObserved`, and `mesh_from_field(step::PerceptionStep, grid; kwargs...)`. Lengths in mm (the
link's FE units). `elem_size` is required: no default fits every specimen.

1. **Seal below the table.** Pad the grid by one NaN voxel on every side. Every voxel with
   z ≤ table_z gets φ = table_z − z (outside). This cuts the NaN region under the material off from
   the NaN region below the table. The padding stays NaN until the fill, so the outside NaN region
   reaches the border.
2. **Fill NaN.** `label_components` on the NaN mask. A component touching the padded border is
   outside: φ = +Δ. Every other component is enclosed by the observed surface: inside, φ = −L with L
   the padded grid's diagonal (below every real value).
3. **Cut at the table.** φ ← max(φ, table_z − z) everywhere, so the bottom cap lies exactly at
   table_z (the inside fill never wins over the plane).
4. **Keep the largest piece of material.** `label_components(φ < 0)`; smaller pieces (spurious blobs)
   become outside, with a warning.
5. **Checks** (see Errors): leak, material at the grid boundary.
6. **Nudge exact zeros** to +1e-6·Δ: a voxel centre exactly on the surface makes marching cubes emit
   zero-area triangles, which open holes in the surface Gmsh receives.
7. **Surface.** MarchingCubes.jl on the padded field at iso 0, with the voxel centres as coordinates
   (world mm). The cap = triangles with all three vertices within 1e-3·Δ of table_z; the rest is the free
   surface. Written as `out_dir/surface.msh` (Gmsh MSH 2.2, shared nodes; elementary entity 1 = cap,
   2 = free).
8. **Volume.** Template `mesh_files/templates/observed_tet.geo`: `Merge "surface.msh"` (resolved next to
   the generated `.geo`, verified for the CLI and the API), `CreateTopology`, `CreateGeometry` (each of
   the two disk-shaped entities is reparametrized whole), `Surface Loop(1) = {1, 2}`, `Volume`, size from
   `elem_size`, `Mesh.SecondOrderLinear = 1`. Physical groups `"Bottom"` = surface 1, `"Free"` =
   surface 2, `"Volume"`. `elem_size` goes through `_generate_mesh_geo`'s `key = value;` substitution.
   `_get_mesh_data` generates only when the `.msh` is missing, so `out_dir/observed.msh` is deleted first.
9. **Mesh struct.** `_get_mesh_data(...; face_groups=["Bottom","Free"],
   node_set_groups=["Bottom","Free"], body_group="Volume")` → `MeshgridObserved`.

## Changed by the prototype (2026-10-31)

- **STL replaced by MSH 2.2, `ClassifySurfaces` and `Surface In BoundingBox` dropped.** Marching cubes
  rounds the 90° bottom edge over several shallow steps, so the angle classification merged the cap
  into the side (2 patches, each spanning the full height) and the Bottom group came out empty. We
  classify the cap ourselves instead. Two STL solids do not share vertices and broke the topology; MSH
  with shared nodes does not.
- **Exact zeros nudged** (step 6): the r = 30 test cylinder on a 1 mm grid has 640 zero-area triangles
  (Pythagorean voxel centres on the surface); dropping them left holes ("Wrong topology of boundary
  mesh for parametrization").
- **`Mesh.SecondOrderLinear = 1`:** P2 mid-edge nodes curved onto the reparametrized surface folded
  52 tets (negative Jacobian); straight edges give minimum det J > 0.
- Measured on the test cylinder (r 30, h 40, Δ 1, `elem_size` 6): volume −1.05 %; free-surface
  vertices within 0.34 mm of the true surface (P2 mid-edge nodes up to 1.5 mm, chords across the
  rounded rims); bottom nodes on the table to 1e-15 mm; the cap is about one voxel smaller in radius
  than the true base (area −5.8 %, the bevel of the bottom edge belongs to the free surface); P1 927
  nodes, P2 6 219, Gmsh ~0.85 s.

## `MeshgridObserved <: AbstractMeshgrid`

Shared fields as every mesh type (`NodeList`, `IEN`, `ID`, `volume_element_shape` = `:Tet`,
`surface_element_shape`, `basis_order`, `nNodes`, `ne`, `initial_state`, `effect_elem_sze`) plus
`IEN_bottom`, `IEN_free`, `bottom_nodes`, `free_nodes`. Keyword constructor with a default for every
field, like the other mesh types in `src/fem/meshes/mesh_types.jl`.

## Errors and warnings

- **Leak:** a hole in the observed surface lets the interior NaN region reach the border, so it is
  filled as outside and the mesh would be a hollow shell. Check: in every horizontal slice above the
  table, the inside set (φ < 0) must not enclose outside voxels (`label_components` on the slice's
  outside set; a component not touching the slice border is a ring's interior). Error ("not closed")
  naming the z-range of the ringed slices, and the likely causes: an unobserved patch, or `table_z`
  too low. Assumes the specimen's horizontal cross-sections have no holes (true for a lump of material
  at t = 0; a ring shape is refused). A lip ring detached from the wall is a separate piece of
  material, so step 4 drops it with a warning before this check.
- `table_z` outside the grid's z range; no material (no φ < 0) above the table; no cap (the material
  does not reach the table): `ArgumentError`.
- Material in the outermost voxel layer of the sides or top: warning, the mesh is cut flat there.
- Gmsh failures pass through; `out_dir` keeps `surface.msh` to inspect.
- `@info`: voxels filled as inside and as outside.

## Tests — `test/mesh_observed_test.jl` (in `Pkg.test()`)

Analytic field: cylinder r = 30 mm, h = 40 mm on a table at table_z = 4.685 (catches a hard-coded 0),
grid Δ = 1 mm, made TSDF-like: NaN beyond a band μ = 3 mm outside, NaN deeper than 7 mm inside, NaN
below the table.

- Closed field has the sign of the true distance wherever that is more than one voxel from 0.
- Surface: closed (every edge in two triangles), volume within 0.5 % of πr²h, cap at table_z.
- Mesh: volume within 2 %; free-surface vertices within one voxel of the true surface; bottom nodes at
  table_z; cap area between π(r − Δ)² and πr²; every tet positively oriented; P1 and P2.
- table_z on a voxel plane (the rig: grid bottom = table = 0); the Gmsh API path (no `gmsh` on PATH);
  `volume_update!` meets its constraints on the result.
- Errors and warnings above: a hole through the side wall, a detached lip ring, a spurious blob,
  material at the grid side, `table_z` outside the grid, an empty field, floating material.

## Ground-truth comparison (experiment, not in `Pkg.test()`: the data is not in the repo)

- Input: step 0 of `reconstruct_bag.py <synthetic Carreau/5 bag> out.h5 --gt <sim data dir>
  --save-field --free-space-veto --steps 0 1` (no saved reconstruction exists on disk now; run it first).
- `experiments/mesh_init/gt_comparison.jl`: `read_perception_steps` → `mesh_from_field(step,
  grid; table_z=4.685)` → writes the mesh's boundary (corner nodes of Bottom + Free faces) as OBJ.
- `smearPerception/tools/analysis/compare_mesh_gt.py`: GT `Blender_obj/play_dough_000.obj` (mm, mesh
  frame) shifted by `cyl_T` (mm), scored with `sdf.load_mesh`/`signed_distance` (Open3D
  `RaycastingScene`) as `--gt` scoring does, both ways: observed boundary vertices → GT signed distance;
  GT vertices → signed distance to the observed surface. Reports signed median, |·| p95 and max, and the
  relative volume difference. GT data: `SMEAR-DataFiles/_keep_carreau5/ground_truth/sim_data/Carreau/5/data/`.

## Files

- New: `src/fem/meshes/mesh_observed.jl`, `mesh_files/templates/observed_tet.geo`,
  `test/mesh_observed_test.jl`, `experiments/mesh_init/gt_comparison.jl`,
  `smearPerception/tools/analysis/compare_mesh_gt.py`, `smearPerception/tests/test_compare_mesh_gt.py`.
- Changed: `src/fem/meshes/mesh_types.jl` (`MeshgridObserved`), `src/smearFEM.jl` (include, export),
  `src/io/perception_link.jl` (`_grid_centres`, shared with `sample_field`), `Project.toml` (deps +
  compat), `test/runtests.jl`, `src/solver/volume_update.jl` (docstring: observed meshes).

## Not decided

- Whether to smooth the marching-cubes surface before Gmsh (start without; the GT comparison decides).
- `elem_size` for paper 2 (required argument; the experiment uses 3 mm).

## Out of scope

Filling unobserved patches; time-averaging the initial field (perception side, question 18); sending
the mesh back to perception; plate squeeze on `MeshgridObserved`.
