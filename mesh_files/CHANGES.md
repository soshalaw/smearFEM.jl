# Mesh Template Changelog

## 2026-06-29 — nz-primary parametrisation + element-count calibration

### Motivation
All 3D templates previously used `elem_size` (a characteristic length in mm) as their
primary input.  This made it hard to reason about how many elements a mesh would contain
and made the mesh resolution dependent on the geometry dimensions.  All 3D templates were
reparametrised so that `nz` (number of elements along the height / z-direction) is the
primary input, targeting a roughly cubic element distribution (`n_total ≈ nz³`).

---

### `cylinder_hex.geo`

**Before**
```
elem_size = 6.0;   // characteristic length in mm (user had to tune manually)
```

**After**
```
nz = 4;            // number of elements along height (primary input)
elem_size = radius * Sqrt(2.65 * Pi) / nz;
```

**Why the factor √(2.65π):**
Frontal-Delaunay on a circular boundary produces approximately 2.65× more quads than the
ideal area formula `n_2d = πR²/h²` predicts.  The corrected formula ensures `n_2d ≈ nz²`,
giving `n_total ≈ nz³`.  Calibrated from one empirical data point (nz=7 → 130 base quads).

`SubdivisionAlgorithm=1` was **removed**.  On a circular boundary Blossom already pairs
nearly all triangle pairs; leaving `SubdivisionAlgorithm=1` active inflated each unpaired
triangle into 3 quads, producing ~6× element count inflation.

---

### `cylinder_tet.geo`

**Before:** `elem_size_2d` (separate 2D variable) + manual `elem_size` for height.

**After**
```
nz = 4;
elem_size = height / nz;
```
All `elem_size_2d` references renamed to `elem_size`.  Explicit top/bottom circles and
lateral surfaces are defined (no extrusion), so the tet mesher fills the volume directly.

---

### `box_hex_sharp.geo` (sharp corners, Transfinite)

**Before:** `elem_size = lz / nz;`

**After**
```
nz = 4;
elem_size_xy = Sqrt(lx * ly) / nz;   // geometric mean of cross-section extents
nx = Ceil(lx / elem_size_xy);
ny = Ceil(ly / elem_size_xy);
```

**Why `√(lx·ly)` instead of `lz/nz`:**
`lz/nz` couples the 2D element size to the height.  For boxes wider than they are tall
(`lx > lz`), this over-refines the XY plane.  Using the geometric mean of the cross-section
extents makes `nx ≈ ny ≈ nz` for any aspect ratio, keeping `n_total ≈ nz³`.
Transfinite meshing gives exact `nx × ny × nz` elements.

---

### `box_hex.geo` (rounded corners, Frontal-Delaunay)

**Before:** `elem_size = lz / nz;`

**After**
```
nz = 4;
elem_size_xy = Sqrt(3.5 * lx * ly) / nz;
```

**Why the factor √3.5:**
`RecombinationAlgorithm=2` (full-quad) runs a second subdivision pass that converts every
existing quad into 4 sub-quads (in addition to converting unpaired triangles to 3 quads).
This inflates the base quad count by approximately 3.5× for a rounded-rectangle surface.
Multiplying the char length by √3.5 pre-compensates for this, targeting `n_2d ≈ nz²`.

`SubdivisionAlgorithm=1` was removed (it was redundant — algorithm 2 already performs
the equivalent subdivision internally).

Calibrated from: nz=6, 50×50mm → 48 base quads → 288 volume Q2 hexes (target nz³=216,
~33% over; acceptable for an unstructured rounded-corner mesh).

---

### `box_tet.geo` (rounded corners)

**Before:** `elem_size = lz / nz;`

**After**
```
nz = 4;
elem_size_xy = Sqrt(lx * ly) / nz;
```
No recombination correction needed (tet meshes do not recombine).

---

### `box_tet_sharp.geo` (sharp corners)

**Before:** `elem_size = lz / nz;`

**After**
```
nz = 4;
elem_size_xy = Sqrt(lx * ly) / nz;
```
Same reasoning as `box_hex_sharp.geo`.  No extrusion — prisms from the layered extrusion
are split into tets automatically by omitting `Recombine`.

---

### 2D / 1D templates (unchanged)

`disk_quad.geo`, `disk_tri.geo`, `square_quad.geo`, `square_tri.geo`, `line.geo` retain
`elem_size` as a direct characteristic length (in mm).  There is no z-direction to anchor
`nz` for these templates.

---

### Julia API (`mesh_gmsh.jl`)

`_get_gmsh_cylinder` and `_get_gmsh_box` were updated to pass `"nz" => round(Int, elem_size)`
instead of a raw characteristic length.  The public `meshgrid_cylinder` / `meshgrid_cuboid`
API is unchanged — callers still pass `elem_size=Float64(nz)`.
