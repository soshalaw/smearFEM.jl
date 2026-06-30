# Mesh Templates — Usage Guide

Generated `.msh` files are cached in subdirectories (e.g. `mesh_files/cylinder/Hex/`).
Delete the cached `.msh` and `.geo` files whenever you modify a template so the next run
regenerates the mesh.

---

## 3D — Cylinder

### Hex mesh (`cylinder_hex.geo`)

```julia
mesh = meshgrid_cylinder(radius, height;
    mesh_type      = :unstructured,
    element_shape  = :Hex,
    basis_order    = 1,      # or 2 for quadratic elements
    elem_size      = 6.0,    # ← this is nz (number of elements along height)
    mesh_path      = "mesh_files",
    ndof           = 3)
```

**Element count:** `n_total ≈ nz³` (cubic element distribution).
The 2D characteristic length is `radius * √(2.65π) / nz` — calibrated so the circular
base has `n_2d ≈ nz²` quads after Blossom recombination.

**Boundary node sets:** `mesh.top_nodes`, `mesh.bottom_nodes`, `mesh.side_nodes`
**Surface IENs:** `mesh.IEN_top`, `mesh.IEN_bottom`, `mesh.IEN_lateral`

---

### Tet mesh (`cylinder_tet.geo`)

```julia
mesh = meshgrid_cylinder(radius, height;
    mesh_type      = :unstructured,
    element_shape  = :Tet,
    basis_order    = 1,
    elem_size      = 6.0,    # nz
    mesh_path      = "mesh_files",
    ndof           = 3)
```

**Element count:** approximately `6 × nz³` (tets are ~6× denser than equivalent hexes).
The characteristic length is `height / nz` (no recombination correction needed).

---

## 3D — Cuboid (box)

### Hex mesh, sharp corners (`box_hex_sharp.geo`)

```julia
mesh = meshgrid_cuboid(lx, ly, lz;
    mesh_type      = :unstructured,
    element_shape  = :Hex,
    basis_order    = 1,
    elem_size      = 6.0,    # nz
    mesh_path      = "mesh_files",
    ndof           = 3)
```

Uses Transfinite (structured) meshing — gives exact `nx × ny × nz` elements where
`nx = Ceil(lx/lz · nz)`, `ny = Ceil(ly/lz · nz)`.
For a cubic box (`lx = ly = lz`): exactly `nz³` elements.

---

### Hex mesh, rounded corners (`box_hex.geo`)

```julia
mesh = meshgrid_cuboid(lx, ly, lz;
    mesh_type      = :unstructured,
    element_shape  = :Hex,
    basis_order    = 1,
    elem_size      = 6.0,    # nz
    mesh_path      = "mesh_files",
    ndof           = 3,
    edge_radius    = 3.0)    # fillet radius; must be < min(lx, ly) / 2
```

Uses Frontal-Delaunay + full-quad recombination.
The 2D characteristic length is `√(3.5 · lx · ly) / nz` — the √3.5 factor compensates
for the ~3.5× inflation caused by `RecombinationAlgorithm=2` subdividing every quad into
4 sub-quads.  Typical element count: `~1.3 × nz³`.

**Constraint:** `0 < edge_radius < min(lx, ly) / 2`.

---

### Tet mesh, sharp corners (`box_tet_sharp.geo`)

```julia
mesh = meshgrid_cuboid(lx, ly, lz;
    mesh_type      = :unstructured,
    element_shape  = :Tet,
    basis_order    = 1,
    elem_size      = 6.0,    # nz
    mesh_path      = "mesh_files",
    ndof           = 3)
```

Extruded prisms are split into tets.  Element count: `~6 × nz³`.

---

### Tet mesh, rounded corners (`box_tet.geo`)

```julia
mesh = meshgrid_cuboid(lx, ly, lz;
    mesh_type      = :unstructured,
    element_shape  = :Tet,
    basis_order    = 1,
    elem_size      = 6.0,    # nz
    mesh_path      = "mesh_files",
    ndof           = 3,
    edge_radius    = 3.0)
```

---

## 2D — Disk

```julia
mesh = meshgrid_disk(radius;
    element_shape  = :Quad,   # or :Tri
    basis_order    = 1,
    elem_size      = 8.0,     # characteristic length in mm (not nz)
    mesh_path      = "mesh_files",
    ndof           = 2)
```

**Boundary nodes:** `mesh.boundary_nodes`

---

## 2D — Square

```julia
mesh = meshgrid_square(lx, ly;
    mesh_type      = :unstructured,
    element_shape  = :Quad,   # or :Tri
    basis_order    = 1,
    elem_size      = 10.0,    # characteristic length in mm (not nz)
    mesh_path      = "mesh_files",
    ndof           = 2)
```

The square is **centered at the origin**: node x-coordinates span `[-lx/2, lx/2]`.
**Boundary node sets:** `mesh.top_nodes`, `mesh.bottom_nodes`, `mesh.side_nodes`

---

## 1D — Line

```julia
mesh = meshgrid_line(length;
    mesh_type      = :unstructured,
    element_shape  = :Line,
    basis_order    = 1,
    elem_size      = 8.0,     # characteristic length in mm (not nz)
    mesh_path      = "mesh_files",
    ndof           = 1)
```

The line is **centered at the origin**: node x-coordinates span `[-length/2, length/2]`.
**Boundary nodes:** `mesh.boundary_nodes`

---

## Choosing `nz`

| nz | Approx. cylinder hex elements | Approx. sharp-box hex elements |
|----|-------------------------------|-------------------------------|
| 3  | ~27                           | 27 (exact)                    |
| 4  | ~64                           | 64 (exact)                    |
| 6  | ~192                          | 216 (exact)                   |
| 8  | ~512                          | 512 (exact)                   |
| 11 | ~1331                         | 1331 (exact)                  |

Rounded-corner hex meshes (`box_hex.geo`) produce ~1.3× the sharp-corner count.
Tet meshes produce ~6× the equivalent hex count.

---

## Mesh caching

Meshes are cached as `<mesh_path>/<geometry>/<ElementShape>/<params>.msh`.
Examples:
```
mesh_files/cylinder/Hex/6.0_1.msh
mesh_files/box/Hex/50.0x50.0x40.0_6.0_1_r3.0.msh
mesh_files/disk/Quad/8.0_1.msh
```

**When to delete the cache:** any time you modify a `.geo` template or change geometry
parameters that were not already part of the filename (e.g. `radius`, `height`, `r`).
The Julia code regenerates the `.msh` automatically if the file is absent.
