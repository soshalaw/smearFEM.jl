# Meshes.jl Refactor Plan

## Goal

Clean up `src/fem/Meshes.jl` (1521 lines) to remove dead code, make struct types consistent across structured and unstructured mesh paths, and split the file into logical modules.

---

## Status

| Step | Description | Status |
|------|-------------|--------|
| 1 | Delete dead structs | ✅ Done |
| 2 | Add `initial_state` to all structs | ✅ Done |
| 3 | Unify `MeshgridLine` as the return type for line meshes | ✅ Done |
| 4 | Add `MeshgridDisk` to replace `Meshgrid2D` for disk meshes | ✅ Done |
| 5 | Standardise boundary storage field names | ✅ Done |
| 6 | Split into separate files | ✅ Done |

---

## Step 1 — Delete dead structs ✅

**Removed:**
- `Meshgrid{T}` — unused wrapper shell
- `Meshgrid1D` — replaced by `MeshgridLine`
- `Meshgrid2D` — replaced by `MeshgridDisk`
- `Meshgrid3D` — never used as a return type anywhere

**Updated downstream:**
- `models.jl`: default sentinel `Meshgrid1D()` → `MeshgridLine()`
- `smearFEM.jl`: exports updated

---

## Step 2 — Add `initial_state` to all structs ⬜

`reset_mesh!` and `update_initial_state!` call `mesh.initial_state`, but only `MeshgridCube` and `MeshgridCylinder` have that field. Calling either function on a `MeshgridLine`, `MeshgridDisk`, or `MeshgridSquare` crashes at runtime.

**Changes needed in `Meshes.jl`:**

Add to `MeshgridLine`:
```julia
initial_state::Matrix{Float64}
```
Constructor: `initial_state=copy(NodeList)`

Add to `MeshgridDisk`:
```julia
initial_state::Matrix{Float64}
```
Constructor: `initial_state=copy(NodeList)`

Add to `MeshgridSquare`:
```julia
initial_state::Matrix{Float64}
```
Constructor: `initial_state=copy(NodeList)`

`MeshgridCube` and `MeshgridCylinder` already have it — no change needed.

**Affected callers:** None — `reset_mesh!` and `update_initial_state!` dispatch on `AbstractMeshgrid` and will work on all types once the field is present.

---

## Step 3 — Unify `MeshgridLine` as the return type ✅

Both `meshgrid_line(:structured)` and `_get_gmsh_line` now return `MeshgridLine`.
`MeshgridLine` gained a `boundary_nodes::Vector{Int}` field.

---

## Step 4 — Add `MeshgridDisk` ✅

New struct with fields: `r`, `NodeList`, `IEN`, `IEN_boundary`, `ID`,
`volume_element_shape`, `surface_element_shape`, `basis_order`, `nNodes`, `ne`, `boundary_nodes`.

`_get_gmsh_disk` and `meshgrid_disk` now return `MeshgridDisk`.

---

## Step 5 — Standardise boundary storage field names ⬜

Boundary storage is inconsistent across structs:

| Struct | Boundary IEN storage | Boundary node storage |
|--------|---------------------|-----------------------|
| `MeshgridLine` | — | `boundary_nodes::Vector{Int}` |
| `MeshgridDisk` | `IEN_boundary::Matrix{Int}` | `boundary_nodes::Vector{Int}` |
| `MeshgridSquare` | `IEN_top`, `IEN_bottom`, `IEN_side` | `top_nodes`, `bottom_nodes`, `side_nodes` |
| `MeshgridCube` | `IEN_top`, `IEN_bottom`, `IEN_sides::Dict{Symbol,Matrix{Int}}` | `top_nodes`, `bottom_nodes`, `side_nodes` |
| `MeshgridCylinder` | `IEN_top`, `IEN_bottom`, `IEN_sides`, `IEN_cp`, `IEN_vis` | `top_nodes`, `bottom_nodes`, `side_nodes` |

**Target convention:**
- Named face IEN: `IEN_top`, `IEN_bottom`, `IEN_side(s)` (Matrix{Int})
- Named node sets: `top_nodes`, `bottom_nodes`, `side_nodes` (Vector{Int})
- For `MeshgridCube`: replace `IEN_sides::Dict` with four named fields (`IEN_front`, `IEN_back`, `IEN_left`, `IEN_right`) or keep the Dict — decide before implementing

**Risk:** `IEN_sides` (Dict) and `IEN_cp`/`IEN_vis` on `MeshgridCylinder` are used in downstream FEM assembly. Audit all call sites before renaming.

---

## Step 6 — Split into separate files ⬜

At ~1420 lines after cleanup, the file is still large. Proposed split:

| File | Contents | Est. lines |
|------|----------|-----------|
| `mesh_types.jl` | All struct definitions | ~220 |
| `mesh_structured.jl` | `_meshgrid_line/square/cube`, `_meshgrid_ring`, `_meshgrid_cylinder`, `_inflate_cylinder` | ~700 |
| `mesh_gmsh.jl` | `_get_gmsh_*` functions | ~230 |
| `mesh_api.jl` | Public `meshgrid_*`, `reset_mesh!`, `update_initial_state!` | ~270 |

**Changes to `smearFEM.jl`:** replace `include("fem/Meshes.jl")` with the four new includes (order matters: types first).

---

## Notes

- Steps 2 and 5 are independent and can be done in either order.
- Step 6 should be done last — it is purely mechanical and does not change behaviour.
- Step 5 carries the most risk due to downstream usage of `IEN_sides` (Dict) in `PostProcess.jl` and `smearOptimize.jl`. Audit those files before starting.
