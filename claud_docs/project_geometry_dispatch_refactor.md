---
name: geometry-dispatch-refactor
description: Planned refactor — geometry wrapper types + squeeze_stokes.jl split into stokes_solver.jl + stokes_setup.jl
metadata:
  type: project
---

Refactor of `squeeze_stokes.jl` using geometry wrapper types for `def_problem`/`set_model` dispatch. **COMPLETE as of 2026-06-20.**

**Why:** Multiple dispatch on raw number args is ambiguous for 2D/1D geometries (Square vs Cylinder both have 2 dim args; Disk vs Segment both have 1). Geometry types solve this cleanly.

**How to apply:** Implement in order below; update call sites last.

---

## Step 1 — New file: `src/fem/geometries.jl`

```julia
abstract type AbstractGeometry end
struct Cylinder{R<:Number,H<:Number} <: AbstractGeometry; r::R; h::H; end
struct Cuboid{X<:Number,Y<:Number,Z<:Number} <: AbstractGeometry; lx::X; ly::Y; lz::Z; end
struct Disk{R<:Number} <: AbstractGeometry; r::R; end
struct Square{X<:Number,Y<:Number} <: AbstractGeometry; lx::X; ly::Y; end
struct Segment{L<:Number} <: AbstractGeometry; l::L; end

ndim(::Union{Cylinder,Cuboid}) = 3
ndim(::Union{Disk,Square})     = 2
ndim(::Segment)                = 1
```

Export: `Cylinder, Cuboid, Disk, Square, Segment, ndim`
Include before `squeeze_stokes.jl` in `smearFEM.jl`.

---

## Step 2 — Split `squeeze_stokes.jl` (1552 lines) into two files

**`stokes_solver.jl`** (lines 1–958):
- `assemble_system_A` / `_dense`
- `assemble_system_B` / `_dense`
- `apply_boundary_conditions` / `_dense`
- `set_boundary_cond` ← needs fix (see Step 5)
- `simulate`

**`stokes_setup.jl`** (lines 959–1552):
- `get_η_power_law` (cylinder-only)
- `set_model` (5 dispatch methods)
- `def_problem` (5 dispatch methods)

Delete `squeeze_stokes.jl`. Update `smearFEM.jl` includes.

---

## Step 3 — `set_model` dispatch (in `stokes_setup.jl`)

Remove explicit `ndim` arg (infer via `ndim(geom)`). Keep `nDof_u`, `nDof_p`, element shapes.

```julia
function set_model(geom::Cylinder, ne::Float64, η, es_u, bo_u, nDof_u, es_p, bo_p, nDof_p, es_x, bo_x; GMESH_MESH=true, filepath_mesh="")
    _dim = ndim(geom)
    _mesh(es, bo, ndof) = meshgrid_cylinder(geom.r, geom.h; mesh_type=GMESH_MESH ? :unstructured : :structured, element_shape=es, basis_order=bo, ne=round(Int,ne), ndof=ndof, elem_size=ne, mesh_path=filepath_mesh)
    ...
    return Stokes(ndim=_dim, ...)
end

# Cuboid: extra edge_radius kwarg
# Disk: meshgrid_disk (unstructured only, no :structured path)
# Square: meshgrid_square
# Segment: meshgrid_line
```

---

## Step 4 — `def_problem` dispatch (in `stokes_setup.jl`)

Power law (`get_η_power_law`) only for `Cylinder`. For other geometries, throw if `viscosity_type == "bulk_viscosity" && viscosity_model == "power_law"`.

```julia
def_problem(geom::Cylinder, ne, η_0, es_u, bo_u, nDof_u, es_p, bo_p, nDof_p, es_x, bo_x, β, cParam, control, viscosity_type, sim_time, t_steps; viscosity_model="power_law", GMESH_MESH=true, mesh_path=...)
def_problem(geom::Cuboid,   ...; edge_radius=nothing, ...)
def_problem(geom::Disk,     ...; ...)
def_problem(geom::Square,   ...; ...)
def_problem(geom::Segment,  ...; ...)
```

---

## Step 5 — Fix `set_boundary_cond` for 2D/1D

In `stokes_solver.jl`, change hardcoded `coord[3]` / `NodeList[3,:]` to use `ndim_cached`.

For `nDof_u > 1` block (currently 3D only):
```julia
# OLD:
z1Bound = maximum(mdl.mesh_u.NodeList[3, :])
coord[3] == z1Bound  →  q_upper[ID_cached[3,nNode]]
# NEW:
z1Bound = maximum(mdl.mesh_u.NodeList[ndim_cached, :])
coord[ndim_cached] == z1Bound  →  q_upper[ID_cached[ndim_cached,nNode]]
```

For `nDof_u == 1` block, add `ndim_cached == 1` branch (uses `coord[1]`).

---

## Step 6 — Update call sites

**`scripts/stokes/single_simulation.jl`:**
- `const_vel(lx, ly, lz, ...)` calls `def_problem(Cuboid(lx, ly, lz), ne, ...)`
- Cylinder version: `def_problem(Cylinder(r, h), ne, ...)`

**`src/examples/run_example.jl`:**
- Cylinder: `def_problem(Cylinder(r, h), ne_gt, ...)`
- Cuboid: `def_problem(Cuboid(lx, ly, lz), ne_gt, ...; edge_radius=edge_radius)`

**`smearOptimize.jl`:** `meshgrid_cuboid` call is direct (no `def_problem`), unchanged.

---

## Step 7 — `smearFEM.jl` updates

```julia
# Exports to add:
export Cylinder, Cuboid, Disk, Square, Segment, ndim

# Includes: replace
include("examples/squeeze_stokes.jl")
# with:
include("fem/geometries.jl")
include("examples/stokes_solver.jl")
include("examples/stokes_setup.jl")
```
