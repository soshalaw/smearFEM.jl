# Refactor: Split Element Shape from Basis Order

## Motivation

`FunctionClass::String` currently conflates two independent concepts into one token:

- **Element topology** (the prefix `Q`/`T`): determines reference element shape, number of
  nodes per element, quadrature rule, and which `.geo` template gmsh uses.
- **Polynomial order** (the trailing digit `1`/`2`/`3`): determines DOF count per element,
  derivative expressions, and the `-order` flag passed to gmsh.

Mixing them into one string makes it impossible to express combinations independently — e.g.
a quadratic tet geometry mesh (`mesh_x`) paired with a linear tet velocity mesh (`mesh_u`)
would need a new string token, when in reality it is just two orthogonal choices.

---

## Proposed Data Model

Replace `FunctionClass::String` with two fields everywhere it appears:

| Field | Type | Values | What it replaces |
|-------|------|--------|-----------------|
| `element_shape::Symbol` | `Symbol` | `:Hex`, `:Tet`, `:Quad`, `:Tri`, `:Line` | prefix of `FunctionClass` |
| `basis_order::Int` | `Int` | `1`, `2`, `3` | trailing digit of `FunctionClass` |

### Mapping from current tokens

| `FunctionClass` | `element_shape` | `basis_order` |
|----------------|-----------------|---------------|
| `"Q1"` | `:Hex` (3D) / `:Quad` (2D) / `:Line` (1D) | `1` |
| `"Q2"` | `:Hex` / `:Quad` / `:Line` | `2` |
| `"Q3"` | `:Hex` | `3` |
| `"T1"` | `:Tet` (3D) / `:Tri` (2D) | `1` |
| `"T2"` | `:Tet` / `:Tri` | `2` |
| `"S1"` | (removed — NURBS support dropped) | — |

### Why `Symbol` for shape, `Int` for order

`element_shape` is a pure categorical label — it is only ever compared, never manipulated
as a string and never passed directly to external APIs. `Symbol` comparison uses `===`
(identity), is O(1), and visually signals "this is an enum-like value".

`basis_order` is a numeric quantity used in arithmetic and passed directly to the gmsh
`-order` CLI flag. Using `Int` means no conversion is needed at the gmsh boundary.

---

## File-by-File Changes

### 1. `src/fem/Meshes.jl`

**All 7 mesh structs** (`Meshgrid1D`, `Meshgrid2D`, `Meshgrid3D`, `MeshgridLine`,
`MeshgridSquare`, `MeshgridCube`, `MeshgridCylinder`):

```julia
# Before
FunctionClass::String

# After
element_shape::Symbol
basis_order::Int
```

**All `meshgrid_*` constructor keyword arguments**:

```julia
# Before
FunctionClass::String = "Q1"

# After
element_shape::Symbol = :Hex,  basis_order::Int = 1
```

**Internal branch points** — topology check becomes shape comparison, order check becomes
integer comparison:

```julia
# Before
if FunctionClass == "Q1" ... elseif FunctionClass == "Q2" ... elseif FunctionClass == "Q3"

# After
if basis_order == 1 ... elseif basis_order == 2 ... elseif basis_order == 3
```

Tet-specific branches (e.g. in `meshgrid_cylinder`):

```julia
# Before
if string(FunctionClass[1]) == "S" ... elseif string(FunctionClass[1]) == "Q"

# After
if element_shape === :Tet ... elseif element_shape === :Hex
```

**`get_gmsh_cylinder`** — the two derived values now come directly from the two fields:

```julia
# Before
template_path = startswith(FunctionClass, "T") ?
    joinpath(mesh_dir, "mesh_tet.geo") : joinpath(mesh_dir, "mesh.geo")
mesh_order = (FunctionClass == "Q2" || FunctionClass == "T2") ? 2 : 1

# After — no derivation needed
template_path = element_shape === :Tet ?
    joinpath(mesh_dir, "mesh_tet.geo") : joinpath(mesh_dir, "mesh.geo")
mesh_order = basis_order   # Int pass-through, same value gmsh expects
```

---

### 2. `src/fem/fem.jl`

**`_basis_1d_components`** — takes `basis_order::Int` instead of `FunctionClass`:

```julia
# Before
function _basis_1d_components(ξ::Float64, FunctionClass::String)
    if FunctionClass == "Q1" || FunctionClass == "S1"
    elseif FunctionClass == "Q2"

# After
function _basis_1d_components(ξ::Float64, basis_order::Int)
    if basis_order == 1
    elseif basis_order == 2
```

**`_format_1d_derivatives`** — same change, `basis_order::Int`:

```julia
# Before
function _format_1d_derivatives(FunctionClass::String, ΔN, Δ2N)
    if FunctionClass == "Q1" || FunctionClass == "S1"

# After
function _format_1d_derivatives(basis_order::Int, ΔN, Δ2N)
    if basis_order == 1
```

**`__basis_tet_3d` and `_basis_tet_2d`** — take `basis_order::Int`:

```julia
# Before
function __basis_tet_3d(ξ, η, ζ, FunctionClass::String)
    if FunctionClass == "T1" ... elseif FunctionClass == "T2"

# After
function __basis_tet_3d(ξ, η, ζ, basis_order::Int)
    if basis_order == 1 ... elseif basis_order == 2
```

**`basis_function` dispatch** — all overloads replace the combined string with the two fields:

```julia
# Before
function basis_function(ξ, η, ζ, FunctionClass::String = "Q2"; second_derivatives=false)
    if FunctionClass == "Q1" || FunctionClass == "S1"
        N, ΔN, Δ2N = _tensor_basis_3d(..., _Q1_3D_NODE_TRIPLES)
    elseif FunctionClass == "Q2"
        N, ΔN, Δ2N = _tensor_basis_3d(..., _Q2_3D_NODE_TRIPLES)

# After
function basis_function(ξ, η, ζ, element_shape::Symbol = :Hex, basis_order::Int = 2; second_derivatives=false)
    if element_shape === :Hex || element_shape === :Quad || element_shape === :Line
        node_triples = basis_order == 1 ? _Q1_3D_NODE_TRIPLES : _Q2_3D_NODE_TRIPLES
        N, ΔN, Δ2N = _tensor_basis_3d(..., node_triples)
    elseif element_shape === :Tet || element_shape === :Tri
        N, ΔN, Δ2N = __basis_tet_3d(ξ, η, ζ, basis_order)
```

**`get_basis_volume_functions`** — reads the two fields from mesh structs and branches on
shape for quadrature selection (simplex quadrature for tet/tri, tensor-product for hex/quad):

```julia
element_shape_u = mdl.mesh_u.element_shape
basis_order_u   = mdl.mesh_u.basis_order

is_tet = element_shape_u === :Tet || element_shape_u === :Tri
if is_tet
    pts, wpoints = _tet_quadrature(basis_order_u)   # to be added
else
    # existing tensor-product Gauss setup
end
```

---

### 3. `src/fem/models.jl`

**`LinearElasticity`** struct replaces `FunctionClass::String` with:

```julia
element_shape::Symbol
basis_order::Int
```

Constructor keyword arguments change correspondingly. `Stokes` needs no struct change —
it holds mesh references and reads shape/order from `mesh_u.element_shape` etc.

---

### 4. `src/io/io.jl`

**`write_vtk`** and **`write_scene`** replace their `FunctionClass::String` parameters with
`element_shape::Symbol` and `basis_order::Int`. VTK cell type selection branches on both:

```julia
# Before
if FunctionClass == "Q1"
    cellType = ndim == 3 ? VTK_HEXAHEDRON : ...
elseif FunctionClass == "Q2"
    cellType = ndim == 3 ? VTK_LAGRANGE_HEXAHEDRON : ...
elseif FunctionClass == "T1"
    cellType = VTK_TETRAHEDRON
elseif FunctionClass == "T2"
    cellType = VTK_LAGRANGE_TETRAHEDRON

# After
if element_shape === :Hex && basis_order == 1
    cellType = ndim == 3 ? VTK_HEXAHEDRON : ...
elseif element_shape === :Hex && basis_order == 2
    cellType = ndim == 3 ? VTK_LAGRANGE_HEXAHEDRON : ...
elseif element_shape === :Tet && basis_order == 1
    cellType = VTK_TETRAHEDRON
elseif element_shape === :Tet && basis_order == 2
    cellType = VTK_LAGRANGE_TETRAHEDRON
```

---

### 5. `src/examples/run_example.jl` and `squeeze_stokes.jl`

**Function signatures** replace `FunctionClass_u::String, FunctionClass_p::String`:

```julia
# Before
function simulate_single_tstep_stokes(..., FunctionClass_u::String, FunctionClass_p::String, ...; FunctionClass_x::String=FunctionClass_u, ...)

# After
function simulate_single_tstep_stokes(..., element_shape::Symbol, basis_order_u::Int, basis_order_p::Int, ...; ...)
```

**`@unpack` calls** from mesh structs:

```julia
# Before
@unpack FunctionClass, IEN, ID, NodeList = mdl.mesh_u
FunctionClass_u_cached::String = FunctionClass

# After
@unpack element_shape, basis_order, IEN, ID, NodeList = mdl.mesh_u
```

**Config parsing boundary** — `run_example.jl` reads from TOML config:

```julia
# Before (reads combined string from config)
FunctionClass_u::String = exp_params["FunctionClass_u"]

# After (reads two keys from config)
element_shape   = Symbol(exp_params["element_shape"])
basis_order_u   = parse(Int, exp_params["basis_order_u"])
basis_order_p   = parse(Int, exp_params["basis_order_p"])
```

This requires updating `config.toml` and `config.toml.example` to use the two separate keys.

---

## What Does NOT Change

- Mutable struct layout and keyword constructor pattern
- `AbstractMeshgrid` / `AbstractModel` type hierarchy
- Caching pattern in `BasisFunctionCache` and `get_basis_volume_functions`
- Assembly loop structure in example files (reads from cache, not directly from shape/order)
- `ΔN`/`Δ2N` array shapes and column ordering conventions
- `IEN`, `NodeList`, `ID` field names and semantics

---

## Implementation Order

Steps should be done in this sequence since each layer depends on the one below it:

1. **`Meshes.jl` structs + constructors** — all other files read mesh fields, so this is the
   source of truth. Update all 7 structs and their `meshgrid_*` constructors including
   `get_gmsh_cylinder`.

2. **`fem.jl` helper functions** — update `_basis_1d_components`, `_format_1d_derivatives`,
   `__basis_tet_3d`, `_basis_tet_2d` to take `basis_order::Int`.

3. **`fem.jl` `basis_function` dispatch** — update all overloads to take
   `element_shape::Symbol, basis_order::Int`.

4. **`fem.jl` `get_basis_volume_functions` + `get_surface_basis_functions`** — update field
   reads and add simplex quadrature branch.

5. **`models.jl`** — update `LinearElasticity` struct and constructor.

6. **`io/io.jl`** — update `write_vtk` and `write_scene` VTK cell type mapping.

7. **Example files** — update `run_example.jl` and `squeeze_stokes.jl` function signatures,
   `@unpack` calls, and config parsing.

8. **`config.toml` / `config.toml.example`** — replace `FunctionClass_*` keys with
   `element_shape`, `basis_order_u`, `basis_order_p`, `basis_order_x`.

---

## Extending gmsh Integration

### Current state and coupling points

The gmsh layer currently has three hard coupling points that prevent it from supporting
anything other than a 3D cylinder:

1. **`run_gmsh` always passes `-3`** — hardcoded to 3D meshing.
2. **`get_mesh_data` hardcodes physical group names** — `"Cylinder"`, `"Top"`, `"Bottom"`,
   `"Lateral"` are cylinder-specific strings embedded directly in the function body.
3. **`generate_mesh_geo` hardcodes substitution placeholders** — only knows about `radius`,
   `height`, and `elem_size_2d`; cannot parameterise a box or 2D geometry.

All three need to be generalised before new geometry/topology combinations can be added.

---

### Layer 1 — Generalise `run_gmsh` (dimension flag)

Add a `dim::Int` keyword argument so the CLI flag matches the mesh dimension:

```julia
function run_gmsh(geo_path, msh_path, mesh_order::Int; dim::Int=3)
    cmd = `$gmsh_path $geo_path -$dim -format msh -order $mesh_order -o $msh_path`
```

Callers pass `dim=3` (volume), `dim=2` (surface/2-D domain), or `dim=1` (curve).

---

### Layer 2 — Generalise `generate_mesh_geo` (parameter dict)

Replace the three hardcoded substitutions with a `Dict{String,String}` of placeholder →
value pairs so each template can declare its own parameters:

```julia
function generate_mesh_geo(params::Dict{String,String}, output_path::String,
                            template_path::String)
    content = read(template_path, String)
    for (placeholder, value) in params
        content = replace(content, placeholder => value)
    end
    write(output_path, content)
end
```

Each `.geo` template keeps the same `variable = default_value;` comment convention so the
function can also verify staleness by checking the generated content.

---

### Layer 3 — Generalise `get_mesh_data` (physical group names)

Replace hardcoded group names and the fixed return structure with keyword arguments.
The function returns a `NamedTuple` keyed by the caller-supplied group names so downstream
code never needs to know the internal gmsh group labels:

```julia
function get_mesh_data(filePath::String;
                       volume_group::String,          # e.g. "Cylinder" or "Volume"
                       boundary_groups::Vector{String}, # e.g. ["Top","Bottom","Lateral"]
                       dim::Int=3,
                       mesh_order::Int=1,
                       ...)
    # extract IEN for volume_group at dimension dim
    # extract IEN + node sets for each name in boundary_groups at dimension dim-1
    # return NamedTuple with keys :NodeList, :IEN, :ne, :nNodes,
    #   :boundary_IEN  => Dict(name => matrix),
    #   :boundary_nodes => Dict(name => vector)
end
```

The IEN remapping block already handles the four node-count cases (4, 8, 10, 27); it needs
two more for 2-D elements:

| `npe` | Element | Permutation needed? |
|-------|---------|---------------------|
| 2 | Line Q1 | none |
| 3 | Tri T1 or Line Q2 | none |
| 4 | Quad Q1 | none |
| 6 | Tri T2 | none |
| 8 | Hex Q1 | none |
| 9 | Quad Q2 | yes — gmsh serendipity → smearFEM ordering |
| 10 | Tet T2 | yes — existing map `[1,2,3,4,5,6,7,8,10,9]` |
| 27 | Hex Q2 | yes — existing 27-entry map |

---

### Layer 4 — Standardise physical group naming across `.geo` templates

Consistent names across all templates let `get_mesh_data` work without callers knowing
geometry internals. Proposed convention:

| Dimension | Group name | Meaning |
|-----------|-----------|---------|
| 3D | `"Volume"` | whole domain |
| 3D cylinder | `"Top"`, `"Bottom"`, `"Lateral"` | face sets |
| 3D box | `"Top"`, `"Bottom"`, `"Left"`, `"Right"`, `"Front"`, `"Back"` | face sets |
| 2D | `"Surface"` | whole domain |
| 2D rectangle | `"Top"`, `"Bottom"`, `"Left"`, `"Right"` | edge sets |
| 2D disk | `"Boundary"` | outer edge |
| 1D | `"Curve"` | whole domain |
| 1D | `"Left"`, `"Right"` | endpoints |

The existing `mesh.geo` and `mesh_tet.geo` should be updated to rename `"Cylinder"` →
`"Volume"` so they are consistent with the new convention.

---

### Layer 5 — New `.geo` templates

One template per geometry × topology family (topology is controlled by `Recombine`):

| File | Geometry | Shape | Notes |
|------|----------|-------|-------|
| `mesh.geo` *(rename to `cylinder_hex.geo`)* | 3D cylinder | Hex | existing, rename "Cylinder"→"Volume" |
| `mesh_tet.geo` *(rename to `cylinder_tet.geo`)* | 3D cylinder | Tet | existing, rename "Cylinder"→"Volume" |
| `box_hex.geo` | 3D box | Hex | structured extrusion + `Recombine` |
| `box_tet.geo` | 3D box | Tet | same geometry, no `Recombine` |
| `square_quad.geo` | 2D rectangle | Quad | `Recombine Surface` |
| `square_tri.geo` | 2D rectangle | Tri | no recombine |
| `disk_quad.geo` | 2D circle | Quad | `Recombine Surface` |
| `disk_tri.geo` | 2D circle | Tri | no recombine |
| `line.geo` | 1D line | Line | single `Curve`, no surface |

Within each template, the parameter placeholders follow the convention
`variable = default_value;` so `generate_mesh_geo` can substitute them.

Common parameters across templates:

| Placeholder in `.geo` | Meaning |
|-----------------------|---------|
| `lx = 1.0;` | length / radius in x |
| `ly = 1.0;` | length in y (or height for cylinder) |
| `lz = 1.0;` | length in z (3D only) |
| `elem_size = 0.1;` | characteristic element size |

The cylinder templates keep `radius` and `height` for backward compatibility with existing
generated `.msh` files.

---

### Layer 5b — `mesh_files` Directory Structure

The current `mesh_files/` root mixes three distinct kinds of files:

- **Template `.geo` files** — reusable geometry definitions parameterised by size
- **Generated instance `.geo` files** — concrete geometry files produced by `generate_mesh_geo`
- **Generated `.msh` files** — meshed output consumed by `get_mesh_data`

Keeping them flat makes it ambiguous which `.geo` files are canonical sources and which are
ephemeral outputs. The reorganisation separates them by role and by shape:

```
mesh_files/
  templates/                    ← canonical source templates (version-controlled)
    cylinder_hex.geo
    cylinder_tet.geo
    box_hex.geo
    box_tet.geo
    square_quad.geo
    disk_tri.geo
    line.geo
  cylinder/                     ← generated cylinder instances
    cylinder_x_9.geo
    cylinder_x_9.msh
    cylinder_x_11.0.geo
    cylinder_p_11.0.geo
    cylinder_p_6.0.msh
    cylinder_p_6.0_pys.msh
  box/                          ← generated box instances (empty until first use)
  square/
  disk/
  line/
```

**Naming convention for templates:** `{shape}_{topology}.geo`

| File | Shape | Topology |
|------|-------|----------|
| `cylinder_hex.geo` | cylinder | hex |
| `cylinder_tet.geo` | cylinder | tet |
| `box_hex.geo` | box | hex |
| `box_tet.geo` | box | tet |
| `square_quad.geo` | rectangle | quad |
| `disk_tri.geo` | disk | tri |
| `line.geo` | line | line (no topology ambiguity) |

**Code path changes required:**

All five `get_gmsh_*` constructors currently resolve templates as
`joinpath(mesh_dir, "mesh_*.geo")`. After the reorganisation they resolve to
`joinpath(mesh_dir, "templates", "{shape}_{topology}.geo")`.

Callers that pass explicit `file_path` arguments (e.g. `set_model` in `squeeze_stokes.jl`
and `initialize_mesh` in `run_example.jl`) must change from
`joinpath(filepath_mesh, "cylinder_x_$(ne).msh")` to
`joinpath(filepath_mesh, "cylinder", "cylinder_x_$(ne).msh")`.

The auto-generation path in `get_mesh_data` places the intermediate `.geo` at
`splitext(filePath)[1] * ".geo"` — this naturally lands in the correct shape subdirectory
as long as the caller passes a `file_path` under `mesh_files/{shape}/`.

**`.gitignore` implication:** the `mesh_files/{shape}/` subdirectories contain generated
artefacts and can be git-ignored; only `mesh_files/templates/` should be tracked.

---

### Layer 6 — High-level mesh constructor functions in `Meshes.jl`

Replace the single `get_gmsh_cylinder` with a family of constructors, one per geometry.
Each wraps `get_mesh_data` with the correct physical group names, template selection, and
`MeshgridCylinder` / `MeshgridCube` / `MeshgridSquare` construction:

```julia
get_gmsh_cylinder(file_path, ndof, r, h, element_shape, basis_order; ...)
get_gmsh_box(file_path, ndof, lx, ly, lz, element_shape, basis_order; ...)
get_gmsh_square(file_path, ndof, lx, ly, element_shape, basis_order; ...)
get_gmsh_disk(file_path, ndof, r, element_shape, basis_order; ...)
get_gmsh_line(file_path, ndof, l, basis_order; ...)
```

Each follows the same internal pattern:

```julia
function get_gmsh_box(file_path, ndof, lx, ly, lz,
                      element_shape::Symbol, basis_order::Int; ...)
    template_path = element_shape === :Tet ?
        joinpath(mesh_dir, "box_tet.geo") : joinpath(mesh_dir, "box_hex.geo")

    params = Dict("lx = 1.0;" => "lx = $lx;",
                  "ly = 1.0;" => "ly = $ly;",
                  "lz = 1.0;" => "lz = $lz;",
                  "elem_size = 0.1;" => "elem_size = $elem_size;")

    data = get_mesh_data(file_path;
        volume_group    = "Volume",
        boundary_groups = ["Top","Bottom","Left","Right","Front","Back"],
        dim             = 3,
        mesh_order      = basis_order,
        params          = params,
        template_path   = template_path)

    # build ID, return MeshgridCube(...)
end
```

---

### Updated implementation order for gmsh work

These steps slot into the main implementation order after step 1 (mesh structs):

1. Rename/update existing `.geo` files (`"Cylinder"` → `"Volume"`).
2. Generalise `run_gmsh` (add `dim` parameter).
3. Generalise `generate_mesh_geo` (parameter dict).
4. Generalise `get_mesh_data` (group names + 2-D IEN remapping).
5. Write new `.geo` templates (box, square, disk, line — both topologies where applicable).
6. Add `get_gmsh_box`, `get_gmsh_square`, `get_gmsh_disk`, `get_gmsh_line` in `Meshes.jl`.
7. Rename `get_gmsh_cylinder` parameters from `FunctionClass` to `element_shape`/`basis_order`.

---

## Open Items

- **Simplex quadrature**: `_tet_quadrature(order::Int)` function needs to be added to
  `fem.jl` before step 4 can be completed. Needs Gauss points on the reference tetrahedron
  `{ξ≥0, η≥0, ζ≥0, ξ+η+ζ≤1}` for orders 1 and 2.

- **`ndim` vs `element_shape` for 2D/3D tet**: Currently `ndim` determines whether the
  3D or 2D tet basis is called. After the refactor, `:Tet` means 3D and `:Tri` means 2D —
  the `ndim` field on the model is the authority, not the shape symbol. The shape symbol
  encodes topology family, not spatial dimension.

- **`PostProcess.jl` `eval_on_cylinder`**: Already updated to use the standard
  `basis_function` signature. Will pick up the new signature automatically in step 3.

---

## `fem.jl` Correctness Bugs and Future-Proofing

### Correctness bugs (fix before using Tri/Tet T1 elements)

- **`_basis_tri` T1: `ΔN = ones(3,3)` is wrong** (`fem.jl:243`)
  For linear triangles `N_i = λ_i`, so derivatives should be `copy(∇λ)` (3×2 matrix). The
  current `ones(3,3)` produces incorrect physical-space gradients. T2 already does this
  correctly via `∂N∂λ * ∇λ`.

- **`_basis_tet` T1: `ΔN = ones(4,4)` is wrong** (`fem.jl:291`)
  Same issue — T1 tet derivatives should be `copy(∇λ)` (4×3), not `ones(4,4)`.

- **`_tri_quadrature` and `_tet_quadrature` bad defaults** (`fem.jl:77,97`)
  Both default to `n_gauss_pts=2` but neither implements that case (tri has 1/3/4, tet has
  1/4/5). Calling with the default returns undefined variables. Fix: change defaults to a
  valid value (3 for tri, 4 for tet) and add `else error(...)` fallthrough.

- **`_basis_tri` dead code: `∇2λ = zeros(3, 3)`** (`fem.jl:239`)
  Assigned but never used anywhere in the function. Remove it.

- **2D `basis_function` ignores `element_shape`** (`fem.jl:368`)
  `basis_function(ξ, η, nothing, :Tri, order)` silently computes tensor-product Quad basis
  instead of triangle. The `:Tri` branch only exists in the 3-arg overload. Add a guard:
  `element_shape ∉ (:Quad, :Line) && error(...)` or add a proper `:Tri` branch that
  dispatches to `_basis_tri`.

### Future-proofing improvements

- **`gaussian_quadrature` missing fallthrough error** (`fem.jl:65`)
  Only handles `n=2` and `n=3`. Any other value silently returns undefined `ξ` and `w`.
  Add `else error("Unsupported n_gauss_pts: $n_gauss_pts")`.

- **Element-shape branch inside the GP loop** (`fem.jl:503`, `fem.jl:572`)
  The `if element_shape == :Tet / :Hex` check in both `get_basis_volume_functions` and
  `get_surface_basis_functions` is a static predicate re-evaluated on every Gauss point.
  Should be factored outside the loop. Superseded if `get_quadrature` is introduced (see below).

- **`ndim_cached` assigned but never used** (`fem.jl:454`, `fem.jl:545`)
  Extracted in both `get_basis_volume_functions` and `get_surface_basis_functions` but
  never read. Remove.

- **`BasisFunctionCache.bf_surface` is `NTuple{3, Any}`** (`fem.jl:16`)
  `get_surface_basis_functions` now computes N_p and N_x surface arrays but returns only
  the u tuple (3 values). When those are wired into the return the tuple size must change
  to `NTuple{7, Any}`. Update before callers multiply.

- **Centralise quadrature: introduce `get_quadrature(element_shape, n_gauss_pts)`**
  Every assembly function that needs Gauss points replicates the same `if :Tet / :Hex`
  branching inline. Adding a new element type currently requires editing every one of these
  sites. A single helper:
  ```julia
  get_quadrature(element_shape::Symbol, n_gauss_pts::Int) -> (points, weights)
  ```
  where `points` is a `Vector{NTuple}` of reference coordinates and `weights` is a
  `Vector{Float64}`, makes adding a new shape a one-place change and removes all
  duplicated branching from `get_basis_volume_functions`, `get_surface_basis_functions`,
  and the dense assemblers in `squeeze_stokes.jl`. This also naturally fixes the bad
  defaults and missing errors in `_tri_quadrature`, `_tet_quadrature`, and
  `gaussian_quadrature` — they become internal, called only through `get_quadrature`.
  **Do this before adding any further element types.**

- **Expand `BasisFunctionCache.bf_surface` from `NTuple{3}` to `NTuple{7}`** (`fem.jl:16`)
  `get_surface_basis_functions` already computes N_p, ΔN_p, N_x, ΔN_x surface arrays but
  currently returns only the 3-tuple `(N_u_surf_gp, ΔN_u_surf_gp, wpoints)`. Callers in
  `squeeze_stokes.jl` unpack exactly 3 values. To expose p/x surface functions for
  pressure/geometry BC assembly the return must expand to the 7-tuple
  `(N_u, ΔN_u, N_p, ΔN_p, N_x, ΔN_x, wpoints)` matching `bf_volume`, and the struct
  field and all unpack sites updated together. **Do before surface p/x integration is
  needed to avoid a scattered breaking change.**
