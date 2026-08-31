# smearFEM.jl — Coding Guide & Architecture Reference

## Overview

smearFEM.jl is a Julia finite element library for soft-body mechanics and Stokes flow simulation, developed for the SMEAR robotics project. It targets **cylinder-shaped soft robots** undergoing squeeze flow, with forward simulation, parameter identification, and optimization built in.

The library is structured around:
- A **mesh abstraction layer** that wraps gmsh-generated or analytically-defined grids
- A **FEM core** that computes basis functions and assembles sparse systems
- **Physics models** (linear elasticity and Stokes flow)
- **I/O and visualization** via VTK, HDF5, CSV, and JSON
- An **optimization layer** implementing Gauss-Newton and Levenberg-Marquardt

---

## Directory Layout

```
smearFEM.jl/
├── src/
│   ├── smearFEM.jl          # Module entry: exports, includes, BLAS config
│   ├── config.jl            # Path resolution (data, mesh, scratch dirs)
│   ├── utils.jl             # Shared math utilities
│   ├── fem/
│   │   ├── Meshes.jl        # Mesh struct definitions + constructors
│   │   ├── models.jl        # AbstractModel, LinearElasticity, Stokes
│   │   ├── fem.jl           # Basis functions, quadrature, assembly
│   │   ├── PostProcess.jl   # Post-processing: fitting, projection, NURBS
│   │   ├── types.jl         # Conditions, EnvConditions
│   │   └── scenarios.jl     # SqueezeFlow scenario type
│   ├── io/
│   │   ├── io.jl            # VTK, HDF5, CSV, JSON read/write
│   │   ├── gmsh_utils.jl    # gmsh mesh loading and generation
│   │   ├── plotting.jl      # Makie-based mesh and field visualization
│   │   └── analysis_plots.jl# Slip/viscosity/height analysis plots
│   ├── optimization/
│   │   └── smearOptimize.jl # GN and LM optimizers, cylinder fitting
│   └── examples/
│       ├── squeeze_stokes.jl# Stokes squeeze flow simulation pipeline
│       └── run_example.jl   # Top-level runner and data management
├── mesh_files/
│   ├── templates/           # Canonical .geo templates (version-controlled)
│   │   ├── cylinder_tet.geo
│   │   ├── cylinder_hex.geo
│   │   ├── box_tet.geo
│   │   ├── box_hex.geo
│   │   ├── square_tri.geo
│   │   ├── square_quad.geo
│   │   ├── disk_tri.geo
│   │   ├── disk_quad.geo
│   │   └── line.geo
│   ├── cylinder/            # Generated cylinder .geo and .msh files
│   └── ...                  # Other shape subdirectories (auto-generated)
└── docs/
```

---

## Type Hierarchy

### Meshes

```
AbstractMeshgrid
├── Meshgrid1D          # Generic 1D mesh (Gmsh-loaded)
├── Meshgrid2D          # Generic 2D mesh (Gmsh-loaded)
├── Meshgrid3D          # Generic 3D mesh (Gmsh-loaded)
├── MeshgridLine        # Analytical 1D line mesh
├── MeshgridSquare      # Analytical/Gmsh square with named faces
├── MeshgridCube        # Analytical/Gmsh box with named faces + initial_state
└── MeshgridCylinder    # Gmsh cylinder: carries precomputed BasisFunctionCache (C_vol, W, etc.)
```

All mesh structs are **mutable** with **keyword constructors** that set safe defaults for every field. The wrapper `Meshgrid{T}` exists but is not widely used — work directly with the concrete type.

Every mesh carries:
- `NodeList::Matrix{Float64}` — `ndim × nNodes` coordinate matrix
- `IEN::Matrix{Int}` — `nNodesPerElem × nElems` connectivity
- `ID::Matrix{Int}` — DOF numbering map
- `element_shape::Symbol` — one of `:Hex`, `:Tet`, `:Quad`, `:Tri`, `:Line`
- `basis_order::Int` — polynomial order (1, 2, or 3)
- `nNodes::Int`, `ne::Int`

Boundary-aware types (Square, Cube, Cylinder) also carry `IEN_top`, `IEN_bottom`, named `IEN_sides`, `top_nodes`, `bottom_nodes`, `side_nodes`.

`MeshgridCylinder` additionally carries precomputed quadrature arrays: `C_vol`, `C_top`, `C_btm`, `C_vol_vis`, `W`. These are expensive to compute, so they are stored once at construction time.

### Models

```
AbstractModel
├── LinearElasticity    # Single-mesh elasticity; holds material params θ1, θ2
└── Stokes              # Three-mesh Stokes flow; holds η::Vector{Float64}
```

`Stokes` holds three meshes for a mixed-element discretization:
- `mesh_x` — geometry / reference coordinates (Q1 or Q2)
- `mesh_u` — velocity field (typically P2 / Q2)
- `mesh_p` — pressure field (typically P1 / Q1)

The wrapped type `Model{T<:AbstractModel}` exists for generic dispatch.

---

## Element Shape and Basis Order

Elements are described by two independent fields everywhere in the codebase:

| Field | Type | Values | Meaning |
|---|---|---|---|
| `element_shape` | `Symbol` | `:Hex`, `:Tet`, `:Quad`, `:Tri`, `:Line` | Geometric topology |
| `basis_order` | `Int` | `1`, `2`, `3` | Polynomial order of shape functions |

Common combinations: Q1 Hex (`:Hex`, 1), Q2 Hex (`:Hex`, 2), T2 Tet (`:Tet`, 2), Q2 Quad (`:Quad`, 2), T1 Tri (`:Tri`, 1).

These replaced the old `FunctionClass::String` field (which concatenated shape and order, e.g. `"Q2"`).

---

## FEM Core (`fem/fem.jl`)

### Basis Functions

`basis_function` is overloaded for 1D, 2D, and 3D:

```julia
basis_function(ξ, element_shape, basis_order)           # 1D: returns (N, ΔN, Δ2N)
basis_function(ξ, η, element_shape, basis_order)         # 2D
basis_function(ξ, η, ζ, element_shape, basis_order)      # 3D
```

Hex/Quad elements use tensor-product construction via internal helpers:
- `_basis_1d_components` — returns 1D Lagrange values + derivatives
- `_tensor_basis_2d` / `_tensor_basis_3d` — outer-product assembly using precomputed `_Q1_2D_NODE_PAIRS`, `_Q2_3D_NODE_TRIPLES`, etc.

Tet/Tri elements use barycentric (simplex) coordinates.

### Quadrature

```julia
gaussian_quadrature(a, b, n_gauss_pts)   # Gauss-Legendre on [a,b], n=2 or 3
_tri_quadrature(n_gauss_pts)             # Barycentric tri quadrature
_tet_quadrature(n_gauss_pts)             # Barycentric tet quadrature
```

### BasisFunctionCache

For performance, basis functions at all Gauss points are precomputed once and stored:

```julia
cache = BasisFunctionCache(mdl::Stokes)
# cache.bf_volume  → (N_u_gp, ΔN_u_gp, N_p_gp, ΔN_p_gp, N_x_gp, ΔN_x_gp, wpoints)
# cache.bf_surface → (N_u_surf_gp, ΔN_u_surf_gp, wpoints)
```

Construct once before the assembly loop; pass the cache into assembly functions. Never recompute inside a time-stepping or optimization loop.

### Sparse Assembly Pattern

System matrices are assembled in COO (coordinate) format:

```julia
E = Int[];  J = Int[];  V = Float64[]
# ... element loop ...
push!(E, row_idx)
push!(J, col_idx)
push!(V, value)
# ...
K = sparse(E, J, V)
```

This avoids repeated sparse indexing (which is O(nnz) per access). The COO vectors are pre-sized when the element count is known.

---

## @unpack + `_cached` Pattern

Parameters.jl `@unpack` extracts mesh struct fields into local variables for JIT-friendly access. The convention is to suffix these locals with `_cached`:

```julia
@unpack element_shape, basis_order, IEN, IEN_cp, ID, NodeList, C_vol, W = mdl.mesh_x
element_shape_x_cached = element_shape
basis_order_x_cached   = basis_order
IEN_cached             = IEN
```

The `_cached` suffix signals to readers that these are snapshot locals, not live struct fields. Mutating `IEN_cached` does not affect `mdl.mesh_x.IEN`.

---

## Mesh Construction

### Gmsh-Based Constructors

High-level constructors handle the full pipeline: template → `.geo` → `.msh` → mesh struct.

```julia
# Auto-generates mesh if .msh file is missing
mesh = get_gmsh_cylinder(file_path, nDof, r, h, element_shape, basis_order, ne)
mesh = get_gmsh_box(file_path, nDof, lx, ly, lz, element_shape, basis_order, ne)
mesh = get_gmsh_square(file_path, nDof, lx, ly, element_shape, basis_order, ne)
mesh = get_gmsh_disk(file_path, nDof, r, element_shape, basis_order, ne)
mesh = get_gmsh_line(file_path, nDof, lx, element_shape, basis_order, ne)
```

These call `get_mesh_data` from `io/gmsh_utils.jl`, which:
1. Checks if the `.msh` file exists; if not, generates it from the template
2. Loads the mesh via the Gmsh Julia API (thread-safe via `GMSH_LOCK`)
3. Remaps node ordering to smearFEM convention (critical for Q2 hex: 27-node map)
4. Returns `(NodeList, IEN, face_IENs, nNodes, node_sets, ne)`

### `get_mesh_data` Keyword Interface

```julia
NodeList, IEN, face_IENs, nNodes, node_sets, ne = get_mesh_data(
    file_path;
    params          = Dict("radius" => 25.0, "height" => 40.0, "elem_size_2d" => 10.0),
    template_path   = joinpath(mesh_dir, "templates", "cylinder_tet.geo"),
    mesh_order      = 2,       # gmsh element order
    mesh_dim        = 3,       # 1, 2, or 3
    body_group      = "Volume",
    body_dim        = 3,
    face_groups     = ["Top", "Bottom", "Lateral"],
    node_set_groups = ["Lateral", "Bottom", "Top"]
)
```

`face_IENs` and `node_sets` are `Vector`s aligned with `face_groups` and `node_set_groups` respectively. Entries are `nothing` if the group was not found in the mesh.

### Template Files

Templates live in `mesh_files/templates/` and follow the naming convention `{shape}_{topology}.geo`:

| Template | Shape | Topology |
|---|---|---|
| `cylinder_tet.geo` | cylinder | tet mesh |
| `cylinder_hex.geo` | cylinder | hex mesh |
| `box_tet.geo` / `box_hex.geo` | box | tet / hex |
| `square_tri.geo` / `square_quad.geo` | square | tri / quad |
| `disk_tri.geo` / `disk_quad.geo` | disk | tri / quad |
| `line.geo` | line | — |

Templates use named physical groups matching the constructor's `body_group` and `face_groups` arguments. Cylinder and box templates use `Physical Volume("Volume")`. Shape-specific names are used for 2D: `"Square"`, `"Disk"`, `"Line"`.

Generated `.msh` files go in `mesh_files/{shape}/` subdirectories and are **not version-controlled** (only templates are).

### Analytical Constructors

For structured meshes without gmsh:

```julia
meshgrid_line(lx, ne; element_shape=:Line, basis_order=1)    # → MeshgridLine
meshgrid_square(lx, ly, ne; ...)                              # → MeshgridSquare
meshgrid_cube(lx, ly, lz, ne; ...)                            # → MeshgridCube
inflate_cylinder(r, h, ne; ...)                               # → MeshgridCylinder
```

---

## I/O

### Configuration and Path Resolution (`config.jl`)

Three-tier fallback for all data paths:

1. **Environment variable** (`SMEAR_DATA_DIR`, `SMEAR_MESH_DIR`, `SMEAR_SCRATCH_DIR`)
2. **`config.toml`**, searched for from the project root upwards — it sits in the shared
   `smear-modules` directory, so smearPerception's `src/config.py` reads the same file
3. **Hard default** (`~/SMEAR-Data`, `~/SMEAR-Data/meshes`, `/tmp/smear`)

Use the resolution functions rather than hardcoding paths:

```julia
get_data_dir()                        # root data directory
get_mesh_dir()                        # mesh subdirectory
get_scratch_dir()                     # temporary/output directory
resolve_data_path("ground_truth/...")
resolve_mesh_path("cylinder/cyl.msh")
create_output_dir("my_experiment")
```

Never hardcode absolute paths in source files.

### Output Formats

| Format | Function | Use case |
|---|---|---|
| VTK | `write_vtk`, `write_scene` | Paraview visualization of fields |
| HDF5 | `write_data`, `read_h5` | Bulk simulation data, large arrays |
| JSON | `write_json`, `read_json` | Metadata, small structured data |
| CSV | `write_csv`, `read_csv` | Tabular data, parameter sweeps |

VTK output uses WriteVTK.jl. `write_scene` writes mesh geometry and optionally field data in a single call.

### Gmsh Thread Safety

All Gmsh API calls must go through `GMSH_LOCK`:

```julia
lock(GMSH_LOCK) do
    gmsh.initialize()
    # ... gmsh calls ...
    gmsh.finalize()
end
```

Do not call `gmsh.*` outside this lock — the Gmsh C library is not thread-safe.

---

## Physics Examples

### Stokes Squeeze Flow (`examples/squeeze_stokes.jl`)

The main simulation pipeline:

```
def_problem(r, h, ne, η_0, ndim, element_shape_u, basis_order_u, ...) → Stokes
set_model(r, h, ne, η, ndim, ...; GMESH_MESH, filepath_mesh) → Stokes
simulate(model, conditions, cache) → (u, p)
apply_boundary_conditions(model, u) → u_constrained
```

`def_problem` defines the physics (boundary conditions, DOF numbering).  
`set_model` constructs three meshes (u, p, x) and wires them into a `Stokes` struct.  
`simulate` runs one or more time steps.

### Mixed P2-P1 Discretization

Stokes flow uses mixed elements for inf-sup stability:
- **Velocity** (`mesh_u`): Q2 hex (27 nodes/element) or T2 tet (10 nodes/element)
- **Pressure** (`mesh_p`): Q1 hex (8 nodes/element) or T1 tet (4 nodes/element)
- **Geometry** (`mesh_x`): matches velocity mesh topology but with `nDof=1`

### Node Ordering Convention

Gmsh and smearFEM use different node orderings for higher-order elements. `get_mesh_data` applies a remapping:
- **Q2 hex (27 nodes)**: explicit 27-entry permutation vector
- **T2 tet (10 nodes)**: `[1,2,3,4,5,6,7,8,10,9]` (swaps last two midpoint nodes)
- All others: identity (no remap)

If you add a new element type, add the remap in the `body_map` block in `get_mesh_data`.

---

## Optimization (`optimization/smearOptimize.jl`)

Implements parameter identification via:
- **Gauss-Newton (GN)**: `closest_point`, `match_points`
- **Levenberg-Marquardt (LM)**: `fit_model`

The optimizer is initialized with `init_cylinder()`, which builds a reference mesh using `meshgrid_cube`.

---

## Coding Style

### Naming

| Thing | Convention | Example |
|---|---|---|
| Types | `PascalCase` | `MeshgridCylinder`, `LinearElasticity` |
| Functions | `snake_case` | `get_mesh_data`, `basis_function` |
| Constants | `UPPER_SNAKE_CASE` | `GMSH_LOCK`, `_Q2_3D_NODE_TRIPLES` |
| Internal helpers | leading `_` | `_basis_1d_components`, `_get_elements_for_physical` |
| Cached locals | `_cached` suffix | `IEN_cached`, `element_shape_x_cached` |
| Keyword args with defaults | always prefer over positional | `function f(; a=1, b=2)` |

### Structs

All structs are `mutable` with keyword-only inner constructors and safe defaults for every field. This lets callers build partial objects and fill fields later:

```julia
mutable struct Meshgrid2D <: AbstractMeshgrid
    NodeList::Matrix{Float64}
    ...
    function Meshgrid2D(;
        NodeList::Matrix{Float64}=Matrix{Float64}(undef, 2, 1),
        ...
    )
        new(NodeList, ...)
    end
end
```

### Docstrings

All exported functions require a docstring. Use the standard Julia format:

```julia
"""
    function_name(arg1, arg2; kwarg=default) -> ReturnType

One-line summary.

# Arguments
- `arg1::Type`: description
- `arg2::Type`: description

# Returns
- `ReturnType`: description

# Example
```julia
result = function_name(x, y)
```
"""
```

Non-exported helpers (`_` prefix) do not need docstrings unless the logic is non-obvious.

### Comments

Only comment the **why**, never the **what**. Well-named functions and variables are self-documenting. Acceptable comments:
- A non-obvious invariant (e.g., node ordering convention)
- A workaround for a specific external library behavior
- A constraint that is not visible from the types alone

Not acceptable: restate what the next line does, reference the PR/issue that introduced the code, or describe callers.

### Logging

Use Julia's `@debug`, `@info`, `@warn`, `@error` macros — never `println` in library code:

```julia
@info "Generated $(basename(output_path)) from $(basename(template_path))"
@debug "get_mesh_data: NodeList size=$(size(NodeList))"
@warn "gmsh executable not found in PATH"
@error "Failed to run gmsh: $e"
```

The module configures a stderr logger at load time so log output does not interfere with stdout data streams.

### Performance

- Pre-allocate inside hot loops; avoid `push!` inside element loops when size is known
- Use `@unpack` + `_cached` locals before any inner loop
- Construct `BasisFunctionCache` once per simulation, not per time step
- Use COO format for sparse assembly; call `sparse(E, J, V)` once at the end
- BLAS threading is set to `Threads.nthreads()` at module load — do not override inside library code

### No Absolute Paths

Never write absolute paths in source files. Use `get_data_dir()`, `get_mesh_dir()`, `resolve_data_path()`, or `joinpath(@__DIR__, ...)` relative to a known anchor.

---

## Adding a New Geometry

1. Write a `.geo` template in `mesh_files/templates/{shape}_{topology}.geo`
   - Name physical groups consistently (`Physical Volume("Volume")` for 3D bodies)
   - Expose all variable parameters at the top of the file

2. Add a constructor `get_gmsh_{shape}` in `fem/Meshes.jl`
   - Follow the pattern of `get_gmsh_cylinder` or `get_gmsh_square`
   - Call `get_mesh_data` with the appropriate `body_group`, `face_groups`, `node_set_groups`
   - Return the appropriate mesh struct type

3. Export the constructor from `smearFEM.jl`

4. Add the template path to the constructor using `joinpath(mesh_dir, "templates", "{name}.geo")`

Generated instances go in `mesh_files/{shape}/` — add that directory to `.gitignore`.

---

## Adding a New Element Type

1. Add quadrature support in `fem.jl`: implement `_xxx_quadrature` and wire it into `gaussian_quadrature` or a new dispatch
2. Add `basis_function` overloads for the new topology
3. Add node-pair/triple constants if it is a tensor-product element
4. Add the node-ordering remap in the `body_map` block of `get_mesh_data` in `gmsh_utils.jl`
5. Update `element_shape` handling in assembly routines

---

## Dependencies

| Package | Purpose |
|---|---|
| `LinearAlgebra` | Dense linear algebra, BLAS |
| `SparseArrays` | COO and CSC sparse matrices |
| `Parameters` | `@unpack` macro |
| `ArgCheck` | Argument validation in constructors |
| `LinearSolve` | Unified sparse direct/iterative solvers |
| `Gmsh` | Mesh loading via Gmsh Julia API |
| `WriteVTK` | VTK output for Paraview |
| `HDF5` | Binary data storage |
| `TOML` | Config file parsing |
| `DataInterpolations` | 1D interpolation in post-processing |
| `ConvexHulls2d` | 2D convex hull for boundary detection |
