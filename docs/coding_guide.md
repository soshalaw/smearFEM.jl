# smearFEM.jl — Coding Guide & Architecture Reference

## Overview

smearFEM.jl is a Julia finite element library developed for the SMEAR robotics project. It
simulates the **squeeze flow** of a soft object compressed between two plates, renders its
deforming silhouette through a camera model, and inverts that pipeline to recover the material
parameters `θ = [η, β]` from observed border contours.

The library is structured around:
- **Geometry wrapper types** that select meshing and problem setup by dispatch
- A **mesh layer** with structured (pure Julia) and unstructured (Gmsh) generators
- A **FEM core** that computes basis functions and quadrature
- A **mixed-element Stokes solver** that carries analytic parameter sensitivities
- An **optimization layer**: Gauss-Newton (optionally Tikhonov-regularized) and
  Levenberg-Marquardt over residual-based contour costs
- **I/O, plotting and analysis** via VTK, HDF5, CSV, JSON and Plots.jl

---

## Directory Layout

```
smearFEM.jl/
├── src/
│   ├── smearFEM.jl              # Module entry: single export list, includes, BLAS config
│   ├── config.jl                # Path resolution (data, mesh, scratch dirs)
│   ├── utils.jl                 # Shared helpers (camera matrix, logging, time log)
│   ├── fem/
│   │   ├── models.jl            # AbstractModel, Model, LinearElasticity, Stokes
│   │   ├── fem.jl               # Basis functions, quadrature, BasisFunctionCache
│   │   ├── geometries.jl        # Cylinder, Cuboid, Disk, Square, Segment
│   │   ├── meshes/
│   │   │   ├── Meshes.jl        # Public meshgrid_* entry points, reset_mesh!
│   │   │   ├── mesh_types.jl    # Meshgrid* structs
│   │   │   ├── mesh_structured.jl  # Pure-Julia structured generators (_meshgrid_*)
│   │   │   └── mesh_gmsh.jl     # Gmsh-backed generators (_get_gmsh_*)
│   │   ├── PostProcess.jl       # Camera projection, border extraction, curve fitting
│   │   ├── types.jl             # EnvConditions, Conditions
│   │   └── scenarios.jl         # AbstractScenario, SqueezeFlow
│   ├── solver/
│   │   ├── stokes_solver.jl     # Assembly, boundary conditions, time stepping (simulate)
│   │   ├── squeeze_system.jl    # Fixed-pattern saddle-point system and lagged solves behind simulate
│   │   ├── sparse_lu.jl         # Sparse factorization: MKL Pardiso, UMFPACK fallback
│   │   ├── stokes_setup.jl      # def_problem / set_model per geometry
│   │   └── run_example.jl       # High-level drivers (write_sim_data, write_gt_data)
│   ├── optimization/
│   │   ├── smearOptimize.jl     # fit_model: GN, GN-Tikhonov, LM; line searches
│   │   ├── cost_functions.jl    # ContourCost types and contour_cost
│   │   └── window_estimation.jl # Moving-window estimation and prediction
│   ├── io/
│   │   ├── io.jl                # VTK, HDF5, CSV, JSON read/write
│   │   ├── gmsh_utils.jl        # Gmsh loading, template generation, node remapping
│   │   ├── plotting.jl          # Plots.jl mesh and field visualization
│   │   └── analysis_plots.jl    # Slip/viscosity/height analysis plots
│   └── analysis/
│       ├── pointcloud_metrics.jl   # Hausdorff / Chamfer / closest-point metrics
│       ├── experiment_tree.jl      # Result-tree traversal (ExpLeaf, ExpGroup)
│       └── statistical_analysis.jl # Cross-replicate statistics
├── experiments/                 # Driver scripts (not part of the package)
├── mesh_files/
│   ├── templates/               # Canonical .geo templates (version-controlled)
│   └── <geometry>/<Shape>/      # Generated .geo/.msh cache (git-ignored)
├── test/
└── docs/
```

`src/smearFEM.jl` is a flat include list: every file shares one namespace and one export list.
When adding a file, add both the `include(...)` and its `export` line there.

---

## Type Hierarchy

### Geometries

```
AbstractGeometry
├── Cylinder(r, h)                         # 3D
├── Cuboid(lx, ly, lz, edge_radius)        # 3D; edge_radius=nothing → sharp vertical edges
├── Disk(r)                                # 2D
├── Square(lx, ly)                         # 2D
└── Segment(l)                             # 1D
```

The geometry type selects `set_model` (meshing) and `def_problem` (full squeeze-flow setup) by
dispatch. `def_problem` currently exists for `Cylinder` and `Cuboid` only. Add a new geometry as
a wrapper type plus a method pair — never by branching inside an existing method.

### Meshes

```
AbstractMeshgrid
├── MeshgridLine        # 1D; carries boundary_nodes instead of face IENs
├── MeshgridDisk        # 2D
├── MeshgridSquare      # 2D
├── MeshgridCuboid      # 3D
└── MeshgridCylinder    # 3D
```

Every mesh carries:
- `NodeList::Matrix{Float64}` — `ndim × nNodes` coordinates
- `IEN::Matrix{Int}` — `nNodesPerElem × ne` connectivity
- `ID::Matrix{Int}` — DOF numbering
- `volume_element_shape::Symbol` — `:Hex`, `:Tet`, `:Quad`, `:Tri` or `:Line`
- `basis_order::Int` — 1 or 2
- `nNodes::Int` (total node count), `ne::Int`
- `initial_state::Matrix{Float64}` — copy of `NodeList` used by `reset_mesh!`
- `effect_elem_sze` — realised element size

2D and 3D meshes also carry `surface_element_shape`, face connectivities (`IEN_top`,
`IEN_bottom`, `IEN_sides`) and node sets (`top_nodes`, `bottom_nodes`, `side_nodes`).

`nNodes` is the **total** node count. Size DOF vectors as `nDof * nNodes`, never
`nDof * nNodes^ndim` (a leftover from the old nodes-per-edge convention).

### Models and scenarios

```
AbstractModel
├── LinearElasticity    # Single-mesh elasticity
└── Stokes              # Three-mesh Stokes flow; η::Vector{Float64}

AbstractScenario
└── SqueezeFlow         # β, Dirichlet operators, control history, viscosity type

EnvConditions
└── Conditions          # Camera, object pose, viewing angles, output flags
```

`Stokes` holds three meshes for a mixed discretization, all built from the same geometry:
- `mesh_u` — velocity (typically quadratic)
- `mesh_p` — pressure (typically linear)
- `mesh_x` — geometry / node tracking; velocity topology with `nDof=1`

### Contour costs

```
ContourCost
├── ClosestPointCost
├── ChamferCost
└── SignedDistanceCost
```

Costs are residual-based: a new cost is a `_residuals` method, not a scalar metric, because the
Gauss-Newton step needs the residual vector and its Jacobian.

---

## Element Shape and Basis Order

| Field | Type | Values | Meaning |
|---|---|---|---|
| `volume_element_shape` | `Symbol` | `:Hex`, `:Tet`, `:Quad`, `:Tri`, `:Line` | Geometric topology |
| `basis_order` | `Int` | `1`, `2` | Polynomial order of shape functions |

Constructors and public functions take the shape as an `element_shape` keyword. Shape and order
are always two independent values, never a combined `"Q2"`-style string.

---

## FEM Core (`fem/fem.jl`)

### Basis Functions

`basis_function` dispatches on the number of reference coordinates:

```julia
basis_function(ξ, nothing, nothing, :Line, order)    # 1D
basis_function(ξ, η, nothing, :Quad, order)          # 2D tensor-product
basis_function(ξ, η, ζ, :Hex, order)                 # 3D tensor-product
basis_function(ξ, η, ζ, θ, :Tet, order)              # simplex (barycentric)
```

Hex/Quad use tensor-product construction (`_basis_1d_components`, `_tensor_basis_2d/3d`, and
the `_Q1_2D_NODE_PAIRS` / `_Q2_3D_NODE_TRIPLES` tables). Tet/Tri use barycentric coordinates.

### Quadrature

```julia
gaussian_quadrature(a, b, n_gauss_pts)   # Gauss-Legendre on [a,b]
_tri_quadrature(n_gauss_pts)             # simplex, weights sum to √3/4
_tet_quadrature(n_gauss_pts)             # simplex, weights sum to 1/(6√2)
```

The simplex gradients `∇λ` are defined on the **equilateral / regular unit-edge** master
element, so simplex weights must sum to its measure (`TRI_MASTER_AREA`, `TET_MASTER_VOLUME`),
not the right-simplex 1/2 and 1/6 of most published tables. `test/verify_meshes_test.jl` guards
this.

### BasisFunctionCache

```julia
cache = BasisFunctionCache(mdl::Stokes)
# cache.bf_volume  → (N_u, ΔN_u, N_p, ΔN_p, N_x, ΔN_x, wpoints) at volume Gauss points
# cache.bf_surface → the same tuple at surface Gauss points
```

Build once per simulation, before any time-stepping or optimization loop, and pass it into
assembly. Never recompute it inside a loop.

### Sparse Assembly Pattern

System matrices are assembled in COO format:

```julia
E = Int[];  J = Int[];  V = Float64[]
# ... element loop ...
push!(E, row_idx); push!(J, col_idx); push!(V, value)
# ...
K = sparse(E, J, V)
```

This avoids incremental sparse indexing (O(nnz) per insertion). `sizehint!` the vectors when the
entry count is known.

---

## @unpack + `_cached` Pattern

Struct fields are pulled into typed locals with Parameters.jl's `@unpack` before a hot loop, and
suffixed `_cached` to mark them as snapshots rather than live struct fields:

```julia
@unpack volume_element_shape, basis_order, IEN, ID, NodeList = mdl.mesh_x
element_shape_x_cached = volume_element_shape
IEN_x_cached::Matrix{Int} = IEN
```

---

## Mesh Construction

### Public entry points

```julia
meshgrid_cylinder(r, h; mesh_type=:structured, element_shape=:Hex, basis_order=1, ne, ...)
meshgrid_cuboid(lx, ly, lz; ...)
meshgrid_square(lx, ly; ...)
meshgrid_disk(r; ...)
meshgrid_line(l; ...)
```

`mesh_type=:structured` calls the pure-Julia `_meshgrid_*` generators (needs `ne`);
`mesh_type=:unstructured` calls the Gmsh-backed `_get_gmsh_*` generators (takes `elem_size`).

### Gmsh pipeline

`_get_gmsh_*` call `_get_mesh_data` (`io/gmsh_utils.jl`), which:
1. Generates the `.msh` from a `.geo` template if it is missing
2. Loads it through the Gmsh API under `GMSH_LOCK`
3. Remaps node ordering to the smearFEM convention
4. Returns `(NodeList, IEN, face_IENs, nNodes, node_sets, ne)`, with `face_IENs` / `node_sets`
   aligned to the `face_groups` / `node_set_groups` keywords (`nothing` for a missing group)

Generated meshes are cached in the repository under `mesh_files/<geometry>/<Shape>/`. The cache
is **not** invalidated automatically — delete it after editing a template. (`get_mesh_dir()` is
resolved for tools sharing `config.toml`; it does not relocate this cache.)

### Templates

Templates live in `mesh_files/templates/` as `{shape}_{topology}.geo`: `cylinder_{hex,tet}`,
`box_{hex,tet}` (plus `_sharp` variants without edge fillets), `square_{quad,tri}`,
`disk_{quad,tri}`, `line`. 3D templates use `Physical Volume("Volume")`; 2D/1D use the shape
name (`"Square"`, `"Disk"`, `"Line"`). See `mesh_files/USAGE.md` for element-count calibration.

### Node ordering

Gmsh and smearFEM order higher-order element nodes differently. `_get_mesh_data` remaps:
- **Q2 hex (27 nodes)**: explicit 27-entry permutation
- **T2 tet (10 nodes)**: `[1,2,3,4,5,6,7,8,10,9]`
- All others: identity

A new element type needs its own entry in the `body_map` block.

### Gmsh thread safety

All Gmsh API calls must go through `GMSH_LOCK`; the Gmsh C library is not thread-safe.

```julia
lock(GMSH_LOCK) do
    gmsh.initialize()
    # ... gmsh calls ...
    gmsh.finalize()
end
```

---

## Squeeze-Flow Pipeline

```julia
model, scene = def_problem(geom, ne, η_0,
                           shape_u, order_u, nDof_u,
                           shape_p, order_p, nDof_p,
                           shape_x, order_x,
                           β, cParam, control, viscosity_type, sim_time, t_steps)
output, gradList, borderPts2D, ... = simulate(model, scene, conditions)
stats = fit_model(model, scene, conditions, obsBorderPts, θ_init; method=:gn)
```

- `def_problem` builds the meshes via `set_model`, the boundary operators, and the
  `SqueezeFlow` scenario.
- `simulate` time-steps under `"force"` or `"velocity"` control and returns the state plus the
  analytic sensitivities `∂q/∂η`, `∂q/∂β`.
- `fit_model` dispatches on `method` (`:gn`, `:gn_tikhonov`, `:lm`) and returns a `Dict` with
  `"η"`, `"β"`, `"ηList"`, `"βList"`, `"cost_list"`, `"iterList"`.

Force is in **kg·mm/s²** (the solver works in mm), so 1 N = 1e3 solver units.

---

## Configuration and Paths (`config.jl`)

Three-tier fallback, highest priority first:

1. **Environment variable** (`SMEAR_DATA_DIR`, `SMEAR_MESH_DIR`, `SMEAR_SCRATCH_DIR`)
2. **`config.toml`**, searched for from the package root upwards — it sits in the shared
   `smear-modules` directory, so smearPerception's `src/config.py` reads the same file
   (`SMEAR_CONFIG_FILE` overrides the search)
3. **Default** (`~/SMEAR-Data`, `<data_dir>/meshes`, `/tmp/smear` or `~/SMEAR-Scratch` on Windows)

```julia
get_data_dir(); get_mesh_dir(); get_scratch_dir()
resolve_data_path("ground_truth/...")
create_output_dir("my_experiment")
```

Never write absolute paths in source files. Use the functions above, or `joinpath(@__DIR__, ...)`
relative to a known anchor.

### Output Formats

| Format | Function | Use case |
|---|---|---|
| VTK | `write_vtk`, `write_scene` | ParaView visualization of fields |
| HDF5 | `write_data`, `read_h5` | Bulk simulation data |
| JSON | `write_json`, `read_json` | Metadata, small structured data |
| CSV | `write_csv`, `read_csv` | Tabular data, parameter sweeps |

---

## Coding Style

### Naming

| Thing | Convention | Example |
|---|---|---|
| Types | `PascalCase` | `MeshgridCylinder`, `SqueezeFlow` |
| Functions | `snake_case` | `basis_function`, `fit_model` |
| Constants | `UPPER_SNAKE_CASE` | `GMSH_LOCK`, `TET_MASTER_VOLUME` |
| Internal helpers | leading `_` | `_basis_1d_components`, `_get_mesh_data` |
| Cached locals | `_cached` suffix | `IEN_cached`, `element_shape_x_cached` |
| Optional parameters | keyword args with defaults over positional | `f(x; tol=1e-6)` |

### Structs

Structs that hold model or run state are `mutable` with keyword-only inner constructors and safe
defaults for every field, so callers can build partial objects and fill fields later:

```julia
mutable struct Conditions <: EnvConditions
    ANIMATE::Bool
    filepath::String
    ...
    function Conditions(; ANIMATE::Bool=false, filepath::String="", ...)
        return new(ANIMATE, filepath, ...)
    end
end
```

Two deliberate exceptions:
- **Geometry wrappers** (`Cylinder`, `Cuboid`, ...) and `SqueezeFlow` take **positional**
  arguments, each with a default — every call site constructs them positionally.
- **Small value types** with no mutable state (`BasisFunctionCache`, the `ContourCost` types,
  `ExpLeaf`, `ExpGroup`) are immutable `struct`s.

### Docstrings

Every exported symbol requires a docstring:

```julia
"""
    function_name(arg1, arg2; kwarg=default) -> ReturnType

One-line summary.

# Arguments
- `arg1::Type`: description
- `kwarg::Type`: description (default: `default`)

# Returns
- `ReturnType`: description

# Example
```julia
result = function_name(x, y)
```
"""
```

Keyword arguments go in `# Arguments`; there is no separate `# Keyword Arguments` section.
`# Example` is optional. `_`-prefixed helpers do not need docstrings unless the logic is
non-obvious. This file is authoritative if `CONTRIBUTING.md` ever disagrees.

### Comments

Comment only the **why**, never the **what**. Acceptable:
- A non-obvious invariant (e.g. node ordering, quadrature master element)
- A workaround for an external library's behavior
- A constraint not visible from the types

Not acceptable: restating the next line, referencing the PR/issue that introduced the code, or
describing callers.

### Logging

Use `@debug`, `@info`, `@warn`, `@error` — never `println` or `printstyled` in library code:

```julia
@info "Iteration $iter"
@debug "NodeList size=$(size(NodeList))"
@warn "η: Below minimum ($η_min), clamped"
```

Enable debug output with `ENV["JULIA_DEBUG"] = "smearFEM"`. The module routes logging to stderr
at load time so it does not collide with stdout data streams.

Exceptions: functions whose sole purpose is formatted console output (`show_config`,
`print_replicate_table`) use `println`, since log prefixes break their layout; and
`println(io, ...)` against an open file handle (`write_time_log`) is file formatting, not
logging.

### Performance

- `@unpack` + typed `_cached` locals before any inner loop
- Construct `BasisFunctionCache` once per simulation, not per time step
- COO assembly with a single `sparse(E, J, V)` at the end
- BLAS threading is set to `Threads.nthreads()` at module load — do not override it in library
  code; run scripts with `--threads=auto`

### File size

Over 600 lines invites a review of whether the file should be split; over 1000 is a strong flag.

---

## Adding a New Geometry

1. Add a wrapper type in `fem/geometries.jl` (positional constructor, default for every field)
2. If meshed with Gmsh, add `mesh_files/templates/{shape}_{topology}.geo` with consistent physical
   group names and all parameters exposed at the top
3. Add a mesh struct in `mesh_types.jl`, structured / Gmsh generators, and a public
   `meshgrid_{shape}` in `Meshes.jl`
4. Add `set_model` (and, for squeeze flow, `def_problem`) methods in `solver/stokes_setup.jl`
5. Export the new symbols from `src/smearFEM.jl`

## Adding a New Element Type

1. Add quadrature in `fem.jl` (`_xxx_quadrature`); weights must sum to the measure of the master
   element the basis gradients are defined on
2. Add `basis_function` support for the new topology
3. Add node-pair/triple tables if it is a tensor-product element
4. Add the node-ordering remap in `_get_mesh_data`'s `body_map` block
5. Update `_surface_element_shape` and the assembly routines

---

## Testing

```julia
using Pkg; Pkg.test()      # full suite, ~98 000 checks, ~20 min
```

Individual files run with `julia --project -e 'include("test/<file>.jl")'`, except `qa.jl`
(Aqua), which needs the test environment. Test dependencies live in the root `Project.toml`'s
`[extras]` / `[targets]` — do not add a `test/Project.toml`, which would shadow them.

---

## Dependencies

All runtime dependencies are listed in `Project.toml`. The ones that shape the code:

| Package | Purpose |
|---|---|
| `LinearAlgebra`, `SparseArrays` | Dense / sparse linear algebra, COO assembly |
| `LinearSolve`, `IterativeSolvers` | Linear solvers |
| `Parameters` | `@unpack` |
| `ArgCheck` | Argument validation in constructors |
| `Gmsh` | Unstructured meshing via the Gmsh API |
| `NearestNeighbors` | KD-trees for contour matching and point-cloud metrics |
| `ConvexHulls2d`, `DataInterpolations` | Border extraction and contour resampling |
| `WriteVTK`, `HDF5`, `JLD2`, `CSV`, `DataFrames`, `JSON`, `TOML` | I/O and configuration |
| `Plots`, `StatsPlots`, `LaTeXStrings`, `Measures`, `Colors` | Plotting |
| `Distributions`, `Statistics` | Noise models and replicate statistics |
| `ProgressMeter` | Progress reporting |
