# Examples

## Squeeze between two plates

```julia
using smearFEM

mesh = meshgrid_cylinder(25.0, 40.0; mesh_type=:unstructured, element_shape=:Hex, basis_order=2,
                         elem_size=6.0, mesh_path=joinpath(pkgdir(smearFEM), "mesh_files"))
λ, μ = lame_parameters(10.0, 0.35)                  # kPa
mdl = LinearElasticity(mesh=mesh, λ=λ, μ=μ)

r = simulate(mdl, fill(-0.05, 100); control="displacement", β=1.0, sensitivities=true)
```

## Volume update

```julia
reset_model!(mdl)
X = mesh.NodeList
side = mesh.side_nodes
normals = vcat(X[1:2, side] ./ sqrt.(sum(abs2, X[1:2, side], dims=1)), zeros(1, length(side)))
u = volume_update!(mdl, side, normals, fill(0.5, length(side)); β_table=1e-3)
```
