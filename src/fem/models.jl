# Abstract type for all models
abstract type AbstractModel end

"""
    Model{T}

Mutable box holding one `AbstractModel`, so a model can be swapped in place without rebuilding
the structures that reference it.

# Arguments
- `mdl::T`: The wrapped model.
"""
mutable struct Model{T<:AbstractModel}
    mdl::T
end

"""
    LinearElasticity(; mesh=MeshgridCylinder(), λ=0.0, μ=0.0)

Isotropic linear elastic material on one isoparametric 3D mesh: small-strain increments
solved on the current geometry. With lengths in mm and forces in kg·mm/s², stresses and
the Lamé parameters are in kPa.

# Arguments
- `mesh::AbstractMeshgrid`: The mesh; its `NodeList` is the current geometry, which `simulate`
  and `volume_update!` move.
- `λ::Real`: First Lamé parameter.
- `μ::Real`: Shear modulus.
"""
mutable struct LinearElasticity <: AbstractModel
    mesh::AbstractMeshgrid
    λ::Float64
    μ::Float64

    LinearElasticity(; mesh::AbstractMeshgrid=MeshgridCylinder(), λ::Real=0.0, μ::Real=0.0) =
        new(mesh, Float64(λ), Float64(μ))
end

"""
    lame_parameters(E, ν)

Lamé parameters of an isotropic material from its Young's modulus and Poisson's ratio.

# Arguments
- `E::Real`: Young's modulus.
- `ν::Real`: Poisson's ratio, in (-1, 0.5).

# Returns
- `λ::Float64`: First Lamé parameter, in the units of `E`.
- `μ::Float64`: Shear modulus, in the units of `E`.
"""
function lame_parameters(E::Real, ν::Real)
    -1 < ν < 0.5 || throw(ArgumentError("Poisson's ratio must lie in (-1, 0.5), got $ν"))
    return E * ν / ((1 + ν) * (1 - 2ν)), E / (2 * (1 + ν))
end

"""
    update_model!(model::LinearElasticity)

Make the current geometry the model's initial state, the one `reset_model!` returns to.

# Arguments
- `model::LinearElasticity`: The model; modified in place.
"""
update_model!(model::LinearElasticity) = update_initial_state!(model.mesh, model.mesh.NodeList)

"""
    reset_model!(model::LinearElasticity)

Move the model's mesh back to its initial state.

# Arguments
- `model::LinearElasticity`: The model; modified in place.
"""
reset_model!(model::LinearElasticity) = reset_mesh!(model.mesh)
