"""
    AbstractScenario

Supertype for the physical experiment a model is driven through. Concrete subtype: `SqueezeFlow`.
"""
abstract type AbstractScenario end

"""
    SqueezeFlow(model, β=Float64[], q=…, C_uc=…, control="force", sim_time=0.0, t_steps=0.0, viscosity_type="constant", cParam=Float64[])

An object compressed between two plates. The constructor takes the three prescribed-velocity
blocks as a positional vector and stores them keyed by surface, so callers index by name rather
than by position.

`model` is the one parameter without a default — it is abstractly typed, so no safe value
exists. Every other parameter defaults to an empty or zero-sized value, except `control` and
`viscosity_type`, which default to the most common configuration.

# Arguments
- `model::AbstractModel`: Model being driven.
- `β::Vector{Float64}`: Boundary slip/friction parameter.
- `q::Vector{SparseMatrixCSC{Float64, Int64}}`: Prescribed velocities, ordered top plate, bottom
  plate, lateral surface; stored as `q_d[:top]`, `q_d[:bottom]`, `q_d[:border]`. The lateral
  entry is all-zero — the side is a free surface, and `set_boundary_cond` never marks it.
- `C_uc::AbstractMatrix`: Constraint matrix mapping unconstrained to full degrees of freedom.
- `control::String`: Driving mode, e.g. constant velocity or constant force.
- `sim_time::Float64`: Total simulated time in seconds.
- `t_steps::Float64`: Time step size in seconds.
- `viscosity_type::String`: Viscosity law, e.g. `"constant"` or `"bulk_viscosity"`.
- `cParam::Vector{Float64}`: Control parameters; a force here is in kg*mm/s^2, not newtons.
"""
mutable struct SqueezeFlow <: AbstractScenario
    model::AbstractModel
    β::Vector{Float64}
    q_d::Dict{Symbol, Matrix{Float64}}
    C_uc::AbstractMatrix
    control::String
    sim_time::Float64
    t_steps::Float64
    viscosity_type::String
    cParam::Vector{Float64}

    function SqueezeFlow(
        model::AbstractModel,
        β::Vector{Float64}=Float64[],
        q::Vector{SparseMatrixCSC{Float64, Int64}}=[spzeros(Float64, 0, 0) for _ in 1:3],
        C_uc::AbstractMatrix=spzeros(Float64, 0, 0),
        control::String="force",
        sim_time::Float64=0.0,
        t_steps::Float64=0.0,
        viscosity_type::String="constant",
        cParam::Vector{Float64}=Float64[]
    )
        q_d = Dict{Symbol, Matrix{Float64}}(:top => q[1], :bottom => q[2], :border => q[3])
        new(model, β, q_d, C_uc, control, sim_time, t_steps, viscosity_type, cParam)
    end
end