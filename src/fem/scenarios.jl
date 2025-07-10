abstract type AbstractScenario end

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
        β::Vector{Float64},
        q::Vector{SparseMatrixCSC{Float64, Int64}},
        C_uc::AbstractMatrix,
        control::String,
        sim_time::Float64,
        t_steps::Float64,
        viscosity_type::String,
        cParam::Vector{Float64}
    )   
        q_d = Dict{Symbol, Matrix{Float64}}(:top => q[1], :bottom => q[2], :border => q[3])
        new(model, β, q_d, C_uc, control, sim_time, t_steps, viscosity_type, cParam)
    end
end