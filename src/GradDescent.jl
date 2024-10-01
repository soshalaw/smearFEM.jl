using LinearAlgebra

mutable struct GradDescent
    η::Float64
    t::Int64
end

function GradDescent(; η::Float64 = 0.01)
    @assert η > 0 "η must be positive"
    GradDescent(η, 0)
end

params(opt::GradDescent) = "η = $(opt.η)"

function update(opt::GradDescent, x, ∇f)
    opt.t += 1
    x -= opt.η * ∇f

    return x
end