using LinearAlgebra

mutable struct Newton
    η::Float64
    t::Int64
end

function Newton(; η::Float64 = 0.01)
    @assert η > 0 "η must be positive"

    Newton(η, 0)
end

params(opt::Newton) = "η = $(opt.η)"

function update(opt::Newton, x, ∇f, ∇²f)
    opt.t += 1
    x -= opt.η * ∇f ./ ∇²f

    return x
end