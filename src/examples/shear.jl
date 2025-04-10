using Plots

F = 3.0
R_0 = 1.2
H_0 = 12.4
η_0 = 20.0
n = 0.5
K = 2.0

function get_η(t::Float64, F::Float64, R_0::Float64, H_0::Float64, η_0::Float64, n::Float64, K::Float64)

    H(t) = H_0*(1+8*H_0^2*F*t/(3*π*η_0*R_0^4))^(-1/4) # height with time
    
    R(t) = R_0*(1+8*H_0^2*F*t/(3*π*η_0*R_0^4))^(1/8) # radius with time

    H_dot(t) = 8/3*(-2*F*H(t)^3/(8*π*η_0*R(t)^4)) # rate of change of height

    γ_dot(t) = H_dot(t)/H(t) # shear rate

    η(t) = K*(abs(γ_dot(t)))^(n-1) # shear viscosity

    return η(t)
end

endTime = 10
steps = 10
tSteps = endTime/steps
time = collect(range(start=tSteps,stop=endTime,step=tSteps))

η = get_η.(time, F, R_0, H_0, η_0, n, K)
