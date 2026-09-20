# Moving-window parameter estimation — Section 3 of "Vision-Informed Short-Horizon Prediction
# for the Squeeze Flow of Soft Materials".
#
# The observation interval [0, T] is partitioned into N non-overlapping windows
# w_k = [t_k, t_{k+1}] (Eq. 20). In each window the parameters θ = [η, β] are calibrated to the
# observed borders (Eq. 25), and the result is propagated forward (step iii.c):
#
#     θ̂_{k+1} ← θ_k*        and        M_{k+1,0} ← M_{k,N_k}
#
# `estimate_window` performs the calibration and the propagation; `predict_window` simulates a
# window forward under a *fixed* estimate from the preceding one — the k → k+1 prediction the
# paper reports against. Both belong here rather than in `experiments/`: they are the method,
# and they depend only on library types and functions.

"""
    set_time_window(time_step_len, data; method="linear", window_size=10.0, global_window_sz=30.0) -> (time_windows, windows, data_ranges, t_windows)

Split `data` (indexed along its first dimension by time) into successive time
windows whose cumulative end time grows according to `method`
(`"linear"`: `window_size*iter`; `"quadratic"`: `window_size*iter^2`;
`"exponential"`: `window_size*exp(3*(iter-1))`), stopping once the window end
reaches the end of `data` or `global_window_sz`.

# Arguments
- `time_step_len::Float64`: time step length of `data`.
- `data::AbstractArray`: array indexed along its first dimension by time.

# Keyword Arguments
- `method::String`: window growth schedule, `"linear"`, `"quadratic"`, or
  `"exponential"` (default: `"linear"`).
- `window_size::Float64`: base window size (default: `10.0`).
- `global_window_sz::Float64`: overall time budget the windows are truncated
  to (default: `30.0`).

# Returns
- `time_windows::Vector{Float64}`: duration of each window.
- `windows::Vector{AbstractArray}`: the slice of `data` in each window.
- `data_ranges::Vector{AbstractArray}`: the index range (relative to `data`)
  of each window.
- `t_windows::Vector{Float64}`: cumulative end time of each window.
"""
function set_time_window(time_step_len::Float64, data::AbstractArray; method::String="linear", window_size::Float64=10.0, global_window_sz::Float64=50.0)
    windows::Vector{AbstractArray} = Vector{AbstractArray}()
    time_windows::Vector{Float64} = Vector{Float64}()
    data_ranges::Vector{AbstractArray} = Vector{AbstractArray}()
    t_windows::Vector{Float64} = Vector{Float64}()

    function get_t_window(window_size::Float64, step_len::Float64, iter::Int, method::String)::Float64
        t_window = 0.0
        if method == "linear"
            t_window = round(window_size*iter, digits=1)
        elseif method == "quadratic"
            t_window = round(window_size*iter^2, digits=1)
        elseif method == "exponential"
            if iter == 1
                t_window = round(window_size*exp(0.5*(iter-1)), digits=1)
            else 
                δt = window_size*round(exp(3*(iter-1)), digits=1)
                if δt < window_size
                    @warn "Computed time window increment $δt is less than minimum Window $window_size; using minimum."
                    δt = window_size
                end
                t_window = round(window_size*exp(3*(iter-1)), digits=1)
            end
        end
        return t_window
    end

    iter::Int = 1
    start_point::Int = 1
    t_window_prev::Float64 = 0.0
    t_window_end::Float64 = get_t_window(window_size, time_step_len, iter, method)
    end_point::Int = round(Int,t_window_end*time_step_len)+1

    END_FLAG::Bool = false
    while true
            # If the computed end_point reaches or exceeds the data length, adjust.
            # Keep the original requested value for clearer messages.
            if end_point >= size(data, 1) || t_window_end >= global_window_sz
                if end_point == size(data, 1)
                    @info "Reached end of data at point $end_point."
                elseif end_point > size(data, 1)
                    requested_end = end_point
                    # adjust the time window end to the last available sample
                    t_window_end = round((size(data, 1)-1)/time_step_len, digits=1)
                    end_point = size(data, 1)
                    @warn "Requested end point $requested_end exceeds data size $(size(data,1)); adjusting to end of data (end_point=$end_point). Adjusted end time to $t_window_end seconds."
                elseif t_window_end >= global_window_sz
                    requested_end = t_window_end
                    # adjust the end point to match the global window size
                    t_window_end = global_window_sz
                    end_point = round(Int,t_window_end*time_step_len)+1
                    @warn "Requested time window end $requested_end seconds exceeds global window size $global_window_sz seconds; adjusting to global window size (end_point=$end_point)."
                end
                END_FLAG = true
            end

        data_range = start_point:end_point
        data_range_ = start_point:(end_point-1)

        @debug "Data frame : $data_range"
        @debug "time windows from : $t_window_prev to $t_window_end"
        @debug "time Window : $(t_window_end - t_window_prev) seconds"
        @debug "data length : $(size(data[data_range], 1))"
        @debug "data range size : $(length(data_range_))"
        @debug "----------"
        t_window_size = round(t_window_end - t_window_prev, digits=1)
        push!(time_windows, t_window_size)
        push!(windows, data[data_range])
        push!(data_ranges, data_range_)
        push!(t_windows, t_window_end)

        if END_FLAG == true
            break
        end
        iter = iter + 1
        start_point = end_point
        t_window_prev = t_window_end
        t_window_end = get_t_window(window_size, time_step_len, iter, method)
        end_point = round(Int,t_window_end*time_step_len)+1
    end
    return time_windows, windows, data_ranges, t_windows
end

"""
    estimate_window(model, scene, conditions, obs, θ; outliers, method, kwargs...)

Calibrate `θ = [η, β]` to one window's border observations — Eq. (25),
`θ_k* = argmin_θ d_k(θ)` — then propagate the result forward.

On return `θ` holds the converged estimate (so it seeds the next window, `θ̂_{k+1} ← θ_k*`) and
`model`'s mesh state has been committed via `update_model!`, so the next window starts from
`M_{k+1,0} ← M_{k,N_k}`. Both mutations are the paper's step (iii.c); callers relying on the
pre-fit state must copy it first.

# Arguments
- `model::Stokes`: model for this window; mutated in place.
- `scene::SqueezeFlow`: scenario supplying the control history and this window's timing.
- `conditions::Conditions`: camera model and output flags.
- `obs`: observed border points for this window, `{B_{k,j}}`.
- `θ::Vector{Float64}`: initial guess, overwritten with the converged estimate.
- `outliers::Vector{Int}`: frame indices excluded from the cost.
- `method::Symbol`: `:gn`, `:lm` or `:gn_tikhonov`.
- `kwargs...`: forwarded to `fit_model`.

# Returns
- `stats::Dict`: as `fit_model` returns — `"η"`, `"β"`, `"cost_list"`, `"iterList"`, and the
  method's own extra keys.
"""
function estimate_window(model::Stokes, scene::SqueezeFlow, conditions::Conditions,
                         obs, θ::Vector{Float64};
                         outliers::Vector{Int}=Int[], method::Symbol=:gn, kwargs...)
    stats = fit_model(model, scene, conditions, obs, θ; outliers=outliers, method=method, kwargs...)
    θ[1] = stats["η"]
    θ[2] = stats["β"]
    update_model!(model)
    return stats
end

"""
    predict_window(model, scene, conditions, θ_prev)

Simulate one window forward under a parameter estimate held fixed — the `k → k+1` prediction:
parameters identified in window `k` used to predict window `k+1`.

The model is reset to its committed reference first, so the prediction starts from the same
mesh state the estimation for this window started from, not from wherever a previous
prediction left it.

# Arguments
- `model::Stokes`: model used for prediction; mutated in place.
- `scene::SqueezeFlow`: scenario; its `sim_time` and `cParam` must already describe this window.
- `conditions::Conditions`: camera model and output flags.
- `θ_prev`: the preceding window's converged `(η, β)` — a 2-tuple or a vector.

# Returns
- the full `simulate` tuple, so callers take only what they need.
"""
function predict_window(model::Stokes, scene::SqueezeFlow, conditions::Conditions,
                        θ_prev::Union{AbstractVector{<:Real},Tuple{Real,Real}})
    reset_model!(model)
    model.η = [float(θ_prev[1])]
    scene.β = [float(θ_prev[2])]
    return simulate(model, scene, conditions)
end
