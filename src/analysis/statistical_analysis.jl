# Cross-replicate statistics for post-processing: given the same normalized quantity measured
# on several nominally identical experiments, how well do they agree?
#
# Only quantities normalized per run belong here. A raw material parameter is a property of
# the individual specimen and legitimately differs between runs, so pooling it across
# replicates describes the specimens rather than the method. Normalizing each run against its
# own reference (h/h_m, MOSD_est/MOSD_gt, …) removes the specimen and leaves the agreement
# between model and measurement, which is the thing replicates can actually speak to.
#
# Sample sizes here are small (n = 5 is typical), so intervals use the t quantile and every
# mean carries its median and range alongside.

using Statistics
using Distributions
using Printf
using DelimitedFiles

"""
    replicate_stats(v) -> NamedTuple

Summary statistics for one sample of replicate values.

At n = 5 the interval must use the t quantile (t₄ = 2.776 against the normal 1.96 — a 42%
wider interval), and the median/range is carried alongside because five points cannot support
a normality assumption. `rc` is the Bland-Altman repeatability coefficient `2.77·SD`: the
spread within which two replicates agree 95% of the time.

# Arguments
- `v::AbstractVector{<:Real}`: one value per replicate.

# Returns
- `NamedTuple` with `n`, `mean`, `sd`, `ci` (half-width), `lo`, `hi`, `cv` (%), `median`,
  `min`, `max`, `rc`.
"""
function replicate_stats(v::AbstractVector{<:Real})
    n = length(v); m = mean(v); sd = n > 1 ? std(v) : 0.0
    half = n > 1 ? quantile(TDist(n - 1), 0.975) * sd / sqrt(n) : 0.0
    return (n=n, mean=m, sd=sd, ci=half, lo=m - half, hi=m + half,
            cv=100 * sd / abs(m), median=median(v), min=minimum(v), max=maximum(v), rc=2.77 * sd)
end

"""
    bias_vs_reference(v, ref) -> (t, p)

Two-sided one-sample t-test of `mean(v) == ref`.

Applied to the per-run mean of a normalized quantity against its ideal value, it separates a
systematic model offset — which more replicates will not average away — from run-to-run
scatter. Only meaningful for a quantity that can fall on either side of `ref`; see
[`normalized_replicate_stats`](@ref).
"""
function bias_vs_reference(v::AbstractVector{<:Real}, ref::Real)
    n = length(v); sd = std(v)
    sd == 0 && return (Inf, 0.0)
    t = (mean(v) - ref) / (sd / sqrt(n))
    return (t, 2 * ccdf(TDist(n - 1), abs(t)))
end

"""
    replicate_ci(sd, n) -> Vector

Half-width of the 95% t interval for each element of a pointwise SD, for plotting a
mean ± CI ribbon across `n` replicates.
"""
replicate_ci(sd, n::Int) = n > 1 ? quantile(TDist(n - 1), 0.975) .* sd ./ sqrt(n) : zero(sd)

"""
    normalized_replicate_stats(series, label; reference=1.0, signed=true) -> NamedTuple

Cross-replicate statistics for one normalized quantity.

`series` holds one already-normalized vector per replicate; runs of differing length are
truncated to their common span. `reference` is the quantity's ideal value (1 for a ratio, 0
for an error).

`signed` marks whether the quantity can fall on either side of `reference`. The bias test
runs only when it can: an absolute error is non-negative, so testing its mean against zero
rejects whenever any error exists and says nothing about whether that error is systematic.
The signed ratio against 1 is what answers that.

# Returns
- `NamedTuple` with the pointwise `mean`/`sd` across replicates (for a ribbon), the raw
  per-run values `per_run_mean`/`per_run_rmse`/`per_run_maxdev`, their summaries
  `mean`/`rmse`/`maxdev` (each a [`replicate_stats`](@ref)), `bias` (`nothing` when
  `signed=false`), and the worst pointwise SD with its index.
"""
function normalized_replicate_stats(series::Vector{Vector{Float64}}, label::String;
                                    reference::Float64=1.0, signed::Bool=true)
    isempty(series) && return nothing
    n = minimum(length.(series))
    M = hcat((v[1:n] for v in series)...)          # time × replicate

    per_run_mean   = vec(mean(M; dims=1))
    per_run_rmse   = vec(sqrt.(mean((M .- reference) .^ 2; dims=1)))
    per_run_maxdev = vec(maximum(abs.(M .- reference); dims=1))
    return (label=label, n_runs=size(M, 2), n_steps=n, reference=reference,
            pointwise_mean=vec(mean(M; dims=2)), pointwise_sd=vec(std(M; dims=2)),
            # The individual per-run values, not just their summaries: at these sample sizes
            # a plot has to be able to show every observation.
            per_run_mean=per_run_mean, per_run_rmse=per_run_rmse, per_run_maxdev=per_run_maxdev,
            mean=replicate_stats(per_run_mean),
            rmse=replicate_stats(per_run_rmse),
            maxdev=replicate_stats(per_run_maxdev),
            bias=signed ? bias_vs_reference(per_run_mean, reference) : nothing,
            worst_sd=maximum(vec(std(M; dims=2))), worst_sd_step=argmax(vec(std(M; dims=2))))
end

"""
    print_replicate_table(stats, title)

Print the replicate statistics for several quantities as one table.
"""
function print_replicate_table(stats::AbstractVector, title::String)
    isempty(stats) && return nothing
    println("\n", "="^96)
    println("$title — n = $(stats[1].n_runs) runs over $(stats[1].n_steps) steps")
    println("="^96)
    @printf("%-22s %-24s %-22s %8s %10s\n", "quantity", "mean ± SD", "95% CI (t)", "CV", "repeat.")
    println("-"^96)
    for st in stats
        for (what, sm) in (("run mean", st.mean), ("RMSE vs ref", st.rmse), ("max deviation", st.maxdev))
            @printf("%-22s %10.4f ± %-11.4f [%8.4f, %8.4f] %7.1f%% %9.4f\n",
                           what == "run mean" ? st.label : "  " * what,
                           sm.mean, sm.sd, sm.lo, sm.hi, sm.cv, sm.rc)
        end
        if isnothing(st.bias)
            @printf("  %-20s worst between-run SD %.4f at step %d\n", "", st.worst_sd, st.worst_sd_step)
        else
            @printf("  %-20s worst between-run SD %.4f at step %d;  bias t(%d) = %.2f, p = %.4f → %s\n",
                           "", st.worst_sd, st.worst_sd_step, st.n_runs - 1, st.bias[1], st.bias[2],
                           st.bias[2] < 0.05 ? "SYSTEMATIC" : "consistent with noise")
        end
    end
    println("="^96)
    return nothing
end

"""
    replicate_report(stats, title, csvpath) -> Bool

Write `stats` to `csvpath`.

Returns `false` without writing when there is nothing to report, so a caller can guard its
plotting on the return value instead of repeating the emptiness check.

Nothing is printed: these run inside a long post-analysis pass where a table per tree buries
the warnings that matter. Call [`print_replicate_table`](@ref) directly to see one. `title`
is kept so a caller names what it is reporting at the call site.
"""
function replicate_report(stats::AbstractVector, title::String, csvpath::String)
    isempty(stats) && return false
    write_replicate_stats(csvpath, stats)
    @info "Wrote $title to $csvpath"
    return true
end

"""
    write_replicate_stats(filepath, stats)

Write the per-run scalar summaries for several quantities to one CSV.
"""
function write_replicate_stats(filepath::String, stats::AbstractVector)
    rows = Any[["quantity" "statistic" "n" "mean" "sd" "ci_lo" "ci_hi" "cv_pct" "median" "min" "max" "repeat_coef"]]
    for st in stats, (what, sm) in (("run_mean", st.mean), ("rmse_vs_ref", st.rmse), ("max_deviation", st.maxdev))
        push!(rows, [st.label what sm.n sm.mean sm.sd sm.lo sm.hi sm.cv sm.median sm.min sm.max sm.rc])
    end
    writedlm(filepath, vcat(rows...), ',')
    return filepath
end
