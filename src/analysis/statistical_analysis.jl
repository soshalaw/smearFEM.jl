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

Returns `(NaN, NaN)` for a single replicate: `std` is undefined there and `TDist(0)` has no
degrees of freedom, so the question cannot be asked. `replicate_stats` and `replicate_ci`
already degrade the same way at `n = 1`, and a method fitted on only one run is a real
configuration — it must not take the whole post-analysis down with it.
"""
function bias_vs_reference(v::AbstractVector{<:Real}, ref::Real)
    n = length(v)
    n > 1 || return (NaN, NaN)
    sd = std(v)
    sd == 0 && return (Inf, 0.0)
    t = (mean(v) - ref) / (sd / sqrt(n))
    return (t, 2 * ccdf(TDist(n - 1), abs(t)))
end

"""
    pooled_stats(series) -> NamedTuple

Grand mean and SD computed in two stages — per-replicate first, then pooled across
replicates — rather than by flattening every raw value into one sample.

`series` holds one time series per replicate, of possibly differing length — no truncation to a
common span is needed here.

**Stage 1, per replicate** (purely descriptive, no independence assumption needed — this is
just characterising *one* trajectory): `per_run_mean[i]` and `per_run_sd[i]`, the time-average
and temporal SD of replicate `i` around its own mean.

**Stage 2, across replicates** (independence matters here, since the `R` replicates — not the
individual time steps — are the actual independent unit):
- `mean` = the unweighted average of `per_run_mean` (matches [`replicate_stats`](@ref)'s
  `.mean.mean` applied to the same vector).
- within-run variance = `Σᵢ(Tᵢ-1)·per_run_sd[i]² / Σᵢ(Tᵢ-1)`, the standard pooled-variance
  formula for combining `R` independent variance estimates, each contributing its own degrees
  of freedom.
- between-run variance = the sample variance of `per_run_mean` (`R-1` denominator) — the same
  quantity `replicate_stats(per_run_mean).sd^2` reports.
- `sd = √(within-run variance + between-run variance)`.

This is deliberately *not* `std(vcat(series...))`, i.e. not the SD of every `(time, replicate)`
value flattened into one sample: that computes a single sum of squared deviations across
`n_obs - 1` "degrees of freedom" as if every time step were an independent draw, when time steps
within one replicate are autocorrelated. Computing each replicate's own SD first, then combining
those `R` independent estimates via the pooled-variance formula, never treats a within-run
timestep as contributing independent information — only the `R` replicates do. (For a balanced
design the two approaches are numerically close but not identical; this one is the
methodologically correct decomposition, not an approximation to the flattened one.)

The CI's standard error still divides by `√R`, not `√n_obs`, for the same reason: only the
replicates are independent.

# Returns
- `NamedTuple` with `n` (= number of replicates, for the CI's degrees of freedom), `n_obs`
  (total pooled point count, for reference), `mean` (grand mean), `sd` (pooled SD, see above),
  `ci` (half-width for the grand mean), `lo`, `hi`, `cv` (%), `median`, `min`, `max` (of every
  raw pooled value — informative for eyeballing range regardless of autocorrelation), `rc`
  (`2.77·sd`, as in [`replicate_stats`](@ref)).
"""
function pooled_stats(series::Vector{Vector{Float64}})
    isempty(series) && return nothing
    r = length(series)
    per_run_mean = mean.(series)
    per_run_sd   = [length(v) > 1 ? std(v) : 0.0 for v in series]
    Ti = length.(series)
    n_obs = sum(Ti)

    m = mean(per_run_mean)

    within_df  = sum(Ti .- 1)
    within_var = within_df > 0 ? sum((Ti .- 1) .* per_run_sd .^ 2) / within_df : 0.0
    between_var = r > 1 ? var(per_run_mean) : 0.0
    sd = sqrt(within_var + between_var)

    half = r > 1 ? quantile(TDist(r - 1), 0.975) * sd / sqrt(r) : 0.0

    pooled = reduce(vcat, series)
    return (n=r, n_obs=n_obs, mean=m, sd=sd, ci=half, lo=m - half, hi=m + half,
            cv=100 * sd / abs(m), median=median(pooled), min=minimum(pooled), max=maximum(pooled),
            rc=2.77 * sd)
end

"""
    pointwise_pooled_stats(pointwise_mean, pointwise_sd, n_runs) -> NamedTuple

Grand mean and SD computed by pooling the *pointwise* (cross-replicate, at each time step)
statistics over time — collapsing across replicates first, at every instant, then combining
those `T` estimates across time. This is the third of three approaches to cross-replicate
spread; see [`pooled_stats`](@ref) for the other (collapse time first, then combine across
replicates) and the module docstring for how they relate.

`pointwise_mean`/`pointwise_sd` are the same arrays [`normalized_replicate_stats`](@ref) already
computes for the ribbon: at each time step `t`, the mean and sample SD across the `n_runs`
replicates.

`sd` is the RMS average of the pointwise SD over time, `√(mean(pointwise_sd.^2))`. Because the
replicates are collapsed *first*, at each fixed `t`, anything identical across every replicate
at that instant — a pattern shared by the whole batch, e.g. a time-varying artifact of the
fitting method rather than genuine specimen-to-specimen noise — cancels out of `pointwise_sd`
exactly, before time is ever touched. [`pooled_stats`](@ref) collapses time first instead, so it
cannot make that distinction: a shared pattern inflates its SD along with genuine noise. The two
agree only when there is no such shared pattern.

The CI's standard error still divides by `√n_runs`, matching [`replicate_stats`](@ref)/
[`pooled_stats`](@ref): the replicates are the only independent unit anywhere in this analysis,
regardless of which axis the SD itself is built from.

# Returns
- `NamedTuple` with the same shape as [`replicate_stats`](@ref)/[`pooled_stats`](@ref) — `n`,
  `mean`, `sd`, `ci`, `lo`, `hi`, `cv`, `median`, `min`, `max`, `rc` — so any of the three can be
  used interchangeably wherever a summary is expected. `median`/`min`/`max` describe
  `pointwise_mean` (the ensemble-mean trajectory over time), the closest analogue here to
  `pooled_stats`'s pooled raw values.
"""
function pointwise_pooled_stats(pointwise_mean::AbstractVector{<:Real},
                                pointwise_sd::AbstractVector{<:Real}, n_runs::Int)
    isempty(pointwise_mean) && return nothing
    m = mean(pointwise_mean)
    sd = sqrt(mean(pointwise_sd .^ 2))
    half = n_runs > 1 ? quantile(TDist(n_runs - 1), 0.975) * sd / sqrt(n_runs) : 0.0
    return (n=n_runs, mean=m, sd=sd, ci=half, lo=m - half, hi=m + half,
            cv=100 * sd / abs(m), median=median(pointwise_mean), min=minimum(pointwise_mean),
            max=maximum(pointwise_mean), rc=2.77 * sd)
end

"""
    approach_stats(st, approach) -> NamedTuple

Pick the [`replicate_stats`](@ref)-shaped cross-replicate summary for one of the three spread
approaches out of a [`normalized_replicate_stats`](@ref) result `st`:

- `:A` — `st.mean`: collapse each replicate to a scalar first, then analyze across replicates.
- `:B` — `st.pointwise_pooled`: collapse across replicates first (at each time step), then pool
  over time. See [`pointwise_pooled_stats`](@ref).
- `:C` — `st.pooled`: collapse each replicate over time first, then pool across replicates. See
  [`pooled_stats`](@ref).
"""
function approach_stats(st, approach::Symbol)
    approach === :A && return st.mean
    approach === :B && return st.pointwise_pooled
    approach === :C && return st.pooled
    error("approach_stats: unknown approach :$approach — use :A, :B or :C")
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
The signed ratio against 1 is what answers that. It is also skipped for a single replicate,
where `bias` comes back `nothing` — one run cannot separate an offset from scatter.

# Returns
- `NamedTuple` with the pointwise `mean`/`sd` across replicates (for a ribbon), the raw
  per-run values `per_run_mean`/`per_run_rmse`/`per_run_maxdev`, their summaries
  `mean`/`rmse`/`maxdev` (each a [`replicate_stats`](@ref)), `pooled` (approach C — a
  [`pooled_stats`](@ref) computed by collapsing each replicate over time first, then pooling
  across replicates), `pointwise_pooled` (approach B — a [`pointwise_pooled_stats`](@ref)
  computed by collapsing across replicates first, then pooling over time), `bias` (`nothing`
  when `signed=false`), and the worst pointwise SD with its index. See [`approach_stats`](@ref)
  to pick one of `mean`/`pointwise_pooled`/`pooled` by `:A`/`:B`/`:C`.
"""
function normalized_replicate_stats(series::Vector{Vector{Float64}}, label::String;
                                    reference::Float64=1.0, signed::Bool=true)
    isempty(series) && return nothing
    n = minimum(length.(series))
    M = hcat((v[1:n] for v in series)...)          # time × replicate

    per_run_mean   = vec(mean(M; dims=1))
    per_run_rmse   = vec(sqrt.(mean((M .- reference) .^ 2; dims=1)))
    per_run_maxdev = vec(maximum(abs.(M .- reference); dims=1))
    _pointwise_mean = vec(mean(M; dims=2))
    _pointwise_sd   = vec(std(M; dims=2))
    return (label=label, n_runs=size(M, 2), n_steps=n, reference=reference,
            pointwise_mean=_pointwise_mean, pointwise_sd=_pointwise_sd,
            # The individual per-run values, not just their summaries: at these sample sizes
            # a plot has to be able to show every observation.
            per_run_mean=per_run_mean, per_run_rmse=per_run_rmse, per_run_maxdev=per_run_maxdev,
            mean=replicate_stats(per_run_mean),
            rmse=replicate_stats(per_run_rmse),
            maxdev=replicate_stats(per_run_maxdev),
            # Pooled over the untruncated series (not `M`'s common-length columns): every raw
            # observation gets to speak, unlike `per_run_mean`, which drops the temporal spread
            # the moment it averages each run down to one number.
            pooled=pooled_stats(series),
            pointwise_pooled=pointwise_pooled_stats(_pointwise_mean, _pointwise_sd, size(M, 2)),
            # `nothing`, not a NaN pair: `print_replicate_table` then omits the verdict
            # instead of reading NaN < 0.05 as "consistent with noise".
            bias=(signed && size(M, 2) > 1) ? bias_vs_reference(per_run_mean, reference) : nothing,
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
        for (what, sm) in (("run mean", st.mean), ("pointwise pooled", st.pointwise_pooled),
                           ("pooled mean", st.pooled),
                           ("RMSE vs ref", st.rmse), ("max deviation", st.maxdev))
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
    for st in stats, (what, sm) in (("run_mean", st.mean), ("pointwise_pooled", st.pointwise_pooled),
                                    ("pooled_mean", st.pooled),
                                    ("rmse_vs_ref", st.rmse), ("max_deviation", st.maxdev))
        push!(rows, [st.label what sm.n sm.mean sm.sd sm.lo sm.hi sm.cv sm.median sm.min sm.max sm.rc])
    end
    writedlm(filepath, vcat(rows...), ',')
    return filepath
end

"""
    write_approach_stats(filepath, stats, approach)

Write mean/SD/95% CI for one cross-replicate spread approach (`:A`/`:B`/`:C`, see
[`approach_stats`](@ref)) across several quantities to one CSV — one row per quantity, unlike
[`write_replicate_stats`](@ref) which writes every approach as separate rows of one combined
file.
"""
function write_approach_stats(filepath::String, stats::AbstractVector, approach::Symbol)
    rows = Any[["quantity" "n" "mean" "sd" "ci_lo" "ci_hi" "cv_pct" "median" "min" "max" "repeat_coef"]]
    for st in stats
        sm = approach_stats(st, approach)
        push!(rows, [st.label sm.n sm.mean sm.sd sm.lo sm.hi sm.cv sm.median sm.min sm.max sm.rc])
    end
    writedlm(filepath, vcat(rows...), ',')
    return filepath
end
