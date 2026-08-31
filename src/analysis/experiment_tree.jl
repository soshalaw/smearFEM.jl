# Experiment result-tree types and walker: the on-disk convention is
# `elem_size/simtime_.../noise_.../dt_.../<opt_method>/view_...`, grouped for
# post-processing. The `<opt_method>` level (`gn`, `lm`, `gn_tikhonov`) was added on
# 2026-08-26 so runs of different optimizers do not overwrite each other; trees written
# before then have `view_...` directly under `dt_...` and are still walked correctly.

"""
    ExpLeaf

One `dt_.../<opt_method>/view_...` result directory (`exp_path`) within an `ExpGroup`.
`method_folder` is the optimizer that produced it, or `""` for a pre-2026-08-26 tree
that has no method level.
"""
struct ExpLeaf
    step_folder::String
    view_folder::String
    exp_path::String
    method_folder::String
end

# Pre-2026-08-26 trees have no method level; default it so old positional construction
# and any external caller keep working.
ExpLeaf(step_folder, view_folder, exp_path) = ExpLeaf(step_folder, view_folder, exp_path, "")

"""
    ExpGroup

All `ExpLeaf`s sharing the same `(elem_size_folder, sim_time_folder,
noise_folder)` combination, as produced by `collect_experiment_groups`.
"""
struct ExpGroup
    elem_size_folder::String
    sim_time_folder::String
    noise_folder::String
    step_path::String
    leaves::Vector{ExpLeaf}
end

"""
    collect_experiment_groups(filepath; method=nothing) -> Vector{ExpGroup}

Walk the `elem_size/simtime/noise/dt/<opt_method>/view` result folder tree under
`filepath` once and group leaves by `(elem_size_folder, sim_time_folder, noise_folder)`.

A `dt_*` directory holding `view_*` directly (no method level) is a pre-2026-08-26 tree
and yields leaves with `method_folder == ""`, so old and migrated trees can be mixed.

# Arguments
- `filepath::String`: root of the result tree to walk.

# Keyword Arguments
- `method`: keep only leaves produced by this optimizer (a `Symbol` or `String`, e.g.
  `:gn_tikhonov`); `nothing` (default) keeps every method. Use it when several methods
  have been run over the same experiment and a caller wants one of them — notably
  `_single_leaf`, which errors when a group holds more than one leaf.

# Returns
- `Vector{ExpGroup}`: one group per `(elem_size_folder, sim_time_folder, noise_folder)` combination.
"""
function collect_experiment_groups(filepath; method=nothing)
    want_method = isnothing(method) ? nothing : String(method)
    groups = ExpGroup[]
    for elem_size_folder in readdir(filepath)
        elem_size_folder == "post_analysis" && continue

        for sim_time_folder in readdir(joinpath(filepath, elem_size_folder))
            startswith(sim_time_folder, "simtime_") || continue
            noise_base = joinpath(filepath, elem_size_folder, sim_time_folder)

            for noise_folder in readdir(noise_base)
                startswith(noise_folder, "noise_") || continue
                step_path = joinpath(noise_base, noise_folder)

                leaves = ExpLeaf[]
                for step_folder in readdir(step_path)
                    startswith(step_folder, "dt_") || continue
                    dt_path = joinpath(step_path, step_folder)

                    for entry in readdir(dt_path)
                        entry_path = joinpath(dt_path, entry)
                        isdir(entry_path) || continue

                        if startswith(entry, "view_")
                            # Pre-2026-08-26 layout: view_* sits directly under dt_*.
                            isnothing(want_method) || continue
                            push!(leaves, ExpLeaf(step_folder, entry, entry_path, ""))
                            continue
                        end

                        isnothing(want_method) || entry == want_method || continue
                        for view_folder in readdir(entry_path)
                            startswith(view_folder, "view_") || continue
                            push!(leaves, ExpLeaf(step_folder, view_folder,
                                                  joinpath(entry_path, view_folder), entry))
                        end
                    end
                end
                push!(groups, ExpGroup(elem_size_folder, sim_time_folder, noise_folder, step_path, leaves))
            end
        end
    end
    return groups
end
