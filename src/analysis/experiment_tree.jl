# Experiment result-tree types and walker: the on-disk convention is
# `elem_size/simtime_.../noise_.../dt_.../view_...`, grouped for post-processing.

"""
    ExpLeaf

One `dt_.../view_...` result directory (`exp_path`) within an `ExpGroup`.
"""
struct ExpLeaf
    step_folder::String
    view_folder::String
    exp_path::String
end

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
    collect_experiment_groups(filepath) -> Vector{ExpGroup}

Walk the `elem_size/simtime/noise/dt/view` result folder tree under `filepath`
once and group leaves by `(elem_size_folder, sim_time_folder, noise_folder)`.

# Arguments
- `filepath::String`: root of the result tree to walk.

# Returns
- `Vector{ExpGroup}`: one group per `(elem_size_folder, sim_time_folder, noise_folder)` combination.
"""
function collect_experiment_groups(filepath)
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

                    for view_folder in readdir(dt_path)
                        startswith(view_folder, "view_") || continue
                        push!(leaves, ExpLeaf(step_folder, view_folder, joinpath(dt_path, view_folder)))
                    end
                end
                push!(groups, ExpGroup(elem_size_folder, sim_time_folder, noise_folder, step_path, leaves))
            end
        end
    end
    return groups
end
