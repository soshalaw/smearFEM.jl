# Point-cloud distance metrics (Hausdorff, Chamfer, closest-point RMSE).
# KDTree-based variants use NearestNeighbors; brute-force variants are kept for
# reference/benchmarking. `mean` (Statistics) and `KDTree`/`nn` (NearestNeighbors)
# are brought into the module scope by other included files.

"""
    _as_points_by_rows(pts) -> AbstractArray

Normalize a point cloud to a `(n_points, n_dims)` layout. Upstream writers are
inconsistent about whether points are stored by row or by column (e.g.
`surface_nodes` vs `3D_points_*`), so orient defensively: the coordinate axis
(2 or 3) is always the smaller of the two dimensions for a real point cloud.

# Arguments
- `pts::AbstractArray`: point cloud, stored as either `(n_points, n_dims)` or
  `(n_dims, n_points)`.

# Returns
- `AbstractArray`: `pts`, transposed if necessary, oriented as `(n_points, n_dims)`.
"""
function _as_points_by_rows(pts::AbstractArray)
    return size(pts, 1) < size(pts, 2) ? permutedims(pts) : pts
end

"""
    _nn_min_dists(query_pts, ref_pts) -> Vector{Float64}

Nearest-neighbor distance from each row of `query_pts` to the closest row of
`ref_pts`, via a `KDTree` instead of brute-force O(N·M) broadcasting.

# Arguments
- `query_pts::AbstractArray`, `ref_pts::AbstractArray`: point clouds (row- or
  column-major; normalized internally via `_as_points_by_rows`).

# Returns
- `dists::Vector{Float64}`: nearest-neighbor distance for each point in `query_pts`.
"""
function _nn_min_dists(query_pts::AbstractArray, ref_pts::AbstractArray)
    tree = KDTree(permutedims(_as_points_by_rows(ref_pts)))
    _, dists = nn(tree, permutedims(_as_points_by_rows(query_pts)))
    return dists
end

"""
    hausdorff_distance_kdtree(pred_pts, gt_pts) -> Float64

KDTree-accelerated Hausdorff distance: the larger of the two directional
max nearest-neighbor distances between `pred_pts` and `gt_pts`.

# Arguments
- `pred_pts::AbstractArray`, `gt_pts::AbstractArray`: point clouds, one row per point.

# Returns
- `Float64`: the Hausdorff distance between the two point clouds.
"""
function hausdorff_distance_kdtree(pred_pts::AbstractArray, gt_pts::AbstractArray)
    d_pred_to_gt = _nn_min_dists(pred_pts, gt_pts)
    d_gt_to_pred = _nn_min_dists(gt_pts, pred_pts)
    return max(maximum(d_pred_to_gt), maximum(d_gt_to_pred))
end

"""
    chamfer_distance_kdtree(pred_pts, gt_pts) -> Float64

KDTree-accelerated Chamfer distance: the sum of the two directional mean
nearest-neighbor distances between `pred_pts` and `gt_pts`.

# Arguments
- `pred_pts::AbstractArray`, `gt_pts::AbstractArray`: point clouds, one row per point.

# Returns
- `Float64`: the Chamfer distance between the two point clouds.
"""
function chamfer_distance_kdtree(pred_pts::AbstractArray, gt_pts::AbstractArray)
    d_pred_to_gt = _nn_min_dists(pred_pts, gt_pts)
    d_gt_to_pred = _nn_min_dists(gt_pts, pred_pts)
    return mean(d_pred_to_gt) + mean(d_gt_to_pred)
end

"""
    closest_point_distance_kdtree(pred_pts, gt_pts) -> Float64

RMSE of the `pred_pts` → `gt_pts` nearest-neighbor distances (one-directional,
unlike the symmetric `chamfer_distance_kdtree`/`hausdorff_distance_kdtree`).

# Arguments
- `pred_pts::AbstractArray`, `gt_pts::AbstractArray`: point clouds, one row per point.

# Returns
- `Float64`: RMSE of the `pred_pts` → `gt_pts` nearest-neighbor distances.
"""
function closest_point_distance_kdtree(pred_pts::AbstractArray, gt_pts::AbstractArray)
    d_pred_to_gt = _nn_min_dists(pred_pts, gt_pts)
    return sqrt(mean(d_pred_to_gt.^2))
end

"""
    hausdorff_distance(pred_pts, gt_pts) -> Float64

Brute-force (O(N·M)) Hausdorff distance between two point clouds: the largest
of the two directional sup-min distances between `pred_pts` and `gt_pts`.
Superseded by `hausdorff_distance_kdtree`; kept for reference/benchmarking.

# Arguments
- `pred_pts::AbstractArray`, `gt_pts::AbstractArray`: point clouds, one row per point.

# Returns
- `hausdorff_dist::Float64`: the Hausdorff distance between the two point clouds.
"""
function hausdorff_distance(pred_pts::AbstractArray, gt_pts::AbstractArray)
    sup_pred = maximum([minimum(sqrt.(sum((gt_pts .- permutedims(pred_pt)).^2, dims=2))) for pred_pt in eachrow(pred_pts)])
    sup_gt = maximum([minimum(sqrt.(sum((pred_pts .- permutedims(gt_pt)).^2, dims=2))) for gt_pt in eachrow(gt_pts)])
    hausdorff_dist = max(sup_pred, sup_gt)
    return hausdorff_dist
end

"""
    chamfer_distance(pred_pts, gt_pts) -> Float64

Brute-force (O(N·M)) Chamfer distance between two point clouds: the sum of
the two directional mean-min distances between `pred_pts` and `gt_pts`.
Superseded by `chamfer_distance_kdtree`; kept for reference/benchmarking.

# Arguments
- `pred_pts::AbstractArray`, `gt_pts::AbstractArray`: point clouds, one row per point.

# Returns
- `chamfer_dist::Float64`: the Chamfer distance between the two point clouds.
"""
function chamfer_distance(pred_pts::AbstractArray, gt_pts::AbstractArray)
    mean_pred = mean([minimum(sqrt.(sum((gt_pts .- permutedims(pred_pt)).^2, dims=2))) for pred_pt in eachrow(pred_pts)])
    mean_gt = mean([minimum(sqrt.(sum((pred_pts .- permutedims(gt_pt)).^2, dims=2))) for gt_pt in eachrow(gt_pts)])
    chamfer_dist = mean_pred + mean_gt
    return chamfer_dist
end

"""
    compare_pt_clouds(pred_pts, gt_pts) -> (hausdorff_distances, chamfer_distances, closest_pt_distances)

Compute, frame by frame, the KD-tree-based Hausdorff, Chamfer, and
(one-directional) closest-point RMSE distances between matching pairs of
predicted and ground-truth point clouds.

# Arguments
- `pred_pts::AbstractArray`, `gt_pts::AbstractArray`: matching iterables of
  point-cloud arrays, one pair per frame.

# Returns
- `hausdorff_distances::Vector{Float64}`, `chamfer_distances::Vector{Float64}`,
  `closest_pt_distances::Vector{Float64}`: one distance per frame.
"""
function compare_pt_clouds(pred_pts::AbstractArray, gt_pts::AbstractArray)
    hausdorff_distances = Float64[]
    chamfer_distances = Float64[]
    closest_pt_distances = Float64[]
    for (sim_pts, gt_pts) in zip(pred_pts, gt_pts)
        hausdorff_dist = hausdorff_distance_kdtree(sim_pts, gt_pts)
        chamfer_dist = chamfer_distance_kdtree(sim_pts, gt_pts)
        closest_pt_dist = closest_point_distance_kdtree(sim_pts, gt_pts)
        push!(hausdorff_distances, hausdorff_dist)
        push!(chamfer_distances, chamfer_dist)
        push!(closest_pt_distances, closest_pt_dist)
    end
    return  hausdorff_distances, chamfer_distances, closest_pt_distances
end
