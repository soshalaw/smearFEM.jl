# Archive: dead code removed from `features/volume_update`

This branch is `features/volume_update` at e01e4c0 (2026-11-04), kept so the code below can be
recovered. It was removed from `features/volume_update` because nothing on that branch calls it
(checked by grep over `src/`, `test/`, `experiments/`, `docs/` and the README).

- `src/fem/PostProcess.jl` — the Stokes-era camera pipeline: border extraction
  (`extract_borders`, convex hull), projection and its gradient (`project_to`, `back_project`,
  `∇π`, `project_to_camera_frame`), contour fitting and resampling (`fit_curve`,
  `upsample_contour`, `filter_points`, `add_noise`), outlier frames
  (`detect_outlier_observations`, `get_pose`), `get_height`, and the NURBS/IGA helpers
  (`eval_on_cylinder`, `get_lagrange_proj`, `get_nurbs_2_lagrange_proj`, `get_lagrange_pts`).
  Its users (the Stokes solver and the fits) were deleted earlier. Two functions stayed:
  `rearrange` (moved to `src/io/io.jl`, which uses it) and the vector method of
  `_set_cam_frame` (moved to `src/io/plotting.jl`, which uses it). Its 4×4-matrix method (a
  measured camera pose) is archived: `animate_fields` only takes a vector `cam_pose`.
- The dependencies only that file used: `ConvexHulls2d`, `DataInterpolations`, `ElasticArrays`
  (removed from `Project.toml`; Aqua flags a dependency the package never loads).
- `src/fem/types.jl` — `EnvConditions` and `Conditions`, the output/camera settings of a Stokes
  run. No function takes them any more.
- `Model{T}` (`src/fem/models.jl`) — a mutable box around a model; never constructed.
- `hausdorff_distance`, `chamfer_distance` (`src/analysis/pointcloud_metrics.jl`) — brute-force
  O(N·M) versions of the KDTree metrics, kept for reference; no caller or test.

To recover a file: `git checkout archive/dead_code -- <path>`.
