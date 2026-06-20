using LinearAlgebra
using DataInterpolations
using ElasticArrays
using ConvexHulls2d
import ConvexHulls2d as ch
using Distributions
using Parameters
using Statistics

"""
    extract_borders(NodeList, camera_matrix, obj_pose; BorderNodesList, GRAD, dqdθ)

Project 3D mesh nodes to 2D, compute the convex-hull border, and optionally propagate
parameter gradients through the projection.

# Arguments
- `NodeList`: `(3, nNodes)` matrix of 3D node coordinates.
- `camera_matrix`: `(3, 3)` camera intrinsic matrix.
- `obj_pose`: `(4, 4)` object-to-camera rigid transformation.

# Keyword Arguments
- `BorderNodesList`: Subset of node indices to project (default: all nodes).
- `GRAD`: Propagate gradients when `true` (default: `false`).
- `dqdθ`: `(3, nNodes, nParams)` gradient tensor (required when `GRAD=true`).
# Returns (GRAD=false)
`border_pts_sorted, surface_pts_2d`

# Returns (GRAD=true)
`border_pts_sorted, ∇border_pts_sorted, surface_pts_2d, ∇surface_pts_2d`
"""
function extract_borders(NodeList::AbstractMatrix{Float64}, camera_matrix::AbstractMatrix{Float64}, obj_pose::Vector{Float64}, height::Float64, rot_angle::Float64=[0.0]; BorderNodesList::Vector{Int64}=zeros(Int64,0), GRAD::Bool=false, dqdθ::AbstractArray{Float64}=zeros(2,2,2))

    if length(BorderNodesList) != 0
        surface_nodes = NodeList[:,BorderNodesList]
    else
        surface_nodes = NodeList
    end

    if GRAD
        ∇surface_nodes = length(BorderNodesList) != 0 ? dqdθ[:,BorderNodesList,:] : dqdθ
        surface_pts_2d, ∇surface_pts_2d = back_project(surface_nodes, camera_matrix, obj_pose, ∇surface_nodes, height, rot_angle)
    else
        surface_pts_2d = back_project(surface_nodes, camera_matrix, obj_pose, height, rot_angle)
    end

    p = [surface_pts_2d[:,i] for i in axes(surface_pts_2d, 2)]

    hull = ch.ConvexHull(p)
    hull_vertex_ids = ch.indices(hull)

    border_pt_ids = get_convex_hull_boundary_points(surface_pts_2d, hull_vertex_ids)
    border_pts = surface_pts_2d[:,border_pt_ids]
    border_pts_sorted, border_pts_sorted_ids = sort_points(border_pts)


    if GRAD
        ∇border_pts = ∇surface_pts_2d[:,border_pt_ids,:]
        ∇border_pts_sorted = ∇border_pts[:,border_pts_sorted_ids,:]
    end

    if GRAD
        return border_pts_sorted, ∇border_pts_sorted, surface_pts_2d, ∇surface_pts_2d
    else
        return border_pts_sorted, surface_pts_2d
    end
end

function _get_2D_data(NodeList_cached::Matrix{Float64}, camera_matrix_cached::AbstractMatrix{Float64}, obj_pose_cached::Vector{Float64}, h_cached::Float64; BorderNodesList::Vector{Int64}=zeros(Int64,0), GRAD::Bool=false, dqdθ::AbstractArray{Float64}=zeros(2,2,2), angles::Vector{Float64}=[0.0])
    BorderPts2D   = Vector{Matrix{Float64}}()
    SurfacePts2D  = Vector{Matrix{Float64}}()
    dudθ       = Vector{Array{Float64,3}}()
    ∇SurfacePts2D = Vector{Array{Float64,3}}()
    obs_border_pts = Vector{Matrix{Float64}}()

    for rot_angle in angles
        if GRAD
            BorderPts2D_θ, dudθ_θ, SurfacePts2D_θ, ∇SurfacePts2D_θ = extract_borders(NodeList_cached, camera_matrix_cached, obj_pose_cached, h_cached, rot_angle; BorderNodesList=BorderNodesList, GRAD=true, dqdθ=dqdθ)
            pi_θ, qi_θ = fit_curve(border=BorderPts2D_θ)
            push!(BorderPts2D,   BorderPts2D_θ)
            push!(dudθ,          dudθ_θ)
            push!(SurfacePts2D,  SurfacePts2D_θ)
            push!(∇SurfacePts2D, ∇SurfacePts2D_θ)
            push!(obs_border_pts, vcat(pi_θ', qi_θ'))
        else
            BorderPts2D_i, SurfacePts2D_i = extract_borders(NodeList_cached, camera_matrix_cached, obj_pose_cached, h_cached, rot_angle; BorderNodesList=BorderNodesList)
            pi_i, qi_i = fit_curve(border=BorderPts2D_i)
            push!(BorderPts2D,  BorderPts2D_i)
            push!(SurfacePts2D, SurfacePts2D_i)
            push!(obs_border_pts, vcat(pi_i', qi_i'))
        end
    end

    if GRAD
        return BorderPts2D, dudθ, SurfacePts2D, ∇SurfacePts2D, obs_border_pts
    else
        return BorderPts2D, SurfacePts2D, obs_border_pts
    end
end


"""
    get_convex_hull_boundary_points(surface_pts_2d, hull_vertex_ids; tol) -> Vector{Int}
-
Return sorted indices of all 2D points that lie on the convex hull boundary, including
points on straight hull edges (not just the vertex extremes).

# Arguments
- `surface_pts_2d`: `(2, nPoints)` matrix of projected 2D coordinates.
- `hull_vertex_ids`: Vertex indices of the convex hull (from `ConvexHulls2d`).

# Keyword Arguments
- `tol`: Perpendicular-distance tolerance for edge membership (default: `1e-6`).
"""
function get_convex_hull_boundary_points(surface_pts_2d::AbstractMatrix{Float64}, hull_vertex_ids::Vector{Int64}; tol::Float64=5.0)
    n_vertices = length(hull_vertex_ids)
    all_boundary_ids = Set(hull_vertex_ids)
    
    # Check each edge of the hull
    for i in 1:n_vertices
        v1_id = hull_vertex_ids[i]
        v2_id = hull_vertex_ids[mod(i, n_vertices) + 1]
        
        v1 = surface_pts_2d[:, v1_id]
        v2 = surface_pts_2d[:, v2_id]
        
        edge_vec = v2 - v1
        edge_len = norm(edge_vec)
        
        if edge_len < tol
            continue
        end
        
        edge_dir = edge_vec / edge_len
        
        # Check all points to see if they lie on this edge
        for j in 1:size(surface_pts_2d, 2)
            if j == v1_id || j == v2_id
                continue
            end
            
            pt = surface_pts_2d[:, j]
            
            # Vector from v1 to point
            to_pt = pt - v1
            
            # Project onto edge
            proj_len = dot(to_pt, edge_dir)
            
            # Check if projection is within edge bounds and point is on the line
            if proj_len >= -tol && proj_len <= edge_len + tol
                # Point perpendicular distance to edge
                proj_pt = v1 + proj_len * edge_dir
                perp_dist = norm(pt - proj_pt)
                
                if perp_dist < tol
                    push!(all_boundary_ids, j)
                end
            end
        end
    end
    
    return sort(collect(all_boundary_ids))
end

"""
    sort_points(Data) -> (sorted_pts, sorted_ids)

Sort a 2D point cloud for contour traversal: left-half ascending in y, right-half
descending in y, producing a closed-loop ordering.

# Arguments
- `Data::Matrix{Float64}`: 2×N matrix of 2D point coordinates.

# Returns
- `sorted_pts::Matrix{Float64}`: 2×N reordered points.
- `sorted_ids::Vector{Int}`: original column indices in the new order.
"""
function sort_points(Data::Matrix{Float64})

    y_sorted_id = sortperm(Data[2,:])
    y_sorted = Data[:,y_sorted_id]
    
    pointsLeft_id = findall(y_sorted[1,:].<=y_sorted[1,1])
    leftElements = y_sorted_id[pointsLeft_id]
    pointsLeft = y_sorted[:,pointsLeft_id]

    pointsRight_id = findall(y_sorted[1,:].>y_sorted[1,1])
    rightElements = y_sorted_id[pointsRight_id]
    pointsRight = y_sorted[:,pointsRight_id]

    return hcat(pointsLeft, reverse(pointsRight,dims=2)), vcat(leftElements, reverse(rightElements))
end

""" 
    back_project(x::AbstractMatrix{Float64}, camera_matrix::AbstractMatrix{Float64}, dxdλ::AbstractMatrix{Float64}=nothing)

Project the 3D mesh to 2D image plane
    
# Arguments:
- `x::Matrix{Float64}{3,nNodes}`: 3D mesh grid
- `camera_matrix::AbstractMatrix{Float64}{3,3}`: Camera matrix

# Returns:
- `x2D::Matrix{Float64}{2,nNodes}`: 2D coordinates of the nodes
- `dπdx::Array{Float64, nNodes}`{nNodes×2×nParams}: ∇π(x)
"""
function back_project(x::AbstractMatrix{Float64}, camera_matrix::AbstractMatrix{Float64}, obj_pose::Vector{Float64}, dxdθ::AbstractArray{Float64}, height::Float64, rot_angle::Float64=0.0)
    
    nx = size(x,2)
    nθ = size(dxdθ,3)
    xNorm = zeros(3,nx)
    xProj = zeros(3,nx)
    dudθ = zeros(3,nx,nθ)

    x_rot = _rotate_round_z(x, rot_angle)  # Rotate the points around the z-axis by the object's yaw angle
    obj_trans = _set_cam_frame(obj_pose, height)

    R = obj_trans[1:3,1:3]    # rotation matrix from object frame to camera frame

    # transform point cloud wrt to camera frame 
    pad = ones(1,size(x,2))
    x_trans_padded = vcat(x_rot, pad)
    x_trans_mat = obj_trans*x_trans_padded

    x_trans = x_trans_mat[1:3,:]

    iter = 1:nx
    iterθ = 1:nθ
    for j in iter
        xNorm[1,j] = x_trans[1,j]/x_trans[3,j]
        xNorm[2,j] = x_trans[2,j]/x_trans[3,j]
        xNorm[3,j] = x_trans[3,j]/x_trans[3,j]

        xProj[:,j] = camera_matrix*xNorm[:,j]

        dπdx = ∇π(x_trans[:,j],camera_matrix)

        # iterate over the model parameters
        for i in iterθ
            dxdθi = R*dxdθ[:,:,i]           # rigid transformation applied to the gradient
            dudθ[:,j,i] = dπdx*dxdθi[:,j]
        end
    end

    # NOTE: Outout should be the 2D matrix removed for testing
    x2D = xProj[1:2,:]            # extract x and y coordinates
    dudθ2D = dudθ[1:2,:,:]

    return x2D, dudθ2D
end 

"""
    back_project(x, camera_matrix, obj_trans) -> x2D

Project 3D points into the 2D image plane without gradient computation.

# Arguments
- `x::AbstractMatrix{Float64}`: 3×N matrix of 3D coordinates (object frame).
- `camera_matrix::AbstractMatrix{Float64}`: 3×3 camera intrinsic matrix.
- `obj_trans::AbstractMatrix{Float64}`: 4×4 object-to-camera rigid transform.

# Returns
- `x2D::Matrix{Float64}`: 2×N matrix of pixel coordinates.
"""
function back_project(x::AbstractMatrix{Float64}, camera_matrix::AbstractMatrix{Float64}, obj_pose::Vector{Float64}, height::Float64, rot_angle::Float64=0.0)

    x_rot = _rotate_round_z(x, rot_angle)
    obj_trans = _set_cam_frame(obj_pose, height)

    pad = ones(1,size(x,2))
    x_trans_padded = vcat(x_rot, pad)
    x_trans_mat = obj_trans*x_trans_padded
    x_trans = x_trans_mat[1:3,:]

    xProj = project_to(x_trans, camera_matrix)
    x2D = xProj[1:2,:]

    return x2D
end

"""
    project_to(x, camera_matrix) -> xProj

Apply the pinhole camera projection to pre-transformed 3D points.

# Arguments
- `x::AbstractMatrix{Float64}`: 3×N matrix of points already in the camera frame.
- `camera_matrix::AbstractMatrix{Float64}`: 3×3 camera intrinsic matrix.

# Returns
- `xProj::Matrix{Float64}`: 3×N homogeneous pixel coordinates.
"""
function project_to(x::AbstractMatrix{Float64}, camera_matrix::AbstractMatrix{Float64})

    nx = size(x,2)
    xNorm = zeros(3,nx)
    xProj = zeros(3,nx)

    iter = 1:nx
    for i in iter
        xNorm[1,i] = x[1,i]/x[3,i]
        xNorm[2,i] = x[2,i]/x[3,i]
        xNorm[3,i] = x[3,i]/x[3,i]

        xProj[:,i] = camera_matrix*xNorm[:,i]
    end

    return xProj
end

"""
    ∇π(x, camera_matrix) -> dπdx

Jacobian of the pinhole projection at a single 3D point.

# Arguments
- `x::Array{Float64}`: length-3 vector `[X, Y, Z]` in the camera frame.
- `camera_matrix::AbstractMatrix{Float64}`: 3×3 camera intrinsic matrix.

# Returns
- `dπdx::Matrix{Float64}`: 3×3 Jacobian matrix `∂π/∂x`.
"""
function ∇π(x::Array{Float64},camera_matrix::AbstractMatrix{Float64})

    dπdx = camera_matrix*[1/x[3] 0 -x[1]/x[3]^2; 0 1/x[3] -x[2]/x[3]^2; 0 0 0]

    return dπdx
end

"""
    closest_point_contour(contour1, contour2) -> (hausdorff, [avg, chamfer])

Compute robust distances between two 2D contours using Hausdorff, average, and
Chamfer (RMS) metrics. If the contours have different numbers of points they are
first resampled to `min(n1, n2)` uniformly-spaced points via linear interpolation.

# Arguments
- `contour1::AbstractMatrix{Float64}`: 2×N₁ matrix of 2D contour points.
- `contour2::AbstractMatrix{Float64}`: 2×N₂ matrix of 2D contour points.

# Returns
- `[hausdorff_dist]::Vector{Float64}`: maximum point-wise distance.
- `[average_dist, chamfer_dist]::Vector{Float64}`: mean and RMS point-wise distances.
"""
function closest_point_contour(contour1::AbstractMatrix{Float64}, contour2::AbstractMatrix{Float64})
    # Ensure contours have same number of points or interpolate
    n1, n2 = size(contour1, 2), size(contour2, 2)
    
    if n1 != n2
        # Interpolate to same number of points
        n_target = min(n1, n2)
        t1 = range(0, 1, length=n1)
        t2 = range(0, 1, length=n2)
        t_target = range(0, 1, length=n_target)
        
        # Simple linear interpolation for each coordinate
        x1_interp = [LinearInterpolation(t1, contour1[1,:])(t) for t in t_target]
        y1_interp = [LinearInterpolation(t1, contour1[2,:])(t) for t in t_target]
        x2_interp = [LinearInterpolation(t2, contour2[1,:])(t) for t in t_target]
        y2_interp = [LinearInterpolation(t2, contour2[2,:])(t) for t in t_target]
        
        contour1_norm = hcat(x1_interp, y1_interp)'
        contour2_norm = hcat(x2_interp, y2_interp)'
    else
        contour1_norm = contour1
        contour2_norm = contour2
    end
    
    # Compute point-wise distances
    distances = [norm(contour1_norm[:,i] - contour2_norm[:,i]) for i in 1:size(contour1_norm,2)]
    
    # Multiple distance metrics
    hausdorff_dist = maximum(distances)
    average_dist = mean(distances)
    chamfer_dist = sqrt(mean(distances.^2))  # RMS distance
    
    return [hausdorff_dist], [average_dist, chamfer_dist]
end

"""
    detect_outlier_observations(contour_list; area_outlier_threshold, centroid_outlier_threshold,
                                perimeter_outlier_threshold, rolling_window, min_valid_neighbors)
        -> (valid_frames, outlier_frames)

Identify outlier frames in a time-series of 2D contours using local rolling-window
statistics on area, centroid position, and perimeter.

# Arguments
- `contour_list::AbstractArray`: vector of 2×Nᵢ contour matrices, one per frame.

# Keyword Arguments
- `area_outlier_threshold`: absolute area deviation threshold (default: `3e5`).
- `centroid_outlier_threshold`: combined centroid z-score threshold (default: `2.5`).
- `perimeter_outlier_threshold`: perimeter z-score threshold (default: `3.0`).
- `rolling_window`: half-width of the local statistics window (default: `20`).
- `min_valid_neighbors`: minimum neighbours required for outlier analysis (default: `3`).

# Returns
- `valid_frames::Vector{Int64}`: indices of frames classified as valid.
- `outlier_frames::Vector{Int64}`: indices of frames classified as outliers.
"""
function detect_outlier_observations(contour_list::AbstractArray;
                                   area_outlier_threshold=3e5,
                                   centroid_outlier_threshold=2.5,
                                   perimeter_outlier_threshold=3.0,
                                   rolling_window=20,
                                   min_valid_neighbors=3)

    frame_len = length(contour_list)
    valid_frames = Vector{Int64}()
    outlier_frames = Vector{Int64}()
    
    # Pre-compute geometric features for all contours
    areas = Float64[]
    perimeters = Float64[]
    centroids = Vector{Vector{Float64}}()
    
    for contour in contour_list
        # Area using shoelace formula
        area = 0.5 * abs(sum(contour[1,i] * contour[2,i+1] - contour[1,i+1] * contour[2,i] 
                            for i in 1:(size(contour,2)-1)))
        push!(areas, area)
        
        # Perimeter
        perimeter = sum(norm([contour[1,i+1] - contour[1,i], contour[2,i+1] - contour[2,i]]) 
                       for i in 1:(size(contour,2)-1))
        push!(perimeters, perimeter)
        
        # Centroid
        centroid = [mean(contour[1,:]), mean(contour[2,:])]
        push!(centroids, centroid)
    end
    
    # Convert centroids to matrix for easier processing
    centroid_x = [c[1] for c in centroids]
    centroid_y = [c[2] for c in centroids]
    
    # Analyze each frame for outliers
    for i in 1:frame_len
        is_outlier = false
        
        # Define analysis window around current frame
        window_start = max(1, i - rolling_window÷2)
        window_end = min(frame_len, i + rolling_window÷2)
        
        # Exclude current frame from statistics to avoid self-influence
        window_indices = [j for j in window_start:window_end if j != i]
        
        if length(window_indices) >= min_valid_neighbors
            # Local statistics for area
            local_areas = areas[window_indices]
            area_mean = mean(local_areas)
            area_std = std(local_areas)
            area_z_score = abs(areas[i] - area_mean)
            
            # Local statistics for centroid
            local_centroid_x = centroid_x[window_indices]
            centroid_x_mean = mean(local_centroid_x)
            centroid_x_std = std(local_centroid_x)
            
            local_centroid_y = centroid_y[window_indices]
            centroid_y_mean = mean(local_centroid_y)
            centroid_y_std = std(local_centroid_y)
            
            # Combined centroid displacement
            centroid_displacement = sqrt((centroid_x[i] - centroid_x_mean)^2 + (centroid_y[i] - centroid_y_mean)^2)
            expected_centroid_std = sqrt(centroid_x_std^2 + centroid_y_std^2)
            centroid_combined_z_score = centroid_displacement / (expected_centroid_std + 1e-10)
            
            # Local statistics for perimeter
            local_perimeters = perimeters[window_indices]
            perimeter_mean = mean(local_perimeters)
            perimeter_std = std(local_perimeters)
            perimeter_z_score = abs(perimeters[i] - perimeter_mean) / (perimeter_std + 1e-10)
            
            # Check for outliers
            if area_z_score > area_outlier_threshold
                is_outlier = true
            end
            
            # if centroid_combined_z_score > centroid_outlier_threshold
            #     is_outlier = true
            # end
            
            # if perimeter_z_score > perimeter_outlier_threshold
            #     is_outlier = true
            # end
            
            # Additional checks for extreme cases
            if areas[i] <= 0 || size(contour_list[i], 2) < 3
                is_outlier = true
            end
        end

        # Classify frame
        if is_outlier
            push!(outlier_frames, i)
        else
            push!(valid_frames, i)
            # p = Plots.Plots.scatter(contour_list[i][1,:], contour_list[i][2,:]; color=:green, label="Valid Frame $i, Area score: $(round(area_z_score,digits=2))")
            # xlabel!(p, "X (pixels)")
            # ylabel!(p, "Y (pixels)")
            # xlims!(p, 0, 2048)
            # ylims!(p, 0, 1536)
            # yflip!(p, true)  # if needed to match image coords
            # display(p)
            # sleep(0.5)  # Pause to ensure plot updates
        end
    end
    
    return valid_frames, outlier_frames
end

"""
    get_pose(frames) -> obj_pose

Estimate the object pose by averaging a stack of 4×4 homogeneous transformation
matrices across all frames.

# Arguments
- `frames::AbstractArray`: 4×4×T array of per-frame pose matrices.

# Returns
- `obj_pose::Matrix{Float64}`: 4×4 mean pose matrix.
"""
function get_pose(frames::AbstractArray)
    @info "Computing specimen pose from $(size(frames,3)) frames"
    pose_mat = zeros(Float64, 4,4)
    frame_len = size(frames,3)

    iter = 1:frame_len
    for i in iter
        pose_mat = pose_mat + frames[:,:,i]
    end
    obj_pose = pose_mat/frame_len
    return obj_pose
end
"""
filter_points(border, centerx)

Select the nodes on the right side of the centerline and sort them

# Arguments:
- `border::Matrix{Float64}{2,nbNodes}`: 2D coordinates of the border nodes
- `centerx::Float64`: x-coordinate of the centerline

# Returns:
- `newBorderxSrt::Vector{Float64}`: x coordinates of the sorted border nodes
- `newBorderySrt::Vector{Float64}`: y coordinates of the sorted border nodes
"""
function filter_points(border, centerx)

    # half_border = ElasticArray{Float64}(undef, 2, size(border,2))
    new_borderx = Array{Float64}(undef, 0)
    new_bordery = Array{Float64}(undef, 0)
    
    iter = 1:size(border,2)
    for i in iter
        if border[1,i] ≥ centerx
            # append!(half_border, border[:,i])
            push!(new_borderx, border[1,i])
            push!(new_bordery, border[2,i])
        end
    end
    ids = sortperm(new_bordery)

    newBorderxSrt = new_borderx[ids]
    newBorderySrt = new_bordery[ids]

    return newBorderxSrt, newBorderySrt
end

"""
    fit_curve(; border, borderx, bordery)

Fit a curve to the border nodes of the 2D mesh

# Arguments:
- `border::Matrix{Float64}{2,nbNodes}`: 2D coordinates of the border nodes
- `borderx::Vector{Float64}`: x coordinates of the border nodes
- `bordery::Vector{Float64}`: y coordinates of the border nodes

# Returns:
- `pi::Vector{Float64}`: x coordinates of the fitted curve
- `qi::Vector{Float64}`: y coordinates of the fitted curve
"""
function fit_curve(;border=nothing, borderx=nothing, bordery=nothing, samples=nothing)
    if isnothing(border)
        x = borderx
        y = bordery

        p = CubicSpline(x,y)

        pi = [p(i) for i in samples]
        
        return pi

    elseif isnothing(borderx) && isnothing(bordery)
        x = push!(border[1,:], border[1,1])
        y = push!(border[2,:], border[2,1])

        len = length(x)
        seq = 1:(len)

        p = CubicSpline(x,seq)
        q = CubicSpline(y,seq)
        pi = [p(i) for i in 1:0.1:len]
        qi = [q(i) for i in 1:0.1:len]

        return pi, qi
    else
        throw(AssertionError("No border provided"))
        return nothing
    end
end

"""
    fit_curve_2D(x,y, n)

Fit a curve to the border nodes of the 2D mesh

# Arguments:
- `x::Vector{Float64}`: x coordinates
- `y::Vector{Float64}`: y coordinates
- `n::Integer`: number of sampled points

# Returns:
- `points::Vector{Float64}`: vector of n sampled points
"""
function fit_curve_2D(x,y, n)
    spl = CubicSpline(x,y)
    points = [spl(i) for i in range(y[1],stop=y[end],length=n)]

    return points
end

"""
    rearrange(q, ne, ndim, IEN, basis_order)

Rearrange the connectivity vector to be visualised in paraview when lagrangian basis functions are used.

# Arguments:
- `ndim::Integer`: number of dimensions
- `IEN::Matrix{Float64}{nElem, nNodes}`: Connectivity matrix

# Returns:
- `IEN_new::Matrix{Int64}`: rearranged connectivity matrix
"""
function rearrange(ndim, IEN)
    if size(IEN, 1) != 27
        return copy(IEN)
    end
    IEN_new = zeros(Int64,size(IEN))
    if ndim == 2
        error("rearrange: 2D case is not yet implemented (only 3D is supported)")
    elseif ndim == 3
        iter = 1:size(IEN,2)
        for n in iter
            IEN_new[1,n] = IEN[1,n]
            IEN_new[2,n] = IEN[2,n]
            IEN_new[3,n] = IEN[3,n]
            IEN_new[4,n] = IEN[4,n]
            IEN_new[5,n] = IEN[5,n]
            IEN_new[6,n] = IEN[6,n]
            IEN_new[7,n] = IEN[7,n]
            IEN_new[8,n] = IEN[8,n]
            IEN_new[9,n] = IEN[9,n]
            IEN_new[10,n] = IEN[10,n]
            IEN_new[11,n] = IEN[11,n]
            IEN_new[12,n] = IEN[12,n]
            IEN_new[13,n] = IEN[13,n]
            IEN_new[14,n] = IEN[14,n]
            IEN_new[15,n] = IEN[15,n]
            IEN_new[16,n] = IEN[16,n]
            IEN_new[17,n] = IEN[17,n]
            IEN_new[18,n] = IEN[18,n]
            IEN_new[19,n] = IEN[20,n]
            IEN_new[20,n] = IEN[19,n]
            IEN_new[21,n] = IEN[24,n]
            IEN_new[22,n] = IEN[22,n]
            IEN_new[23,n] = IEN[21,n]
            IEN_new[24,n] = IEN[23,n]
            IEN_new[25,n] = IEN[25,n]
            IEN_new[26,n] = IEN[26,n]
            IEN_new[27,n] = IEN[27,n]
        end
    end
    return IEN_new
end

"""
    add_noise(obsScene; nFactor) -> (nScene, nSplinex, nSpliney, pd)

Add zero-mean Gaussian noise to a sequence of 2D observation contours.

# Arguments
- `obsScene`: vector of 2×Nᵢ contour matrices.

# Keyword Arguments
- `nFactor`: standard deviation of the Gaussian noise (default: `0` — no noise).

# Returns
- `nScene::Vector`: noisy contour matrices.
- `nSplinex::Vector`: x-coordinate vectors of the noisy contours.
- `nSpliney::Vector`: y-coordinate vectors of the noisy contours.
- `pd::Normal`: the noise distribution used.
"""
function add_noise(obsScene; nFactor=0)

    nScene = AbstractArray[]      
    nSplinex = AbstractArray[]
    nSpliney = AbstractArray[] 
    pd = Normal(0,nFactor)

    iter=1:size(obsScene,1)                           
    for i in iter
        Data = obsScene[i]
        w = rand(pd, size(Data))
        Data = Data + w
        push!(nScene,Data)
        push!(nSplinex, Data[1,:])
        push!(nSpliney, Data[2,:])
    end
    return nScene, nSplinex, nSpliney, pd
end

"""
    get_height(μ_tp, H_0) -> h

Compute cumulative heights of mesh layers from per-step displacements.

# Arguments
- `μ_tp::Vector{Float64}`: vector of incremental top-plate displacements.
- `H_0::Float64`: initial height of the specimen.

# Returns
- `h::Vector{Float64}`: length `length(μ_tp)+1` vector of cumulative heights,
  where `h[1] == H_0`.
"""
function get_height(μ_tp::Vector{Float64}, H_0::Float64)
    # Get the height of the mesh
    h = H_0*ones(Float64, (length(μ_tp)+1))
    for i in eachindex(μ_tp)
        h[i+1] = h[i] + μ_tp[i]
    end
    return h
end

"""
    eval_on_cylinder(mdl, nsub, sol_u) -> (NodeList_, IEN_list, plot_u)

Evaluate a displacement solution on a refined visualization mesh by subdividing
each element `2^nsub` times per axis using the NURBS/IGA basis functions.

# Arguments
- `mdl::AbstractModel`: model containing `mesh_x` (geometry mesh) and `mesh_u`
  (displacement mesh) fields.
- `nsub::Int64`: number of subdivision levels; each element is split into
  `2^(3*nsub)` sub-elements.
- `sol_u`: displacement solution array (3×nNodes).

# Returns
- `NodeList_::Matrix{Float64}`: 3×M refined node coordinates.
- `IEN_list::Matrix{Int64}`: connectivity array for the refined mesh.
- `plot_u::Matrix{Float64}`: 3×M interpolated displacement at refined nodes.
"""
function eval_on_cylinder(mdl::AbstractModel, nsub::Int64, sol_u)

    @unpack NodeList, IEN, volume_element_shape, basis_order = mdl.mesh_x

    NodeList_x_cached = NodeList
    IEN_x_cached = IEN
    element_shape_x_cached, basis_order_x_cached = volume_element_shape, basis_order

    @unpack IEN, volume_element_shape, basis_order = mdl.mesh_u

    IEN_u_cached = IEN
    element_shape_u_cached, basis_order_u_cached = volume_element_shape, basis_order

    NodeList_ = zeros(Float64, 3, (2^(nsub*3))*size(IEN_x_cached,2)*size(IEN_x_cached,1))
    IEN_list = zeros(Int64, size(IEN_x_cached,1), (2^(nsub*3))*size(IEN_x_cached,2))
    plot_u = zeros(Float64, 3, (2^(nsub*3))*size(IEN_u_cached,2)*size(IEN_u_cached,1))
    # get parent element
    lPoints = zeros(27,3)
    lPoints[1,:] = [-1 , -1, -1]
    lPoints[2,:] = [1 , -1, -1]
    lPoints[3,:] = [1 , 1, -1]
    lPoints[4,:] = [-1 , 1, -1]

    lPoints[5,:] = [-1 , -1, 1]
    lPoints[6,:] = [1 , -1, 1]
    lPoints[7,:] = [1 , 1, 1]
    lPoints[8,:] = [-1 , 1, 1]
    
    lPoints[9,:] = [0 , -1, -1]
    lPoints[10,:] = [1 , 0, -1]
    lPoints[11,:] = [0 , 1, -1]
    lPoints[12,:] = [-1 , 0, -1]

    lPoints[13,:] = [0 , -1, 1]
    lPoints[14,:] = [1 , 0, 1]
    lPoints[15,:] = [0 , 1, 1]
    lPoints[16,:] = [-1 , 0, 1]

    lPoints[17,:] = [-1 , -1, 0]
    lPoints[18,:] = [1 , -1, 0]
    lPoints[19,:] = [1 , 1, 0]
    lPoints[20,:] = [-1 , 1, 0]

    lPoints[21,:] = [0 , -1, 0]
    lPoints[22,:] = [1 , 0, 0]
    lPoints[23,:] = [0 , 1, 0]
    lPoints[24,:] = [-1 , 0, 0]

    lPoints[25,:] = [0, 0, -1]
    lPoints[26,:] = [0, 0, 1]
    lPoints[27,:] = [0, 0, 0]

    scale = 1/(2^nsub)

    eiter = 1:size(IEN_x_cached,2)
    nodeiter = 1:size(lPoints,1)
    for e in eiter
        for xsub in 1:2^nsub
            for ysub in 1:2^nsub
                for zsub in 1:2^nsub
                    # counter within the elements
                    cnte = (e-1)*(2^(nsub*3)) + (xsub-1)*(2^(nsub*2)) + (ysub-1)*(2^(nsub)) + zsub
                    offset = [(xsub-1)*scale, (ysub-1)*scale, (zsub-1)*scale]
                    for j in nodeiter
                        scaledPoint = offset .+ scale*(0.5.+0.5*lPoints[j,:])
                        scaledPoint = -1 .+ 2*scaledPoint
                        
                        Re, ΔRe = basis_function(scaledPoint[1],scaledPoint[2],scaledPoint[3], element_shape_x_cached, basis_order_x_cached)
                        NodeList_[:,(cnte-1)*size(IEN_x_cached,1)+j] = Re'*NodeList_x_cached[:,IEN_x_cached[:,e]]'
                        IEN_list[j,cnte] = (cnte-1)*size(IEN_x_cached,1)+j

                        Re_u, _ = basis_function(scaledPoint[1],scaledPoint[2],scaledPoint[3], element_shape_u_cached, basis_order_u_cached)
                        plot_u[:,(cnte-1)*size(IEN_x_cached,1)+j] = Re_u'*sol_u[:,IEN_u_cached[:,e]]'
                        
                        # println(size(plot_u))
                    end
                end
            end
        end
    end
    return NodeList_, IEN_list, plot_u
end

"""
    get_lagrange_proj(IEN_cp, IEN_l, C_vol, X_cp) -> C

Assemble the global Lagrange extraction operator from per-element Bézier extraction
matrices `C_vol`.

# Arguments
- `IEN_cp`: control-point connectivity array (nCp_per_elem × nElem).
- `IEN_l`: Lagrange-node connectivity array (nL_per_elem × nElem).
- `C_vol`: per-element extraction operators (nCp_per_elem × nL_per_elem × nElem).
- `X_cp`: control-point coordinate matrix (3 × nCp), used only for size inference.

# Returns
- `C::Matrix{Float64}`: nCp × nLagrange global extraction matrix.
"""
function get_lagrange_proj(IEN_cp, IEN_l, C_vol, X_cp)
    n_elem = size(IEN_cp,2)
    @assert size(IEN_cp,2) == size(IEN_l,2) "Number of elements in IEN_cp and IEN_l must be the same"

    n_cp = size(X_cp,2)
    n_lagrange = maximum(IEN_l)
    
    @debug "get_lagrange_proj: n_cp=$n_cp, n_lagrange=$n_lagrange"

    C = zeros(Float64, n_cp, n_lagrange)

    e_iter = 1:n_elem
    for e in e_iter
        C[IEN_cp[:,e],IEN_l[:,e]] = C_vol[:,:,e]
    end
    return C
end

"""
    get_nurbs_2_lagrange_proj(IEN_cp, IEN_l, C_vol, X_cp, W_cp) -> M

Compute the NURBS-to-Lagrange projection matrix weighted by the rational weights.

# Arguments
- `IEN_cp`: control-point connectivity (nCp_per_elem × nElem).
- `IEN_l`: Lagrange connectivity (nL_per_elem × nElem).
- `C_vol`: per-element Bézier extraction operators.
- `X_cp`: control-point coordinates (3 × nCp).
- `W_cp`: control-point weights (nCp,).

# Returns
- `M::Matrix{Float64}`: nCp × nLagrange weighted projection matrix.
"""
function get_nurbs_2_lagrange_proj(IEN_cp, IEN_l, C_vol, X_cp, W_cp)

    P = get_lagrange_proj(IEN_cp, IEN_l, C_vol, X_cp)
    n_l = size(P,2)
    n_cp = size(P,1)
    M = zeros(Float64, n_cp, n_l)
    W_l = P'*W_cp

    n_iter = 1:n_l
    for i in n_iter        
        M[:,i] = (diagm(W_cp)/W_l[i])*P[:,i]
    end
    return M
end

"""
    get_lagrange_pts(IEN_cp, IEN_l, C_vol, X_cp, W_cp) -> X_l

Compute the physical coordinates of Lagrange interpolation points from a NURBS
control mesh using the weighted extraction operator.

# Arguments
- `IEN_cp`: control-point connectivity (nCp_per_elem × nElem).
- `IEN_l`: Lagrange connectivity (nL_per_elem × nElem).
- `C_vol`: per-element Bézier extraction operators.
- `X_cp`: control-point coordinates (3 × nCp).
- `W_cp`: control-point weights (nCp,).

# Returns
- `X_l::Matrix{Float64}`: 3 × nLagrange matrix of Lagrange point coordinates.
"""
function get_lagrange_pts(IEN_cp, IEN_l, C_vol, X_cp, W_cp)

    P = get_lagrange_proj(IEN_cp, IEN_l, C_vol, X_cp)
    n_l = size(P,2)
    n_cp = size(P,1)
    X_l = zeros(Float64, 3, n_l)
    W_l = P'*W_cp

    n_iter = 1:n_l
    for i in n_iter
        X_l[:,i] = (X_cp*diagm(W_cp)/W_l[i])*P[:,i]
    end
    return X_l
end

function _set_cam_frame(cam_pose, height)

    cam_pose_ = cam_pose - [0.0, 0.0, height/2]  # translate to camera frame
    z = [0.0, 0.0, 1.0]
    z_ = -cam_pose_ / norm(cam_pose_)
    x_ = -cross(z, z_) / norm(cross(z, z_))
    y_ = cross(z_, x_) / norm(cross(z_, x_))

    T_o_c = [x_' -dot(x_, cam_pose); y_' -dot(y_, cam_pose); z_' -dot(z_, cam_pose); 0 0 0 1]

    return T_o_c
end

function _rotate_round_z(node_list, d_angle::Float64=0.0)

    angle = d_angle/180*pi # convert to radians

    # Create rotation matrix
    R = [cos(angle) -sin(angle) 0;
         sin(angle)  cos(angle) 0;
         0           0          1]
    
    # Apply rotation to each node in the mesh
    rotated_nodes = R * node_list
    
    # Return a new mesh with rotated nodes
    return rotated_nodes
end