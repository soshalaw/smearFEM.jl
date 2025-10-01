using LinearAlgebra
using DataInterpolations
using ElasticArrays
using ConvexHulls2d
import ConvexHulls2d as ch
using Distributions
using Parameters
using Statistics
using StatsPlots

""" 
    Extract_borders(NodeList::Matrix{Float64}, camera_matrix::AbstractMatrix{Float64}, BorderNodesList::Vector{Vector{Int64}}, nNodes::Int64, GRAD::Bool=false, dqdθ::AbstractMatrix{Float64}=zeros(2,2), SIDES::Bool=false)

Project the 3D mesh to 2D image plane and extract the border nodes (left and right)
    
# Arguments:
- `NodeList::Matrix{Float64}{ndim,nNodes}` : coordinates of the nodes
- `camera_matrix::AbstractMatrix{Float64}{3,3}` : Camera matrix
- `BorderNodesList::Vector{Vector{Any}{4,N}`: List of border nodes
- `nNodes::Int64`: number of nodes
- `GRAD::Bool=false`: gradient flag
- `dqdθ::AbstractMatrix{Float64}{2,2}`: gradient of the solution
- `SIDES::Bool=false`: sides flag

# Returns:
- `NodeList::Matrix{Float64}{ndim,nbNodes}`: 2D coordinates of the border nodes
- `BorderNodes::Vector{Int}`: Indexes of the border nodes
"""
function extract_borders(NodeList::AbstractMatrix{Float64}, camera_matrix::AbstractMatrix{Float64}, camera_pose::AbstractMatrix{Float64}, n_nodes::Int64; BorderNodesList::Vector{Int64}=zeros(Int64,0),GRAD::Bool=false, dqdθ::AbstractArray{Float64}=zeros(2,2,2), SIDES::Bool=false)

    if length(BorderNodesList) != 0
        surface_nodes = NodeList[:,BorderNodesList]  # extract the border nodes from the NodeList
    else        
        surface_nodes = NodeList
    end
    
    if GRAD
        ∇surface_nodes = dqdθ[:,BorderNodesList,:] 
        surface_pts_2d, ∇surface_pts_2d = back_project(surface_nodes, camera_matrix, camera_pose, ∇surface_nodes) 
    else
        surface_pts_2d = back_project(surface_nodes, camera_matrix, camera_pose) 
    end

    left_border_pts = zeros(Float64,2,(n_nodes))        # vector to store indexes of the border nodes
    right_border_pts = zeros(Float64,2,(n_nodes))       # vector to store indexes of the border nodes

    left_border_nodes = Vector{Int64}(undef, 0)        # vector to store indexes of the border nodes
    right_border_nodes = Vector{Int64}(undef, 0)       # vector to store indexes of the border nodes
    # TopBorderNodes = Vector{Int64}(undef, 0)         # vector to store indexes of the border nodes
    # BottomBorderNodes = Vector{Int64}(undef, 0)      # vector to store indexes of the border nodes

    sz_side = size(surface_pts_2d,2)÷(n_nodes)           # size of each layer

    p = Array{Vector{Float64}}(undef,0)

    iter = 1:size(surface_pts_2d,2)
    for i in iter
        push!(p, surface_pts_2d[:,i])
    end

    hull = ch.ConvexHull(p)

    top_bottom_nodes = ch.indices(hull)
    # border_pts = surface_pts_2d[:,border_pt_ids]
    # border_pts_sorted, border_pts_sorted_ids = sort_points(border_pts)

    for layers in 1:n_nodes                                        # loop through each layer
        nodes = surface_pts_2d[:,(layers-1)*sz_side+1:layers*sz_side]
        min_node = (layers-1)*sz_side + argmin(nodes[1,:])
        max_node = (layers-1)*sz_side + argmax(nodes[1,:])
        push!(left_border_nodes, min_node)
        push!(right_border_nodes, max_node)
        left_border_pts[:,layers] = surface_pts_2d[:,min_node]         # left border nodes
        right_border_pts[:,layers] = surface_pts_2d[:,max_node]        # right border nodes
        # if Layers == nNodes
        #     nodeIdi = 1:size(nodes,2)
        #     for nodeId in nodeIdi
        #         if nodes[2,nodeId] > surface_pts_2d[2,min_node]    
        #             push!(TopBorderNodes, (Layers-1)*szSide+nodeId)
        #         end
        #     end 
        # elseif Layers == 1
        #     nodeIdi = 1:size(nodes,2)
        #     for nodeId in nodeIdi
        #         if nodes[2,nodeId] < surface_pts_2d[2,min_node] 
        #             push!(BottomBorderNodes, nodeId)
        #         end
        #     end 
        # end
    end  

    border_pt_ids = vcat(left_border_nodes, top_bottom_nodes, right_border_nodes)
    border_pts = surface_pts_2d[:,border_pt_ids]
    border_pts_sorted, border_pts_sorted_ids = sort_points(border_pts)

    if GRAD
        ∇border_pts = ∇surface_pts_2d[:,border_pt_ids,:]
        ∇border_pts_sorted = ∇border_pts[:,border_pts_sorted_ids,:]
        if SIDES
            side_pts, index = get_sides(border_pts_sorted)

            side_pts_ids = border_pts_sorted_ids[index]
            ∇side_pts = ∇border_pts_sorted[:,side_pts_ids,:]

            return side_pts, ∇side_pts, surface_pts_2d
        else
            return border_pts_sorted, ∇border_pts_sorted, surface_pts_2d, ∇surface_pts_2d
        end
    else
        if SIDES
            side_pts, index = get_sides(border_pts_sorted)
    
            side_pts = border_pts_sorted[index]
            side_pts_ids = border_pts_sorted_ids[index]
    
            return side_pts, surface_pts_2d
        else
            return border_pts_sorted, surface_pts_2d
        end
    end 
end

function extract_borders(NodeList::AbstractMatrix{Float64}, camera_matrix::AbstractMatrix{Float64}, camera_pose::AbstractMatrix{Float64}; BorderNodesList::Vector{Int64}=zeros(Int64,0), GRAD::Bool=false, dqdθ::AbstractArray{Float64}=zeros(2,2,2), SIDES::Bool=false)
    
    if length(BorderNodesList) != 0
        surface_nodes = NodeList[:,BorderNodesList]  # extract the border nodes from the NodeList
    else        
        surface_nodes = NodeList
    end

    if GRAD
        ∇surface_nodes = dqdθ[:,BorderNodesList,:] 
        surface_pts_2d, ∇surface_pts_2d = back_project(surface_nodes, camera_matrix, camera_pose, ∇surface_nodes) 
    else
        surface_pts_2d = back_project(surface_nodes, camera_matrix, camera_pose) 
    end

    p = Array{Vector{Float64}}(undef,0)

    iter = 1:size(surface_pts_2d,2)
    for i in iter
        push!(p, surface_pts_2d[:,i])
    end

    hull = ch.ConvexHull(p)

    border_pt_ids = ch.indices(hull)
    border_pts = surface_pts_2d[:,border_pt_ids]
    border_pts_sorted, border_pts_sorted_ids = sort_points(border_pts)

    if GRAD
        ∇border_pts = ∇surface_pts_2d[:,border_pt_ids,:]
        ∇border_pts_sorted = ∇border_pts[:,border_pts_sorted_ids,:]
        if SIDES
            side_pts, index = get_sides(border_pts_sorted)

            side_pts_ids = border_pts_sorted_ids[index]
            ∇side_pts = ∇border_pts_sorted[:,side_pts_ids,:]

            return side_pts, ∇side_pts, surface_pts_2d
        else
            return border_pts_sorted, ∇border_pts_sorted, surface_pts_2d, ∇surface_pts_2d
        end
    else
        if SIDES
            side_pts, index = get_sides(border_pts_sorted)
    
            side_pts = border_pts_sorted[index]
            side_pts_ids = border_pts_sorted_ids[index]
    
            return side_pts, surface_pts_2d
        else
            return border_pts_sorted, surface_pts_2d
        end
    end
end

function get_sides(Data::Matrix{Float64})

    indexes = []
    dataIter = 1:(size(Data,2)-1)
    for i in dataIter
        r = Data[:,i+1] - Data[:,i]
        sθ = dot(r,[0,1])/norm(r)
        if abs(sθ) >= sin(π/4)
            push!(indexes,i)
        elseif i != 1
            r_ = Data[:,i] - Data[:,i-1]
            sθ_ = dot(r_,[0,1])/norm(r_)
            if abs(sθ_) >= sin(π/4)
                push!(indexes,i)
            end
        end
    end
    sides = Data[:,indexes]
    return sides, indexes
end

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
function back_project(x::AbstractMatrix{Float64}, camera_matrix::AbstractMatrix{Float64}, obj_trans::AbstractMatrix{Float64}, dxdθ::AbstractArray{Float64})
    
    nx = size(x,2)
    nθ = size(dxdθ,3)
    xNorm = zeros(3,nx)
    xProj = zeros(3,nx)
    dudθ = zeros(3,nx,nθ)
    R = obj_trans[1:3,1:3]    # rotation matrix from object frame to camera frame

    # transform point cloud wrt to camera frame 
    pad = ones(1,size(x,2))
    x_trans_padded = vcat(x, pad)
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

function back_project(x::AbstractMatrix{Float64}, camera_matrix::AbstractMatrix{Float64}, obj_trans::AbstractMatrix{Float64})

    pad = ones(1,size(x,2))
    x_trans_padded = vcat(x, pad)
    x_trans_mat = obj_trans*x_trans_padded

    x_trans = x_trans_mat[1:3,:]

    p = Plots.scatter3d(x_trans[1,:], x_trans[2,:], x_trans[3,:]; markersize=1, label="Transformed Points")
    xlabel!(p, "X (mm)")
    ylabel!(p, "Y (mm)")
    zlabel!(p, "Z (mm)")
    display(p)
    
    xProj = project_to(x_trans, camera_matrix)

    x2D = xProj[1:2,:]            # extract x and y coordinates

    return x2D
end 

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

function ∇π(x::Array{Float64},camera_matrix::AbstractMatrix{Float64})

    dπdx = camera_matrix*[1/x[3] 0 -x[1]/x[3]^2; 0 1/x[3] -x[2]/x[3]^2; 0 0 0]

    return dπdx
end

function closest_point_contour(contour1::AbstractMatrix{Float64}, contour2::AbstractMatrix{Float64})
    """
    Compute robust distance between two contours using multiple metrics
    Returns: (hausdorff_distance, average_distance, chamfer_distance)
    """
    
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

function filter_frames(contour_list::AbstractArray; 
                     distance_threshold=0.05,
                     area_change_threshold=0.15,
                     perimeter_change_threshold=0.10,
                     centroid_shift_threshold=0.08,
                     shape_deformation_threshold=0.12,
                     consistency_window=3,
                     neighbor_window=4,
                     neighbor_weight=7.0,
                     display_plots=true)

    frame_len = length(contour_list)
    pos_calib_frames = Vector{Int64}()
    compression_frames = Vector{Int64}()

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

    # First pass: compute raw scores for all frames
    raw_scores = Float64[]
    frame_iter = 4:frame_len
    
    for i in frame_iter
        contour_curr = contour_list[i]
        contour_prev = contour_list[i-3]

        # Feature 1: Point-wise distance (your original metric)
        diff_t, _ = closest_point([contour_curr], [contour_prev])
        distance_score = diff_t[1] > distance_threshold ? 1.0 : 0.0

        # Feature 2: Area change
        area_change = abs(areas[i] - areas[i-3]) / max(areas[i-3], 1e-10)
        area_score = area_change > area_change_threshold ? 1.0 : 0.0

        # Feature 3: Perimeter change
        perimeter_change = abs(perimeters[i] - perimeters[i-3]) / max(perimeters[i-3], 1e-10)
        perimeter_score = perimeter_change > perimeter_change_threshold ? 1.0 : 0.0

        # Feature 4: Centroid shift
        centroid_shift = norm(centroids[i] - centroids[i-3])
        centroid_score = centroid_shift > centroid_shift_threshold ? 1.0 : 0.0

        # Feature 5: Shape deformation (Hausdorff-like distance)
        max_displacement = maximum([norm([contour_curr[1,j] - contour_prev[1,j], 
                                         contour_curr[2,j] - contour_prev[2,j]]) 
                                   for j in 1:min(size(contour_curr,2), size(contour_prev,2))])
        shape_score = max_displacement > shape_deformation_threshold ? 1.0 : 0.0

        # Raw compression score (before neighbor influence)
        raw_compression_score = 0.3 * distance_score + 
                               0.25 * area_score + 
                               0.2 * perimeter_score + 
                               0.15 * centroid_score + 
                               0.1 * shape_score
        
        push!(raw_scores, raw_compression_score)
    end

    # Second pass: apply neighbor-aware decision making
    decisions = String[]
    final_scores = Float64[]
    
    for (idx, i) in enumerate(frame_iter)
        # Get current frame's raw score
        current_raw_score = raw_scores[idx]
        
        # Collect neighbor scores within the window
        neighbor_scores = Float64[]
        
        # Look at previous frames
        for offset in 1:neighbor_window
            neighbor_frame = i - offset
            if neighbor_frame >= 4  # Ensure we have valid frame indices
                neighbor_idx = neighbor_frame - 3  # Convert to raw_scores index
                if neighbor_idx >= 1 && neighbor_idx <= length(raw_scores)
                    push!(neighbor_scores, raw_scores[neighbor_idx])
                end
            end
        end
        
        # Look at future frames
        for offset in 1:neighbor_window
            neighbor_frame = i + offset
            if neighbor_frame <= frame_len
                neighbor_idx = neighbor_frame - 3  # Convert to raw_scores index
                if neighbor_idx >= 1 && neighbor_idx <= length(raw_scores)
                    push!(neighbor_scores, raw_scores[neighbor_idx])
                end
            end
        end
        
        # Calculate neighbor influence
        neighbor_influence = 0.0
        if !isempty(neighbor_scores)
            # Strong neighbor influence: if most neighbors are compression, bias toward compression
            neighbor_compression_ratio = sum(s -> s > 0.5, neighbor_scores) / length(neighbor_scores)
            
            if neighbor_compression_ratio >= 0.7  # 70% of neighbors are compression
                neighbor_influence = neighbor_weight
            elseif neighbor_compression_ratio <= 0.3  # 70% of neighbors are position calibration
                neighbor_influence = -neighbor_weight
            else
                # Moderate influence based on average neighbor score
                avg_neighbor_score = mean(neighbor_scores)
                neighbor_influence = neighbor_weight * (avg_neighbor_score - 0.5)
            end
        end
        
        # Apply neighbor influence to get final score
        final_score = current_raw_score + neighbor_influence
        final_score = clamp(final_score, 0.0, 1.0)  # Keep within [0,1] bounds
        push!(final_scores, final_score)
        
        # Make decision based on final score
        current_decision = final_score > 0.5 ? "compression" : "position_calibration"
        push!(decisions, current_decision)
        
        # Apply consistency filter (optional additional smoothing)
        if length(decisions) >= consistency_window
            recent_decisions = decisions[end-consistency_window+1:end]
            compression_votes = count(d -> d == "compression", recent_decisions)
            
            # Majority vote for final classification
            if compression_votes > consistency_window ÷ 2
                status = "compression"
                push!(compression_frames, i)
            else
                status = "position_calibration"
                push!(pos_calib_frames, i)
            end
        else
            # For initial frames, use direct decision
            status = current_decision
            if status == "compression"
                push!(compression_frames, i)
            else
                push!(pos_calib_frames, i)
            end
        end

        # Optional visualization
        if display_plots
            contour_curr = contour_list[i]
            plt = Plots.plot(contour_curr[1,:], contour_curr[2,:], 
                           title="Frame $i: $status\nRaw=$(round(current_raw_score, digits=2)) Final=$(round(final_score, digits=2))", 
                           legend=false)
            xlims!(plt, 0, 2048)
            ylims!(plt, 0, 1536)
            yflip!(plt, true)
            display(plt)
            sleep(0.5)
        end
        
        # Enhanced logging with neighbor information
        neighbor_info = isempty(neighbor_scores) ? "no_neighbors" : 
                       "$(round(mean(neighbor_scores), digits=2))($(length(neighbor_scores)))"
        
        println("Frame $i: $status " *
                "(raw=$(round(current_raw_score, digits=3)), " *
                "neighbors=$neighbor_info, " *
                "influence=$(round(neighbor_influence, digits=3)), " *
                "final=$(round(final_score, digits=3)))")
    end

    return pos_calib_frames, compression_frames
end

function get_pose(frames::AbstractArray)
    pose_mat = zeros(Float64, 4,4)
    frame_len = size(frames,3)

    println("Number of frames: ", frame_len)
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

    out = vcat(newBorderxSrt', newBorderySrt')
    
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
        AssertionError("No border provided")
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
    rearrange(q, ne, ndim, IEN, FunctionClass)

Rearrange the connectivity vector to be visualised in paraview when langrangian basis functions are used.

# Arguments:
- `ndim::Integer`: number of dimensions
- `IEN::Matrix{Float64}{nElem, nNodes}`: Connectivity matrix

# Returns:
- `IEN_new::Matrix{Int64}`: rearranged connectivity matrix
"""
function rearrange(ndim, IEN) 
    IEN_new = zeros(Int64,size(IEN))
    if ndim == 2
        # TODO: Implement the 2D case
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

function get_height(μ_tp::Vector{Float64}, H_0::Float64)
    # Get the height of the mesh
    h = H_0*ones(Float64, (size(μ_tp,1)+1))
    iter = 1:length(h)
    for i::Int in iter
        if i == 1
            h[i] = h[i]
        else
            h[i] = h[i-1] + μ_tp[i-1]
        end
    end
    return h
end

function eval_on_cylinder(mdl::AbstractModel, nsub::Int64, sol_u)

    @unpack NodeList, IEN, C_vol, W, FunctionClass = mdl.mesh_x

    NodeList_x_cached = NodeList
    C_vol_x_cached = C_vol
    IEN_x_cached = IEN
    W_x_cached = W
    FunctionClass_x_cached = FunctionClass

    @unpack NodeList, IEN, C_vol, W, FunctionClass = mdl.mesh_u

    C_vol_u_cached = C_vol
    IEN_u_cached = IEN
    W_u_cached = W
    FunctionClass_u_cached = FunctionClass

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
                        
                        # get the NURBs and the gradients
                        Re, ΔRe = basis_function(scaledPoint[1],scaledPoint[2],scaledPoint[3], C_vol_x_cached[:,:,e], W_x_cached[IEN_x_cached[:,e]], FunctionClass_x_cached)
                        NodeList_[:,(cnte-1)*size(IEN_x_cached,1)+j] = Re'*NodeList_x_cached[:,IEN_x_cached[:,e]]'
                        IEN_list[j,cnte] = (cnte-1)*size(IEN_x_cached,1)+j

                        Re_u, _ = basis_function(scaledPoint[1],scaledPoint[2],scaledPoint[3], FunctionClass_u_cached) #TODO generalise to both cases
                        plot_u[:,(cnte-1)*size(IEN_x_cached,1)+j] = Re_u'*sol_u[:,IEN_u_cached[:,e]]'
                        
                        # println(size(plot_u))
                    end
                end
            end
        end
    end
    return NodeList_, IEN_list, plot_u
end

function get_lagrange_proj(IEN_cp, IEN_l, C_vol, X_cp)
    n_elem = size(IEN_cp,2)
    @assert size(IEN_cp,2) == size(IEN_l,2) "Number of elements in IEN_cp and IEN_l must be the same"

    n_cp = size(X_cp,2)
    n_lagrange = maximum(IEN_l)
    
    println("n_cp: ", n_cp)
    println("n_lagrange: ", n_lagrange)

    C = zeros(Float64, n_cp, n_lagrange)

    e_iter = 1:n_elem
    for e in e_iter
        C[IEN_cp[:,e],IEN_l[:,e]] = C_vol[:,:,e]
    end
    return C
end

"""
    function get_lagrange_pts(IEN_cp, IEN_l, C_vol, X_cp, W_cp)

    Funtion to obtain the lagrange points given the control mesh, extracton oparator and the weights.
    
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