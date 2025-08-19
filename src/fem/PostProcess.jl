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
function extract_borders(NodeList::Matrix{Float64}, camera_matrix::AbstractMatrix{Float64}, camera_pose::AbstractMatrix{Float64}, BorderNodesList::Vector{Int64}, nNodes::Int64; GRAD::Bool=false, dqdθ::AbstractArray{Float64}=zeros(2,2,2), SIDES::Bool=false)

    surface_nodes = NodeList[:,BorderNodesList]  # extract the border nodes from the NodeList
    
    if GRAD
        ∇surface_nodes = dqdθ[:,BorderNodesList,:] 
        SurfacePts2D, ∇SurfacePts2D = back_project(surface_nodes, camera_matrix, camera_pose, ∇surface_nodes) 
    else
        SurfacePts2D = back_project(surface_nodes, camera_matrix, camera_pose) 
    end

    LeftborderPts = zeros(Float64,2,(nNodes))        # vector to store indexes of the border nodes
    RightborderPts = zeros(Float64,2,(nNodes))       # vector to store indexes of the border nodes

    LeftborderNodes = Vector{Int64}(undef, 0)        # vector to store indexes of the border nodes
    RightborderNodes = Vector{Int64}(undef, 0)       # vector to store indexes of the border nodes
    TopBorderNodes = Vector{Int64}(undef, 0)         # vector to store indexes of the border nodes
    BottomBorderNodes = Vector{Int64}(undef, 0)      # vector to store indexes of the border nodes

    szSide = size(SurfacePts2D,2)÷(nNodes)           # size of each layer

    for Layers in 1:nNodes                                        # loop through each layer
        nodes = SurfacePts2D[:,(Layers-1)*szSide+1:Layers*szSide]
        minNode = (Layers-1)*szSide + argmin(nodes[1,:])
        maxNode = (Layers-1)*szSide + argmax(nodes[1,:])
        push!(LeftborderNodes, minNode)
        push!(RightborderNodes, maxNode)
        LeftborderPts[:,Layers] = SurfacePts2D[:,minNode]         # left border nodes
        RightborderPts[:,Layers] = SurfacePts2D[:,maxNode]        # right border nodes
        if Layers == nNodes
            nodeIdi = 1:size(nodes,2)
            for nodeId in nodeIdi
                if nodes[2,nodeId] > SurfacePts2D[2,minNode]    
                    push!(TopBorderNodes, (Layers-1)*szSide+nodeId)
                end
            end 
        elseif Layers == 1
            nodeIdi = 1:size(nodes,2)
            for nodeId in nodeIdi
                if nodes[2,nodeId] < SurfacePts2D[2,minNode] 
                    push!(BottomBorderNodes, nodeId)
                end
            end 
        end
    end  

    BorderPtIds = vcat(LeftborderNodes, TopBorderNodes, RightborderNodes, BottomBorderNodes)
    BorderPts = SurfacePts2D[:,BorderPtIds]
    BorderPtsSorted, BorderPtsSortedIds = sort_points(BorderPts)

    if GRAD
        ∇BorderPts = ∇SurfacePts2D[:,BorderPtIds,:]
        ∇BorderPtsSorted = ∇BorderPts[:,BorderPtsSortedIds,:]
        if SIDES
            sidePts, index = get_sides(BorderPtsSorted)

            sidePtsIds = BorderPtsSortedIds[index]
            ∇sidePts = ∇BorderPtsSorted[:,sidePtsIds,:]

            return sidePts, ∇sidePts, SurfacePts2D
        else
            return BorderPtsSorted, ∇BorderPtsSorted, SurfacePts2D, ∇SurfacePts2D
        end
    else
        if SIDES
            sidePts, index = get_sides(BorderPtsSorted)
    
            sidePts = BorderPtsSorted[index]
            sidePtsIds = BorderPtsSortedIds[index]
    
            return sidePts, SurfacePts2D
        else
            return BorderPtsSorted, SurfacePts2D
        end
    end 
end

function extract_borders(NodeList::Matrix{Float64}, camera_matrix::AbstractMatrix{Float64}, camera_pose::AbstractMatrix{Float64}, BorderNodesList::Vector{Int64}; GRAD::Bool=false, dqdθ::AbstractArray{Float64}=zeros(2,2,2), SIDES::Bool=false)
    
    surface_nodes = NodeList[:,BorderNodesList]  # extract the border nodes from the NodeList

    if GRAD
        ∇surface_nodes = dqdθ[:,BorderNodesList,:] 
        SurfacePts2D, ∇SurfacePts2D = back_project(surface_nodes, camera_matrix, camera_pose, ∇surface_nodes) 
    else
        SurfacePts2D = back_project(surface_nodes, camera_matrix, camera_pose) 
    end

    p = Array{Vector{Float64}}(undef,0)

    iter = 1:size(SurfacePts2D,2)
    for i in iter
        push!(p, SurfacePts2D[:,i])
    end

    hull = ch.ConvexHull(p)

    BorderPtIds = ch.indices(hull)
    BorderPts = SurfacePts2D[:,BorderPtIds]
    BorderPtsSorted, BorderPtsSortedIds = sort_points(BorderPts)

    if GRAD
        ∇BorderPts = ∇SurfacePts2D[:,BorderPtIds,:]
        ∇BorderPtsSorted = ∇BorderPts[:,BorderPtsSortedIds,:]
        if SIDES
            sidePts, index = get_sides(BorderPtsSorted)

            sidePtsIds = BorderPtsSortedIds[index]
            ∇sidePts = ∇BorderPtsSorted[:,sidePtsIds,:]

            return sidePts, ∇sidePts, SurfacePts2D
        else
            return BorderPtsSorted, ∇BorderPtsSorted, SurfacePts2D, ∇SurfacePts2D
        end
    else
        if SIDES
            sidePts, index = get_sides(BorderPtsSorted)
    
            sidePts = BorderPtsSorted[index]
            sidePtsIds = BorderPtsSortedIds[index]
    
            return sidePts, SurfacePts2D
        else
            return BorderPtsSorted, SurfacePts2D
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
function back_project(x::AbstractMatrix{Float64}, camera_matrix::AbstractMatrix{Float64}, camera_pose::AbstractMatrix{Float64}, dxdθ::AbstractArray{Float64})
    
    nx = size(x,2)
    nθ = size(dxdθ,3)
    xNorm = zeros(3,nx)
    xProj = zeros(3,nx)
    dudθ = zeros(3,nx,nθ)

    # transform point cloud wrt to camera frame 
    R = [1 0 0; 0 0 1; 0 -1 0]     # rotation matrix

    xTrans = R*x .+ camera_pose

    iter = 1:nx
    iterθ = 1:nθ
    for j in iter
        xNorm[1,j] = xTrans[1,j]/xTrans[3,j]
        xNorm[2,j] = xTrans[2,j]/xTrans[3,j]
        xNorm[3,j] = xTrans[3,j]/xTrans[3,j]

        xProj[:,j] = camera_matrix*xNorm[:,j]

        dπdx = ∇π(xTrans[:,j],camera_matrix)

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

function back_project(x::AbstractMatrix{Float64}, camera_matrix::AbstractMatrix{Float64}, camera_pose::AbstractMatrix{Float64})
                
    # transform point cloud wrt to camera frame 
    R = [1 0 0; 0 0 1; 0 -1 0]     # rotation matrix

    xTrans = R*x .+ camera_pose
    
    xProj = project_to(xTrans, camera_matrix)

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

"""
    noramlize(q, IEN)

Function normalize the solution vector for plotting
    
# Arguments:
- `q`: solution vector
- `IEN::Matrix{Float64}{nElem, nNodes}`: IEN array

# Returns:
- `qList`: normalized list of solutions 
"""
function noramlize(q, IEN)

    qList = zeros(size(IEN))
    max = maximum(q)
    min = minimum(q)
    iter = 1:size(IEN,1)
    for e in iter
        for n in 1:4
            qList[e,n] = (q[IEN[e,n]] - min) / (max - min)
        end
    end
    return qList
end

"""
    truncate_colormap(minval=0.0, maxval=1.0, n=100)

Function to truncate a colormap
    
# Arguments:
- `minval::Integer`: minimum value of the colormap
- `maxval::Integer`: maximum value of the colormap
- `n::Integer`: number of colors

# Returns:
- `new_cmap`: truncated colormap
"""
function truncate_colormap(minval=0.0, maxval=1.0, n=100)
    new_cmap = matplotlib.colors.LinearSegmentedColormap.from_list("mycmap", get_cmap("jet")(collect(range(maxval, minval, n))))
    return new_cmap
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

function plot_covariance(η_list::Vector{Float64}, β_list::Vector{Float64}, filepath::String)
    set_file(filepath)

    mean_η = mean(η_list)
    mean_β = mean(β_list)

    cov_η = cov(η_list)
    cov_β = cov(β_list)
    cov_ηβ = cov(η_list, β_list)

    cov_mat = [cov_η cov_ηβ; cov_ηβ cov_β]
    mean_vec = [mean_η; mean_β]

    # Plot the covariance matrix
    StatsPlots.covellipse(mean_vec, cov_mat, label="Covariance", color=:red, alpha=0.5, linewidth=2)
    # scatter!(η_list, β_list, label="Data", color=:blue, alpha=0.5)
    xlabel!("η")
    ylabel!("β")    
    title!("Covariance")
    Plots.savefig(string(filepath,"covariance.svg"))
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
                        plot_u[:,(cnte-1)*size(IEN_u_cached,1)+j] = Re_u'*sol_u[:,IEN_u_cached[:,e]]'
                        
                        # println(size(plot_u))
                    end
                end
            end
        end
    end
    return NodeList_, IEN_list, plot_u
end