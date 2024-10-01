using LinearAlgebra
using DataInterpolations
using ElasticArrays
using ConvexHulls2d
import ConvexHulls2d as ch


function greet_pp()
    println("Hello, I am the PostProcess module")
end

"""
    inflate_sphere(NodeList, x0, x1, y0, y1)

Inflate the sphere to a unit sphere 
    
# Arguments:
- `NodeList::Matrix` : {[ndims,nNodes], Matrix{Float64}} : The coordinates of the nodes
- `x0` : Float64 : The lower bound of the x direction
- `x1` : Float64 : The upper bound of the x direction
- `y0` : Float64 : The lower bound of the y direction
- `y1` : Float64 : The upper bound of the y direction
    
# Returns:
- `NodeList::Matrix{Float64}{ndims,nNodes} `:  The coordinates of the nodes after inflation
"""
function inflate_sphere(NodeList, x0, x1, y0, y1)

    x_center = [0.5*(x0 + x1), 0.5*(y0 + y1)]

    iter = 1:size(NodeList,2)
    for i in iter
        scale = maximum(abs.(NodeList[1:2,i] - x_center))
        if scale ≈ 0.
            NodeList[1:2,i] = [0 , 0]
        else
            r = sqrt((NodeList[1,i] - x_center[1])^2 + (NodeList[2,i] - x_center[2])^2)
            NodeList[1:2,i] = scale*(NodeList[1:2,i] - x_center)/r
        end
    end
    return NodeList
end

""" 
    Extract_borders(NodeList, CameraMatrix, BorderNodesList, state, ne = nothing)

Project the 3D mesh to 2D image plane and extract the border nodes (left and right)
    
# Arguments:
- `NodeList::Matrix{Float64}{ndim,nNodes}` : coordinates of the nodes
- `CameraMatrix::Matrix{Float64}{3,3}` : Camera matrix
- `BorderNodesList::Vector{Vector{Any}{4,N}`:  : List of border nodes
- `state::String` : State of the function (init:During the initialization of the mesh or update: when the mesh is updated)
- `ne::Integer`: Number of elements in each direction
# Returns:
- `NodeList::Matrix{Float64}{ndim,nbNodes}`: 2D coordinates of the border nodes
- `BorderNodes::Vector{Int}`: Indexes of the border nodes
"""
function extract_borders(NodeList, CameraMatrix, BorderNodesList, state, ne = nothing)

    SideNodes = NodeList[:,BorderNodesList[1]]  # extract the border nodes 

    SideNodes2D = back_project(SideNodes, CameraMatrix) 

    # project the nodes to the image plane and extract the border nodes as an ordered list
    if state == "init"
        @assert !isnothing(ne) "Number of elements must be provided"

        LeftborderPts = zeros(2,(ne+1))                  # vector to store indexes of the border nodes
        RightborderPts = zeros(2,(ne+1))                 # vector to store indexes of the border nodes
        LeftborderNodes = Vector{Int64}(undef, 0)        # vector to store indexes of the border nodes
        RightborderNodes = Vector{Int64}(undef, 0)   
        TopLayerList = []                                # vector to store indexes of the border nodes
        BottomLayerList = []                             # vector to store indexes of the border nodes
        szSide = size(SideNodes2D,2)÷(ne+1)                            # size of each layer
        BorderNodes = Vector{Int64}(undef, 0)            # vector to store indexes of the border nodes
        for Layers in 1:ne+1                                        # loop through each layer
            nodes = SideNodes2D[:,(Layers-1)*szSide+1:Layers*szSide]
            minNode = (Layers-1)*szSide + argmin(nodes[1,:])
            maxNode = (Layers-1)*szSide + argmax(nodes[1,:])
            push!(LeftborderNodes, minNode)
            push!(RightborderNodes, maxNode)
            LeftborderPts[:,Layers] = SideNodes2D[:,minNode]         # left border nodes
            RightborderPts[:,Layers] = SideNodes2D[:,maxNode]        # right border nodes
            if Layers == ne + 1
                nodeIdi = 1:size(nodes,2)
                for nodeId in nodeIdi
                    if nodes[2,nodeId] > SideNodes2D[2,minNode]    
                        push!(TopLayerList, (Layers-1)*szSide+nodeId)
                    end
                end 
            elseif Layers == 1
                nodeIdi = 1:size(nodes,2)
                for nodeId in nodeIdi
                    if nodes[2,nodeId] < SideNodes2D[2,minNode] 
                        push!(BottomLayerList, nodeId)
                    end
                end 
            end
        end  
        TopLayer = sortslices(SideNodes2D[:,TopLayerList],dims=2)                     # top layer nodes
        BottomLayer = sortslices(SideNodes2D[:,BottomLayerList],dims=2)               # bottom layer nodes
        
        # topNodeList = sortperm(SideNodes2D[2,TopLayerList])
        # TopLayer = SideNodes2D[:,topNodeList]
        # bottomNodeList = sortperm(SideNodes2D[2,BottomLayerList])
        # BottomLayer = SideNodes2D[:,bottomNodeList]

        BorderPoints = hcat(LeftborderPts, TopLayer, reverse(RightborderPts,dims=2), reverse(BottomLayer,dims=2))         # concatenate the left and right border nodes
        BorderNodes = vcat(LeftborderNodes, TopLayerList, reverse(RightborderNodes), reverse(BottomLayerList))            # concatenate the left and right border nodes 
    elseif state == "update"
        p = Array{Vector{Float64}}(undef,0)

        iter = 1:size(SideNodes2D,2)
        for i in iter
            push!(p, SideNodes2D[:,i])
        end

        hull = ch.ConvexHull(p)
        points = ch.vertices(hull)
        sz = length(points)
        BorderPoints = zeros(2,sz)

        iter = 1:sz
        for i in iter
            BorderPoints[:,i] = points[i]
        end

        BorderNodes = ch.vertices(hull)
    end

    return BorderPoints, BorderNodes, SideNodes2D# return the border nodes
end

""" 
    back_project(NodeList, CameraMatrix)

Project the 3D mesh to 2D image plane
    
# Arguments:
- `NodeList::Matrix{Float64}{3,nbNodes}`: 3D mesh grid
- `CameraMatrix::Matrix{Float64}{3,3}`: Camera matrix

# Returns:
- `NodeList2D::Matrix{Float64}{2,nbNodes}`: 2D coordinates of the nodes
"""
function back_project(NodeList, CameraMatrix)
                
    # transform point cloud wrt to camera frame 
    R = [1 0 0; 0 0 1; 0 -1 0]     # rotation matrix
    t = [0; -0.5; 4]               # translation vector

    NodeListTrans = R*NodeList .+ t
    
    NodeListNorm = zeros(3,size(NodeListTrans,2))

    iter = 1:size(NodeListNorm,2)
    for i in iter
        NodeListNorm[1,i] = NodeListTrans[1,i]/NodeListTrans[3,i]
        NodeListNorm[2,i] = NodeListTrans[2,i]/NodeListTrans[3,i]
        NodeListNorm[3,i] = NodeListTrans[3,i]/NodeListTrans[3,i]
    end
    
    NodeListProj = CameraMatrix'*NodeListNorm   # project to image plane

    NodeList2D = NodeListProj[1:2,:]            # extract x and y coordinates

    return NodeList2D
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

        p = CubicSpline(x,y,extrapolate=true)

        pi = [p(i) for i in samples]
        
        return pi

    elseif isnothing(borderx) && isnothing(bordery)
        x = push!(border[1,:], border[1,1])
        y = push!(border[2,:], border[2,1])

        len = length(x)
        seq = 1:(len)
        
        p = CubicSpline(x,seq)
        q = CubicSpline(y,seq)
        pi = [p(i) for i in 1:0.5:len]
        qi = [q(i) for i in 1:0.5:len]

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

function average_pts(border, centerx)
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