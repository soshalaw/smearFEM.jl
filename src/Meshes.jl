"""
    meshgrid_square(x0,x1,y0,y1,ne,ndim;FunctionClass="Q1")

Set up the mesh grid for a 2D square

# Arguments:
- `x0::Float64` : x-coordinate of the lower left corner of the domain
- `x1::Float64` : x-coordinate of the upper right corner of the domain
- `y0::Float64` : y-coordinate of the lower left corner of the domain
- `y1::Float64` : y-coordinate of the upper right corner of the domain
- `ne::Int` : number of elements
- `ndim::Int` : number of dimensions
- `FunctionClass::String` : type of basis function

# Returns:
- `NodeList::Matrix{Float64}{nNodes,ndim}` : array of nodes
- `IEN::Matrix{Int64}{ne^ndim,2^ndim}` : array of elements
- `ID::Matrix{Int64}{nNodes,ndim}` : array of node IDs
- `IEN_top::Matrix{Int64}{ne^(ndim-1),2^(ndim-1)}` : array of elements on the top surface
- `IEN_btm::Matrix{Int64}{ne^(ndim-1),2^(ndim-1)}` : array of elements on the bottom surface
- `BorderNodes::Vector{Int64}` : array of nodes on the boundaries
- `BottomBorderNodes::Vector{Int64}` : array of nodes on the bottom boundary
- `TopBorderNodes::Vector{Int64}` : array of nodes on the top boundary
"""
function meshgrid_square(x0,x1,y0,y1,ne,ndim;FunctionClass="Q1")
            
    BorderNodes = []
    BottomBorderNodes = []
    TopBorderNodes = []

    if FunctionClass == "Q1"
        nNodes = ne+1 # number of nodes in each direction
        NodeList = zeros(ndim,(nNodes)^ndim)
        IEN = zeros(Int64,ne^ndim,2^ndim) # IEN for the 3D mesh
        IEN_top = zeros(Int64,ne^(ndim-1),2^(ndim-1)) # IEN for the top surface
        IEN_btm = zeros(Int64,ne^(ndim-1),2^(ndim-1)) # IEN for the bottom surface
        IEN_side = zeros(Int64,ne^(ndim-1),2^(ndim-1)) # IEN for the side surfaces
        ID = zeros(Int64,(nNodes)^ndim,ndim)

        x = collect(range(x0, x1, length=nNodes))
        y = collect(range(y0, y1, length=nNodes))
    
        m = 1
        for j in 1:nNodes # y direction
            for i in 1:nNodes # x direction
                NodeList[1,m] = x[i]
                NodeList[2,m] = y[j]
                for l in 1:ndim
                    ID[m,l] = ndim*(m-1) + l
                end
                if i == 1 || i == nNodes # populate the BorderNodes with the nodes on the left and right boundaries
                    push!(BorderNodes,m)
                end
                m = m + 1
            end
        end 
        
        n = 1
        for j in 1:ne # y direction
            for i in 1:ne # x direction
                IEN[n,1] = (j-1)*(nNodes) + i
                IEN[n,2] = (j-1)*(nNodes) + i + 1
                IEN[n,3] = j*(nNodes) + i + 1
                IEN[n,4] = j*(nNodes) + i
                if j == 1 # populate the IEN for the bottom surface
                    IEN_btm[i,1] = IEN[n,1]
                    IEN_btm[i,2] = IEN[n,2]
                elseif j == ne # populate the IEN for the top surface
                    IEN_top[i,1] = IEN[n,4]
                    IEN_top[i,2] = IEN[n,3]
                end
                n = n + 1
            end
        end
    elseif FunctionClass == "Q2"
        nNodes = 2*ne+1
        NodeList = zeros(ndim,(nNodes)^ndim)
        IEN = zeros(Int64,ne^ndim,3^ndim) # IEN for the 3D mesh
        IEN_top = zeros(Int64,ne^(ndim-1),3^(ndim-1)) # IEN for the top surface
        IEN_btm = zeros(Int64,ne^(ndim-1),3^(ndim-1)) # IEN for the bottom surface
        ID = zeros(Int64,(nNodes)^ndim,ndim)

        x = collect(range(x0, x1, length=nNodes))
        y = collect(range(y0, y1, length=nNodes))
        
        m = 1
        for j in 1:nNodes # y direction
            for i in 1:nNodes # x direction 
                NodeList[1,m] = x[i]
                NodeList[2,m] = y[j]
                m = m + 1
            end
        end

        n = 1
        for j in 1:ne # y direction
            for i in 1:ne # x direction
                IEN[n,1] = 2*(j-1)*(nNodes) + 2*i - 1
                IEN[n,2] = 2*(j-1)*(nNodes) + 2*i + 1
                IEN[n,3] = 2*j*(nNodes) + 2*i + 1
                IEN[n,4] = 2*j*(nNodes) + 2*i - 1
                IEN[n,5] = 2*(j-1)*(nNodes) + 2*i
                IEN[n,6] = (2*j-1)*(nNodes) + 2*i + 1
                IEN[n,7] = 2*j*(nNodes) + 2*i
                IEN[n,8] = (2*j-1)*(nNodes) + 2*i - 1
                IEN[n,9] = (2*j-1)*(nNodes) + 2*i
                n = n + 1
            end
        end 
    end
    return NodeList, IEN, ID, IEN_top, IEN_btm, [BorderNodes, BottomBorderNodes, TopBorderNodes]
end

"""
    meshgrid_cube(x0,x1,y0,y1,z0,z1,ne,ndim;FunctionClass="Q1")

Set up the mesh grid for a 3D cube

# Arguments:
- `x0::Float64` : x-coordinate of the lower left corner of the domain
- `x1::Float64` : x-coordinate of the upper right corner of the domain
- `y0::Float64` : y-coordinate of the lower left corner of the domain
- `y1::Float64` : y-coordinate of the upper right corner of the domain
- `z0::Float64` : z-coordinate of the lower left corner of the domain
- `z1::Float64` : z-coordinate of the upper right corner of the domain
- `ne::Int` : number of elements
- `ndim::Int` : number of dimensions
- `FunctionClass::String` : type of basis function

# Returns:
- `NodeList::Matrix{Float64}{nNodes,ndim}` : array of nodes
- `IEN::Matrix{Int64}{ne^ndim,2^ndim}` : array of elements
- `ID::Matrix{Int64}{nNodes,ndim}` : array of node IDs
- `IEN_top::Matrix{Int64}{ne^(ndim-1),2^(ndim-1)}` : array of elements on the top surface
- `IEN_btm::Matrix{Int64}{ne^(ndim-1),2^(ndim-1)}` : array of elements on the bottom surface
- `BorderNodes::Vector{Int64}` : array of nodes on the boundaries
- `BottomBorderNodes::Vector{Int64}` : array of nodes on the bottom boundary
- `TopBorderNodes::Vector{Int64}` : array of nodes on the top boundary
"""
function meshgrid_cube(x0,x1,y0,y1,z0,z1,ne,ndim;FunctionClass=FunctionClass)
    BorderNodes = []
    BottomBorderNodes = []
    TopBorderNodes = []

    if FunctionClass == "Q1"
        nNodes = ne+1 # number of nodes in each direction
        NodeList = zeros(ndim,(nNodes)^ndim)
        IEN = zeros(Int64,ne^ndim,2^ndim) # IEN for the 3D mesh
        IEN_top = zeros(Int64,ne^(ndim-1),2^(ndim-1)) # IEN for the top surface
        IEN_btm = zeros(Int64,ne^(ndim-1),2^(ndim-1)) # IEN for the bottom surface
        IEN_side = zeros(Int64,ne^(ndim-1),2^(ndim-1)) # IEN for the side surfaces
        ID = zeros(Int64,(nNodes)^ndim,ndim)

        x = collect(range(x0, x1, length=nNodes))
        y = collect(range(y0, y1, length=nNodes))
        z = collect(range(z0, z1, length=nNodes))
        
        m = 1
        for k in 1:nNodes # z direction
            for j in 1:nNodes # y direction
                for i in 1:nNodes # x direction
                    NodeList[1,m] = x[i]
                    NodeList[2,m] = y[j]
                    NodeList[3,m] = z[k]
                    for l in 1:ndim
                        ID[m,l] = ndim*(m-1) + l
                    end
                    if (i == 1 || i == nNodes || j == 1 || j == nNodes)
                        push!(BorderNodes,m) # populate the BorderNodes with the nodes on the boundaries (excluding the top and bottom surfaces)
                    elseif k == 1
                        push!(BottomBorderNodes,m) # populate the BorderNodes with the nodes on the top and bottom surfaces
                    elseif k == nNodes
                        push!(TopBorderNodes,m) # populate the BorderNodes with the nodes on the top and bottom surfaces
                    end
                    m = m + 1
                end
            end
        end
        
        n = 1           # element number on 3D mesh
        nt = 1          # element number on 2D mesh of top surface
        nb = 1          # element number on 2D mesh of bottom surface
        for k in 1:ne            # z direction
            for j in 1:ne        # y direction
                for i in 1:ne    # x direction
                    IEN[n,1] = (k-1)*(nNodes)^2 + (j-1)*(nNodes) + i
                    IEN[n,2] = (k-1)*(nNodes)^2 + (j-1)*(nNodes) + i + 1
                    IEN[n,3] = (k-1)*(nNodes)^2 + j*(nNodes) + i + 1
                    IEN[n,4] = (k-1)*(nNodes)^2 + j*(nNodes) + i
                    IEN[n,5] = k*(nNodes)^2 + (j-1)*(nNodes) + i
                    IEN[n,6] = k*(nNodes)^2 + (j-1)*(nNodes) + i + 1
                    IEN[n,7] = k*(nNodes)^2 + j*(nNodes) + i + 1
                    IEN[n,8] = k*(nNodes)^2 + j*(nNodes) + i
                    if k == 1 # populate the IEN for the bottom surface
                        IEN_btm[nb,1] = IEN[n,1]
                        IEN_btm[nb,2] = IEN[n,2]
                        IEN_btm[nb,3] = IEN[n,3]
                        IEN_btm[nb,4] = IEN[n,4]
                        nb = nb + 1
                    elseif k == ne # populate the IEN for the top surface
                        IEN_top[nt,1] = IEN[n,5]
                        IEN_top[nt,2] = IEN[n,6]
                        IEN_top[nt,3] = IEN[n,7]
                        IEN_top[nt,4] = IEN[n,8]
                        nt = nt + 1
                    # elseif j == 1 || j == ne || i == 1 || i == ne
                    #     IEN_side[n,1] = IEN[n,1]
                    #     IEN_side[n,2] = IEN[n,2]
                    #     IEN_side[n,3] = IEN[n,3]
                    #     IEN_side[n,4] = IEN[n,4]
                    #     IEN_side[n,5] = IEN[n,5]
                    #     IEN_side[n,6] = IEN[n,6]
                    #     IEN_side[n,7] = IEN[n,7]
                    #     IEN_side[n,8] = IEN[n,8]
                    end
                    n = n + 1
                end
            end
        end
    elseif FunctionClass == "Q2"
        nNodes = 2*ne+1
        NodeList = zeros(ndim,(nNodes)^ndim)
        IEN = zeros(Int64,ne^ndim,3^ndim) # IEN for the 3D mesh
        IEN_top = zeros(Int64,ne^(ndim-1),3^(ndim-1)) # IEN for the top surface
        IEN_btm = zeros(Int64,ne^(ndim-1),3^(ndim-1)) # IEN for the bottom surface
        ID = zeros(Int64,(nNodes)^ndim,ndim)

        x = collect(range(x0, x1, length=nNodes))
        y = collect(range(y0, y1, length=nNodes))
        z = collect(range(z0, z1, length=nNodes))

        m = 1
        for k in 1:nNodes
            for j in 1:nNodes
                for i in 1:nNodes
                    NodeList[1,m] = x[i]
                    NodeList[2,m] = y[j]
                    NodeList[3,m] = z[k]
                    for l in 1:ndim
                        ID[m,l] = ndim*(m-1) + l
                    end
                    if (i == 1 || i == nNodes || j == 1 || j == nNodes)
                        push!(BorderNodes,m) # populate the BorderNodes with the nodes on the boundaries (excluding the top and bottom surfaces)
                    elseif k == 1
                        push!(BottomBorderNodes,m) # populate the BorderNodes with the nodes on the top and bottom surfaces
                    elseif k == nNodes
                        push!(TopBorderNodes,m) # populate the BorderNodes with the nodes on the top and bottom surfaces
                    end
                    m = m + 1
                end
            end
        end
        
        n = 1       # element number on 3D mesh
        nt = 1      # element number on 2D mesh of top surface
        nb = 1      # element number on 2D mesh of bottom surface
        for k in 1:ne
            for j in 1:ne
                for i in 1:ne
                    IEN[n,1] = 2*(k-1)*(nNodes)^2 + 2*(j-1)*(nNodes) + 2*i - 1
                    IEN[n,2] = 2*(k-1)*(nNodes)^2 + 2*(j-1)*(nNodes) + 2*i + 1
                    IEN[n,3] = 2*(k-1)*(nNodes)^2 + 2*j*(nNodes) + 2*i + 1
                    IEN[n,4] = 2*(k-1)*(nNodes)^2 + 2*j*(nNodes) + 2*i - 1
                    IEN[n,5] = 2*k*(nNodes)^2 + 2*(j-1)*(nNodes) + 2*i -1
                    IEN[n,6] = 2*k*(nNodes)^2 + 2*(j-1)*(nNodes) + 2*i + 1
                    IEN[n,7] = 2*k*(nNodes)^2 + 2*j*(nNodes) + 2*i + 1
                    IEN[n,8] = 2*k*(nNodes)^2 + 2*j*(nNodes) + 2*i - 1
                    IEN[n,9] = 2*(k-1)*(nNodes)^2 + 2*(j-1)*(nNodes) + 2*i 
                    IEN[n,10] = 2*(k-1)*(nNodes)^2 + (2*j-1)*(nNodes) + 2*i + 1
                    IEN[n,11] = 2*(k-1)*(nNodes)^2 + 2*j*(nNodes) + 2*i
                    IEN[n,12] = 2*(k-1)*(nNodes)^2 + (2*j-1)*(nNodes) + 2*i - 1
                    IEN[n,13] = 2*k*(nNodes)^2 + 2*(j-1)*(nNodes) + 2*i
                    IEN[n,14] = 2*k*(nNodes)^2 + (2*j-1)*(nNodes) + 2*i + 1
                    IEN[n,15] = 2*k*(nNodes)^2 + 2*j*(nNodes) + 2*i
                    IEN[n,16] = 2*k*(nNodes)^2 + (2*j-1)*(nNodes) + 2*i - 1
                    IEN[n,17] = (2*k-1)*(nNodes)^2 + 2*(j-1)*(nNodes) + 2*i - 1
                    IEN[n,18] = (2*k-1)*(nNodes)^2 + 2*(j-1)*(nNodes) + 2*i + 1
                    IEN[n,19] = (2*k-1)*(nNodes)^2 + 2*j*(nNodes) + 2*i + 1
                    IEN[n,20] = (2*k-1)*(nNodes)^2 + 2*j*(nNodes) + 2*i - 1
                    IEN[n,21] = (2*k-1)*(nNodes)^2 + 2*(j-1)*(nNodes) + 2*i
                    IEN[n,22] = (2*k-1)*(nNodes)^2 + (2*j-1)*(nNodes) + 2*i + 1
                    IEN[n,23] = (2*k-1)*(nNodes)^2 + 2*j*(nNodes) + 2*i
                    IEN[n,24] = (2*k-1)*(nNodes)^2 + (2*j-1)*(nNodes) + 2*i - 1
                    IEN[n,25] = 2*(k-1)*(nNodes)^2 + (2*j-1)*(nNodes) + 2*i
                    IEN[n,26] = 2*k*(nNodes)^2 + (2*j-1)*(nNodes) + 2*i
                    IEN[n,27] = (2*k-1)*(nNodes)^2 + (2*j-1)*(nNodes) + 2*i
                    if k == 1
                        IEN_btm[nb,1] = IEN[n,1]
                        IEN_btm[nb,2] = IEN[n,2]
                        IEN_btm[nb,3] = IEN[n,3]
                        IEN_btm[nb,4] = IEN[n,4]
                        IEN_btm[nb,5] = IEN[n,9]
                        IEN_btm[nb,6] = IEN[n,10]
                        IEN_btm[nb,7] = IEN[n,11]
                        IEN_btm[nb,8] = IEN[n,12]
                        IEN_btm[nb,9] = IEN[n,25]
                        nb = nb + 1
                    elseif k == ne
                        IEN_top[nt,1] = IEN[n,5]
                        IEN_top[nt,2] = IEN[n,6]
                        IEN_top[nt,3] = IEN[n,7]
                        IEN_top[nt,4] = IEN[n,8]
                        IEN_top[nt,5] = IEN[n,13]
                        IEN_top[nt,6] = IEN[n,14]
                        IEN_top[nt,7] = IEN[n,15]
                        IEN_top[nt,8] = IEN[n,16]
                        IEN_top[nt,9] = IEN[n,26]
                        nt = nt + 1
                    end
                    n = n + 1
                end
            end
        end
    end
    return NodeList, IEN, ID, IEN_top, IEN_btm, [BorderNodes, BottomBorderNodes, TopBorderNodes]
end

"""
    meshgrid_ring(r1,r2, theta1, theta2, ne)

Set up the mesh grid for a 2D annulus ring

# Arguments:
- `r1::Float64` : inner radius
- `r2::Float64` : outer radius
- `theta1::Float64` : start angle
- `theta2::Float64` : end angle
- `ne::Int` : number of elements

# Returns:
- `NodeList::Matrix{Float64}{nNodes,ndim}` : array of nodes
- `IEN::Matrix{Int64}{ne^ndim,2^ndim}` : connectivity matrix
"""
function meshgrid_ring(r1,r2, theta1, theta2, ne)
    
    r = collect(range(r1, r2, length=ne+1))
    theta = collect(range(theta1, theta2, length=ne+1))

    NodeList = zeros(2,(ne+1)*(ne+1))

    k = 1
    for i in 1:ne+1
        for j in 1:ne+1
            NodeList[1,k] = r[i]*sin(theta[j])
            NodeList[2,k] = r[i]*cos(theta[j])
            k = k + 1
        end
    end

    IEN = zeros(Int64,ne*ne,4)
    
    l = 1
    for i in 1:ne
        for j in 1:ne
            IEN[l,1] = (i-1)*(ne+1) + j
            IEN[l,2] = (i-1)*(ne+1) + j + 1
            IEN[l,3] = i*(ne+1) + j + 1
            IEN[l,4] = i*(ne+1) + j
            l = l + 1
        end
    end

    return NodeList, IEN
end

"""
    inflate_cylinder(NodeList, x0, x1, y0, y1)

Inflate the sphere to a cylinder of unit radius and height

# Arguments:
- `NodeList::Matrix{Float64}{nNodes,ndim}` : array of nodes
- `x0::Float64` : x-coordinate of the lower left corner of the domain
- `x1::Float64` : x-coordinate of the upper right corner of the domain
- `y0::Float64` : y-coordinate of the lower left corner of the domain
- `y1::Float64` : y-coordinate of the upper right corner of the domain

# Returns:
- `NodeList::Matrix{Float64}{nNodes,ndim}` : array of nodes
"""
function inflate_cylinder(NodeList, x0, x1, y0, y1)

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