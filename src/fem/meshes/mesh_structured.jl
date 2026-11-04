"""
    _meshgrid_line(l,ne;element_shape=:Line,basis_order=1)

Set up the mesh grid for a 1D line

# Arguments:
- `l::Number` : length of the line
- `ne::Int` : number of elements
- `element_shape::Symbol` : element shape (`:Line`)
- `basis_order::Int` : polynomial order (1 or 2)

# Returns:
- `NodeList::Matrix{Number}{nNodes,ndim}` : array of nodes
- `IEN::Matrix{Int64}{2^ndim,ne^ndim}` : array of elements
"""
function _meshgrid_line(l::T, ne::Int64; element_shape::Symbol=:Line, basis_order::Int=1) where {T<:Number}
    BorderNodes = Int64[]
    if basis_order == 1
        nNodes = ne+1 # number of nodes in each direction
        NodeList = zeros(Float64, 1,nNodes)
        IEN = zeros(Int64,2,ne) # IEN for the 1D mesh

        x = collect(Float64, range(-l/2, l/2, length=nNodes))

        m = 1
        for i in 1:nNodes # x direction
            NodeList[1,m] = x[i]
            if i == 1 || i == nNodes # the left and right boundaries
                push!(BorderNodes,m)
            end
            m = m + 1
        end

        n = 1
        for i in 1:ne # x direction
            IEN[1,n] = i
            IEN[2,n] = i + 1
            n = n + 1
        end
    elseif basis_order == 2
        nNodes = 2*ne+1
        NodeList = zeros(Float64,1,nNodes)
        IEN = zeros(Int64,3,ne) # IEN for the 1D mesh

        x = collect(Float64, range(-l/2, l/2, length=nNodes))

        m = 1
        for i in 1:nNodes # x direction
            NodeList[1,m] = x[i]
            if i == 1 || i == nNodes # the left and right boundaries
                push!(BorderNodes,m)
            end
            m = m + 1
        end

        n = 1
        for i in 1:ne # x direction
            IEN[1,n] = 2*i - 1
            IEN[2,n] = 2*i + 1
            IEN[3,n] = 2*i
            n = n + 1
        end
    else
        throw(ArgumentError("Unknown basis_order $basis_order for element_shape $element_shape"))
    end
    return NodeList, IEN, BorderNodes
end

"""
    _meshgrid_square(lx,ly,ne;element_shape=:Quad,basis_order=1)

Set up the mesh grid for a 2D square

# Arguments:
- `x0::Number` : x-coordinate of the lower left corner of the domain
- `x1::Number` : x-coordinate of the upper right corner of the domain
- `y0::Number` : y-coordinate of the lower left corner of the domain
- `y1::Number` : y-coordinate of the upper right corner of the domain
- `ne::Int` : number of elements
- `ndim::Int` : number of dimensions
- `element_shape::Symbol` : element shape (`:Quad`)
- `basis_order::Int` : polynomial order (1 or 2)

# Returns:
- `NodeList::Matrix{Number}{nNodes,ndim}` : array of nodes
- `IEN::Matrix{Int64}{2^ndim,ne^ndim}` : array of elements
- `ID::Matrix{Int64}{nNodes,ndim}` : array of node IDs
- `IEN_top::Matrix{Int64}{2^(ndim-1),ne^(ndim-1)}` : array of elements on the top surface
- `IEN_bottom::Matrix{Int64}{2^(ndim-1),ne^(ndim-1)}` : array of elements on the bottom surface
- `BorderNodes::Vector{Int64}` : array of nodes on the boundaries
- `BottomBorderNodes::Vector{Int64}` : array of nodes on the bottom boundary
- `TopBorderNodes::Vector{Int64}` : array of nodes on the top boundary
"""
function _meshgrid_square(lx::X, ly::Y, ne::Int64; element_shape::Symbol=:Quad,
                          basis_order::Int=1) where {X<:Number,Y<:Number}

    BorderNodes = Int64[]
    BottomBorderNodes = Int64[]
    TopBorderNodes = Int64[]
    ndim = 2
    if basis_order == 1
        nNodes = ne+1 # number of nodes in each direction
        NodeList = zeros(Float64,ndim,(nNodes)^ndim)
        IEN = zeros(Int64,2^ndim,ne^ndim)              # IEN for the 3D mesh
        IEN_top = zeros(Int64,2^(ndim-1),ne^(ndim-1))  # IEN for the top surface
        IEN_bottom = zeros(Int64,2^(ndim-1),ne^(ndim-1))  # IEN for the bottom surface
        ID = zeros(Int64,ndim,(nNodes)^ndim)

        x = collect(Float64, range(-lx/2, lx/2, length=nNodes))
        y = collect(Float64, range(-ly/2, ly/2, length=nNodes))

        m = 1
        for j in 1:nNodes      # y direction
            for i in 1:nNodes  # x direction
                NodeList[1,m] = x[i]
                NodeList[2,m] = y[j]
                for l in 1:ndim
                    ID[l,m] = ndim*(m-1) + l
                end
                if i == 1 || i == nNodes # the left and right boundaries
                    push!(BorderNodes,m)
                end
                m = m + 1
            end
        end

        n = 1
        for j in 1:ne # y direction
            for i in 1:ne # x direction
                IEN[1,n] = (j-1)*(nNodes) + i
                IEN[2,n] = (j-1)*(nNodes) + i + 1
                IEN[3,n] = j*(nNodes) + i + 1
                IEN[4,n] = j*(nNodes) + i
                if j == 1 # populate the IEN for the bottom surface
                    IEN_bottom[1,i] = IEN[1,n]
                    IEN_bottom[2,i] = IEN[2,n]
                elseif j == ne # populate the IEN for the top surface
                    IEN_top[1,i] = IEN[4,n]
                    IEN_top[2,i] = IEN[3,n]
                end
                n = n + 1
            end
        end
    elseif basis_order == 2
        nNodes = 2*ne+1
        NodeList = zeros(Float64,ndim,(nNodes)^ndim)
        IEN = zeros(Int64,3^ndim,ne^ndim) # IEN for the 3D mesh
        IEN_top = zeros(Int64,3^(ndim-1),ne^(ndim-1)) # IEN for the top surface
        IEN_bottom = zeros(Int64,3^(ndim-1),ne^(ndim-1)) # IEN for the bottom surface
        ID = zeros(Int64,ndim,(nNodes)^ndim)

        x = collect(Float64, range(-lx/2, lx/2, length=nNodes))
        y = collect(Float64, range(-ly/2, ly/2, length=nNodes))

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
                IEN[1,n] = 2*(j-1)*(nNodes) + 2*i - 1
                IEN[2,n] = 2*(j-1)*(nNodes) + 2*i + 1
                IEN[3,n] = 2*j*(nNodes) + 2*i + 1
                IEN[4,n] = 2*j*(nNodes) + 2*i - 1
                IEN[5,n] = 2*(j-1)*(nNodes) + 2*i
                IEN[6,n] = (2*j-1)*(nNodes) + 2*i + 1
                IEN[7,n] = 2*j*(nNodes) + 2*i
                IEN[8,n] = (2*j-1)*(nNodes) + 2*i - 1
                IEN[9,n] = (2*j-1)*(nNodes) + 2*i
                n = n + 1
            end
        end
    else
        throw(ArgumentError("Unknown basis_order $basis_order for element_shape $element_shape"))
    end
    return NodeList, IEN, ID, IEN_top, IEN_bottom, [BorderNodes, BottomBorderNodes, TopBorderNodes]
end

"""
    _meshgrid_cuboid(lx, ly, lz, ne; element_shape=:Hex, basis_order=1)

Structured hex mesh of the box [-lx/2, lx/2] × [-ly/2, ly/2] × [0, lz] with `ne` elements along
each edge. Only hexes are generated: `element_shape` is not used beyond the error message.

# Arguments:
- `lx::Number`, `ly::Number`, `lz::Number` : edge lengths
- `ne::Int` : number of elements along each edge
- `element_shape::Symbol` : element shape (`:Hex`)
- `basis_order::Int` : polynomial order (1, 2, or 3)

# Returns:
- `NodeList::Matrix{Float64}` : (3, nNodes) node coordinates, x fastest
- `IEN::Matrix{Int64}` : ((p+1)^3, ne^3) element connectivity
- `ID::Matrix{Int64}` : (3, nNodes) DOF numbering
- `IEN_top::Matrix{Int64}`, `IEN_bottom::Matrix{Int64}` : ((p+1)^2, ne^2) faces on the top and bottom
- `IEN_side::Matrix{Int}` : empty
- `nNodes::Int` : nodes along each edge, `p*ne + 1`
- `[BorderNodes, BottomBorderNodes, TopBorderNodes]` : nodes on the four sides, and on the bottom
  and top but not on a side
"""
function _meshgrid_cuboid(lx::X, ly::Y, lz::Z, ne::Int64; element_shape::Symbol=:Hex,
                          basis_order::Int=1) where {X<:Number,Z<:Number,Y<:Number}
    haskey(_HEX_LOCAL_NODES, basis_order) ||
        throw(ArgumentError("Unknown basis_order $basis_order for element_shape $element_shape"))
    nNodes = basis_order*ne + 1         # nodes along each edge
    NodeList, BorderNodes, BottomBorderNodes, TopBorderNodes = _cuboid_nodes(lx, ly, lz, nNodes)
    IEN = _cuboid_connectivity(ne, basis_order)
    IEN_top, IEN_bottom = _cuboid_faces(IEN, ne, basis_order)
    ID = _dof_ids(3, nNodes^3)
    IEN_side = Matrix{Int}(undef, 0, 0)
    return NodeList, IEN, ID, IEN_top, IEN_bottom, IEN_side, nNodes,
           [BorderNodes, BottomBorderNodes, TopBorderNodes]
end

# Position of each local node of a hex of order p in its element, as steps (0..p) along x, y, z
# of the node grid, in smearFEM's local node order (Q1/Q2: that of `_Q1_3D_NODE_TRIPLES` and
# `_Q2_3D_NODE_TRIPLES` in fem.jl).
const _HEX_LOCAL_NODES = Dict(
    1 => (
        (0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0),
        (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1),
    ),
    2 => (
        (0, 0, 0), (2, 0, 0), (2, 2, 0), (0, 2, 0),
        (0, 0, 2), (2, 0, 2), (2, 2, 2), (0, 2, 2),
        (1, 0, 0), (2, 1, 0), (1, 2, 0), (0, 1, 0),
        (1, 0, 2), (2, 1, 2), (1, 2, 2), (0, 1, 2),
        (0, 0, 1), (2, 0, 1), (2, 2, 1), (0, 2, 1),
        (1, 0, 1), (2, 1, 1), (1, 2, 1), (0, 1, 1),
        (1, 1, 0), (1, 1, 2), (1, 1, 1),
    ),
    3 => (
        (0, 0, 0), (3, 0, 0), (3, 3, 0), (0, 3, 0),
        (0, 0, 3), (3, 0, 3), (3, 3, 3), (0, 3, 3),
        (0, 1, 0), (3, 1, 0), (3, 2, 0), (0, 2, 0),
        (0, 1, 3), (3, 1, 3), (3, 2, 3), (0, 2, 3),
        (1, 0, 0), (2, 0, 0), (2, 3, 0), (1, 3, 0),
        (1, 0, 3), (2, 0, 3), (2, 3, 3), (1, 3, 3),
        (0, 0, 1), (3, 0, 1), (3, 3, 1), (0, 3, 1),
        (0, 0, 2), (3, 0, 2), (3, 3, 2), (0, 3, 2),
        (0, 1, 1), (3, 1, 1), (3, 2, 1), (0, 2, 1),
        (0, 1, 2), (3, 1, 2), (3, 2, 2), (0, 2, 2),
        (1, 0, 1), (2, 0, 1), (2, 3, 1), (1, 3, 1),
        (1, 0, 2), (2, 0, 2), (2, 3, 2), (1, 3, 2),
        (1, 1, 0), (2, 1, 0), (2, 2, 0), (1, 2, 0),
        (1, 1, 3), (2, 1, 3), (2, 2, 3), (1, 2, 3),
        (1, 1, 1), (2, 1, 1), (2, 2, 1), (1, 2, 1),
        (1, 1, 2), (2, 1, 2), (2, 2, 2), (1, 2, 2),
    ),
)

# Local nodes of an element's bottom (z = 0) and top faces, in face node order. Q3 lists 33:36 and
# 37:40, which lie inside the element, where the face's own inner nodes are 49:52 and 53:56.
const _HEX_BOTTOM_FACE = Dict(1 => [1:4;], 2 => [1:4; 9:12; 25], 3 => [1:4; 9:12; 17:20; 33:36])
const _HEX_TOP_FACE = Dict(1 => [5:8;], 2 => [5:8; 13:16; 26], 3 => [5:8; 13:16; 21:24; 37:40])

# Nodes of an nNodes³ grid on [-lx/2, lx/2] × [-ly/2, ly/2] × [0, lz], x fastest, and the nodes on
# the four sides, on the bottom (not on a side) and on the top (not on a side).
function _cuboid_nodes(lx, ly, lz, nNodes::Int)
    x = collect(Float64, range(-lx/2, lx/2, length=nNodes))
    y = collect(Float64, range(-ly/2, ly/2, length=nNodes))
    z = collect(Float64, range(0, lz, length=nNodes))
    NodeList = zeros(Float64, 3, nNodes^3)
    side, bottom, top = Int64[], Int64[], Int64[]
    m = 1
    for k in 1:nNodes, j in 1:nNodes, i in 1:nNodes
        NodeList[1, m] = x[i]
        NodeList[2, m] = y[j]
        NodeList[3, m] = z[k]
        if i == 1 || i == nNodes || j == 1 || j == nNodes
            push!(side, m)
        elseif k == 1
            push!(bottom, m)
        elseif k == nNodes
            push!(top, m)
        end
        m += 1
    end
    return NodeList, side, bottom, top
end

# Connectivity of the ne³ hexes of order p on the grid of `_cuboid_nodes`, elements x fastest.
function _cuboid_connectivity(ne::Int, p::Int)
    nNodes = p*ne + 1
    local_nodes = _HEX_LOCAL_NODES[p]
    IEN = zeros(Int64, length(local_nodes), ne^3)
    n = 1
    for k in 1:ne, j in 1:ne, i in 1:ne
        for (l, (a, b, c)) in enumerate(local_nodes)
            IEN[l, n] = (p*(k-1) + c)*nNodes^2 + (p*(j-1) + b)*nNodes + p*(i-1) + a + 1
        end
        n += 1
    end
    return IEN
end

# Top faces of the last element layer and bottom faces of the first. With one layer (ne = 1) the
# top stays zero: the original generator filled only the bottom then.
function _cuboid_faces(IEN::Matrix{Int64}, ne::Int, p::Int)
    nfaces = ne^2
    IEN_bottom = IEN[_HEX_BOTTOM_FACE[p], 1:nfaces]
    IEN_top = ne > 1 ? IEN[_HEX_TOP_FACE[p], end-nfaces+1:end] : zeros(Int64, (p+1)^2, nfaces)
    return IEN_top, IEN_bottom
end

"""
    _meshgrid_ring(r1,r2, theta1, theta2, ne)

Set up the mesh grid for a 2D annulus ring

# Arguments:
- `r1::Number` : inner radius
- `r2::Number` : outer radius
- `theta1::Number` : start angle
- `theta2::Number` : end angle
- `ne::Int` : number of elements

# Returns:
- `NodeList::Matrix{Number}{nNodes,ndim}` : array of nodes
- `IEN::Matrix{Int64}{ne^ndim,2^ndim}` : connectivity matrix
"""
function _meshgrid_ring(r1::T, r2::U, theta1::V, theta2::W, ne::Int64,
                        ndim::Int64) where {T<:Number,U<:Number,V<:Number,W<:Number}

    r = collect(Float64, range(r1, r2, length=ne+1))
    theta = collect(Float64, range(theta1, theta2, length=ne+1))

    NodeList = zeros(Float64,2,(ne+1)*(ne+1))

    k = 1
    for i in 1:ne+1
        for j in 1:ne+1
            NodeList[1,k] = r[i]*sin(theta[j])
            NodeList[2,k] = r[i]*cos(theta[j])
            k = k + 1
        end
    end

    IEN = zeros(Int64,4,ne^ndim)

    l = 1
    for i in 1:ne
        for j in 1:ne
            IEN[1,l] = (i-1)*(ne+1) + j
            IEN[2,l] = (i-1)*(ne+1) + j + 1
            IEN[3,l] = i*(ne+1) + j + 1
            IEN[4,l] = i*(ne+1) + j
            l = l + 1
        end
    end

    return NodeList, IEN
end

"""
    _meshgrid_cylinder(r,h,ne;element_shape=:Hex,basis_order=1)
Set up the mesh grid for a 3D cylinder
# Arguments:
- `r::Number` : radius of the cylinder
- `h::Number` : height of the cylinder
- `ne::Int` : number of elements
- `element_shape::Symbol` : element shape (`:Hex` or `:Tet`)
- `basis_order::Int` : polynomial order (1 or 2)
# Returns:
- `NodeList::Matrix{Number}{nNodes,ndim}` : array of nodes
- `IEN::Matrix{Int64}{ne^ndim,2^ndim}` : array of elements
- `ID::Matrix{Int64}{nNodes,ndim}` : array of node IDs
- `IEN_top::Matrix{Int64}{ne^(ndim-1),2^(ndim-1)}` : array of elements on the top surface
- `IEN_bottom::Matrix{Int64}{ne^(ndim-1),2^(ndim-1)}` : array of elements on the bottom surface
- `BorderNodes::Vector{Int64}` : array of nodes on the boundaries

"""
function _meshgrid_cylinder(r::T, h::U, ne::Int64; element_shape::Symbol=:Hex,
                            basis_order::Int=1) where {T<:Number,U<:Number}
    @info "Representing fields / geometry using $element_shape order-$basis_order basis functions"
    NodeList, IEN, ID, IEN_top, IEN_bottom, _, _, BorderNodes =
        _meshgrid_cuboid(1, 1, 1, ne; element_shape=element_shape, basis_order=basis_order)
    NodeList = _inflate_cylinder(NodeList, -0.5, 0.5, -0.5, 0.5, r, h)

    mesh = MeshgridCylinder(r=r, h=h, NodeList=NodeList, IEN=IEN, IEN_top=IEN_top, IEN_bottom=IEN_bottom,
        ID=ID, volume_element_shape=element_shape, basis_order=basis_order,
        nNodes=size(NodeList, 2), ne=ne^3,
        side_nodes=BorderNodes[1], top_nodes=BorderNodes[3], bottom_nodes=BorderNodes[2])
    return mesh
end

function _inflate_cylinder(NodeListCube::Matrix{Float64}, x0::T, x1::U, y0::V, y1::W, r::Y, h::Z=1.0;
                          GRAD::Bool=false) where {T<:Number,U<:Number,V<:Number,W<:Number,
                                                   Y<:Number,Z<:Number}
    NodeListCyl = copy(NodeListCube)
    center = [0.5*(x0 + x1), 0.5*(y0 + y1)]
    iter = 1:size(NodeListCyl,2)

    if GRAD == true
        ∇NodeListCyl = zeros(Float64,size(NodeListCube,1),size(NodeListCube,2),(size(r,1)+1))
        for i in iter
            scale = maximum(abs.(NodeListCube[1:2,i] - center))
            if scale ≈ 0.
                NodeListCyl[1:2,i] = [0 , 0]
                NodeListCyl[3,i] = NodeListCube[3,i]*h
            else
                r_ = sqrt((NodeListCube[1,i] - center[1])^2 + (NodeListCube[2,i] - center[2])^2)
                NodeListCyl[1:2,i] = 2*r*scale*(NodeListCube[1:2,i] - center)/r_
                NodeListCyl[3,i] = NodeListCube[3,i]*h

                # ∂x/∂r = K_1 and ∂y/∂r = K_2
                ∇NodeListCyl[1:2,i,1] = 2*scale*(NodeListCube[1:2,i] - center)/r_
                ∇NodeListCyl[3,i,2] = NodeListCube[3,i] # ∂z/∂h = λ_i
            end
        end
        return NodeListCyl, ∇NodeListCyl
    else
        for i in iter
            scale = maximum(abs.(NodeListCube[1:2,i] - center))
            if scale ≈ 0.
                NodeListCyl[1:2,i] = [0 , 0]
                NodeListCyl[3,i] = NodeListCube[3,i]*h
            else
                r_ = sqrt((NodeListCube[1,i] - center[1])^2 + (NodeListCube[2,i] - center[2])^2)
                NodeListCyl[1:2,i] = 2*r*scale*(NodeListCube[1:2,i] - center)/r_
                NodeListCyl[3,i] = NodeListCube[3,i]*h
            end
        end
        return NodeListCyl
    end
end
