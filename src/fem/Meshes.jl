"""
    meshgrid_line(l,ne;FunctionClass="Q1")

Set up the mesh grid for a 1D line

# Arguments:
- `l::Number` : length of the line
- `ne::Int` : number of elements
- `FunctionClass::String` : type of basis function (Q1 or Q2)

# Returns:
- `NodeList::Matrix{Number}{nNodes,ndim}` : array of nodes
- `IEN::Matrix{Int64}{2^ndim,ne^ndim}` : array of elements
"""
function meshgrid_line(l::Number, ne::Int64; FunctionClass::String="Q1")
    BorderNodes = Int64[]
    if FunctionClass == "Q1"
        nNodes = ne+1 # number of nodes in each direction
        NodeList = zeros(Float64, 1,nNodes)
        IEN = zeros(Int64,2,ne) # IEN for the 1D mesh

        x = collect(Float64, range(-l/2, l/2, length=nNodes))

        m = 1
        for i in 1:nNodes # x direction
            NodeList[1,m] = x[i]
            if i == 1 || i == nNodes # populate the BorderNodes with the nodes on the left and right boundaries
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
    elseif FunctionClass == "Q2"
        nNodes = 2*ne+1
        NodeList = zeros(Float64,1,nNodes)
        IEN = zeros(Int64,3,ne) # IEN for the 1D mesh

        x = collect(Float64, range(x0, x1, length=nNodes))

        m = 1
        for i in 1:nNodes # x direction
            NodeList[1,m] = x[i]
            if i == 1 || i == nNodes # populate the BorderNodes with the nodes on the left and right boundaries
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
        throw(ArgumentError("Basis function type $FunctionClass is unknown"))
    end
    return NodeList, IEN, BorderNodes
end

"""
    meshgrid_square(lx,ly,ne,ndim;FunctionClass="Q1")

Set up the mesh grid for a 2D square

# Arguments:
- `x0::Number` : x-coordinate of the lower left corner of the domain
- `x1::Number` : x-coordinate of the upper right corner of the domain
- `y0::Number` : y-coordinate of the lower left corner of the domain
- `y1::Number` : y-coordinate of the upper right corner of the domain
- `ne::Int` : number of elements
- `ndim::Int` : number of dimensions
- `FunctionClass::String` : type of basis function

# Returns:
- `NodeList::Matrix{Number}{nNodes,ndim}` : array of nodes
- `IEN::Matrix{Int64}{2^ndim,ne^ndim}` : array of elements
- `ID::Matrix{Int64}{nNodes,ndim}` : array of node IDs
- `IEN_top::Matrix{Int64}{2^(ndim-1),ne^(ndim-1)}` : array of elements on the top surface
- `IEN_btm::Matrix{Int64}{2^(ndim-1),ne^(ndim-1)}` : array of elements on the bottom surface
- `BorderNodes::Vector{Int64}` : array of nodes on the boundaries
- `BottomBorderNodes::Vector{Int64}` : array of nodes on the bottom boundary
- `TopBorderNodes::Vector{Int64}` : array of nodes on the top boundary
"""
function meshgrid_square(lx::Number, ly::Number, ne::Int64; FunctionClass::String="Q1")
            
    BorderNodes = Int64[]
    BottomBorderNodes = Int64[]
    TopBorderNodes = Int64[]

    if FunctionClass == "Q1"
        nNodes = ne+1 # number of nodes in each direction
        NodeList = zeros(Float64,ndim,(nNodes)^ndim)
        IEN = zeros(Int64,2^ndim,ne^ndim)              # IEN for the 3D mesh
        IEN_top = zeros(Int64,2^(ndim-1),ne^(ndim-1))  # IEN for the top surface
        IEN_btm = zeros(Int64,2^(ndim-1),ne^(ndim-1))  # IEN for the bottom surface
        IEN_side = zeros(Int64,2^(ndim-1),ne^(ndim-1)) # IEN for the side surfaces
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
                if i == 1 || i == nNodes # populate the BorderNodes with the nodes on the left and right boundaries
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
                    IEN_btm[1,i] = IEN[1,n]
                    IEN_btm[2,i] = IEN[2,n]
                elseif j == ne # populate the IEN for the top surface
                    IEN_top[1,i] = IEN[4,n]
                    IEN_top[2,i] = IEN[3,n]
                end
                n = n + 1
            end
        end
    elseif FunctionClass == "Q2"
        nNodes = 2*ne+1
        NodeList = zeros(Float64,ndim,(nNodes)^ndim)
        IEN = zeros(Int64,3^ndim,ne^ndim) # IEN for the 3D mesh
        IEN_top = zeros(Int64,3^(ndim-1),ne^(ndim-1)) # IEN for the top surface
        IEN_btm = zeros(Int64,3^(ndim-1),ne^(ndim-1)) # IEN for the bottom surface
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
        throw(ArgumentError("Basis function type $FunctionClass is unknown"))
    end
    return NodeList, IEN, ID, IEN_top, IEN_btm, [BorderNodes, BottomBorderNodes, TopBorderNodes]
end

"""
    meshgrid_cube(lx,ly,lz,ne,ndim;FunctionClass="Q1")

Set up the mesh grid for a 3D cube

# Arguments:
- `x0::Number` : x-coordinate of the lower left corner of the domain
- `x1::Number` : x-coordinate of the upper right corner of the domain
- `y0::Number` : y-coordinate of the lower left corner of the domain
- `y1::Number` : y-coordinate of the upper right corner of the domain
- `z0::Number` : z-coordinate of the lower left corner of the domain
- `z1::Number` : z-coordinate of the upper right corner of the domain
- `ne::Int` : number of elements
- `ndim::Int` : number of dimensions
- `FunctionClass::String` : type of basis function

# Returns:
- `NodeList::Matrix{Number}{nNodes,ndim}` : array of nodes
- `IEN::Matrix{Int64}{ne^ndim,2^ndim}` : array of elements
- `ID::Matrix{Int64}{nNodes,ndim}` : array of node IDs
- `IEN_top::Matrix{Int64}{ne^(ndim-1),2^(ndim-1)}` : array of elements on the top surface
- `IEN_btm::Matrix{Int64}{ne^(ndim-1),2^(ndim-1)}` : array of elements on the bottom surface
- `BorderNodes::Vector{Int64}` : array of nodes on the boundaries
- `BottomBorderNodes::Vector{Int64}` : array of nodes on the bottom boundary
- `TopBorderNodes::Vector{Int64}` : array of nodes on the top boundary
"""
function meshgrid_cube(lx::Number, ly::Number, lz::Number, ne::Int64; FunctionClass::String="Q1")
    BorderNodes = Int64[]
    BottomBorderNodes = Int64[]
    TopBorderNodes = Int64[]
    ndim = 3

    if FunctionClass == "Q1"
        nNodes = ne+1 # number of nodes in each direction
        NodeList = zeros(Float64,ndim,(nNodes)^ndim)
        IEN = zeros(Int64,2^ndim,ne^ndim) # IEN for the 3D mesh
        IEN_top = zeros(Int64,2^(ndim-1),ne^(ndim-1)) # IEN for the top surface
        IEN_btm = zeros(Int64,2^(ndim-1),ne^(ndim-1)) # IEN for the bottom surface
        IEN_side = zeros(Int64,2^(ndim-1),ne^(ndim-1)) # IEN for the side surfaces
        ID = zeros(Int64,ndim,(nNodes)^ndim)

        x = collect(Float64, range(-lx/2, lx/2, length=nNodes))
        y = collect(Float64, range(-ly/2, ly/2, length=nNodes))
        z = collect(Float64, range( 0, lz, length=nNodes))
        
        m = 1
        for k in 1:nNodes # z direction
            for j in 1:nNodes # y direction
                for i in 1:nNodes # x direction
                    NodeList[1,m] = x[i]
                    NodeList[2,m] = y[j]
                    NodeList[3,m] = z[k]
                    for l in 1:ndim
                        ID[l,m] = ndim*(m-1) + l
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
                    IEN[1,n] = (k-1)*(nNodes)^2 + (j-1)*(nNodes) + i
                    IEN[2,n] = (k-1)*(nNodes)^2 + (j-1)*(nNodes) + i + 1
                    IEN[3,n] = (k-1)*(nNodes)^2 + j*(nNodes) + i + 1
                    IEN[4,n] = (k-1)*(nNodes)^2 + j*(nNodes) + i
                    IEN[5,n] = k*(nNodes)^2 + (j-1)*(nNodes) + i
                    IEN[6,n] = k*(nNodes)^2 + (j-1)*(nNodes) + i + 1
                    IEN[7,n] = k*(nNodes)^2 + j*(nNodes) + i + 1
                    IEN[8,n] = k*(nNodes)^2 + j*(nNodes) + i
                    if k == 1 # populate the IEN for the bottom surface
                        IEN_btm[1,nb] = IEN[1,n]
                        IEN_btm[2,nb] = IEN[2,n]
                        IEN_btm[3,nb] = IEN[3,n]
                        IEN_btm[4,nb] = IEN[4,n]
                        nb = nb + 1
                    elseif k == ne # populate the IEN for the top surface
                        IEN_top[1,nt] = IEN[5,n]
                        IEN_top[2,nt] = IEN[6,n]
                        IEN_top[3,nt] = IEN[7,n]
                        IEN_top[4,nt] = IEN[8,n]
                        nt = nt + 1
                    # elseif j == 1 || j == ne || i == 1 || i == ne
                    #     IEN_side[1,n] = IEN[1,n]
                    #     IEN_side[2,n] = IEN[2,n]
                    #     IEN_side[3,n] = IEN[3,n]
                    #     IEN_side[4,n] = IEN[4,n]
                    #     IEN_side[5,n] = IEN[5,n]
                    #     IEN_side[6,n] = IEN[6,n]
                    #     IEN_side[7,n] = IEN[7,n]
                    #     IEN_side[8,n] = IEN[8,n]
                    end
                    n = n + 1
                end
            end
        end
    elseif FunctionClass == "Q2"
        nNodes = 2*ne+1
        NodeList = zeros(Float64,ndim,(nNodes)^ndim)
        IEN = zeros(Int64,3^ndim,ne^ndim) # IEN for the 3D mesh
        IEN_top = zeros(Int64,3^(ndim-1),ne^(ndim-1)) # IEN for the top surface
        IEN_btm = zeros(Int64,3^(ndim-1),ne^(ndim-1)) # IEN for the bottom surface
        ID = zeros(Int64,ndim,(nNodes)^ndim)

        x = collect(Float64, range(-lx/2, lx/2, length=nNodes))
        y = collect(Float64, range(-ly/2, ly/2, length=nNodes))
        z = collect(Float64, range( 0, lz, length=nNodes))

        m = 1
        for k in 1:nNodes
            for j in 1:nNodes
                for i in 1:nNodes
                    NodeList[1,m] = x[i]
                    NodeList[2,m] = y[j]
                    NodeList[3,m] = z[k]
                    for l in 1:ndim
                        ID[l,m] = ndim*(m-1) + l
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
                    IEN[1,n] = 2*(k-1)*(nNodes)^2 + 2*(j-1)*(nNodes) + 2*i - 1
                    IEN[2,n] = 2*(k-1)*(nNodes)^2 + 2*(j-1)*(nNodes) + 2*i + 1
                    IEN[3,n] = 2*(k-1)*(nNodes)^2 + 2*j*(nNodes) + 2*i + 1
                    IEN[4,n] = 2*(k-1)*(nNodes)^2 + 2*j*(nNodes) + 2*i - 1
                    IEN[5,n] = 2*k*(nNodes)^2 + 2*(j-1)*(nNodes) + 2*i -1
                    IEN[6,n] = 2*k*(nNodes)^2 + 2*(j-1)*(nNodes) + 2*i + 1
                    IEN[7,n] = 2*k*(nNodes)^2 + 2*j*(nNodes) + 2*i + 1
                    IEN[8,n] = 2*k*(nNodes)^2 + 2*j*(nNodes) + 2*i - 1
                    IEN[9,n] = 2*(k-1)*(nNodes)^2 + 2*(j-1)*(nNodes) + 2*i 
                    IEN[10,n] = 2*(k-1)*(nNodes)^2 + (2*j-1)*(nNodes) + 2*i + 1
                    IEN[11,n] = 2*(k-1)*(nNodes)^2 + 2*j*(nNodes) + 2*i
                    IEN[12,n] = 2*(k-1)*(nNodes)^2 + (2*j-1)*(nNodes) + 2*i - 1
                    IEN[13,n] = 2*k*(nNodes)^2 + 2*(j-1)*(nNodes) + 2*i
                    IEN[14,n] = 2*k*(nNodes)^2 + (2*j-1)*(nNodes) + 2*i + 1
                    IEN[15,n] = 2*k*(nNodes)^2 + 2*j*(nNodes) + 2*i
                    IEN[16,n] = 2*k*(nNodes)^2 + (2*j-1)*(nNodes) + 2*i - 1
                    IEN[17,n] = (2*k-1)*(nNodes)^2 + 2*(j-1)*(nNodes) + 2*i - 1
                    IEN[18,n] = (2*k-1)*(nNodes)^2 + 2*(j-1)*(nNodes) + 2*i + 1
                    IEN[19,n] = (2*k-1)*(nNodes)^2 + 2*j*(nNodes) + 2*i + 1
                    IEN[20,n] = (2*k-1)*(nNodes)^2 + 2*j*(nNodes) + 2*i - 1
                    IEN[21,n] = (2*k-1)*(nNodes)^2 + 2*(j-1)*(nNodes) + 2*i
                    IEN[22,n] = (2*k-1)*(nNodes)^2 + (2*j-1)*(nNodes) + 2*i + 1
                    IEN[23,n] = (2*k-1)*(nNodes)^2 + 2*j*(nNodes) + 2*i
                    IEN[24,n] = (2*k-1)*(nNodes)^2 + (2*j-1)*(nNodes) + 2*i - 1
                    IEN[25,n] = 2*(k-1)*(nNodes)^2 + (2*j-1)*(nNodes) + 2*i
                    IEN[26,n] = 2*k*(nNodes)^2 + (2*j-1)*(nNodes) + 2*i
                    IEN[27,n] = (2*k-1)*(nNodes)^2 + (2*j-1)*(nNodes) + 2*i
                    if k == 1
                        IEN_btm[1,nb] = IEN[1,n]
                        IEN_btm[2,nb] = IEN[2,n]
                        IEN_btm[3,nb] = IEN[3,n]
                        IEN_btm[4,nb] = IEN[4,n]
                        IEN_btm[5,nb] = IEN[9,n]
                        IEN_btm[6,nb] = IEN[10,n]
                        IEN_btm[7,nb] = IEN[11,n]
                        IEN_btm[8,nb] = IEN[12,n]
                        IEN_btm[9,nb] = IEN[25,n]
                        nb = nb + 1
                    elseif k == ne
                        IEN_top[1,nt] = IEN[5,n]
                        IEN_top[2,nt] = IEN[6,n]
                        IEN_top[3,nt] = IEN[7,n]
                        IEN_top[4,nt] = IEN[8,n]
                        IEN_top[5,nt] = IEN[13,n]
                        IEN_top[6,nt] = IEN[14,n]
                        IEN_top[7,nt] = IEN[15,n]
                        IEN_top[8,nt] = IEN[16,n]
                        IEN_top[9,nt] = IEN[26,n]
                        nt = nt + 1
                    end
                    n = n + 1
                end
            end
        end
    else
        throw(ArgumentError("Basis function type $FunctionClass is unknown"))
    end

    return NodeList, IEN, ID, IEN_top, IEN_btm, [BorderNodes, BottomBorderNodes, TopBorderNodes]
end

"""
    meshgrid_ring(r1,r2, theta1, theta2, ne)

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
function meshgrid_ring(r1::Number, r2::Number, theta1::Number, theta2::Number, ne::Int64, ndim::Int64)
    
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
    meshgrid_cylinder(r,h,ne,ndim;FunctionClass="Q1")
Set up the mesh grid for a 3D cylinder
# Arguments:
- `r::Number` : radius of the cylinder
- `h::Number` : height of the cylinder
- `ne::Int` : number of elements
- `FunctionClass::String` : type of basis function
# Returns:
- `NodeList::Matrix{Number}{nNodes,ndim}` : array of nodes
- `IEN::Matrix{Int64}{ne^ndim,2^ndim}` : array of elements
- `ID::Matrix{Int64}{nNodes,ndim}` : array of node IDs
- `IEN_top::Matrix{Int64}{ne^(ndim-1),2^(ndim-1)}` : array of elements on the top surface
- `IEN_btm::Matrix{Int64}{ne^(ndim-1),2^(ndim-1)}` : array of elements on the bottom surface
- `BorderNodes::Vector{Int64}` : array of nodes on the boundaries
"""
function meshgrid_cylinder(r::Number, h::Number, ne::Int64; FunctionClass::String="Q1")

    NodeList, IEN, ID, IEN_top, IEN_btm, BorderNodes = meshgrid_cube(1, 1, h, ne, FunctionClass=FunctionClass)
    NodeList = inflate_cylinder(NodeList, -0.5, 0.5, -0.5, 0.5, r)

    return NodeList, IEN, ID, IEN_top, IEN_btm, BorderNodes
end

"""
    inflate_cylinder(NodeList, x0, x1, y0, y1)

Inflate the sphere to a cylinder of unit radius and height

# Arguments:
- `NodeList::Matrix{Number}{nNodes,ndim}` : array of nodes
- `x0::Number` : x-coordinate of the lower left corner of the domain
- `x1::Number` : x-coordinate of the upper right corner of the domain
- `y0::Number` : y-coordinate of the lower left corner of the domain
- `y1::Number` : y-coordinate of the upper right corner of the domain

# Returns:
- `NodeList::Matrix{Number}{nNodes,ndim}` : array of nodes
"""
function inflate_cylinder(NodeListCube::Matrix{Float64}, x0::Number, x1::Number, y0::Number, y1::Number, rad::Number=1.0)
    NodeListCyl = copy(NodeListCube)
    center = [0.5*(x0 + x1), 0.5*(y0 + y1)]

    iter = 1:size(NodeListCyl,2)
    for i in iter
        scale = maximum(abs.(NodeListCyl[1:2,i] - center))
        if scale ≈ 0.
            NodeListCyl[1:2,i] = [0 , 0]
        else
            r = sqrt((NodeListCyl[1,i] - center[1])^2 + (NodeListCyl[2,i] - center[2])^2)
            NodeListCyl[1:2,i] = 2*rad*scale*(NodeListCyl[1:2,i] - center)/r
        end
    end
    return NodeListCyl
end
