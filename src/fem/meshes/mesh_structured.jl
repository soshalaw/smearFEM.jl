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
    elseif basis_order == 2
        nNodes = 2*ne+1
        NodeList = zeros(Float64,1,nNodes)
        IEN = zeros(Int64,3,ne) # IEN for the 1D mesh

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
function _meshgrid_square(lx::X, ly::Y, ne::Int64; element_shape::Symbol=:Quad, basis_order::Int=1) where {X<:Number,Y<:Number}

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
    _meshgrid_cube(lx,ly,lz,ne;element_shape=:Hex,basis_order=1)

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
- `element_shape::Symbol` : element shape (`:Hex`)
- `basis_order::Int` : polynomial order (1, 2, or 3)

# Returns:
- `NodeList::Matrix{Number}{nNodes,ndim}` : array of nodes
- `IEN::Matrix{Int64}{ne^ndim,2^ndim}` : array of elements
- `ID::Matrix{Int64}{nNodes,ndim}` : array of node IDs
- `IEN_top::Matrix{Int64}{ne^(ndim-1),2^(ndim-1)}` : array of elements on the top surface
- `IEN_bottom::Matrix{Int64}{ne^(ndim-1),2^(ndim-1)}` : array of elements on the bottom surface
- `BorderNodes::Vector{Int64}` : array of nodes on the boundaries
- `BottomBorderNodes::Vector{Int64}` : array of nodes on the bottom boundary
- `TopBorderNodes::Vector{Int64}` : array of nodes on the top boundary
"""
function _meshgrid_cube(lx::X, ly::Y, lz::Z, ne::Int64; element_shape::Symbol=:Hex, basis_order::Int=1) where {X<:Number,Z<:Number,Y<:Number}
    BorderNodes = Int64[]
    BottomBorderNodes = Int64[]
    TopBorderNodes = Int64[]
    ndim = 3
    IEN_side = Matrix{Int}(undef, 0, 0)

    if basis_order == 1
        nNodes = ne+1 # number of nodes in each direction
        NodeList = zeros(Float64,ndim,(nNodes)^ndim)
        IEN = zeros(Int64,2^ndim,ne^ndim) # IEN for the 3D mesh
        IEN_top = zeros(Int64,2^(ndim-1),ne^(ndim-1)) # IEN for the top surface
        IEN_bottom = zeros(Int64,2^(ndim-1),ne^(ndim-1)) # IEN for the bottom surface
        ID = zeros(Int64,ndim,(nNodes)^ndim)

        x = collect(Float64, range(-lx/2, lx/2, length=nNodes))
        y = collect(Float64, range(-ly/2, ly/2, length=nNodes))
        z = collect(Float64, range(0, lz, length=nNodes))

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
                        IEN_bottom[1,nb] = IEN[1,n]
                        IEN_bottom[2,nb] = IEN[2,n]
                        IEN_bottom[3,nb] = IEN[3,n]
                        IEN_bottom[4,nb] = IEN[4,n]
                        nb = nb + 1
                    elseif k == ne # populate the IEN for the top surface
                        IEN_top[1,nt] = IEN[5,n]
                        IEN_top[2,nt] = IEN[6,n]
                        IEN_top[3,nt] = IEN[7,n]
                        IEN_top[4,nt] = IEN[8,n]
                        nt = nt + 1
                    end
                    n = n + 1
                end
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
        # z = collect(Float64, range(-lz/2, lz/2, length=nNodes))
        z = collect(Float64, range(0, lz, length=nNodes))

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
                        IEN_bottom[1,nb] = IEN[1,n]
                        IEN_bottom[2,nb] = IEN[2,n]
                        IEN_bottom[3,nb] = IEN[3,n]
                        IEN_bottom[4,nb] = IEN[4,n]
                        IEN_bottom[5,nb] = IEN[9,n]
                        IEN_bottom[6,nb] = IEN[10,n]
                        IEN_bottom[7,nb] = IEN[11,n]
                        IEN_bottom[8,nb] = IEN[12,n]
                        IEN_bottom[9,nb] = IEN[25,n]
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
    elseif basis_order == 3
        nNodes = 3*ne+1
        NodeList = zeros(Float64,ndim,(nNodes)^ndim)
        IEN = zeros(Int64,4^ndim,ne^ndim) # IEN for the 3D mesh
        IEN_top = zeros(Int64,4^(ndim-1),ne^(ndim-1)) # IEN for the top surface
        IEN_bottom = zeros(Int64,4^(ndim-1),ne^(ndim-1)) # IEN for the bottom surface
        ID = zeros(Int64,ndim,(nNodes)^ndim)

        x = collect(Float64, range(-lx/2, lx/2, length=nNodes))
        y = collect(Float64, range(-ly/2, ly/2, length=nNodes))
        z = collect(Float64, range(0, lz, length=nNodes))

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
                    IEN[1,n] = 3*(k-1)*(nNodes)^2 + 3*(j-1)*(nNodes) + 3*i - 2
                    IEN[2,n] = 3*(k-1)*(nNodes)^2 + 3*(j-1)*(nNodes) + 3*i + 1
                    IEN[3,n] = 3*(k-1)*(nNodes)^2 + 3*j*(nNodes) + 3*i + 1
                    IEN[4,n] = 3*(k-1)*(nNodes)^2 + 3*j*(nNodes) + 3*i - 2

                    IEN[5,n] = 3*k*(nNodes)^2 + 3*(j-1)*(nNodes) + 3*i - 2
                    IEN[6,n] = 3*k*(nNodes)^2 + 3*(j-1)*(nNodes) + 3*i + 1
                    IEN[7,n] = 3*k*(nNodes)^2 + 3*j*(nNodes) + 3*i + 1
                    IEN[8,n] = 3*k*(nNodes)^2 + 3*j*(nNodes) + 3*i - 2

                    IEN[9,n] = 3*(k-1)*(nNodes)^2 + (3*j-2)*(nNodes) + 3*i - 2
                    IEN[10,n] = 3*(k-1)*(nNodes)^2 + (3*j-2)*(nNodes) + 3*i + 1
                    IEN[11,n] = 3*(k-1)*(nNodes)^2 + (3*j-1)*(nNodes) + 3*i + 1
                    IEN[12,n] = 3*(k-1)*(nNodes)^2 + (3*j-1)*(nNodes) + 3*i - 2

                    IEN[13,n] = 3*k*(nNodes)^2 + (3*j-2)*(nNodes) + 3*i - 2
                    IEN[14,n] = 3*k*(nNodes)^2 + (3*j-2)*(nNodes) + 3*i + 1
                    IEN[15,n] = 3*k*(nNodes)^2 + (3*j-1)*(nNodes) + 3*i + 1
                    IEN[16,n] = 3*k*(nNodes)^2 + (3*j-1)*(nNodes) + 3*i - 2

                    IEN[17,n] = 3*(k-1)*(nNodes)^2 + 3*(j-1)*(nNodes) + 3*i - 1
                    IEN[18,n] = 3*(k-1)*(nNodes)^2 + 3*(j-1)*(nNodes) + 3*i
                    IEN[19,n] = 3*(k-1)*(nNodes)^2 + 3*j*(nNodes) + 3*i
                    IEN[20,n] = 3*(k-1)*(nNodes)^2 + 3*j*(nNodes) + 3*i - 1

                    IEN[21,n] = 3*k*(nNodes)^2 + 3*(j-1)*(nNodes) + 3*i - 1
                    IEN[22,n] = 3*k*(nNodes)^2 + 3*(j-1)*(nNodes) + 3*i
                    IEN[23,n] = 3*k*(nNodes)^2 + 3*j*(nNodes) + 3*i
                    IEN[24,n] = 3*k*(nNodes)^2 + 3*j*(nNodes) + 3*i - 1

                    IEN[25,n] = (3*k-2)*(nNodes)^2 + 3*(j-1)*(nNodes) + 3*i - 2
                    IEN[26,n] = (3*k-2)*(nNodes)^2 + 3*(j-1)*(nNodes) + 3*i + 1
                    IEN[27,n] = (3*k-2)*(nNodes)^2 + 3*j*(nNodes) + 3*i + 1
                    IEN[28,n] = (3*k-2)*(nNodes)^2 + 3*j*(nNodes) + 3*i - 2

                    IEN[29,n] = (3*k-1)*(nNodes)^2 + 3*(j-1)*(nNodes) + 3*i - 2
                    IEN[30,n] = (3*k-1)*(nNodes)^2 + 3*(j-1)*(nNodes) + 3*i + 1
                    IEN[31,n] = (3*k-1)*(nNodes)^2 + 3*j*(nNodes) + 3*i + 1
                    IEN[32,n] = (3*k-1)*(nNodes)^2 + 3*j*(nNodes) + 3*i - 2

                    IEN[33,n] = (3*k-2)*(nNodes)^2 + (3*j-2)*(nNodes) + 3*i - 2
                    IEN[34,n] = (3*k-2)*(nNodes)^2 + (3*j-2)*(nNodes) + 3*i + 1
                    IEN[35,n] = (3*k-2)*(nNodes)^2 + (3*j-1)*(nNodes) + 3*i + 1
                    IEN[36,n] = (3*k-2)*(nNodes)^2 + (3*j-1)*(nNodes) + 3*i - 2

                    IEN[37,n] = (3*k-1)*(nNodes)^2 + (3*j-2)*(nNodes) + 3*i - 2
                    IEN[38,n] = (3*k-1)*(nNodes)^2 + (3*j-2)*(nNodes) + 3*i + 1
                    IEN[39,n] = (3*k-1)*(nNodes)^2 + (3*j-1)*(nNodes) + 3*i + 1
                    IEN[40,n] = (3*k-1)*(nNodes)^2 + (3*j-1)*(nNodes) + 3*i - 2

                    IEN[41,n] = (3*k-2)*(nNodes)^2 + 3*(j-1)*(nNodes) + 3*i - 1
                    IEN[42,n] = (3*k-2)*(nNodes)^2 + 3*(j-1)*(nNodes) + 3*i
                    IEN[43,n] = (3*k-2)*(nNodes)^2 + 3*j*(nNodes) + 3*i
                    IEN[44,n] = (3*k-2)*(nNodes)^2 + 3*j*(nNodes) + 3*i - 1

                    IEN[45,n] = (3*k-1)*(nNodes)^2 + 3*(j-1)*(nNodes) + 3*i - 1
                    IEN[46,n] = (3*k-1)*(nNodes)^2 + 3*(j-1)*(nNodes) + 3*i
                    IEN[47,n] = (3*k-1)*(nNodes)^2 + 3*j*(nNodes) + 3*i
                    IEN[48,n] = (3*k-1)*(nNodes)^2 + 3*j*(nNodes) + 3*i - 1

                    IEN[49,n] = 3*(k-1)*(nNodes)^2 + (3*j-2)*(nNodes) + 3*i - 1
                    IEN[50,n] = 3*(k-1)*(nNodes)^2 + (3*j-2)*(nNodes) + 3*i
                    IEN[51,n] = 3*(k-1)*(nNodes)^2 + (3*j-1)*(nNodes) + 3*i
                    IEN[52,n] = 3*(k-1)*(nNodes)^2 + (3*j-1)*(nNodes) + 3*i - 1

                    IEN[53,n] = 3*k*(nNodes)^2 + (3*j-2)*(nNodes) + 3*i - 1
                    IEN[54,n] = 3*k*(nNodes)^2 + (3*j-2)*(nNodes) + 3*i
                    IEN[55,n] = 3*k*(nNodes)^2 + (3*j-1)*(nNodes) + 3*i
                    IEN[56,n] = 3*k*(nNodes)^2 + (3*j-1)*(nNodes) + 3*i - 1

                    IEN[57,n] = (3*k-2)*(nNodes)^2 + (3*j-2)*(nNodes) + 3*i - 1
                    IEN[58,n] = (3*k-2)*(nNodes)^2 + (3*j-2)*(nNodes) + 3*i
                    IEN[59,n] = (3*k-2)*(nNodes)^2 + (3*j-1)*(nNodes) + 3*i
                    IEN[60,n] = (3*k-2)*(nNodes)^2 + (3*j-1)*(nNodes) + 3*i - 1

                    IEN[61,n] = (3*k-1)*(nNodes)^2 + (3*j-2)*(nNodes) + 3*i - 1
                    IEN[62,n] = (3*k-1)*(nNodes)^2 + (3*j-2)*(nNodes) + 3*i
                    IEN[63,n] = (3*k-1)*(nNodes)^2 + (3*j-1)*(nNodes) + 3*i
                    IEN[64,n] = (3*k-1)*(nNodes)^2 + (3*j-1)*(nNodes) + 3*i - 1

                    if k == 1
                        IEN_bottom[1,nb] = IEN[1,n]
                        IEN_bottom[2,nb] = IEN[2,n]
                        IEN_bottom[3,nb] = IEN[3,n]
                        IEN_bottom[4,nb] = IEN[4,n]
                        IEN_bottom[5,nb] = IEN[9,n]
                        IEN_bottom[6,nb] = IEN[10,n]
                        IEN_bottom[7,nb] = IEN[11,n]
                        IEN_bottom[8,nb] = IEN[12,n]
                        IEN_bottom[9,nb] = IEN[17,n]
                        IEN_bottom[10,nb] = IEN[18,n]
                        IEN_bottom[11,nb] = IEN[19,n]
                        IEN_bottom[12,nb] = IEN[20,n]
                        IEN_bottom[13,nb] = IEN[33,n]
                        IEN_bottom[14,nb] = IEN[34,n]
                        IEN_bottom[15,nb] = IEN[35,n]
                        IEN_bottom[16,nb] = IEN[36,n]

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
                        IEN_top[9,nt] = IEN[21,n]
                        IEN_top[10,nt] = IEN[22,n]
                        IEN_top[11,nt] = IEN[23,n]
                        IEN_top[12,nt] = IEN[24,n]
                        IEN_top[13,nt] = IEN[37,n]
                        IEN_top[14,nt] = IEN[38,n]
                        IEN_top[15,nt] = IEN[39,n]
                        IEN_top[16,nt] = IEN[40,n]

                        nt = nt + 1
                    end
                    n = n + 1
                end
            end
        end
    else
        throw(ArgumentError("Unknown basis_order $basis_order for element_shape $element_shape"))
    end

    # display(ID)
    return NodeList, IEN, ID, IEN_top, IEN_bottom, IEN_side, nNodes, [BorderNodes, BottomBorderNodes, TopBorderNodes]
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
function _meshgrid_ring(r1::T, r2::U, theta1::V, theta2::W, ne::Int64, ndim::Int64) where {T<:Number,U<:Number,V<:Number,W<:Number}

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
function _meshgrid_cylinder(r::T, h::U, ne::Int64; element_shape::Symbol=:Hex, basis_order::Int=1) where {T<:Number,U<:Number}
    @info "Representing fields / geometry using $element_shape order-$basis_order basis functions"
    NodeList, IEN, ID, IEN_top, IEN_bottom, _, nNodes, BorderNodes = _meshgrid_cube(1, 1, 1, ne; element_shape=element_shape, basis_order=basis_order)
    NodeList = _inflate_cylinder(NodeList, -0.5, 0.5, -0.5, 0.5, r, h)

    mesh = MeshgridCylinder(r=r, h=h, NodeList=NodeList, IEN=IEN, IEN_top=IEN_top, IEN_bottom=IEN_bottom, ID=ID,
        volume_element_shape=element_shape, basis_order=basis_order,
        nNodes=nNodes^3, ne=ne^3, side_nodes=BorderNodes[1], top_nodes=BorderNodes[3], bottom_nodes=BorderNodes[2])
    return mesh
end

function _inflate_cylinder(NodeListCube::Matrix{Float64}, x0::T, x1::U, y0::V, y1::W, r::Y, h::Z=1.0; GRAD::Bool=false) where {T<:Number,U<:Number,V<:Number,W<:Number,Y<:Number,Z<:Number}
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

                ∇NodeListCyl[1:2,i,1] = 2*scale*(NodeListCube[1:2,i] - center)/r_ # ∂x/∂r = K_1 and ∂y/∂r = K_2
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
