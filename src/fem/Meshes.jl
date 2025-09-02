# abstract type Mesh end
using ArgCheck

struct Meshgrid{T <: AbstractMeshgrid} 
    msh::T
end

mutable struct Meshgrid1D <: AbstractMeshgrid
    # Define the properties of the Meshgrid1D struct
    NodeList::Matrix{Float64}
    IEN::Matrix{Int}
    ID::Matrix{Int}
    FunctionClass::String
    nNodes::Int
    ne::Int
    boundary_nodes::Vector{Int}

    function Meshgrid1D(;
        NodeList::Matrix{Float64}=Matrix{Float64}(undef, 1, 1),
        IEN::Matrix{Int}=Matrix{Int}(undef, 2, 1),
        ID::Matrix{Int}=Matrix{Int}(undef, 1, 1),
        FunctionClass::String="Q1",
        nNodes::Int=0,
        ne::Int=0,
        boundary_nodes::Vector{Int}=Vector{Int}()
    )
        # Constructor for Meshgrid1D
        new(NodeList, IEN, ID, FunctionClass, nNodes, ne)
        
    end
end
Int
mutable struct Meshgrid2D <: AbstractMeshgrid
    # Define the properties of the Meshgrid2D struct
    NodeList::Matrix{Float64}
    IEN::Matrix{Int}
    IEN_boundaries::Vector{Matrix{Int}}
    ID::Matrix{Int}
    FunctionClass::String
    nNodes::Int
    ne::Int
    boundary_nodes::Vector{Vector{Int}}

    function Meshgrid2D(;
        NodeList::Matrix{Float64}=Matrix{Float64}(undef, 2, 1),
        IEN::Matrix{Int}=Matrix{Int}(undef, 4, 1),
        IEN_boundaries::Vector{Matrix{Int}}=[Matrix{Int}(undef, 4, 1)],
        ID::Matrix{Int}=Matrix{Int}(undef, 2, 1),
        FunctionClass::String="Q1",
        nNodes::Int=0,
        ne::Int=0,
        boundary_nodes::Vector{Vector{Int}}=[Vector{Int}()]
    )
        # Constructor for Meshgrid2D
        new(NodeList, IEN, IEN_boundaries, ID, FunctionClass, nNodes, ne, boundary_nodes)
        
    end
end

mutable struct Meshgrid3D <: AbstractMeshgrid
    # Define the properties of the Meshgrid3D struct
    NodeList::Matrix{Float64}
    IEN::Matrix{Int}
    IEN_boundaries::Vector{Matrix{Int}}
    ID::Matrix{Int}
    FunctionClass::String
    nNodes::Int
    ne::Int
    boundary_nodes::Vector{Vector{Int}}

    function Meshgrid3D(;
        NodeList::Matrix{Float64}=Matrix{Float64}(undef, 3, 1),
        IEN::Matrix{Int}=Matrix{Int}(undef, 8, 1),
        IEN_boundaries::Vector{Matrix{Int}}=[Matrix{Int}(undef, 4, 1)],
        ID::Matrix{Int}=Matrix{Int}(undef, 3, 1),
        FunctionClass::String="Q1",
        nNodes::Int=0,
        ne::Int=0,
        boundary_nodes::Vector{Vector{Int}}=[Vector{Int}()]
    )
        # Constructor for Meshgrid3D
        new(NodeList, IEN, IEN_boundaries, ID, FunctionClass, nNodes, ne, boundary_nodes)
        
    end
end
mutable struct MeshgridLine <: AbstractMeshgrid
    # Define the properties of the MeshgridLine struct
    lx::Number
    NodeList::Matrix{Float64}
    IEN::Matrix{Int}
    ID::Matrix{Int}
    FunctionClass::String
    nNodes::Int
    ne::Int

    function MeshgridLine(;
        NodeList::Matrix{Float64}=Matrix{Float64}(undef, 1, 1),
        IEN::Matrix{Int}=Matrix{Int}(undef, 2, 1),
        ID::Matrix{Int}=Matrix{Int}(undef, 1, 1),
        FunctionClass::String="Q1",
        nNodes::Int=0,
        ne::Int=0,
        lx::Number=0.0
    )
        # Constructor for MeshgridLine
        new(NodeList, IEN, ID, FunctionClass, nNodes, ne, lx)
        
    end
end

mutable struct MeshgridSquare <: AbstractMeshgrid
    # Define the properties of the MeshgridSquare struct
    lx::Number
    ly::Number
    NodeList::Matrix{Float64}
    IEN::Matrix{Int}
    IEN_top::Matrix{Int}
    IEN_bottom::Matrix{Int}
    IEN_side::Matrix{Int}
    ID::Matrix{Int}
    FunctionClass::String
    nNodes::Int
    ne::Int
    top_nodes::Vector{Int}
    bottom_nodes::Vector{Int}
    side_nodes::Vector{Int}

    function MeshgridSquare(;
        lx::Number=0.0,
        ly::Number=0.0,
        NodeList::Matrix{Float64}=Matrix{Float64}(undef, 2, 1),
        IEN::Matrix{Int}=Matrix{Int}(undef, 4, 1),
        IEN_top::Matrix{Int}=Matrix{Int}(undef, 2, 1),
        IEN_bottom::Matrix{Int}=Matrix{Int}(undef, 2, 1),
        IEN_side::Matrix{Int}=Matrix{Int}(undef, 2, 1),
        ID::Matrix{Int}=Matrix{Int}(undef, 2, 1),
        FunctionClass::String="Q1",
        nNodes::Int=0,
        ne::Int=0,
        top_nodes::Vector{Int}=Vector{Int}(),
        bottom_nodes::Vector{Int}=Vector{Int}(),
        side_nodes::Vector{Int}=Vector{Int}()
    )
        # Constructor for MeshgridSquare
        new(lx, ly, NodeList, IEN, IEN_top, IEN_bottom, IEN_side, ID,
            FunctionClass, nNodes, ne, top_nodes, bottom_nodes, side_nodes)
        
    end
end

mutable struct MeshgridCube <: AbstractMeshgrid
    # Define the properties of the MeshgridCube struct
    lx::Number
    ly::Number
    lz::Number
    NodeList::Matrix{Float64}
    IEN::Matrix{Int}
    IEN_top::Matrix{Int}
    IEN_bottom::Matrix{Int}
    IEN_sides::Dict{Symbol, Matrix{Int}}
    ID::Matrix{Int}
    FunctionClass::String
    nNodes::Int
    ne::Int
    top_nodes::Vector{Int}
    bottom_nodes::Vector{Int}
    side_nodes::Vector{Int}
    initial_state::Matrix{Float64}

    function MeshgridCube(;
        lx::Number=0.0,
        ly::Number=0.0,
        lz::Number=0.0,
        NodeList::Matrix{Float64}=Matrix{Float64}(undef, 3, 1),
        IEN::Matrix{Int}=Matrix{Int}(undef, 8, 1),
        IEN_top::Matrix{Int}=Matrix{Int}(undef, 4, 1),
        IEN_bottom::Matrix{Int}=Matrix{Int}(undef, 4, 1),
        IEN_sides_::Vector{Matrix{Int}}=[Matrix{Int}(undef, 4, 1), Matrix{Int}(undef, 4, 1), Matrix{Int}(undef, 4, 1), Matrix{Int}(undef, 4, 1)],
        ID::Matrix{Int}=Matrix{Int}(undef, 3, 1),
        FunctionClass::String="Q1",
        nNodes::Int=0,
        ne::Int=0,
        top_nodes::Vector{Int}=Vector{Int}(),
        bottom_nodes::Vector{Int}=Vector{Int}(),
        side_nodes::Vector{Int}=Vector{Int}()
    )
        # Constructor for MeshgridCube
        IEN_sides = Dict(:front => IEN_sides_[1], :back => IEN_sides_[2], :left => IEN_sides_[3], :right => IEN_sides_[4])

        new(lx, ly, lz, NodeList, IEN, IEN_top, IEN_bottom, IEN_sides,
            ID, FunctionClass, nNodes, ne, top_nodes, bottom_nodes, side_nodes, NodeList)
        
    end
end

mutable struct MeshgridCylinder <: AbstractMeshgrid
    # Define the properties of the MeshgridCylinder struct
    r::Number
    h::Number
    NodeList::Matrix{Float64}
    IEN::Matrix{Int}
    IEN_cp::Matrix{Int}
    IEN_top::Matrix{Int}
    IEN_bottom::Matrix{Int}
    IEN_sides::Matrix{Int}
    ID::Matrix{Int}
    FunctionClass::String
    nNodes::Int
    ne::Int
    top_nodes::Vector{Int}
    bottom_nodes::Vector{Int}
    side_nodes::Vector{Int}
    initial_state::Matrix{Float64}
    C_vol::Array{Float64}
    C_top::Array{Float64}
    C_btm::Array{Float64}
    W::Vector{Float64}
    
    function MeshgridCylinder(;
        r::Number=0.0,
        h::Number=0.0,
        NodeList::Matrix{Float64}=Matrix{Float64}(undef, 3, 1),
        IEN::Matrix{Int}=Matrix{Int}(undef, 8, 1),
        IEN_cp::Matrix{Int}=Matrix{Int}(undef, 8, 1),
        IEN_top::Matrix{Int}=Matrix{Int}(undef, 4, 1),
        IEN_bottom::Matrix{Int}=Matrix{Int}(undef, 4, 1),
        IEN_sides::Matrix{Int}=Matrix{Int}(undef, 4, 1),
        ID::Matrix{Int}=Matrix{Int}(undef, 3, 1),
        FunctionClass::String="Q1",
        nNodes::Int=0,
        ne::Int=0,
        top_nodes::Vector{Int}=Vector{Int}(),
        bottom_nodes::Vector{Int}=Vector{Int}(),
        side_nodes::Vector{Int}=Vector{Int}(),
        C_vol::Array{Float64, 3} = Array{Float64, 3}(undef, 1, 1, 1),
        C_top::Array{Float64, 3} = Array{Float64, 3}(undef, 1, 1, 1),
        C_btm::Array{Float64, 3} = Array{Float64, 3}(undef, 1, 1, 1),
        W::Vector{Float64} = Vector{Float64}(undef, 1)
    )
        # Constructor for MeshgridCylinder
        new(r, h, NodeList, IEN, IEN_cp, IEN_top, IEN_bottom, IEN_sides,
            ID, FunctionClass, nNodes, ne, top_nodes, bottom_nodes, side_nodes, NodeList, C_vol, C_top, C_btm, W)
        
    end
end

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
function meshgrid_line(l::T, ne::Int64; FunctionClass::String="Q1") where {T<:Number}
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
- `IEN_bottom::Matrix{Int64}{2^(ndim-1),ne^(ndim-1)}` : array of elements on the bottom surface
- `BorderNodes::Vector{Int64}` : array of nodes on the boundaries
- `BottomBorderNodes::Vector{Int64}` : array of nodes on the bottom boundary
- `TopBorderNodes::Vector{Int64}` : array of nodes on the top boundary
"""
function meshgrid_square(lx::X, ly::Y, ne::Int64; FunctionClass::String="Q1") where {X<:Number,Y<:Number}
            
    BorderNodes = Int64[]
    BottomBorderNodes = Int64[]
    TopBorderNodes = Int64[]
    ndim = 2
    if FunctionClass == "Q1"
        nNodes = ne+1 # number of nodes in each direction
        NodeList = zeros(Float64,ndim,(nNodes)^ndim)
        IEN = zeros(Int64,2^ndim,ne^ndim)              # IEN for the 3D mesh
        IEN_top = zeros(Int64,2^(ndim-1),ne^(ndim-1))  # IEN for the top surface
        IEN_bottom = zeros(Int64,2^(ndim-1),ne^(ndim-1))  # IEN for the bottom surface
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
                    IEN_bottom[1,i] = IEN[1,n]
                    IEN_bottom[2,i] = IEN[2,n]
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
        throw(ArgumentError("Basis function type $FunctionClass is unknown"))
    end
    return NodeList, IEN, ID, IEN_top, IEN_bottom, [BorderNodes, BottomBorderNodes, TopBorderNodes]
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
- `IEN_bottom::Matrix{Int64}{ne^(ndim-1),2^(ndim-1)}` : array of elements on the bottom surface
- `BorderNodes::Vector{Int64}` : array of nodes on the boundaries
- `BottomBorderNodes::Vector{Int64}` : array of nodes on the bottom boundary
- `TopBorderNodes::Vector{Int64}` : array of nodes on the top boundary
"""
function meshgrid_cube(lx::X, ly::Y, lz::Z, ne::Int64; FunctionClass::String="Q1") where {X<:Number,Z<:Number,Y<:Number}
    BorderNodes = Int64[]
    BottomBorderNodes = Int64[]
    TopBorderNodes = Int64[]
    ndim = 3

    if FunctionClass == "Q1"
        nNodes = ne+1 # number of nodes in each direction
        NodeList = zeros(Float64,ndim,(nNodes)^ndim)
        IEN = zeros(Int64,2^ndim,ne^ndim) # IEN for the 3D mesh
        IEN_top = zeros(Int64,2^(ndim-1),ne^(ndim-1)) # IEN for the top surface
        IEN_bottom = zeros(Int64,2^(ndim-1),ne^(ndim-1)) # IEN for the bottom surface
        IEN_side = zeros(Int64,2^(ndim-1),ne^(ndim-1)) # IEN for the side surfaces
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
                    # elseif j == 1 || j == ne || i == 1 || i == ne # populate the IEN for the side surfaces
                    #     IEN_side[1,n] = IEN[1,n]
                    #     IEN_side[2,n] = IEN[2,n]
                    #     IEN_side[3,n] = IEN[3,n]
                    #     IEN_side[4,n] = IEN[4,n]
                    #     IEN_side[5,n] = IEN[5,n]
                    #     IEN_side[6,n] = IEN[6,n]
                    #     IEN_side[7,n] = IEN[7,n]
                    #     IEN_side[8,n] = IEN[8,n]
                    #     IEN_side[9,n] = IEN[9,n]
                    #     IEN_side[10,n] = IEN[10,n]
                    #     IEN_side[11,n] = IEN[11,n]
                    #     IEN_side[12,n] = IEN[12,n]
                    #     IEN_side[13,n] = IEN[13,n]
                    #     IEN_side[14,n] = IEN[14,n]
                    #     IEN_side[15,n] = IEN[15,n]
                    #     IEN_side[16,n] = IEN[16,n]
                    #     IEN_side[17,n] = IEN[17,n]
                    #     IEN_side[18,n] = IEN[18,n]
                    #     IEN_side[19,n] = IEN[19,n]
                    #     IEN_side[20,n] = IEN[20,n]
                    #     IEN_side[21,n] = IEN[21,n]
                    #     IEN_side[22,n] = IEN[22,n]
                    #     IEN_side[23,n] = IEN[23,n]
                    #     IEN_side[24,n] = IEN[24,n]
                    #     IEN_side[25,n] = IEN[25,n]
                    #     IEN_side[26,n] = IEN[26,n]
                    #     IEN_side[27,n] = IEN[27,n]
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
        IEN_bottom = zeros(Int64,3^(ndim-1),ne^(ndim-1)) # IEN for the bottom surface
        IEN_side = zeros(Int64,3^(ndim-1),ne^(ndim-1)) # IEN for the side surfaces
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
                    # elseif j == 1 || j == ne || i == 1 || i == ne # populate the IEN for the side surfaces
                    #     IEN_side[1,n] = IEN[1,n]
                    #     IEN_side[2,n] = IEN[2,n]
                    #     IEN_side[3,n] = IEN[3,n]
                    #     IEN_side[4,n] = IEN[4,n]
                    #     IEN_side[5,n] = IEN[5,n]
                    #     IEN_side[6,n] = IEN[6,n]
                    #     IEN_side[7,n] = IEN[7,n]
                    #     IEN_side[8,n] = IEN[8,n]
                    #     IEN_side[9,n] = IEN[9,n]
                    #     IEN_side[10,n] = IEN[10,n]
                    #     IEN_side[11,n] = IEN[11,n]
                    #     IEN_side[12,n] = IEN[12,n]
                    #     IEN_side[13,n] = IEN[13,n]
                    #     IEN_side[14,n] = IEN[14,n]
                    #     IEN_side[15,n] = IEN[15,n]
                    #     IEN_side[16,n] = IEN[16,n]
                    #     IEN_side[17,n] = IEN[17,n]
                    #     IEN_side[18,n] = IEN[18,n]
                    #     IEN_side[19,n] = IEN[19,n]
                    #     IEN_side[20,n] = IEN[20,n]
                    #     IEN_side[21,n] = IEN[21,n]
                    #     IEN_side[22,n] = IEN[22,n]
                    #     IEN_side[23,n] = IEN[23,n]
                    #     IEN_side[24,n] = IEN[24,n]
                    #     IEN_side[25,n] = IEN[25,n]
                    #     IEN_side[26,n] = IEN[26,n]
                    #     IEN_side[27,n] = IEN[27,n]
                    end
                    n = n + 1
                end
            end
        end
    else
        throw(ArgumentError("Basis function type $FunctionClass is unknown"))
    end

    return NodeList, IEN, ID, IEN_top, IEN_bottom, IEN_side, nNodes, [BorderNodes, BottomBorderNodes, TopBorderNodes] 
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
function meshgrid_ring(r1::T, r2::U, theta1::V, theta2::W, ne::Int64, ndim::Int64) where {T<:Number,U<:Number,V<:Number,W<:Number}
    
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
- `IEN_bottom::Matrix{Int64}{ne^(ndim-1),2^(ndim-1)}` : array of elements on the bottom surface
- `BorderNodes::Vector{Int64}` : array of nodes on the boundaries
"""
function meshgrid_cylinder(r::T, h::U, ne::Int64; FunctionClass::String="Q1", filePath::String=" ") where {T<:Number,U<:Number}
    if string(FunctionClass[1]) == "S"
        @info "Representing fields / geometry using $FunctionClass NURBS"
        
        if ne == 2
            CPointList, W, C, IEN, IEN_cp, IEN_top, C_top, IEN_btm, C_btm = read_h5(string(filePath,"/cylinder_1.h5"),"sim")
        elseif ne == 4
            CPointList, W, C, IEN, IEN_cp, IEN_top, C_top, IEN_btm, C_btm = read_h5(string(filePath,"/cylinder_2.h5"),"sim")
        elseif ne == 8
            CPointList, W, C, IEN, IEN_cp, IEN_top, C_top, IEN_btm, C_btm = read_h5(string(filePath,"/cylinder_3.h5"),"sim")
        elseif ne == 16
            CPointList, W, C, IEN, IEN_cp, IEN_top, C_top, IEN_btm, C_btm = read_h5(string(filePath,"/cylinder_4.h5"),"sim")
        elseif ne == 16
            CPointList, W, C, IEN, IEN_cp, IEN_top, C_top, IEN_btm, C_btm = read_h5(string(filePath,"/cylinder_5.h5"),"sim")
        else
            throw(ArgumentError("No mesh file found for ne=$ne"))
        end

        ndim = size(CPointList,1)
        @argcheck ne^ndim == size(IEN,2) "Loaded mesh does not have the correct number of elements" 
        nNodes = ne + 2
        @argcheck nNodes^ndim == size(CPointList,2) "Loaded mesh does not have the correct number of nodes"

        ID = zeros(Int64, ndim, size(CPointList,2))
        cpiter = 1:nNodes^ndim
        for m in cpiter
            for l in 1:ndim
                ID[l,m] = ndim*(m-1) + l
            end
        end 
        
        mesh = MeshgridCylinder(r=r, h=h, NodeList=CPointList, IEN=IEN, IEN_cp=IEN_cp, IEN_top=IEN_top, IEN_bottom=IEN_btm, ID=ID, FunctionClass=FunctionClass,
            C_vol=C, C_top=C_top, C_btm=C_btm, W=W, nNodes=nNodes, ne=ne)
    elseif string(FunctionClass[1]) == "Q"  
        @info "Representing fields / geometry using $FunctionClass basis functions"         
        NodeList, IEN, ID, IEN_top, IEN_bottom, IEN_side, nNodes, BorderNodes = meshgrid_cube(1, 1, 1, ne, FunctionClass=FunctionClass)
        NodeList = inflate_cylinder(NodeList, -0.5, 0.5, -0.5, 0.5, r, h)

        mesh = MeshgridCylinder(r=r, h=h, NodeList=NodeList, IEN=IEN, IEN_top=IEN_top, IEN_bottom=IEN_bottom, ID=ID, FunctionClass=FunctionClass,
            nNodes=nNodes, ne=ne, side_nodes=BorderNodes[1], top_nodes=BorderNodes[3], bottom_nodes=BorderNodes[2])
    else
        ArgumentError("FunctionClass type unknown")
    end

    return mesh
end

# Functions to manipulate the meshgrid

"""
    inflate_cylinder(NodeList, x0, x1, y0, y1)

Inflate the cube into a cylinder.

# Arguments:
- `NodeListCube::Matrix{Float64}{nNodes,ndim}` : array of nodes
- `x0::Number` : x-coordinate of the lower left corner of the domain
- `x1::Number` : x-coordinate of the upper right corner of the domain
- `y0::Number` : y-coordinate of the lower left corner of the domain
- `y1::Number` : y-coordinate of the upper right corner of the domain
- `r::Number`  : Radius of the output cylindrical mesh

# Returns:
- `NodeListCyl::Matrix{Float64}{nNodes,ndim}` : array of nodes of the cylindrical mesh
"""
# function inflate_cylinder(NodeListCube::Matrix{Float64}, x0::T, x1::U, y0::V, y1::W, r::Array{Y}, h::Z=1.0; GRAD::Bool=false) where {T<:Number,U<:Number,V<:Number,W<:Number,Y<:Number,Z<:Number}
#     NodeListCyl = copy(NodeListCube)
#     ∇NodeListCyl = zeros(Float64,size(NodeListCube,1),size(NodeListCube,2),2)

#     nNodes = 1
#     center = [0.5*(x0 + x1), 0.5*(y0 + y1)]
#     sz_layer = size(NodeListCyl,2)/nNodes
    
#     iter = 1:sz_layer
#     layer_iter = 1:nNodes

#     for layer in layer_iter
#         if GRAD == true
#             for i in iter
#                 scale = maximum(abs.(NodeListCube[1:2,((layer-1)+i)] - center))
#                 if scale ≈ 0.
#                     NodeListCyl[1:2,((layer-1)+i)] = [0 , 0]
#                 else
#                     r_ = sqrt((NodeListCube[1,((layer-1)+i)] - center[1])^2 + (NodeListCube[2,((layer-1)+i)] - center[2])^2)
#                     NodeListCyl[1:2,i] = 2*r[layer]*scale*(NodeListCube[1:2,((layer-1)+i)] - center)/r_
#                     NodeListCyl[3,i] = NodeListCube[3,((layer-1)+i)]*h

#                     ∇NodeListCyl[1:2,i,1] = 2*scale*(NodeListCube[1:2,((layer-1)+i)] - center)/r_ # ∂x/∂r = K_1 and ∂y/∂r = K_2
#                     ∇NodeListCyl[3,i,2] = NodeListCube[3,((layer-1)+i)] # ∂z/∂h = λ_i 
#                 end
#             end
#             return NodeListCyl, ∇NodeListCyl
#         else
#             for i in iter
#                 scale = maximum(abs.(NodeListCube[1:2,((layer-1)+i)] - center))
#                 if scale ≈ 0.
#                     NodeListCyl[1:2,((layer-1)+i)] = [0 , 0]
#                 else
#                     r_ = sqrt((NodeListCube[1,((layer-1)+i)] - center[1])^2 + (NodeListCube[2,((layer-1)+i)] - center[2])^2)
#                     NodeListCyl[1:2,i] = 2*r[layer]*scale*(NodeListCube[1:2,((layer-1)+i)] - center)/r_
#                     NodeListCyl[3,i] = NodeListCube[3,((layer-1)+i)]*h
#                 end
#             end
#             return NodeListCyl
#         end
#     end
# end

function inflate_cylinder(NodeListCube::Matrix{Float64}, x0::T, x1::U, y0::V, y1::W, r::Y, h::Z=1.0; GRAD::Bool=false) where {T<:Number,U<:Number,V<:Number,W<:Number,Y<:Number,Z<:Number}
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

"""
    reset_mesh!(mesh::AbstractMeshgrid)

Reset the mesh to its initial state by updating the `NodeList` to the initial state of the mesh.

# Arguments:
- `mesh::AbstractMeshgrid`: The mesh object to be reset.
"""
function reset_mesh!(mesh::T) where {T<:AbstractMeshgrid}
    mesh.NodeList = mesh.initial_state
end

"""
    update_initial_state!(mesh::AbstractMeshgrid, new_state::Matrix{Float64})

Update the initial state of the mesh and reset the current state of the mesh to the initial state.

# Arguments:
- `mesh::AbstractMeshgrid`: The mesh object whose state is to be updated.
- `new_state::Matrix{Float64}`: The new state to be assigned to the mesh.

"""
function update_initial_state!(mesh::T, new_state::Matrix{Float64}) where {T<:AbstractMeshgrid}
    mesh.initial_state = copy(new_state)  # Use `copy` to avoid unintended mutations
    reset_mesh!(mesh)  # Reset the mesh to the initial state
end

function init_cylinder(camera_matrix, camera_pose,side_node_list, nNodes_u)

    NodeList, IEN, ID, IEN_top, IEN_bottom, IEN_side, nNodes, BorderNodes = meshgrid_cube(1, 1, h, ne, FunctionClass=FunctionClass)
    r_0 = 1
    h_0 = 1
    NodeList = inflate_cylinder(NodeList, -0.5, 0.5, -0.5, 0.5, r_0, GRAD=true)

    BorderPts2D, SurfacePts2D = extract_borders(NodeList, camera_matrix, camera_pose, side_node_list, nNodes_u)

end

