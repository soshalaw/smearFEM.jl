mutable struct model
    ne::Int64
    NodeList::Matrix{Float64}
    IEN::Matrix{Int}
    IEN_top::Matrix{Int}
    IEN_btm::Matrix{Int}
    IEN_border::Matrix{Int}
    ndim::Int64
    nDof::Int64
    FunctionClass::String
    ID::Matrix{Int}
    Young::Float64
    ν::Float64
    cMat::Matrix{Float64}
    dcMatdλ::Matrix{Float64}
    dcMatdμ::Matrix{Float64}
end

mutable struct Node
    3DCoordinates::Vector{Float64}
    2DCoordinates::Vector{Float64}
    ID::Int64
end

mutable struct Element
    Nodes::Vector{Node}
    ID::Int64
end