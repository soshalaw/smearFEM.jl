
abstract type model end

function def_model(mdl::String; ne::Int64 = 1, 
    NodeList::Matrix{Float64} = [0.0 1.0 1.1 0.1], 
    IEN::Matrix{Int} = [1 2 3 4], 
    ndim::Int64 = 2, 
    nDof::Int64 = 1, 
    FunctionClass::String = "Q1", 
    Young = 1.0, 
    ν = 0.3, 
    ID::Matrix{Int} = [1 2 3 4],
    cMat::Matrix{Float64} = [0.0 1.0 1.1 0.1],
    IEN_top::Matrix{Int} = [1 2 3 4],
    IEN_btm::Matrix{Int} = [1 2 3 4],
    IEN_border::Matrix{Int} = [1 2 3 4],
    dcMatdλ::Matrix{Float64} = [0.0 1.0 1.1 0.1],
    dcMatdμ::Matrix{Float64} = [0.0 1.0 1.1 0.1],
    IEN_u::Matrix{Int} = [1 2 3 4],
    IEN_u_top::Matrix{Int} = [1 2 3 4],
    IEN_u_btm::Matrix{Int} = [1 2 3 4],
    IEN_u_border::Matrix{Int} = [1 2 3 4],
    IEN_p::Matrix{Int} = [1 2 3 4],
    IEN_p_top::Matrix{Int} = [1 2 3 4],
    IEN_p_btm::Matrix{Int} = [1 2 3 4],
    IEN_p_border::Matrix{Int} = [1 2 3 4]
    )

    if mdl == "linear_elasticity"
        return linearElasticity(ne, NodeList, IEN, IEN_top, IEN_btm, IEN_border, ndim, nDof, FunctionClass, ID, Young, ν, cMat, dcMatdλ, dcMatdμ)
    elseif mdl == "stokes"
        return stokes(ne, NodeList, IEN_u, IEN_u_top, IEN_u_btm, IEN_u_border, IEN_p, IEN_p_top, IEN_p_btm, IEN_p_border, ndim, nDof, FunctionClass, ID)
    end
end

mutable struct linearElasticity <: model
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

mutable struct stokes <: model
    ne::Int64
    NodeList::Matrix{Float64}
    IEN_u::Matrix{Int}
    IEN_u_top::Matrix{Int}
    IEN_u_btm::Matrix{Int}
    IEN_u_border::Matrix{Int}
    IEN_p::Matrix{Int}
    IEN_p_top::Matrix{Int}
    IEN_p_btm::Matrix{Int}
    IEN_p_border::Matrix{Int}
    ndim::Int64
    nDof::Int64
    FunctionClass::String
    ID::Matrix{Int}
end

# mutable struct Node
#     Coordinates::Vector{Float64}
#     ID::Int64
# end
# mutable struct Element
#     Nodes::Vector{Node}
#     ID::Int64
# end