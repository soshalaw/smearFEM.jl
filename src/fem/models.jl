
abstract type model end

function def_model(
    
    mdl::String; 
    ne::Int64 = 1, 

    NodeList::Matrix{Float64} = [0.0 1.0 1.0 0.0; 0.0 0.0 1.0 1.0],
    IEN::Matrix{Int} = [1 2 3 4],
    IEN_top::Matrix{Int} = [1 2],
    IEN_btm::Matrix{Int} = [3 4],
    IEN_border::Matrix{Int} = [1 2 3 4],
    ID::Matrix{Int} = [1 2 3 4],
    ndim::Int64 = 2,
    nDof::Int64 = 2,
    FunctionClass::String = "Q1",
    C::Array{Float64, 3} = zeros(2,2,2),
    W::Vector{Float64}= zeros(1),  

    NodeList_2::Matrix{Float64} = [0.0 1.0 1.0 0.0; 0.0 0.0 1.0 1.0],
    IEN_2::Matrix{Int} = [1 2 3 4],
    IEN_2_top::Matrix{Int} = [1 2],
    IEN_2_btm::Matrix{Int} = [3 4],
    IEN_2_border::Matrix{Int} = [1 2 3 4],
    ID_2::Matrix{Int} = [1 2 3 4],
    ndim_2::Int64 = 2,
    nDof_2::Int64 = 2,
    FunctionClass_2::String = "Q1",
    C_2::Array{Float64, 3} = zeros(2,2,2),
    W_2::Vector{Float64}= zeros(1), 

    cMat::Matrix{Float64} = [1.0 0.3 0.3 0.0; 0.3 1.0 0.3 0.0; 0.3 0.3 1.0 0.0; 0.0 0.0 0.0 0.5],
    Young = 1.0, 
    ν = 1, 
    dcMatdλ::Matrix{Float64} = [0.0 1.0 1.1 0.1],
    dcMatdμ::Matrix{Float64} = [0.0 1.0 1.1 0.1])

    if mdl == "linear_elasticity" && FunctionClass == "S2"
        return linearElasticity(ne, ndim, NodeList, IEN, IEN_top, IEN_btm, IEN_border, ID, nDof, FunctionClass, C, W, Young, ν, cMat, dcMatdλ, dcMatdμ)
    elseif mdl == "linear_elasticity"
        return linearElasticity(ne, ndim, NodeList, IEN, IEN_top, IEN_btm, IEN_border, ID, nDof, FunctionClass, C, W, Young, ν, cMat, dcMatdλ, dcMatdμ)
    elseif mdl == "linear_elasticity" && FunctionClass == "S2"
        return stokes(ne, ndim, NodeList, IEN, IEN_top, IEN_btm, IEN_border, ID, nDof, FunctionClass,
                            NodeList_2, IEN_2, IEN_2_top, IEN_2_btm, IEN_2_border, ID_2, nDof_2, FunctionClass_2, ν)
    elseif mdl == "stokes"
        return stokes(ne, ndim, NodeList, IEN, IEN_top, IEN_btm, IEN_border, ID, nDof, FunctionClass,
                            NodeList_2, IEN_2, IEN_2_top, IEN_2_btm, IEN_2_border, ID_2, nDof_2, FunctionClass_2, ν)
    end
end

mutable struct linearElasticity <: model
    ne::Int64
    ndim::Int64

    NodeList::Matrix{Float64}
    IEN::Matrix{Int}
    IEN_top::Matrix{Int}
    IEN_btm::Matrix{Int}
    IEN_border::Matrix{Int}
    ID::Matrix{Int}
    nDof::Int64
    FunctionClass::String
    C::Array{Float64}
    W::Vector{Float64}

    Young::Float64
    ν::Float64
    cMat::Matrix{Float64}
    dcMatdλ::Matrix{Float64}
    dcMatdμ::Matrix{Float64}
end

mutable struct stokes <: model
    ne::Int64
    ndim::Int64

    NodeList_u::Matrix{Float64}
    IEN_u::Matrix{Int}
    IEN_u_top::Matrix{Int}
    IEN_u_btm::Matrix{Int}
    IEN_u_border::Matrix{Int}
    ID_u::Matrix{Int}
    nDof_u::Int64
    FunctionClass_u::String

    NodeList_p::Matrix{Float64}
    IEN_p::Matrix{Int}
    IEN_p_top::Matrix{Int}
    IEN_p_btm::Matrix{Int}
    IEN_p_border::Matrix{Int}
    ID_p::Matrix{Int}
    nDof_p::Int64
    FunctionClass_p::String

    η::Float64
end

# mutable struct Node
#     Coordinates::Vector{Float64}
#     ID::Int64
# end
# mutable struct Element
#     Nodes::Vector{Node}
#     ID::Int64
# end