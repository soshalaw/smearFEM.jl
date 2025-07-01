# Abstract type for all models
abstract type AbstractModel end

mutable struct Model{T<:AbstractModel}
    mdl::T
end

"""
    LinearElasticity

Represents a linear elasticity model with associated parameters.
"""
mutable struct LinearElasticity <: AbstractModel
    ne::Int
    ndim::Int

    NodeList::Matrix{Float64}
    IEN::Matrix{Int}
    IEN_top::Matrix{Int}
    IEN_btm::Matrix{Int}
    IEN_border::Matrix{Int}
    ID::Matrix{Int}
    nDof::Int
    FunctionClass::String

    C::Array{Float64, 3}
    C_top::Array{Float64, 3}
    C_btm::Array{Float64, 3}
    W::Vector{Float64}

    θ1::Float64
    θ2::Float64
    cMat::Matrix{Float64}
    dcMatdθ1::Matrix{Float64}
    dcMatdθ2::Matrix{Float64}

    function LinearElasticity(;
        ne::Int=1, 
        ndim::Int=1,
        NodeList::Matrix{Float64}=Matrix{Float64}(undef, 1, 1),
        IEN::Matrix{Int}=Matrix{Int}(undef, 1, 1),
        IEN_top::Matrix{Int} = Matrix{Int}(undef, 1, 1),
        IEN_btm::Matrix{Int} = Matrix{Int}(undef, 1, 1), 
        IEN_border::Matrix{Int} = Matrix{Int}(undef, 1, 1),
        ID::Matrix{Int} = Matrix{Int}(undef, 1, 1),
        nDof::Int = 0,
        FunctionClass::String = "None",
        C::Array{Float64, 3} = Array{Float64, 3}(undef, 1, 1, 1),
        C_top::Array{Float64, 3} = Array{Float64, 3}(undef, 1, 1, 1),
        C_btm::Array{Float64, 3} = Array{Float64, 3}(undef, 1, 1, 1),
        W::Vector{Float64} = Vector{Float64}(undef, 1),  
        θ1::Float64 = 0.0, 
        θ2::Float64 = 0.0, 
        cMat::Matrix{Float64} = Matrix{Float64}(undef, 1, 1),
        dcMatdθ1::Matrix{Float64} = Matrix{Float64}(undef, 1, 1),
        dcMatdθ2::Matrix{Float64} = Matrix{Float64}(undef, 1, 1)
    )
        # Initialize the linear elasticity model with the provided parameters
        new(ne, ndim, NodeList, IEN, IEN_top, IEN_btm, IEN_border, ID, nDof, FunctionClass,
            C, C_top, C_btm, W, θ1, θ2, cMat, dcMatdθ1, dcMatdθ2)
    end
end

"""
    Stokes

Represents a Stokes flow model with associated parameters.
"""
mutable struct Stokes <: AbstractModel
    # Model dimensions
    ne::Int
    ndim::Int

    # Mesh and degrees of freedom for velocity
    mesh_u::AbstractMeshgrid
    nDof_u::Int

    # Mesh and degrees of freedom for pressure
    mesh_p::AbstractMeshgrid
    nDof_p::Int

    # Viscosity
    η::Vector{Float64}
    tick::Int

    """
        Stokes(; ndim::Int=1, mesh_u::Mesh=Meshgrid1D(), mesh_p::Mesh=Meshgrid1D(),
               nDof_u::Int=0, nDof_p::Int=0, η::Number=0.0)

        Initialize the Stokes model with the provided parameters.
    """
    function Stokes(;
        ndim::Int=1, 
        mesh_u::AbstractMeshgrid=Meshgrid1D(),
        mesh_p::AbstractMeshgrid=Meshgrid1D(),
        nDof_u::Int=0,
        nDof_p::Int=0,
        η::Vector{Float64}=Vector{Float64}(undef, 1),
        tick::Int=0
    )   
        ne = mesh_u.ne
        # Initialize the Stokes model with the provided parameters
        new(ne, ndim, mesh_u, nDof_u, mesh_p, nDof_p, η, tick)
    end
end

"""
    update_model(model::Stokes)

Updates the initial states of the nodal meshes of the model and updates the currect states of the nodes in the mesh to the initial state.

# Arguments
- `model::AbstractModel`: The model containing the meshes.

# Notes
This function modifies the meshes in place.
"""
function update_model!(model::Stokes)
    mesh_u = model.mesh_u
    mesh_p = model.mesh_p
    # Update the nodal coordinates of the meshes
    update_initial_state!(mesh_u,mesh_u.NodeList)
    update_initial_state!(mesh_p,mesh_p.NodeList)
end

function update_model!(model::LinearElasticity)
    mesh = model.mesh
    # Update the nodal coordinates of the meshes
    update_initial_state!(mesh,mesh.NodeList)
end

"""
    reset_mesh(model::Stokes)

Resets the nodal coordinates of the meshes of the model.

# Arguments
- `model::AbstractModel`: The model containing the meshes.

# Notes
This function modifies the meshes in place.
"""
function reset_model!(model::Stokes)
    mesh_u = model.mesh_u
    mesh_p = model.mesh_p
    # Reset the nodal coordinates of the meshes
    reset_mesh!(mesh_u)
    reset_mesh!(mesh_p)
end

function reset_model!(model::LinearElasticity)
    mesh = model.mesh
    # Reset the nodal coordinates of the meshes
    reset_mesh!(mesh)
end
