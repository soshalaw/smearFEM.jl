"""
    AbstractGeometry

Supertype for the shape of the object being squeezed. A geometry is a thin wrapper carrying only
the dimensions; it exists so that `def_problem` and `set_model` dispatch on it. Adding a shape
means adding a subtype plus a `def_problem`/`set_model` method pair, not branching inside the
existing methods.

Concrete subtypes: `Cylinder`, `Cuboid` (3D), `Disk`, `Square` (2D), `Segment` (1D).
"""
abstract type AbstractGeometry end

"""
    Cylinder(r=1.0, h=1.0)

3D cylinder — the default squeeze-flow geometry.

# Arguments
- `r::Number`: Radius, in mm.
- `h::Number`: Height, in mm.
"""
mutable struct Cylinder{R<:Number,H<:Number} <: AbstractGeometry
    r::R
    h::H

    Cylinder(r::R=1.0, h::H=1.0) where {R<:Number,H<:Number} = new{R,H}(r, h)
end

"""
    Cuboid(lx=1.0, ly=1.0, lz=1.0, edge_radius=nothing)

3D box, compressed along `lz`.

# Arguments
- `lx::Number`, `ly::Number`: Cross-section dimensions, in mm.
- `lz::Number`: Height, in mm.
- `edge_radius::Union{Float64,Nothing}`: Fillet radius on the vertical edges, `nothing` for
  sharp ones. Forwarded to the mesher, which needs it at `.geo` template level.
"""
mutable struct Cuboid{X<:Number,Y<:Number,Z<:Number} <: AbstractGeometry
    lx::X
    ly::Y
    lz::Z
    edge_radius::Union{Float64,Nothing}

    Cuboid(lx::X=1.0, ly::Y=1.0, lz::Z=1.0,
           edge_radius::Union{Float64,Nothing}=nothing) where {X<:Number,Y<:Number,Z<:Number} =
        new{X,Y,Z}(lx, ly, lz, edge_radius)
end

"""
    Disk(r=1.0)

2D disk — the axisymmetric cross-section standing in for a cylinder.

# Arguments
- `r::Number`: Radius, in mm.
"""
mutable struct Disk{R<:Number} <: AbstractGeometry
    r::R

    Disk(r::R=1.0) where {R<:Number} = new{R}(r)
end

"""
    Square(lx=1.0, ly=1.0)

2D rectangle — the plane-strain cross-section standing in for a box.

# Arguments
- `lx::Number`, `ly::Number`: Dimensions, in mm.
"""
mutable struct Square{X<:Number,Y<:Number} <: AbstractGeometry
    lx::X
    ly::Y

    Square(lx::X=1.0, ly::Y=1.0) where {X<:Number,Y<:Number} = new{X,Y}(lx, ly)
end

"""
    Segment(l=1.0)

1D segment, used for the one-dimensional verification cases.

# Arguments
- `l::Number`: Length, in mm.
"""
mutable struct Segment{L<:Number} <: AbstractGeometry
    l::L

    Segment(l::L=1.0) where {L<:Number} = new{L}(l)
end

"""
    ndim(geom)

Spatial dimension of a geometry: 3 for `Cylinder` and `Cuboid`, 2 for `Disk` and `Square`, 1 for
`Segment`.

# Arguments
- `geom::AbstractGeometry`: The geometry.

# Returns
- `::Int`: Number of spatial dimensions.
"""
ndim(::Union{Cylinder,Cuboid}) = 3
ndim(::Union{Disk,Square})     = 2
ndim(::Segment)                = 1
