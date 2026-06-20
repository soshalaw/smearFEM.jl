abstract type AbstractGeometry end

struct Cylinder{R<:Number,H<:Number} <: AbstractGeometry
    r::R
    h::H
end

struct Cuboid{X<:Number,Y<:Number,Z<:Number} <: AbstractGeometry
    lx::X
    ly::Y
    lz::Z
end

struct Disk{R<:Number} <: AbstractGeometry
    r::R
end

struct Square{X<:Number,Y<:Number} <: AbstractGeometry
    lx::X
    ly::Y
end

struct Segment{L<:Number} <: AbstractGeometry
    l::L
end

ndim(::Union{Cylinder,Cuboid}) = 3
ndim(::Union{Disk,Square})     = 2
ndim(::Segment)                = 1
