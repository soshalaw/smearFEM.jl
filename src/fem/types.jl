abstract type EnvConditions end

mutable struct Conditions <: EnvConditions
    # Define the properties of the Conditions struct
    ANIMATE::Bool
    WRITEVTK::Bool
    WRITECONTOUR::Bool
    RENDER::Bool
    filepath::String
    camera_matrix::AbstractMatrix{Float64}
    obj_pose::Vector{Float64}
    viewing_angles::Vector{Float64}

    # Constructor with keyword arguments and default values
    function Conditions(;
        ANIMATE::Bool = false,
        WRITEVTK::Bool = false,
        WRITECONTOUR::Bool = false,
        RENDER::Bool = false,
        filepath::String = "",
        camera_matrix::AbstractMatrix{Float64} = Matrix{Float64}(undef, 4, 4),
        obj_pose::Vector{Float64} = zeros(Float64, 3),
        viewing_angles::Vector{Float64} = zeros(Float64, 1)
    )
        return new(ANIMATE, WRITEVTK, WRITECONTOUR, RENDER, filepath, camera_matrix, obj_pose, viewing_angles)
    end
end