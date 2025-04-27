abstract type EnvConditions end

mutable struct Conditions <: EnvConditions
    # Define the properties of the Conditions struct
    ANIMATE::Bool
    WRITEVTK::Bool
    SIDES::Bool
    WRITECONTOUR::Bool
    RENDER::Bool
    filepath::String
    CameraMatrix::AbstractMatrix{Float64}

    # Constructor with keyword arguments and default values
    function Conditions(; 
        ANIMATE::Bool = false,
        WRITEVTK::Bool = false,
        WRITECONTOUR::Bool = false,
        RENDER::Bool = false,
        SIDES::Bool = false,
        filepath::String = "",
        CameraMatrix::AbstractMatrix{Float64} = Matrix{Float64}(undef, 4, 4)
    )
        return new(ANIMATE, WRITEVTK, SIDES, WRITECONTOUR, RENDER, filepath, CameraMatrix)
    end
end