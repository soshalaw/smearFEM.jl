"""
    EnvConditions

Supertype for the non-physical settings of a run: what is recorded, rendered and written where.
Concrete subtype: `Conditions`.
"""
abstract type EnvConditions end

"""
    Conditions

Output and camera settings for a simulation run. Every field has a safe default, so a bare
`Conditions()` records nothing and renders nothing.

# Arguments
- `ANIMATE::Bool`: Write animations of the deforming mesh.
- `WRITEVTK::Bool`: Write VTK files for ParaView.
- `WRITECONTOUR::Bool`: Write the projected 2D border contours.
- `RENDER::Bool`: Project the deforming object through the camera model.
- `filepath::String`: Output directory for everything written.
- `camera_matrix::AbstractMatrix{Float64}`: Camera intrinsics.
- `obj_pose::AbstractArray{Float64}`: Either a camera position (3-vector), whose orientation is
  synthesised by `_set_cam_frame`, or a measured 4x4 object-to-camera transform (physical data).
- `viewing_angles::Vector{Float64}`: Angles, in radians, of the views to render.
"""
mutable struct Conditions <: EnvConditions
    # Define the properties of the Conditions struct
    ANIMATE::Bool
    WRITEVTK::Bool
    WRITECONTOUR::Bool
    RENDER::Bool
    filepath::String
    camera_matrix::AbstractMatrix{Float64}
    # Either a camera position (3-vector), whose orientation is synthesised by
    # `_set_cam_frame`, or a measured 4×4 object-to-camera transform (physical data).
    obj_pose::AbstractArray{Float64}
    viewing_angles::Vector{Float64}

    # Constructor with keyword arguments and default values
    function Conditions(;
        ANIMATE::Bool = false,
        WRITEVTK::Bool = false,
        WRITECONTOUR::Bool = false,
        RENDER::Bool = false,
        filepath::String = "",
        camera_matrix::AbstractMatrix{Float64} = Matrix{Float64}(undef, 4, 4),
        obj_pose::AbstractArray{Float64} = zeros(Float64, 3),
        viewing_angles::Vector{Float64} = zeros(Float64, 1)
    )
        return new(ANIMATE, WRITEVTK, WRITECONTOUR, RENDER, filepath, camera_matrix, obj_pose, viewing_angles)
    end
end