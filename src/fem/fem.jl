using LinearAlgebra
using SparseArrays
using Parameters

"""
    BasisFunctionCache

Structure to store pre-computed basis functions for volume and surface integration.

# Fields:
- `bf_volume` : Tuple of (N_u_gp, ΔN_u_gp, N_p_gp, ΔN_p_gp, N_x_gp, ΔN_x_gp, wpoints) for volume integration
- `bf_surface` : Tuple of (N_u_surf_gp, ΔN_u_surf_gp, wpoints) for surface integration
"""
struct BasisFunctionCache
    bf_volume::NTuple{7, Any}
    bf_surface::NTuple{3, Any}
    
    function BasisFunctionCache(bf_vol::NTuple{7, Any}, bf_surf::NTuple{3, Any})
        new(bf_vol, bf_surf)
    end
end

"""
    BasisFunctionCache(mdl::Stokes)

Construct a BasisFunctionCache by pre-computing basis functions for the given Stokes model.

# Arguments:
- `mdl::Stokes` : Material model

# Returns:
- `BasisFunctionCache` : Cache containing pre-computed basis functions
"""
function BasisFunctionCache(mdl::Stokes)
    bf_volume = get_basis_volume_functions(mdl)
    bf_surface = get_surface_basis_functions(mdl)
    return BasisFunctionCache(bf_volume, bf_surface)
end

# Node index pairs for tensor-product basis functions
const _Q1_2D_NODE_PAIRS = ((1, 1), (2, 1), (2, 2), (1, 2))
const _Q2_2D_NODE_PAIRS = ((1, 1), (2, 1), (2, 2), (1, 2), (3, 1), (2, 3), (3, 2), (1, 3), (3, 3))

const _Q1_3D_NODE_TRIPLES = ((1, 1, 1), (2, 1, 1), (2, 2, 1), (1, 2, 1), (1, 1, 2), (2, 1, 2), (2, 2, 2), (1, 2, 2))
const _Q2_3D_NODE_TRIPLES = ((1, 1, 1), (2, 1, 1), (2, 2, 1), (1, 2, 1), (1, 1, 2), (2, 1, 2), (2, 2, 2), (1, 2, 2),
    (3, 1, 1), (2, 3, 1), (3, 2, 1), (1, 3, 1),
    (3, 1, 2), (2, 3, 2), (3, 2, 2), (1, 3, 2),
    (1, 1, 3), (2, 1, 3), (2, 2, 3), (1, 2, 3),
    (3, 1, 3), (2, 3, 3), (3, 2, 3), (1, 3, 3),
    (3, 3, 1), (3, 3, 2), (3, 3, 3))

""" 
    smearFEM.gaussian_quadrature(a,b,n_gauss_pts)

Compute the nodes and weights for the Gaussian quadrature of order 2
    
# Arguments:    
- `a,b::Integer` : the limits of the integration interval
- `n_gauss_pts::Integer` : number of Gauss points to be considered (2 or 3)

# Returns:    
- `ξ::Vector{Float64}{,n_gauss_pts}`: nodes.
- `w::Vector{Float64}{,n_gauss_pts}`: weights 
"""
function gaussian_quadrature(a::Int64,b::Int64,n_gauss_pts::Int64)
  
    if n_gauss_pts == 2
        ξ = [-(b-a)/(2*sqrt(3))+(b+a)/2, (b-a)/(2*sqrt(3))+(b+a)/2]
        w = [(b-a)/2, (b-a)/2]
    elseif n_gauss_pts == 3
        ξ = [-(b-a)/(2*sqrt(5/3))+(b+a)/2, 0, (b-a)/(2*sqrt(5/3))+(b+a)/2]
        w = [(b-a)/2*5/9, (b-a)/2*8/9, (b-a)/2*5/9]
    end
    return ξ, w
end

function gaussian_quadrature(a::Int64, b::Int64; n_gauss_pts::Int64=2, nGaussPoints::Union{Nothing,Int64}=nothing)
    return gaussian_quadrature(a, b, isnothing(nGaussPoints) ? n_gauss_pts : nGaussPoints)
end

"""
    _basis_1d_components(ξ, FunctionClass)

Return 1D basis values, first derivatives, and second derivatives for the
requested polynomial basis family.
"""
function _basis_1d_components(ξ::Float64, FunctionClass::String)
    if FunctionClass == "Q1" || FunctionClass == "S1"
        N = [0.5 - 0.5 * ξ, 0.5 + 0.5 * ξ]
        ΔN = [-0.5, 0.5]
        Δ2N = [0.0, 0.0]
    elseif FunctionClass == "Q2"
        N = [-(1 - ξ) * ξ / 2,
              ξ * (1 + ξ) / 2,
             (1 - ξ) * (1 + ξ)]

        ΔN = [-(1 - 2 * ξ) / 2,
               (1 + 2 * ξ) / 2,
              -2 * ξ]

        Δ2N = [1.0, 1.0, -2.0]
    else
        error("Unsupported FunctionClass: $FunctionClass")
    end

    return N, ΔN, Δ2N
end

"""
    _format_1d_derivatives(FunctionClass, ΔN, Δ2N)

Normalize derivative return shapes for 1D basis functions so legacy callers
keep receiving matrix output for `Q1`/`S1` and vector output for `Q2`.
"""
function _format_1d_derivatives(FunctionClass::String, ΔN, Δ2N)
    if FunctionClass == "Q1" || FunctionClass == "S1"
        return reshape(ΔN, 1, :), reshape(Δ2N, 1, :)
    end

    return ΔN, Δ2N
end

"""
    _tensor_basis_2d(Nξ, ΔNξ, Δ2Nξ, Nη, ΔNη, Δ2Nη, node_pairs)

Build 2D tensor-product basis values and derivatives from 1D components.

The returned `Δ2N` columns follow `[ξξ, ξη, ηη]`.
"""
function _tensor_basis_2d(Nξ, ΔNξ, Δ2Nξ, Nη, ΔNη, Δ2Nη, node_pairs)
    nBasis = length(node_pairs)
    N = Vector{Float64}(undef, nBasis)
    ΔN = Matrix{Float64}(undef, nBasis, 2)
    Δ2N = Matrix{Float64}(undef, nBasis, 3)

    for (i, (iξ, iη)) in enumerate(node_pairs)
        N[i] = Nξ[iξ] * Nη[iη]
        ΔN[i, 1] = ΔNξ[iξ] * Nη[iη]
        ΔN[i, 2] = Nξ[iξ] * ΔNη[iη]
        Δ2N[i, 1] = Δ2Nξ[iξ] * Nη[iη]
        Δ2N[i, 2] = ΔNξ[iξ] * ΔNη[iη]
        Δ2N[i, 3] = Nξ[iξ] * Δ2Nη[iη]
    end

    return N, ΔN, Δ2N
end

"""
    _tensor_basis_3d(Nξ, ΔNξ, Δ2Nξ, Nη, ΔNη, Δ2Nη, Nζ, ΔNζ, Δ2Nζ, node_triples)

Build 3D tensor-product basis values and derivatives from 1D components.

The returned `Δ2N` columns follow `[ξξ, ξη, ξζ, ηη, ηζ, ζζ]`.
"""
function _tensor_basis_3d(Nξ, ΔNξ, Δ2Nξ, Nη, ΔNη, Δ2Nη, Nζ, ΔNζ, Δ2Nζ, node_triples)
    nBasis = length(node_triples)
    N = Vector{Float64}(undef, nBasis)
    ΔN = Matrix{Float64}(undef, nBasis, 3)
    Δ2N = Matrix{Float64}(undef, nBasis, 6)

    for (i, (iξ, iη, iζ)) in enumerate(node_triples)
        N[i] = Nξ[iξ] * Nη[iη] * Nζ[iζ]
        ΔN[i, 1] = ΔNξ[iξ] * Nη[iη] * Nζ[iζ]
        ΔN[i, 2] = Nξ[iξ] * ΔNη[iη] * Nζ[iζ]
        ΔN[i, 3] = Nξ[iξ] * Nη[iη] * ΔNζ[iζ]
        Δ2N[i, 1] = Δ2Nξ[iξ] * Nη[iη] * Nζ[iζ]
        Δ2N[i, 2] = ΔNξ[iξ] * ΔNη[iη] * Nζ[iζ]
        Δ2N[i, 3] = ΔNξ[iξ] * Nη[iη] * ΔNζ[iζ]
        Δ2N[i, 4] = Nξ[iξ] * Δ2Nη[iη] * Nζ[iζ]
        Δ2N[i, 5] = Nξ[iξ] * ΔNη[iη] * ΔNζ[iζ]
        Δ2N[i, 6] = Nξ[iξ] * Nη[iη] * Δ2Nζ[iζ]
    end

    return N, ΔN, Δ2N
end

"""
    _second_derivative_pairs(ndim)

Return the `(i, j)` index pairs used to map second-derivative columns for a
problem dimension.
"""
function _second_derivative_pairs(ndim::Int)
    if ndim == 1
        return ((1, 1),)
    elseif ndim == 2
        return ((1, 1), (1, 2), (2, 2))
    elseif ndim == 3
        return ((1, 1), (1, 2), (1, 3), (2, 2), (2, 3), (3, 3))
    end

    error("Unsupported derivative dimension: $ndim")
end

"""
    _rational_basis(Be, ΔBe, Δ2Be, We)

Convert B-spline basis values and derivatives to rational (NURBS) basis values,
including optional second derivatives.
"""
function _rational_basis(Be, ΔBe, Δ2Be, We)
    we = Be' * We
    Δwe = ΔBe' * We

    Re = (We .* Be) / we
    ΔRe = (We .* ΔBe) / we - ((We .* Be) * Δwe') / we^2

    if Δ2Be isa AbstractVector
        Δ2we = Δ2Be' * We
        Δ2Re = (We .* Δ2Be) / we - 2 * (We .* ΔBe) * Δwe / we^2 - ((We .* Be) * Δ2we) / we^2 + 2 * (We .* Be) * (Δwe^2) / we^3
        return Re, ΔRe, Δ2Re
    end

    pairs = _second_derivative_pairs(size(ΔBe, 2))
    Δ2Re = Matrix{Float64}(undef, length(Be), size(Δ2Be, 2))

    for (i, (a, b)) in enumerate(pairs)
        Δ2we = Δ2Be[:, i]' * We
        Δ2Re[:, i] = (We .* Δ2Be[:, i]) / we -
            (((We .* ΔBe[:, a]) * Δwe[b]) .+ ((We .* ΔBe[:, b]) * Δwe[a]) .+ ((We .* Be) * Δ2we)) / we^2 +
            2 .* (We .* Be) * (Δwe[a] * Δwe[b]) / we^3
    end

    return Re, ΔRe, Δ2Re
end

"""
    basis_function(ξ, η=nothing, ζ=nothing, FunctionClass="Q2"; second_derivatives=false)

Evaluate reference-element basis functions.

# Arguments
- `ξ`, `η`, `ζ`: parent-space coordinates (typically in `[-1, 1]`).
- `FunctionClass::String`: basis family (`"Q1"`, `"Q2"`, `"S1"`).
- `second_derivatives::Bool`: when `true`, include second derivatives in the return tuple.

# Returns
- `N`: basis function values.
- `ΔN`: first derivatives.
- `Δ2N` (optional): second derivatives.

Second-derivative ordering:
- 1D: `[∂²/∂ξ²]`
- 2D: `[∂²/∂ξ², ∂²/(∂ξ∂η), ∂²/∂η²]`
- 3D: `[∂²/∂ξ², ∂²/(∂ξ∂η), ∂²/(∂ξ∂ζ), ∂²/∂η², ∂²/(∂η∂ζ), ∂²/∂ζ²]`
"""
function basis_function(ξ::Float64, η=nothing, ζ=nothing, FunctionClass::String = "Q2"; second_derivatives::Bool=false)
    N, ΔN, Δ2N = _basis_1d_components(ξ, FunctionClass)
    ΔN_out, Δ2N_out = _format_1d_derivatives(FunctionClass, ΔN, Δ2N)

    if second_derivatives
        return N, ΔN_out, Δ2N_out
    end

    return N, ΔN_out
end

function basis_function(ξ::Float64, η::Float64, ζ=nothing, FunctionClass::String = "Q2"; second_derivatives::Bool=false)
    Nξ, ΔNξ, Δ2Nξ = _basis_1d_components(ξ, FunctionClass)
    Nη, ΔNη, Δ2Nη = _basis_1d_components(η, FunctionClass)

    if FunctionClass == "Q1" || FunctionClass == "S1"
        N, ΔN, Δ2N = _tensor_basis_2d(Nξ, ΔNξ, Δ2Nξ, Nη, ΔNη, Δ2Nη, _Q1_2D_NODE_PAIRS)
    elseif FunctionClass == "Q2"
        N, ΔN, Δ2N = _tensor_basis_2d(Nξ, ΔNξ, Δ2Nξ, Nη, ΔNη, Δ2Nη, _Q2_2D_NODE_PAIRS)
    else
        error("Unsupported FunctionClass: $FunctionClass")
    end

    if second_derivatives
        return N, ΔN, Δ2N
    end

    return N, ΔN
end

function basis_function(ξ::Float64, η::Float64, ζ::Float64, FunctionClass::String = "Q2"; second_derivatives::Bool=false)
    Nξ, ΔNξ, Δ2Nξ = _basis_1d_components(ξ, FunctionClass)
    Nη, ΔNη, Δ2Nη = _basis_1d_components(η, FunctionClass)
    Nζ, ΔNζ, Δ2Nζ = _basis_1d_components(ζ, FunctionClass)

    if FunctionClass == "Q1" || FunctionClass == "S1"
        N, ΔN, Δ2N = _tensor_basis_3d(Nξ, ΔNξ, Δ2Nξ, Nη, ΔNη, Δ2Nη, Nζ, ΔNζ, Δ2Nζ, _Q1_3D_NODE_TRIPLES)
    elseif FunctionClass == "Q2"
        N, ΔN, Δ2N = _tensor_basis_3d(Nξ, ΔNξ, Δ2Nξ, Nη, ΔNη, Δ2Nη, Nζ, ΔNζ, Δ2Nζ, _Q2_3D_NODE_TRIPLES)
    else
        error("Unsupported FunctionClass: $FunctionClass")
    end

    if second_derivatives
        return N, ΔN, Δ2N
    end

    return N, ΔN
end

"""
    basis_function(ξ, Ce, We, FunctionClass="S2"; second_derivatives=false)

Evaluate 1D NURBS basis functions and derivatives from extraction operator `Ce`
and weights `We`.
"""
function basis_function(ξ::Float64, Ce::Matrix{Float64}, We::Vector{Float64}, FunctionClass::String = "S2"; second_derivatives::Bool=false)
    N, ΔN, Δ2N = _basis_1d_components(ξ, "Q2")
    Be = Ce * N
    ΔBe = Ce * ΔN
    Δ2Be = Ce * Δ2N

    Re, ΔRe, Δ2Re = _rational_basis(Be, ΔBe, Δ2Be, We)

    if second_derivatives
        return Re, ΔRe, Δ2Re
    end

    return Re, ΔRe
end

"""
    basis_function(ξ, η, Ce, We, FunctionClass="S2"; second_derivatives=false)

Evaluate 2D NURBS basis functions and derivatives from extraction operator `Ce`
and weights `We`.
"""
function basis_function(ξ::Float64, η::Float64, Ce::Matrix{Float64}, We::Vector{Float64}, FunctionClass::String = "S2"; second_derivatives::Bool=false)
    if second_derivatives
        N, ΔN, Δ2N = basis_function(ξ, η, nothing, "Q2"; second_derivatives=true)
        Be = Ce * N
        ΔBe = Ce * ΔN
        Δ2Be = Ce * Δ2N
        Re, ΔRe, Δ2Re = _rational_basis(Be, ΔBe, Δ2Be, We)
        return Re, ΔRe, Δ2Re
    end

    N, ΔN = basis_function(ξ, η, nothing, "Q2")

    # Compute the Bspline basis
    Be = Ce * N
    ΔBe = Ce * ΔN

    we = Be' * We
    Δwe = ΔBe' * We

    # Compute the NURBS basis
    Re = (We .* Be) / we
    ΔRe = (We .* ΔBe) / we - ((We .* Be) * Δwe') / we^2

    return Re, ΔRe
end

"""
    basis_function(ξ, η, ζ, Ce, We, FunctionClass="S2"; second_derivatives=false)

Evaluate 3D NURBS basis functions and derivatives from extraction operator `Ce`
and weights `We`.
"""
function basis_function(ξ::Float64, η::Float64, ζ::Float64, Ce::Matrix{Float64}, We::Vector{Float64}, FunctionClass::String = "S2"; second_derivatives::Bool=false)
    if second_derivatives
        N, ΔN, Δ2N = basis_function(ξ, η, ζ, "Q2"; second_derivatives=true)
        Be = Ce * N
        ΔBe = Ce * ΔN
        Δ2Be = Ce * Δ2N
        Re, ΔRe, Δ2Re = _rational_basis(Be, ΔBe, Δ2Be, We)
        return Re, ΔRe, Δ2Re
    end

    N, ΔN = basis_function(ξ, η, ζ, "Q2")

    # Compute the Bspline basis
    Be = Ce * N
    ΔBe = Ce * ΔN

    we = Be' * We
    Δwe = ΔBe' * We

    # Compute the NURBS basis
    Re = (We .* Be) / we
    ΔRe = (We .* ΔBe) / we - ((We .* Be) * Δwe') / we^2

    # return Be, ΔBe
    return Re, ΔRe
end

"""
    get_basis_volume_functions(mdl::Stokes)

Pre-compute basis functions at all Gauss points in volume elements for assembly.

# Arguments:
- `mdl::Stokes` : Material model

# Returns:
- `N_u_gp` : Velocity basis function values at each volume GP
- `ΔN_u_gp` : Velocity basis function derivatives at each volume GP
- `N_p_gp` : Pressure basis function values at each volume GP
- `ΔN_p_gp` : Pressure basis function derivatives at each volume GP
- `N_x_gp` : Geometry basis function values at each volume GP
- `ΔN_x_gp` : Geometry basis function derivatives at each volume GP
- `wpoints` : Quadrature weights for volume integration
"""
function get_basis_volume_functions(mdl::Stokes)
    
    @unpack ndim, ne = mdl
    FunctionClass_u = mdl.mesh_u.FunctionClass
    FunctionClass_p = mdl.mesh_p.FunctionClass
    FunctionClass_x = mdl.mesh_x.FunctionClass

    ndim_cached = ndim
    FunctionClass_u_cached = FunctionClass_u
    FunctionClass_p_cached = FunctionClass_p
    FunctionClass_x_cached = FunctionClass_x

    if ndim_cached == 1
        # gaussian quadrature points for the element [-1,1] 
        ξ, w_ξ = gaussian_quadrature(-1,1,3)
    
        wpoints = [w_ξ[1], w_ξ[2]]
        
        x = [ξ[1], ξ[2]]
    elseif ndim_cached == 2
        # gaussian quadrature points for the element [-1,1]x[-1,1] 
        ξ, w_ξ = gaussian_quadrature(-1,1,3)
        η, w_η = gaussian_quadrature(-1,1,3)

        # Pre-allocate arrays for efficiency (faster than push!)
        npts = size(ξ,1) * size(η,1)
        x = Vector{Float64}(undef, npts)
        y = Vector{Float64}(undef, npts)
        wpoints = Vector{Float64}(undef, npts)
        
        idx = 1
        size_η = size(η,1)
        size_ξ = size(ξ,1)
        for j in 1:size_η
            for i in 1:size_ξ
                x[idx] = ξ[i]
                y[idx] = η[j]
                wpoints[idx] = w_ξ[i]*w_η[j]
                idx += 1
            end
        end

    elseif ndim_cached == 3
        # gaussian quadrature points for the element [-1,1]x[-1,1]x[-1,1] 
        ξ, w_ξ = gaussian_quadrature(-1,1,3)
        η, w_η = gaussian_quadrature(-1,1,3)
        ζ, w_ζ = gaussian_quadrature(-1,1,3)

        # Pre-allocate arrays for efficiency (faster than push!)
        npts = size(ξ,1) * size(η,1) * size(ζ,1)
        x = Vector{Float64}(undef, npts)
        y = Vector{Float64}(undef, npts)
        z = Vector{Float64}(undef, npts)
        wpoints = Vector{Float64}(undef, npts)
        
        idx = 1
        size_ζ = size(ζ,1)
        size_η = size(η,1)
        size_ξ = size(ξ,1)
        for k::Int in 1:size_ζ
            for j::Int in 1:size_η
                for i::Int in 1:size_ξ
                    x[idx] = ξ[i]
                    y[idx] = η[j]
                    z[idx] = ζ[k]
                    wpoints[idx] = w_ξ[i]*w_η[j]*w_ζ[k]
                    idx += 1
                end
            end
        end
    end

    gpiter = 1:length(wpoints)  # iterator for integration loop
    
    # Pre-allocate vectors for storing basis functions at all Gauss points
    N_u_gp = Vector{Vector{Float64}}(undef, length(wpoints))
    ΔN_u_gp = Vector{Matrix{Float64}}(undef, length(wpoints))
    N_p_gp = Vector{Vector{Float64}}(undef, length(wpoints))
    ΔN_p_gp = Vector{Matrix{Float64}}(undef, length(wpoints))
    N_x_gp = Vector{Vector{Float64}}(undef, length(wpoints))
    ΔN_x_gp = Vector{Matrix{Float64}}(undef, length(wpoints))
    
    # Compute and store basis functions at all Gauss points
    for gp::Int in gpiter
        if ndim_cached == 1
            N_u_gp[gp], ΔN_u_gp[gp] = basis_function(x[gp], nothing, nothing, FunctionClass_u_cached)
            N_p_gp[gp], ΔN_p_gp[gp] = basis_function(x[gp], nothing, nothing, FunctionClass_p_cached)
            N_x_gp[gp], ΔN_x_gp[gp] = basis_function(x[gp], nothing, nothing, FunctionClass_x_cached)
        elseif ndim_cached == 2
            N_u_gp[gp], ΔN_u_gp[gp] = basis_function(x[gp], y[gp], nothing, FunctionClass_u_cached)
            N_p_gp[gp], ΔN_p_gp[gp] = basis_function(x[gp], y[gp], nothing, FunctionClass_p_cached)
            N_x_gp[gp], ΔN_x_gp[gp] = basis_function(x[gp], y[gp], nothing, FunctionClass_x_cached)
        elseif ndim_cached == 3
            N_u_gp[gp], ΔN_u_gp[gp] = basis_function(x[gp], y[gp], z[gp], FunctionClass_u_cached)
            N_p_gp[gp], ΔN_p_gp[gp] = basis_function(x[gp], y[gp], z[gp], FunctionClass_p_cached)
            N_x_gp[gp], ΔN_x_gp[gp] = basis_function(x[gp], y[gp], z[gp], FunctionClass_x_cached)
        end
    end
    
    return N_u_gp, ΔN_u_gp, N_p_gp, ΔN_p_gp, N_x_gp, ΔN_x_gp, wpoints
end

"""
    get_surface_basis_functions(mdl::Stokes)

Pre-compute basis functions at all Gauss points on surfaces for boundary condition assembly.
For 2D problems: computes 1D basis functions for edge integration.
For 3D problems: computes 2D basis functions for surface integration.

# Arguments:
- `mdl::Stokes` : Material model

# Returns:
- `N_u_surf_gp` : Velocity basis function values at each surface GP
- `ΔN_u_surf_gp` : Velocity basis function derivatives at each surface GP
- `wpoints` : Quadrature weights for surface integration
"""
function get_surface_basis_functions(mdl::Stokes)
    
    @unpack ndim = mdl
    FunctionClass_u = mdl.mesh_u.FunctionClass

    ndim_cached = ndim
    FunctionClass_u_cached = FunctionClass_u

    if ndim_cached == 2
        # 1D Gauss quadrature for edge integration
        ξ, w_ξ = gaussian_quadrature(-1, 1)
        x = [ξ[1], ξ[2]]
        wpoints = [w_ξ[1], w_ξ[2]]
        
    elseif ndim_cached == 3
        # 2D Gauss quadrature for surface integration
        ξ, w_ξ = gaussian_quadrature(-1, 1)
        η, w_η = gaussian_quadrature(-1, 1)
        
        x = [ξ[1], ξ[2], ξ[2], ξ[1]]
        y = [η[1], η[1], η[2], η[2]]
        wpoints = [w_ξ[1]*w_η[1], w_ξ[2]*w_η[1], w_ξ[2]*w_η[2], w_ξ[1]*w_η[2]]
    else
        error("Surface basis functions only supported for ndim >= 2")
    end
    
    # Pre-allocate and compute basis functions at all surface Gauss points
    N_u_surf_gp = Vector{Vector{Float64}}(undef, length(wpoints))
    ΔN_u_surf_gp = Vector{Matrix{Float64}}(undef, length(wpoints))
    
    for gp in 1:length(wpoints)
        if ndim_cached == 2
            N_u_surf_gp[gp], ΔN_u_surf_gp[gp] = basis_function(x[gp], nothing, nothing, FunctionClass_u_cached)
        elseif ndim_cached == 3
            N_u_surf_gp[gp], ΔN_u_surf_gp[gp] = basis_function(x[gp], y[gp], nothing, FunctionClass_u_cached)
        end
    end
    
    return N_u_surf_gp, ΔN_u_surf_gp, wpoints
end