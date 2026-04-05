using LinearAlgebra
using SparseArrays


""" 
    smearFEM.gaussian_quadrature(a,b,nGaussPoints)

Compute the nodes and weights for the Gaussian quadrature of order 2
    
# Arguments:    
- `a,b::Integer` : the limits of the integration interval
- `nGaussPoints::Integer` : number of Gauss points to be considered (2 or 3)

# Returns:    
- `ξ::Vector{Float64}{,nGaussPoints}`: nodes.
- `w::Vector{Float64}{,nGaussPoints}`: weights 
"""
function gaussian_quadrature(a::Int64,b::Int64;nGaussPoints::Int64=2)
  
    if nGaussPoints == 2
        ξ = [-(b-a)/(2*sqrt(3))+(b+a)/2, (b-a)/(2*sqrt(3))+(b+a)/2]
        w = [(b-a)/2, (b-a)/2]
    elseif nGaussPoints == 3
        ξ = [-(b-a)/(2*sqrt(5/3))+(b+a)/2, 0, (b-a)/(2*sqrt(5/3))+(b+a)/2]
        w = [(b-a)/2*5/9, (b-a)/2*8/9, (b-a)/2*5/9]
    end
    return ξ, w
end

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