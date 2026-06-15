using LinearAlgebra
using SparseArrays
using Parameters

"""
    BasisFunctionCache

Structure to store pre-computed basis functions for volume and surface integration.

# Fields:
- `bf_volume`  : `(N_u_gp, ΔN_u_gp, N_p_gp, ΔN_p_gp, N_x_gp, ΔN_x_gp, wpoints)` — volume
- `bf_surface` : `(N_u_surf_gp, ΔN_u_surf_gp, N_p_surf_gp, ΔN_p_surf_gp, N_x_surf_gp, ΔN_x_surf_gp, wpoints)` — surface
"""
struct BasisFunctionCache
    bf_volume::NTuple{7, Any}
    bf_surface::NTuple{7, Any}

    function BasisFunctionCache(bf_vol::NTuple{7, Any}, bf_surf::NTuple{7, Any})
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
    bf_volume = get_volume_basis_functions(mdl)
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
function gaussian_quadrature(a::Int64,b::Int64,n_gauss_pts::Int64=2)

    if n_gauss_pts == 2
        ξ = [-(b-a)/(2*sqrt(3))+(b+a)/2, (b-a)/(2*sqrt(3))+(b+a)/2]
        w = [(b-a)/2, (b-a)/2]
    elseif n_gauss_pts == 3
        ξ = [-(b-a)/(2*sqrt(5/3))+(b+a)/2, 0, (b-a)/(2*sqrt(5/3))+(b+a)/2]
        w = [(b-a)/2*5/9, (b-a)/2*8/9, (b-a)/2*5/9]
    else
        error("Unsupported n_gauss_pts=$n_gauss_pts for gaussian_quadrature (supported: 2, 3)")
    end
    return ξ, w
end

function _tri_quadrature(n_gauss_pts::Int64=3)
    if n_gauss_pts == 1
        λ1 = [1/3]; λ2 = [1/3]; λ3 = [1/3]; w = [0.5]
    elseif n_gauss_pts == 3
        λ1 = [0.5, 0.0, 0.5]
        λ2 = [0.5, 0.5, 0.0]
        λ3 = [0.0, 0.5, 0.5]
        w  = [1/6, 1/6, 1/6]
    elseif n_gauss_pts == 4
        λ1 = [3/5, 0.2, 0.2, 1/3]
        λ2 = [0.2, 3/5, 0.2, 1/3]
        λ3 = [0.2, 0.2, 3/5, 1/3]
        w  = [25/48, 25/48, 25/48, -27/48]
    else
        error("Unsupported n_gauss_pts=$n_gauss_pts for _tri_quadrature (supported: 1, 3, 4)")
    end
    return λ1, λ2, λ3, w
end

function _tet_quadrature(n_gauss_pts::Int64=4)
    if n_gauss_pts == 1
        λ1 = [1/4]; λ2 = [1/4]; λ3 = [1/4]; λ4 = [1/4]; w = [1.0]
    elseif n_gauss_pts == 4
        λ1 = [0.58541020168919, 0.13819660, 0.13819660, 0.13819660]
        λ2 = [0.13819660, 0.58541020168919, 0.13819660, 0.13819660]
        λ3 = [0.13819660, 0.13819660, 0.58541020168919, 0.13819660]
        λ4 = [0.13819660, 0.13819660, 0.13819660, 0.58541020168919]
        w  = [0.25, 0.25, 0.25, 0.25]
    elseif n_gauss_pts == 5
        λ1 = [1/2, 1/6, 1/6, 1/6, 1/4]
        λ2 = [1/6, 1/2, 1/6, 1/6, 1/4]
        λ3 = [1/6, 1/6, 1/2, 1/6, 1/4]
        λ4 = [1/6, 1/6, 1/6, 1/2, 1/4]
        w  = [9/20, 9/20, 9/20, 9/20, -4/5]
    else
        error("Unsupported n_gauss_pts=$n_gauss_pts for _tet_quadrature (supported: 1, 4, 5)")
    end
    return λ1, λ2, λ3, λ4, w
end

"""
    get_quadrature(element_shape, n_gauss_pts) -> (points, weights)

Return reference-element quadrature points and weights for any supported element topology.

`points` is a `Vector` of coordinate tuples — one tuple per Gauss point. Splat a tuple
directly into `basis_function`: `basis_function(points[gp]..., element_shape, order)`.
`weights` is a `Vector{Float64}`.

Adding a new element type requires only a new branch here; no other function needs changing.

# Arguments
- `element_shape::Symbol`: `:Hex`, `:Tet`, `:Quad`, `:Tri`, `:Line`
- `n_gauss_pts::Int`: number of Gauss points (shape-dependent; see individual quadrature functions)

# Supported `n_gauss_pts` per shape
| Shape  | Supported values |
|--------|-----------------|
| `:Hex` | 2, 3 (per direction; total = n³) |
| `:Quad`| 2, 3 (per direction; total = n²) |
| `:Line`| 2, 3 |
| `:Tet` | 1, 4, 5 |
| `:Tri` | 1, 3, 4 |
"""
function get_quadrature(element_shape::Symbol, n_gauss_pts::Int)
    if element_shape == :Tet
        λ1, λ2, λ3, λ4, w = _tet_quadrature(n_gauss_pts)
        pts = NTuple{4,Float64}[(λ1[i], λ2[i], λ3[i], λ4[i]) for i in eachindex(w)]
        return pts, Vector{Float64}(w)
    elseif element_shape == :Tri
        λ1, λ2, λ3, w = _tri_quadrature(n_gauss_pts)
        pts = NTuple{3,Float64}[(λ1[i], λ2[i], λ3[i]) for i in eachindex(w)]
        return pts, Vector{Float64}(w)
    elseif element_shape == :Hex
        ξ, wξ = gaussian_quadrature(-1, 1, n_gauss_pts)
        η, wη = gaussian_quadrature(-1, 1, n_gauss_pts)
        ζ, wζ = gaussian_quadrature(-1, 1, n_gauss_pts)
        pts = NTuple{3,Float64}[]; wts = Float64[]
        for k in eachindex(ζ), j in eachindex(η), i in eachindex(ξ)
            push!(pts, (ξ[i], η[j], ζ[k]))
            push!(wts, wξ[i] * wη[j] * wζ[k])
        end
        return pts, wts
    elseif element_shape == :Quad
        ξ, wξ = gaussian_quadrature(-1, 1, n_gauss_pts)
        η, wη = gaussian_quadrature(-1, 1, n_gauss_pts)
        pts = NTuple{2,Float64}[]; wts = Float64[]
        for j in eachindex(η), i in eachindex(ξ)
            push!(pts, (ξ[i], η[j]))
            push!(wts, wξ[i] * wη[j])
        end
        return pts, wts
    elseif element_shape == :Line
        ξ, w = gaussian_quadrature(-1, 1, n_gauss_pts)
        return NTuple{1,Float64}[(ξ[i],) for i in eachindex(w)], Vector{Float64}(w)
    else
        error("Unsupported element_shape=$element_shape for get_quadrature")
    end
end

"""
    _basis_1d_components(ξ, basis_order)

Return 1D basis values, first derivatives, and second derivatives for the
requested polynomial order.
"""
function _basis_1d_components(ξ::Float64, basis_order::Int)
    if basis_order == 1
        N = [0.5 - 0.5 * ξ, 0.5 + 0.5 * ξ]
        ΔN = [-0.5, 0.5]
        Δ2N = [0.0, 0.0]
    elseif basis_order == 2
        N = [-(1 - ξ) * ξ / 2,
              ξ * (1 + ξ) / 2,
             (1 - ξ) * (1 + ξ)]

        ΔN = [-(1 - 2 * ξ) / 2,
               (1 + 2 * ξ) / 2,
              -2 * ξ]

        Δ2N = [1.0, 1.0, -2.0]
    else
        error("Unsupported basis_order: $basis_order")
    end

    return N, ΔN, Δ2N
end

"""
    _format_1d_derivatives(basis_order, ΔN, Δ2N)

Normalize derivative return shapes for 1D basis functions so legacy callers
keep receiving matrix output for order 1 and vector output for order 2.
"""
function _format_1d_derivatives(basis_order::Int, ΔN, Δ2N)
    if basis_order == 1
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

function _basis_tri(λ1::Float64, λ2::Float64, λ3::Float64, basis_order::Int)

    ∇λ = [-1.0 -1/sqrt(3);
            1.0 -1/sqrt(3);
            0.0 2/sqrt(3)]

    if basis_order == 1
        N    = [λ1, λ2, λ3]
        ΔN   = copy(∇λ)        # (3×2): ∂N_i/∂x_a = ∂λ_i/∂x_a
        Δ2N  = zeros(3, 3)     # second derivatives vanish for linear basis
    elseif basis_order == 2
        N = [λ1*(2*λ1-1), λ2*(2*λ2-1), λ3*(2*λ3-1), 4*λ1*λ2, 4*λ2*λ3, 4*λ3*λ1]
        # ∂Nᵢ/∂λₖ: corner Nᵢ=λᵢ(2λᵢ-1) → 4λᵢ-1; edge Nᵢⱼ=4λᵢλⱼ → product rule
        ∂N∂λ = [4*λ1-1  0         0        ;
                0         4*λ2-1  0        ;
                0         0         4*λ3-1 ;
                4*λ2    4*λ1    0        ;
                0         4*λ3    4*λ2   ;
                4*λ3    0         4*λ1   ]
        ΔN = ∂N∂λ * ∇λ

        # Second derivatives (constant since λ are linear):
        #   corner i:    ∂²Nᵢ/∂xₐ∂xᵦ = 4*(∂λᵢ/∂xₐ)*(∂λᵢ/∂xᵦ)
        #   edge (i,j):  ∂²Nᵢⱼ/∂xₐ∂xᵦ = 4*((∂λᵢ/∂xₐ)(∂λⱼ/∂xᵦ) + (∂λⱼ/∂xₐ)(∂λᵢ/∂xᵦ))
        # Columns: [ξξ, ξη, ηη]
        pairs_ab   = ((1,1),(1,2),(2,2))
        edge_pairs = ((1,2),(2,3),(3,1))
        Δ2N = Matrix{Float64}(undef, 6, 3)
        for i in 1:3
            g = ∇λ[i,:]
            for (col,(a,b)) in enumerate(pairs_ab)
                Δ2N[i,col] = 4*g[a]*g[b]
            end
        end
        for (k,(i,j)) in enumerate(edge_pairs)
            gi, gj = ∇λ[i,:], ∇λ[j,:]
            for (col,(a,b)) in enumerate(pairs_ab)
                Δ2N[3+k,col] = 4*(gi[a]*gj[b] + gj[a]*gi[b])
            end
        end
    else
        error("Unsupported basis_order: $basis_order")
    end

    return N, ΔN, Δ2N
end

function _basis_tet(λ1::Float64, λ2::Float64, λ3::Float64, λ4::Float64, basis_order::Int)
    
    ∇λ = [-1.0 -1/sqrt(3) -1/sqrt(6);
             1.0 -1/sqrt(3) -1/sqrt(6);
             0.0 2/sqrt(3) -1/sqrt(6);
             0.0 0.0 sqrt(3/2)]

    if basis_order == 1
        N    = [λ1, λ2, λ3, λ4]
        ΔN   = copy(∇λ)        # (4×3): ∂N_i/∂x_a = ∂λ_i/∂x_a
        Δ2N  = zeros(4, 6)     # second derivatives vanish for linear basis
    elseif basis_order == 2
        N = [λ1*(2*λ1-1), λ2*(2*λ2-1), λ3*(2*λ3-1), λ4*(2*λ4-1),
             4*λ1*λ2, 4*λ2*λ3, 4*λ3*λ1, 4*λ1*λ4, 4*λ2*λ4, 4*λ3*λ4]

        # ∂Nᵢ/∂λₖ: corner Nᵢ=λᵢ(2λᵢ-1) → 4λᵢ-1; edge Nᵢⱼ=4λᵢλⱼ → product rule
        ∂N∂λ = [4*λ1-1  0         0         0        ;
                0         4*λ2-1  0         0        ;
                0         0         4*λ3-1  0        ;
                0         0         0         4*λ4-1 ;
                4*λ2    4*λ1    0         0        ;
                0         4*λ3    4*λ2    0        ;
                4*λ3    0         4*λ1    0        ;
                4*λ4    0         0         4*λ1   ;
                0         4*λ4    0         4*λ2   ;
                0         0         4*λ4    4*λ3   ]
        ΔN = ∂N∂λ * ∇λ

        # Second derivatives (constant since λ are linear):
        #   corner i:    ∂²Nᵢ/∂xₐ∂xᵦ = 4*(∂λᵢ/∂xₐ)*(∂λᵢ/∂xᵦ)
        #   edge (i,j):  ∂²Nᵢⱼ/∂xₐ∂xᵦ = 4*((∂λᵢ/∂xₐ)(∂λⱼ/∂xᵦ) + (∂λⱼ/∂xₐ)(∂λᵢ/∂xᵦ))
        # Columns: [ξξ, ξη, ξζ, ηη, ηζ, ζζ]
        pairs_ab   = ((1,1),(1,2),(1,3),(2,2),(2,3),(3,3))
        edge_pairs = ((1,2),(2,3),(3,1),(1,4),(2,4),(3,4))
        Δ2N = Matrix{Float64}(undef, 10, 6)
        for i in 1:4
            g = ∇λ[i,:]
            for (col,(a,b)) in enumerate(pairs_ab)
                Δ2N[i,col] = 4*g[a]*g[b]
            end
        end
        for (k,(i,j)) in enumerate(edge_pairs)
            gi, gj = ∇λ[i,:], ∇λ[j,:]
            for (col,(a,b)) in enumerate(pairs_ab)
                Δ2N[4+k,col] = 4*(gi[a]*gj[b] + gj[a]*gi[b])
            end
        end
    else
        error("Unsupported basis_order: $basis_order")
    end

    return N, ΔN, Δ2N
end

"""
    basis_function(ξ, η=nothing, ζ=nothing, element_shape=:Hex, basis_order=1; second_derivatives=false)

Evaluate reference-element basis functions.

# Arguments
- `ξ`, `η`, `ζ`: parent-space coordinates (typically in `[-1, 1]`).
- `element_shape::Symbol`: element topology (`:Hex`, `:Quad`, `:Line`, `:Tet`).
- `basis_order::Int`: polynomial order (1, 2, or 3).
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
function basis_function(ξ::Float64, η=nothing, ζ=nothing, element_shape::Symbol=:Line, basis_order::Int=1; second_derivatives::Bool=false)
    N, ΔN, Δ2N = _basis_1d_components(ξ, basis_order)
    ΔN_out, Δ2N_out = _format_1d_derivatives(basis_order, ΔN, Δ2N)

    if second_derivatives
        return N, ΔN_out, Δ2N_out
    end

    return N, ΔN_out
end

function basis_function(ξ::Float64, η::Float64, ζ=nothing, element_shape::Symbol=:Quad, basis_order::Int=1; second_derivatives::Bool=false)

    element_shape == :Tri && error(":Tri requires all 3 barycentric coords — use the 3-argument overload: basis_function(λ1, λ2, λ3, :Tri, order)")

    Nξ, ΔNξ, Δ2Nξ = _basis_1d_components(ξ, basis_order)
    Nη, ΔNη, Δ2Nη = _basis_1d_components(η, basis_order)
    if basis_order == 1
        N, ΔN, Δ2N = _tensor_basis_2d(Nξ, ΔNξ, Δ2Nξ, Nη, ΔNη, Δ2Nη, _Q1_2D_NODE_PAIRS)
    elseif basis_order == 2
        N, ΔN, Δ2N = _tensor_basis_2d(Nξ, ΔNξ, Δ2Nξ, Nη, ΔNη, Δ2Nη, _Q2_2D_NODE_PAIRS)
    else
        error("Unsupported basis_order: $basis_order for element_shape $element_shape")
    end

    if second_derivatives
        return N, ΔN, Δ2N
    end

    return N, ΔN
end

function basis_function(ξ::Float64, η::Float64, ζ::Float64, element_shape::Symbol=:Hex, basis_order::Int=1; second_derivatives::Bool=false)
    if element_shape == :Tri
        N, ΔN, Δ2N = _basis_tri(ξ, η, ζ, basis_order)
    elseif element_shape == :Hex
        Nξ, ΔNξ, Δ2Nξ = _basis_1d_components(ξ, basis_order)
        Nη, ΔNη, Δ2Nη = _basis_1d_components(η, basis_order)
        Nζ, ΔNζ, Δ2Nζ = _basis_1d_components(ζ, basis_order)
        if basis_order == 1
            N, ΔN, Δ2N = _tensor_basis_3d(Nξ, ΔNξ, Δ2Nξ, Nη, ΔNη, Δ2Nη, Nζ, ΔNζ, Δ2Nζ, _Q1_3D_NODE_TRIPLES)
        elseif basis_order == 2
            N, ΔN, Δ2N = _tensor_basis_3d(Nξ, ΔNξ, Δ2Nξ, Nη, ΔNη, Δ2Nη, Nζ, ΔNζ, Δ2Nζ, _Q2_3D_NODE_TRIPLES)
        else
            error("Unsupported basis_order: $basis_order for element_shape $element_shape")
        end
    else
        error("Unsupported element_shape: $element_shape — use the 4-argument overload for :Tet")
    end

    if second_derivatives
        return N, ΔN, Δ2N
    end

    return N, ΔN
end

function basis_function(ξ::Float64, η::Float64, ζ::Float64, θ::Float64, element_shape::Symbol=:Tet, basis_order::Int=1; second_derivatives::Bool=false)

    if element_shape == :Tet
        N, ΔN, Δ2N = _basis_tet(ξ, η, ζ, θ, basis_order)
    else
        error("Unsupported basis_order: $basis_order for element_shape $element_shape")
    end

    if second_derivatives
        return N, ΔN, Δ2N
    end

    return N, ΔN
end

"""
    get_volume_basis_functions(mdl::Stokes)

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
function get_volume_basis_functions(mdl::Stokes)

    @unpack ne = mdl
    element_shape_u = mdl.mesh_u.volume_element_shape
    basis_order_u   = mdl.mesh_u.basis_order
    element_shape_p = mdl.mesh_p.volume_element_shape
    basis_order_p   = mdl.mesh_p.basis_order
    element_shape_x = mdl.mesh_x.volume_element_shape
    basis_order_x   = mdl.mesh_x.basis_order

    element_shape_u_cached = element_shape_u
    element_shape_p_cached = element_shape_p
    element_shape_x_cached = element_shape_x

    # Quadrature: driven by u mesh topology; p and x must share the same family
    n_gp = element_shape_u_cached == :Tet ? 4 : 3
    points, wpoints = get_quadrature(element_shape_u_cached, n_gp)

    N_u_gp  = Vector{Vector{Float64}}(undef, length(wpoints))
    ΔN_u_gp = Vector{Matrix{Float64}}(undef, length(wpoints))
    N_p_gp  = Vector{Vector{Float64}}(undef, length(wpoints))
    ΔN_p_gp = Vector{Matrix{Float64}}(undef, length(wpoints))
    N_x_gp  = Vector{Vector{Float64}}(undef, length(wpoints))
    ΔN_x_gp = Vector{Matrix{Float64}}(undef, length(wpoints))

    for gp::Int in eachindex(wpoints)
        pt = points[gp]
        N_u_gp[gp],  ΔN_u_gp[gp]  = basis_function(pt..., element_shape_u_cached, basis_order_u)
        N_p_gp[gp],  ΔN_p_gp[gp]  = basis_function(pt..., element_shape_p_cached, basis_order_p)
        N_x_gp[gp],  ΔN_x_gp[gp]  = basis_function(pt..., element_shape_x_cached, basis_order_x)
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

    surf_shape_u = mdl.mesh_u.surface_element_shape
    basis_order_u = mdl.mesh_u.basis_order
    surf_shape_p = mdl.mesh_p.surface_element_shape
    basis_order_p = mdl.mesh_p.basis_order
    surf_shape_x = mdl.mesh_x.surface_element_shape
    basis_order_x = mdl.mesh_x.basis_order

    n_gp = surf_shape_u == :Tri ? 3 : 2
    points, wpoints = get_quadrature(surf_shape_u, n_gp)

    N_u_surf_gp  = Vector{Vector{Float64}}(undef, length(wpoints))
    ΔN_u_surf_gp = Vector{Matrix{Float64}}(undef, length(wpoints))
    N_p_surf_gp  = Vector{Vector{Float64}}(undef, length(wpoints))
    ΔN_p_surf_gp = Vector{Matrix{Float64}}(undef, length(wpoints))
    N_x_surf_gp  = Vector{Vector{Float64}}(undef, length(wpoints))
    ΔN_x_surf_gp = Vector{Matrix{Float64}}(undef, length(wpoints))

    for gp in eachindex(wpoints)
        pt = points[gp]
        # Pad to (ξ, η_or_nothing, ζ_or_nothing) so every basis_function call
        # gets exactly 3 coordinate args regardless of surface dimensionality.
        coords = length(pt) == 1 ? (pt[1], nothing, nothing) :
                 length(pt) == 2 ? (pt[1], pt[2], nothing) :
                                   (pt[1], pt[2], pt[3])
        N_u_surf_gp[gp],  ΔN_u_surf_gp[gp]  = basis_function(coords..., surf_shape_u, basis_order_u)
        N_p_surf_gp[gp],  ΔN_p_surf_gp[gp]  = basis_function(coords..., surf_shape_p, basis_order_p)
        N_x_surf_gp[gp],  ΔN_x_surf_gp[gp]  = basis_function(coords..., surf_shape_x, basis_order_x)
    end

    return N_u_surf_gp, ΔN_u_surf_gp, N_p_surf_gp, ΔN_p_surf_gp, N_x_surf_gp, ΔN_x_surf_gp, wpoints
end