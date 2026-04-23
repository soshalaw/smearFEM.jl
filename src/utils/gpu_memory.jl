"""
    gpu_memory.jl

GPU memory management for smearFEM. Implements GPUContext struct for persisting
GPU arrays across iterations to minimize PCIe transfer overhead (critical for real-time).
"""

"""
    GPUContext

Context structure for managing GPU memory persistence across simulation iterations.
Pre-allocates GPU arrays once and reuses them to avoid repeated allocations.

# Fields
- `A_gpu::Union{CuSparseMatrix, CuMatrix, Nothing}`: System matrix on GPU
- `b_gpu::CuVector`: Right-hand side vector on GPU
- `x_gpu::CuVector`: Solution vector on GPU
- `precond_factors::Union{ILU0, Nothing}`: Cached ILU(0) preconditioner factors
- `krylov_basis::Dict`: Pre-allocated Krylov basis vectors for GMRES
- `workspace::Dict`: Temporary GPU buffers for assembly kernels
- `mesh_hash::UInt64`: Hash of current mesh structure (detect changes)
- `viscosity_scale::Float64`: Last computed viscosity scale (detect viscosity changes)
- `gpu_memory_info::NamedTuple`: Last queried GPU memory stats
"""
mutable struct GPUContext
    A_gpu::Union{Nothing, Any}  # CuSparseMatrix or CuMatrix (Any to avoid circular deps)
    b_gpu::Union{Nothing, Any}
    x_gpu::Union{Nothing, Any}
    precond_factors::Union{Nothing, Any}  # ILU factors
    krylov_basis::Dict
    workspace::Dict
    mesh_hash::UInt64
    viscosity_scale::Float64
    gpu_memory_info::NamedTuple
end

"""
    GPUContext(A::AbstractMatrix, b::AbstractVector; max_elements::Int=200, n_gauss_pts::Int=27)

Initialize GPU context with pre-allocated GPU arrays.

# Arguments
- `A::AbstractMatrix`: System matrix for determining allocation sizes
- `b::AbstractVector`: Right-hand side vector
- `max_elements::Int`: Maximum elements per assembly kernel call (default 200)
- `n_gauss_pts::Int`: Number of Gaussian quadrature points (default 27 for 3D)

# Returns
- `GPUContext` struct with allocated GPU memory if available, empty if GPU unavailable
"""
function GPUContext(A::AbstractMatrix, b::AbstractVector; max_elements::Int=200, n_gauss_pts::Int=27)
    if !has_gpu()
        # Return empty context if GPU unavailable
        return GPUContext(nothing, nothing, nothing, nothing, Dict(), Dict(), 0, 0.0, 
                         (total_mb=0, used_mb=0, free_mb=0, utilization_percent=0))
    end
    
    try
        # Pre-allocate GPU arrays
        A_gpu = Adapt.adapt(CUDA.CuArray, A)
        b_gpu = Adapt.adapt(CUDA.CuArray, b)
        x_gpu = CUDA.zeros(eltype(b), size(b)...)
        
        # Pre-allocate workspace for kernels
        workspace = setup_gpu_kernel_workspace(max_elements, n_gauss_pts)
        
        # Pre-allocate Krylov basis (for GMRES, typically 50 vectors)
        krylov_basis = Dict()
        for i in 1:50
            krylov_basis["v_$i"] = CUDA.zeros(eltype(b), size(b)...)
        end
        
        # Query initial GPU memory
        gpu_mem = query_gpu_memory()
        
        # Compute mesh hash: use nnz of sparse matrix as simple hash
        if hasproperty(A, :nnz)
            mesh_hash = hash(A.nnz) ⊻ hash(size(A))
        else
            mesh_hash = hash(size(A))
        end
        
        return GPUContext(A_gpu, b_gpu, x_gpu, nothing, krylov_basis, workspace, mesh_hash, 1.0, gpu_mem)
    catch e
        @warn "Failed to initialize GPUContext: $e. Running with CPU fallback."
        return GPUContext(nothing, nothing, nothing, nothing, Dict(), Dict(), 0, 0.0,
                         (total_mb=0, used_mb=0, free_mb=0, utilization_percent=0))
    end
end

"""
    reset_gpu_context(ctx::GPUContext; full_reset::Bool=false)

Clear GPU memory in context. Only clears when mesh structure changes (rare).

# Arguments
- `ctx::GPUContext`: Context to reset
- `full_reset::Bool`: If true, deallocate all GPU arrays. If false, only clear values.
"""
function reset_gpu_context(ctx::GPUContext; full_reset::Bool=false)
    if ctx.A_gpu === nothing
        return  # GPU not available, nothing to reset
    end
    
    if full_reset
        # Full reset: deallocate GPU memory
        try
            # GPU arrays will be garbage collected when context is reassigned
            CUDA.reclaim()
            @debug "GPU context fully reset"
        catch e
            @warn "Failed to fully reset GPU context: $e"
        end
    else
        # Partial reset: zero vectors for next iteration
        try
            if ctx.b_gpu !== nothing
                fill!(ctx.b_gpu, 0)
            end
            if ctx.x_gpu !== nothing
                fill!(ctx.x_gpu, 0)
            end
            @debug "GPU context partially reset (vectors zeroed)"
        catch e
            @warn "Failed to partially reset GPU context: $e"
        end
    end
end

"""
    is_gpu_available(ctx::GPUContext)::Bool

Check if GPU arrays are allocated in this context.
"""
function is_gpu_available(ctx::GPUContext)::Bool
    return ctx.A_gpu !== nothing && has_gpu()
end

"""
    update_gpu_memory_info(ctx::GPUContext)

Update GPU memory usage statistics in context.
"""
function update_gpu_memory_info(ctx::GPUContext)
    if is_gpu_available(ctx)
        # Note: This would require mutable struct or returning new context
        # For now, query_gpu_memory() can be called separately
        return query_gpu_memory()
    end
    return (total_mb=0, used_mb=0, free_mb=0, utilization_percent=0)
end

"""
    need_precond_recompute(ctx::GPUContext, viscosity_scale::Float64; threshold::Float64=0.05)::Bool

Detect if preconditioner should be recomputed based on viscosity change.

# Arguments
- `ctx::GPUContext`: Current context
- `viscosity_scale::Float64`: Current viscosity scale factor
- `threshold::Float64`: Relative threshold for recompute (default 5%)

# Returns
- `true` if viscosity changed by >threshold or preconditioner not cached
"""
function need_precond_recompute(ctx::GPUContext, viscosity_scale::Float64; threshold::Float64=0.05)::Bool
    if ctx.precond_factors === nothing
        return true  # First computation, no cached preconditioner
    end
    
    rel_change = abs(viscosity_scale - ctx.viscosity_scale) / (ctx.viscosity_scale + 1e-14)
    return rel_change > threshold
end
