"""
    gpu_utils.jl

GPU acceleration utilities for smearFEM. Provides device detection, array allocation,
and GPU-aware data transfer strategies for real-time performance.

Note: GPU support is optional. CUDA.jl and related packages are loaded dynamically.
"""

# Conditional CUDA import - try to use it, but don't require it
const CUDA_AVAILABLE = try
    import CUDA
    CUDA.functional()
catch
    false
end

# Conditional Adapt import
const ADAPT_AVAILABLE = try
    import Adapt
    true
catch
    false
end

"""
    has_gpu()::Bool

Detect if CUDA GPU is available and functional.

# Returns
- `true` if CUDA is available and at least one GPU device is present
- `false` otherwise
"""
function has_gpu()::Bool
    return CUDA_AVAILABLE
end

"""
    alloc_gpu_storage(A::AbstractMatrix, b::AbstractVector)

Pre-allocate GPU arrays for matrix A and vector b if GPU is available.
Critical for real-time performance: allocate once, reuse across iterations.

# Arguments
- `A::AbstractMatrix`: Sparse or dense system matrix
- `b::AbstractVector`: Right-hand side vector

# Returns
- Tuple `(A_gpu, b_gpu)` where arrays are on GPU if available, CPU otherwise
- Return `(A, b)` if GPU not available
"""
function alloc_gpu_storage(A::AbstractMatrix, b::AbstractVector)
    if !has_gpu()
        return A, b
    end
    
    try
        # Transfer arrays to GPU using CUDA (already imported in CUDA_AVAILABLE check)
        A_gpu = CUDA.CuArray(A)
        b_gpu = CUDA.CuArray(b)
        return A_gpu, b_gpu
    catch e
        @warn "Failed to allocate GPU storage: $e. Falling back to CPU."
        return A, b
    end
end

"""
    keep_on_gpu_strategy(use_gpu::Bool)::Bool

Determine if matrices should persist on GPU between iterations.
This is CRITICAL for real-time performance to minimize PCIe transfers.

# Arguments
- `use_gpu::Bool`: User preference for GPU acceleration

# Returns
- `true` if GPU is available and user enables it (minimize transfers)
- `false` otherwise (transfer data each call)
"""
function keep_on_gpu_strategy(use_gpu::Bool)::Bool
    return use_gpu && has_gpu()
end

"""
    setup_gpu_kernel_workspace(max_elements::Int, n_gauss_pts::Int)

Pre-allocate temporary GPU buffers for assembly kernels.
Reduces allocation overhead during assembly loops.

# Arguments
- `max_elements::Int`: Maximum number of elements to process per kernel call
- `n_gauss_pts::Int`: Number of Gaussian quadrature points

# Returns
- Dictionary of GPU arrays for temporary storage (or empty dict if GPU unavailable)
"""
function setup_gpu_kernel_workspace(max_elements::Int, n_gauss_pts::Int)
    workspace = Dict()
    
    if !has_gpu()
        return workspace
    end
    
    try
        # Pre-allocate temporary buffers for assembly
        # These will be reused, not reallocated in tight loops
        workspace["basis_vals"] = CUDA.zeros(Float64, n_gauss_pts, max_elements)
        workspace["jacobian"] = CUDA.zeros(Float64, 3, 3, max_elements)
        workspace["element_stiffness"] = CUDA.zeros(Float64, 24, 24, max_elements)  # For Q2 Stokes
        
        return workspace
    catch e
        @warn "Failed to setup GPU kernel workspace: $e"
        return workspace
    end
end

"""
    query_gpu_memory()::NamedTuple

Query GPU memory status if available.

# Returns
- NamedTuple with fields `(total_mb, used_mb, free_mb, utilization_percent)`
- `(total_mb=0, used_mb=0, free_mb=0, utilization_percent=0)` if GPU unavailable
"""
function query_gpu_memory()
    if !has_gpu()
        return (total_mb=0, used_mb=0, free_mb=0, utilization_percent=0)
    end
    
    try
        total = CUDA.total_memory()
        free = CUDA.available_memory()
        used = total - free
        utilization = 100 * (used / total)
        
        return (
            total_mb = total / 1e6,
            used_mb = used / 1e6,
            free_mb = free / 1e6,
            utilization_percent = utilization
        )
    catch
        return (total_mb=0, used_mb=0, free_mb=0, utilization_percent=0)
    end
end

"""
    log_gpu_status()

Print GPU availability and memory status to stdout.
Called at module initialization for user awareness.
"""
function log_gpu_status()
    if has_gpu()
        mem = query_gpu_memory()
        println("[GPU] CUDA available. Device: $(CUDA.device()). Memory: $(round(mem.used_mb; digits=1))MB / $(round(mem.total_mb; digits=1))MB ($(round(mem.utilization_percent; digits=1))%)")
    else
        println("[GPU] CUDA not available. Using CPU.")
    end
end

"""
    reset_gpu_memory()

Clear GPU memory cache. Use after major simulation phases.
"""
function reset_gpu_memory()
    if has_gpu()
        try
            GC.gc(false)  # Force Julia garbage collection
            CUDA.reclaim()  # Reclaim GPU memory
        catch
            # Silently ignore if CUDA reclaim fails
        end
    end
end
