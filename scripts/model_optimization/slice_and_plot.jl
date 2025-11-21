#!/usr/bin/env julia
"""
Slice a 2D cost matrix along steepest and flattest descent directions and (optionally) plot.

This script provides functions to compute gradients, get directions, sample along a line
with bilinear interpolation, and a small test that runs without plotting. To enable plotting
install `Plots` and run with `--plot` flag.
"""

using LinearAlgebra
using smearFEM

function compute_gradients(A::AbstractMatrix; dx::Real=1.0, dy::Real=1.0)
    ny, nx = size(A)
    fx = zeros(eltype(A), ny, nx)
    fy = zeros(eltype(A), ny, nx)

    for j in 2:ny-1, i in 2:nx-1
        fx[j,i] = (A[j, i+1] - A[j, i-1]) / (2dx)
        fy[j,i] = (A[j+1, i] - A[j-1, i]) / (2dy)
    end

    # boundaries: forward/backward differences
    for j in 1:ny
        fx[j,1]   = (A[j,2] - A[j,1]) / dx
        fx[j,nx]  = (A[j,nx] - A[j,nx-1]) / dx
    end
    for i in 1:nx
        fy[1,i]   = (A[2,i] - A[1,i]) / dy
        fy[ny,i]  = (A[ny,i] - A[ny-1,i]) / dy
    end

    return fx, fy
end

function directions_at(A::AbstractMatrix, row::Integer, col::Integer; dx::Real=1.0, dy::Real=1.0, fx_all=nothing, fy_all=nothing)
    if fx_all === nothing || fy_all === nothing
        fx_all, fy_all = compute_gradients(A; dx=dx, dy=dy)
    end
    gx = fx_all[row, col]
    gy = fy_all[row, col]
    gnorm = hypot(gx, gy)
    if gnorm == 0
        return (grad=(gx,gy), steepest_unit=(0.0,0.0), flat_unit=(0.0,0.0), grad_norm=gnorm)
    end
    s = (-gx/gnorm, -gy/gnorm)   # steepest descent unit vector (physical x,y)
    f = (-gy/gnorm, gx/gnorm)    # orthogonal (flat to first order)
    return (grad=(gx,gy), steepest_unit=s, flat_unit=f, grad_norm=gnorm)
end

"""
Bilinear interpolation of matrix `A` at continuous coordinates (x_idx, y_idx), where
`x_idx` corresponds to column index (1..nx) and `y_idx` to row index (1..ny).
If (x_idx,y_idx) is outside the domain the function returns `NaN`.
"""
function bilinear_interp(A::AbstractMatrix, x_idx::Real, y_idx::Real)
    ny, nx = size(A)
    # if outside domain return NaN
    if x_idx < 1 || x_idx > nx || y_idx < 1 || y_idx > ny
        return NaN
    end
    # corners
    i0 = clamp(floor(Int, x_idx), 1, nx)
    j0 = clamp(floor(Int, y_idx), 1, ny)
    i1 = min(i0 + 1, nx)
    j1 = min(j0 + 1, ny)
    # local weights
    tx = x_idx - i0
    ty = y_idx - j0
    v00 = A[j0, i0]
    v10 = A[j0, i1]
    v01 = A[j1, i0]
    v11 = A[j1, i1]
    return (1-tx)*(1-ty)*v00 + tx*(1-ty)*v10 + (1-tx)*ty*v01 + tx*ty*v11
end

"""
Bilinear interpolation of the matrix `A` at physical coordinates (x_coord, y_coord)
given the coordinate vectors `xcoords` (length nx, for columns) and `ycoords` (length ny, for rows).
If (x_coord,y_coord) is outside the provided coordinate ranges the function returns `NaN`.
"""
function bilinear_interp_coords(A::AbstractMatrix, xcoords::AbstractVector, ycoords::AbstractVector, x_coord::Real, y_coord::Real)
    ny, nx = size(A)
    # quick out-of-range check
    if x_coord < minimum(xcoords) || x_coord > maximum(xcoords) || y_coord < minimum(ycoords) || y_coord > maximum(ycoords)
        return NaN
    end
    # find surrounding column indices (i0,i1) such that xcoords[i0] <= x_coord <= xcoords[i1]
    i1 = searchsortedfirst(xcoords, x_coord)
    i0 = max(i1-1, 1)
    if i1 > length(xcoords)
        i1 = length(xcoords)
        i0 = max(i1-1,1)
    end
    j1 = searchsortedfirst(ycoords, y_coord)
    j0 = max(j1-1, 1)
    if j1 > length(ycoords)
        j1 = length(ycoords)
        j0 = max(j1-1,1)
    end
    x0 = xcoords[i0]; x1 = xcoords[i1]
    y0 = ycoords[j0]; y1 = ycoords[j1]
    # if coordinates coincide with grid points
    if x1 == x0
        tx = 0.0
    else
        tx = (x_coord - x0) / (x1 - x0)
    end
    if y1 == y0
        ty = 0.0
    else
        ty = (y_coord - y0) / (y1 - y0)
    end
    v00 = A[j0, i0]
    v10 = A[j0, i1]
    v01 = A[j1, i0]
    v11 = A[j1, i1]
    return (1-tx)*(1-ty)*v00 + tx*(1-ty)*v10 + (1-tx)*ty*v01 + tx*ty*v11
end

"""
Sample the matrix `A` along a line going through grid index (row, col) in the given
unit direction `dir` = (ux, uy) (physical x,y). `half_length` is the distance to sample
on each side (in same physical units as dx/dy). Returns (distances, values, idx_coords)
where distances are measured from the center (0 at center), values are interpolated A values,
and idx_coords are tuples of (row_idx, col_idx) used for interpolation.
"""
function sample_line_along(A::AbstractMatrix, row::Integer, col::Integer, dir::Tuple{<:Real,<:Real};
                           half_length::Real=20.0, npoints::Int=201, dx::Real=1.0, dy::Real=1.0, clamp::Bool=false)
    # dir is unit vector in physical x (cols) and y (rows)
    ux, uy = dir
    # distances along the line (centered)
    s = range(-half_length, stop=half_length, length=npoints)
    vals = Vector{Float64}(undef, length(s))
    idx_coords = Vector{Tuple{Float64,Float64}}(undef, length(s))
    for (k, dist) in enumerate(s)
        # displacement in physical units -> convert to delta indices
        delta_cols = (dist * ux) / dx
        delta_rows = (dist * uy) / dy
        col_pos = col + delta_cols
        row_pos = row + delta_rows
        idx_coords[k] = (row_pos, col_pos)
        vals[k] = bilinear_interp(A, col_pos, row_pos)
    end
    if clamp
        # keep only in-domain samples (non-NaN values)
        inds = findall(x -> !isnan(x), vals)
        return collect(s)[inds], vals[inds], idx_coords[inds]
    else
        return collect(s), vals, idx_coords
    end
end

"""
Estimate the Hessian (2x2 matrix) of A at index (row,col) using central differences
where possible and one-sided differences on boundaries. Returns a 2x2 symmetric matrix
H = [fxx fxy; fxy fyy]. dx/dy are grid spacings in x (cols) and y (rows).
"""
function compute_hessian_at(A::AbstractMatrix, row::Integer, col::Integer; dx::Real=1.0, dy::Real=1.0)
    ny, nx = size(A)
    # fxx
    if 1 < col < nx
        fxx = (A[row, col+1] - 2*A[row, col] + A[row, col-1]) / (dx^2)
    else
        if col == 1 && nx >= 3
            fxx = (A[row, col] - 2*A[row, col+1] + A[row, col+2]) / (dx^2)
        elseif col == nx && nx >= 3
            fxx = (A[row, col-2] - 2*A[row, col-1] + A[row, col]) / (dx^2)
        else
            fxx = 0.0
        end
    end
    # fyy
    if 1 < row < ny
        fyy = (A[row+1, col] - 2*A[row, col] + A[row-1, col]) / (dy^2)
    else
        if row == 1 && ny >= 3
            fyy = (A[row, col] - 2*A[row+1, col] + A[row+2, col]) / (dy^2)
        elseif row == ny && ny >= 3
            fyy = (A[row-2, col] - 2*A[row-1, col] + A[row, col]) / (dy^2)
        else
            fyy = 0.0
        end
    end
    # fxy (mixed derivative) using central cross difference when possible
    if 1 < row < ny && 1 < col < nx
        fxy = (A[row+1, col+1] - A[row+1, col-1] - A[row-1, col+1] + A[row-1, col-1]) / (4dx*dy)
    else
        # fallback: finite difference approximations using nearest available neighbors
        fxy = 0.0
    end
    return [fxx fxy; fxy fyy]
end

"""
Return eigenvalues and eigenvectors of 2x2 symmetric Hessian H.
Eigenvalues are returned in ascending order; eigenvectors are columns.
"""
function principal_directions(H::AbstractMatrix{T}) where T
    evals, evecs = eigen(Symmetric(H))
    return evals, evecs
end

"""
Compute the principal curvature direction at (row,col): returns the eigenvector
corresponding to the smallest eigenvalue (flattest curvature). Also returns
the Hessian and eigenvalues for inspection.
"""
function curvature_direction_at(A::AbstractMatrix, row::Integer, col::Integer; dx::Real=1.0, dy::Real=1.0)
    H = compute_hessian_at(A, row, col; dx=dx, dy=dy)
    evals, evecs = principal_directions(H)
    # eigenvector for smallest eigenvalue is evecs[:,1]
    v = evecs[:,1]
    # convert Hessian eigenvector to (ux,uy) in physical x (cols), y (rows)
    # ensure unit length
    v = v / norm(v)
    # return as tuple (ux, uy)
    return (hessian=H, eigvals=evals, eigvec=v)
end


# Try to load Plots at top-level so we don't put `using` inside a function
const __PLOTS_AVAILABLE__ = (try
    @eval using Plots
    true
catch
    false
end)

"""
Plot the matrix `A` as a heatmap and overlay sampled lines.

Arguments:
- A: matrix (rows, cols)
- row, col: center grid index
- dirs: Vector of tuples (ux,uy) unit directions to overlay
- labels: optional vector of labels for the directions
- half_length, npoints, dx, dy: passed to sampling

This uses Plots.jl; if it's not available the function will warn.
"""
function overlay_sampled_lines(A::AbstractMatrix, row::Integer, col::Integer, dirs::AbstractVector; labels=nothing, half_length::Real=20.0, npoints::Int=201, dx::Real=1.0, dy::Real=1.0, clamp::Bool=true, show_profiles::Bool=false, savepath::Union{Nothing,String}=nothing, xcoords::Union{Nothing,AbstractVector}=nothing, ycoords::Union{Nothing,AbstractVector}=nothing, xlabel::AbstractString="col", ylabel::AbstractString="row")
    if !__PLOTS_AVAILABLE__
        @warn "Plots.jl not available; cannot plot overlay"
        return nothing
    end
    # sample each direction
    samples = Dict{String,Any}()
    labellist = labels === nothing ? ["dir$(i)" for i in 1:length(dirs)] : labels
    lines = Vector{Tuple{Vector{Float64},Vector{Float64}}}()
    # helper: interpolate a coordinate vector at a fractional 1-based index
    interp_from_index = function(coords::AbstractVector, idx::Real)
        n = length(coords)
        if idx <= 1
            return coords[1]
        elseif idx >= n
            return coords[end]
        else
            i0 = Base.clamp(floor(Int, idx), 1, n-1)
            i1 = i0 + 1
            t = idx - i0
            return (1-t)*coords[i0] + t*coords[i1]
        end
    end

    for (i, d) in enumerate(dirs)
        s, vals, idx_coords = sample_line_along(A, row, col, d; half_length=half_length, npoints=npoints, dx=dx, dy=dy, clamp=clamp)
        # idx_coords entries are (row_pos, col_pos) in 1-based index space
        cols_idx = [c for (r,c) in idx_coords]
        rows_idx = [r for (r,c) in idx_coords]
        # convert index positions to plotting coordinates when xcoords/ycoords available
        if xcoords !== nothing && length(xcoords) >= 2
            xs = [interp_from_index(xcoords, c) for c in cols_idx]
        else
            xs = cols_idx
        end
        if ycoords !== nothing && length(ycoords) >= 2
            ys = [interp_from_index(ycoords, r) for r in rows_idx]
        else
            ys = rows_idx
        end
        push!(lines, (xs, ys))
        samples[labellist[i]] = (s=s, vals=vals, xs=xs, ys=ys)
    end

    # base plotting: optionally show heatmap + profiles or heatmap alone
    colors = [:red, :blue, :green, :yellow, :magenta, :cyan]
    if show_profiles
        plt = Plots.plot(layout = (1,2), size=(1000,400))
        # plot with provided coordinate axes if available
        if xcoords !== nothing && ycoords !== nothing
            Plots.contourf!(plt[1], xcoords, ycoords, A; aspect_ratio=:equal, colorbar=true, xlabel=xlabel, ylabel=ylabel, title="Matrix with sampled lines")
        else
            Plots.contourf!(plt[1], A; aspect_ratio=:equal, colorbar=true, xlabel=xlabel, ylabel=ylabel, title="Matrix with sampled lines")
        end
        for (i, (xs, ys)) in enumerate(lines)
            Plots.plot!(plt[1], xs, ys, color=colors[mod1(i,length(colors))], lw=2, label=labellist[i])
        end
        # center in plotting coords
        if xcoords !== nothing && ycoords !== nothing
            center_x = interp_from_index(xcoords, col)
            center_y = interp_from_index(ycoords, row)
        else
            center_x = col
            center_y = row
        end
        Plots.scatter!(plt[1], [center_x], [center_y], color=:white, markerstrokecolor=:black, markersize=6, label="center")
        # quiver
        arrow_len = half_length*0.15
        uxv = [d[1]*arrow_len for d in dirs]
        uyv = [d[2]*arrow_len for d in dirs]
        Plots.quiver!(plt[1], [center_x for _ in dirs], [center_y for _ in dirs], quiver=(uxv, uyv), color=colors[1:length(dirs)], label="directions")

        # plot 1D profiles in second panel
        for (i, lbl) in enumerate(labellist)
            entry = samples[lbl]
            Plots.plot!(plt[2], entry.s, entry.vals, color=colors[mod1(i,length(colors))], lw=2, label=lbl)
        end
        Plots.xlabel!(plt[2], "distance from center")
        Plots.ylabel!(plt[2], "A value")
        Plots.title!(plt[2], "1D slices")
        display(plt)
        if savepath !== nothing
            try
                savefig(plt, savepath)
                @info "Saved overlay plot to $savepath"
            catch err
                @warn "Failed to save plot" error=(err, catch_backtrace())
            end
        end
        return samples, plt
    else
    # base plot: support plotting with x/y coordinates if provided
    if xcoords !== nothing && ycoords !== nothing
        p = Plots.contourf(xcoords, ycoords, A; aspect_ratio=:equal, colorbar=true, xlabel=xlabel, ylabel=ylabel, title="Matrix with sampled lines")
        center_x = interp_from_index(xcoords, col)
        center_y = interp_from_index(ycoords, row)
    else
        p = Plots.contourf(A; aspect_ratio=:equal, colorbar=true, xlabel=xlabel, ylabel=ylabel, title="Matrix with sampled lines")
        center_x = col
        center_y = row
    end
        for (i, (xs, ys)) in enumerate(lines)
            Plots.plot!(p, xs, ys, color=colors[mod1(i,length(colors))], lw=2, label=labellist[i])
        end
        Plots.scatter!(p, [center_x], [center_y], color=:white, markerstrokecolor=:black, markersize=6, label="center")
        arrow_len = half_length*0.15
        uxv = [d[1]*arrow_len for d in dirs]
        uyv = [d[2]*arrow_len for d in dirs]
        Plots.quiver!(p, [center_x for _ in dirs], [center_y for _ in dirs], quiver=(uxv, uyv), color=colors[1:length(dirs)], label="directions")
        display(p)
        if savepath !== nothing
            try
                savefig(p, savepath)
                @info "Saved overlay plot to $savepath"
            catch err
                @warn "Failed to save plot" error=(err, catch_backtrace())
            end
        end
        return samples, p
    end
end

function run_test(plotflag::Bool=false)

    row , col = 10, 10
    contour_plot_params = read_json("/home/soshala/SMEAR-PhD/SMEAR-DataFiles/Data/experiments/syn_data/optimization/Stokes/force/constant/Q2_16/1/Q2_6/simtime_10.0/Results/data/contour_plot_params.json")
    mat = contour_plot_params["cost_mat"]
    A = hcat(mat)
    A = reshape(A, 10, 10)'

    ηList = haskey(contour_plot_params, "η_list") ? contour_plot_params["η_list"] : (haskey(contour_plot_params, "eta_list") ? contour_plot_params["eta_list"] : nothing)
    βList = haskey(contour_plot_params, "β_list") ? contour_plot_params["β_list"] : (haskey(contour_plot_params, "beta_list") ? contour_plot_params["beta_list"] : nothing)

    # try to obtain x/y coordinates or ranges from the JSON file; fall back to index coords
    ny, nx = size(A)
    function _get_coords(jdict, base::String, n::Int)
        # common keys to check
        for k in (base, string(base, "_coords"), string(base, "s"), string(base, "_vals"), string(base, "_values"))
            if haskey(jdict, k)
                v = jdict[k]
                if isa(v, AbstractVector) && length(v) == n
                    return collect(v)
                end
            end
        end
        # two-element range
        rk = string(base, "_range")
        if haskey(jdict, rk)
            r = jdict[rk]
            if isa(r, AbstractVector) && length(r) >= 2
                return collect(range(r[1], stop=r[2], length=n))
            end
        end
        # explicit min/max
        kmin = string(base, "_min")
        kmax = string(base, "_max")
        if haskey(jdict, kmin) && haskey(jdict, kmax)
            return collect(range(jdict[kmin], stop=jdict[kmax], length=n))
        end
        return collect(1:n)
    end

    # prefer ηList/βList if present (these correspond to x and y axes in the JSON)
    if ηList !== nothing && isa(ηList, AbstractVector) && length(ηList) == nx
        xcoords = collect(ηList)
    else
        xcoords = _get_coords(contour_plot_params, "x", nx)
    end
    if βList !== nothing && isa(βList, AbstractVector) && length(βList) == ny
        ycoords = collect(βList)
    else
        ycoords = _get_coords(contour_plot_params, "y", ny)
    end
    
    fx, fy = compute_gradients(A; dx=1.0, dy=1.0)

    info = directions_at(A, row, col; dx=1.0, dy=1.0, fx_all=fx, fy_all=fy)
    println("Gradient at ($row,$col): ", info.grad, " norm=", info.grad_norm)
    println("Steepest unit: ", info.steepest_unit, "  flat unit: ", info.flat_unit)

    # clamp samples (drop out-of-domain points) and use smaller half_length for a 10x10 grid
    s1, vals1, _ = sample_line_along(A, row, col, info.steepest_unit; half_length=6.0, npoints=201, clamp=true)
    s2, vals2, _ = sample_line_along(A, row, col, info.flat_unit; half_length=6.0, npoints=201, clamp=true)

    if isempty(vals1)
        println("Steepest slice: no in-domain samples (empty)")
            fx, fy = compute_gradients(A; dx=1.0, dy=1.0)
        println("Steepest slice: min=$(minimum(vals1)), max=$(maximum(vals1))")
            # Automatic orientation check: if ηList/βList appear swapped relative to A dims,
            # transpose A and swap x/y coords. This handles differences in flattening order.
            # If both ηList and βList are present and their lengths match swapped dims, we transpose.
            try
                need_transpose = false
                if ηList !== nothing && βList !== nothing
                    if isa(ηList, AbstractVector) && isa(βList, AbstractVector)
                        if length(ηList) == ny && length(βList) == nx
                            need_transpose = true
                        end
                    end
                else
                    # if only one list present and its length matches the other dimension, transpose
                    if ηList !== nothing && isa(ηList, AbstractVector) && length(ηList) == ny && length(ηList) != nx
                        need_transpose = true
                    elseif βList !== nothing && isa(βList, AbstractVector) && length(βList) == nx && length(βList) != ny
                        need_transpose = true
                    end
                end
                if need_transpose
                    @info "Transpose detected: swapping matrix orientation to match η/β ordering"
                    A = A'
                    # swap sizes
                    ny, nx = size(A)
                    # swap coordinate arrays
                    xcoords, ycoords = ycoords, xcoords
                    # swap η/β lists too so downstream logic still refers to correct vectors
                    ηList, βList = βList, ηList
                    # swap center indices
                    row, col = col, row
                end
            catch err
                @warn "Orientation check failed" error=(err, catch_backtrace())
            end
    end
    if isempty(vals2)
        println("Flat slice: no in-domain samples (empty)")
    else
        println("Flat slice:     min=$(minimum(vals2)), max=$(maximum(vals2))")
    end

    # --- Hessian / principal curvature (flattest curvature) ---
    cur = curvature_direction_at(A, row, col; dx=1.0, dy=1.0)
    println("Hessian at point:\n", cur.hessian)
    println("Hessian eigenvalues: ", cur.eigvals)
    println("Curvature flattest eigenvector (col, row) = ", cur.eigvec)
    s3, vals3, _ = sample_line_along(A, row, col, (cur.eigvec[1], cur.eigvec[2]); half_length=6.0, npoints=201, clamp=true)
    if isempty(vals3)
        println("Curvature-flat slice: no in-domain samples (empty)")
    else
        println("Curvature-flat slice: min=$(minimum(vals3)), max=$(maximum(vals3))")
    end

    # plotting of slices handled by overlay_sampled_lines

    if plotflag
        if __PLOTS_AVAILABLE__
            # Use overlay helper that plots the matrix and overlays sampled lines + quivers
            dirs = [info.steepest_unit, info.flat_unit, (cur.eigvec[1], cur.eigvec[2])]
            labels = ["steepest", "flat", "curvature-flat"]
            # clamp out-of-domain samples and show 1D profiles to avoid NaNs in small grids
            overlay_sampled_lines(A, row, col, dirs; labels=labels, half_length=10.0, npoints=201, dx=1.0, dy=1.0, clamp=true, show_profiles=true, savepath="slice_plot.png", xcoords=xcoords, ycoords=ycoords, xlabel="η", ylabel="β")
        else
            @warn "Plots.jl not available; run `import Pkg; Pkg.add(\"Plots\")` to enable plotting"
        end
    end
    println("slice_and_plot.jl done")
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    plotflag = any(x->x=="--plot", ARGS)
    run_test(plotflag)
end
