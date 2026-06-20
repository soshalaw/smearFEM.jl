using Plots, Plots.PlotMeasures
using LaTeXStrings
using StatsPlots
using ProgressMeter
using ArgCheck
using Base.Filesystem


"""
    plot_meshgrid(NodeList, IEN)

Function to plot the mesh

# Arguments:
- `NodeList::Matrix{Float64}{nNodes,ndim}`: array of nodes.
- `IEN::Matrix{nElem,nNodes}`: IEN array.
"""
function plot_mesh(NodeList, IEN)

    sz = size(NodeList,1)
    if sz == 2
        Plots.scatter(NodeList[1,:], NodeList[2,:], markersize=2, label="", dpi=:400)
        iter = 1:size(IEN,1)
        for i in iter
            x = NodeList[1,IEN[i,:]]
            y = NodeList[2,IEN[i,:]]
            Plots.plot!(x, y, marker=1.5, lw=0.5, label="")
        end

    elseif sz == 3
        Plots.scatter3d(NodeList[1,:], NodeList[2,:], NodeList[3,:], markersize=2, label="", dpi=:400, lw=1)
        iter = 1:size(IEN,1)
        for i in iter
            x = NodeList[1,IEN[i,:]]
            y = NodeList[2,IEN[i,:]]
            z = NodeList[3,IEN[i,:]]
            Plots.plot3d!(x, y, z,marker=1.5, lw=0.5, label="", dpi=:400)
        end

        Plots.xlabel!(L"x")
        Plots.ylabel!(L"y")
        Plots.zlabel!(L"z")
        Plots.title!("3D Grid")
    end
end

"""
    animate_fields(;filepath=nothing, fields=nothing , IEN=nothing, border_nodes_2d=nothing, sim_pts_2d=nothing, p=nothing, q=nothing)

Function to animate the fields as a gif

# Arguments:
- `fields::Vector{Vector{Float64}}`: solution vector
- `sim_pts_2d::Vector{Vector{Float64}}`: 2D projection of the solution vector
- `border_nodes_2d::Vector{Vector{Float64}}`: 2D coordinates of the border nodes of the mesh
- `IEN::Matrix{Float64}{nElem, nNodes}`: IEN array
- `p::Vector{Float64}`: x coordinates of the extracted convex hull
- `q::Vector{Float64}`: y coordinates of the extracted convex hull
"""
function animate_fields(; filepath::String="None", Nodes=nothing , SurfaceNodes3D = nothing, IEN=nothing, border_nodes_2d=nothing, sim_pts_2d=nothing, p=nothing, q=nothing, pObs=nothing, qObs=nothing, pgt=nothing, qgt=nothing, cam_pose::Union{Nothing,Vector{Float64}}=nothing, h::Union{Nothing,Float64}=nothing)
    set_file(filepath) # create the directory to store the VTK files
    if isnothing(Nodes) && isnothing(sim_pts_2d) && isnothing(border_nodes_2d) && isnothing(IEN) && isnothing(p) && isnothing(q) && isnothing(pObs) && isnothing(qObs)
        throw(AssertionError("No fields provided"))
    elseif isnothing(Nodes) && isnothing(SurfaceNodes3D)
        animate2D(border_nodes_2d=border_nodes_2d, sim_pts_2d=sim_pts_2d, p=p, q=q, pObs=pObs, qObs=qObs, pgt=pgt, qgt=qgt, filepath=filepath)
        return
    elseif isnothing(sim_pts_2d) && isnothing(border_nodes_2d) && isnothing(p) && isnothing(q)
        animate3D(fields=Nodes, IEN=IEN, filepath=filepath, cam_pose=cam_pose, height=h)
        return
    else 
        if !isnothing(SurfaceNodes3D)
            animate3D(surface_pts=SurfaceNodes3D, filepath=filepath, cam_pose=cam_pose, height=h)
        else
            animate3D(fields=Nodes, IEN=IEN, filepath=filepath, cam_pose=cam_pose, height=h)
        end
        animate2D(border_nodes_2d=border_nodes_2d, sim_pts_2d=sim_pts_2d, p=p, q=q, pObs=pObs, qObs=qObs, filepath=filepath)
    end
end 

"""
    animate2D(;border_nodes_2d=nothing, sim_pts_2d=nothing, p=nothing, q=nothing)

Function to animate the 2D fields as a gif

# Arguments:
- `border_nodes_2d::Vector{Vector{Float64}}`: 2D coordinates of the border nodes of the mesh
- `sim_pts_2d::Vector{Vector{Float64}}`: 2D projection of the solution vector
- `p::Vector{Float64}`: x coordinates of the extracted convex hull
- `q::Vector{Float64}`: y coordinates of the extracted convex hull
"""
function animate2D(;border_nodes_2d=nothing, sim_pts_2d=nothing, p=nothing, q=nothing, pObs=nothing, qObs=nothing, pgt=nothing, qgt=nothing, filepath="images/2D_grid.gif")

    if isnothing(border_nodes_2d) && isnothing(sim_pts_2d) && isnothing(p) && isnothing(pObs)
        throw(AssertionError("No fields provided"))
    else
        if !isnothing(border_nodes_2d)
            sz = length(border_nodes_2d)
        elseif !isnothing(p)
            sz = length(p)
        elseif !isnothing(sim_pts_2d)
            sz = length(sim_pts_2d)
        elseif !isnothing(pObs)
            sz = length(pObs)
        elseif !isnothing(pgt)
            sz = length(pgt)
        end

        # Preprocess [t][a] → [a][t] for both border_nodes_2d and sim_pts_2d
        _multi_angle(x) = !isnothing(x) && !isempty(x) && x[1] isa AbstractVector
        n_angles = _multi_angle(border_nodes_2d) ? length(border_nodes_2d[1]) :
                   _multi_angle(sim_pts_2d)      ? length(sim_pts_2d[1])      : 1

                   
        borders_θ = _multi_angle(border_nodes_2d) ?
            [[border_nodes_2d[t][a] for t in 1:sz] for a in 1:n_angles] : nothing
        sim_pts_θ = _multi_angle(sim_pts_2d) ?
            [[sim_pts_2d[t][a] for t in 1:sz] for a in 1:n_angles] :
            [sim_pts_2d for _ in 1:n_angles]

        for a in 1:n_angles
            pr = Progress(sz; desc="Animating 2D fields (angle $a of $n_angles)...", showspeed=true)
            borders_a  = !isnothing(borders_θ) ? borders_θ[a] : nothing
            sim_pts_2d_a = sim_pts_θ[a]
            animation = @animate for i in 1:sz
                plt = set_plot(12, frame=true, bottom_margin=0mm)
                if !isnothing(pObs) && !isnothing(p) && !isnothing(pgt)
                    Plots.plot!(plt, p[i],q[i], labels=L"\mathrm{Lagrange basis}", aspect_ratio=:equal, dpi=:400, lw=2, color=:darkorange)
                    Plots.plot!(plt, pObs[i],qObs[i], labels=L"\mathrm{NURBS basis}", aspect_ratio=:equal, dpi=:400, lw=2, color=:cyan4)
                    Plots.plot!(plt, pgt[i],qgt[i], labels=L"\mathrm{Observation}", aspect_ratio=:equal, dpi=:400, lw=2, color=:green4)
                else
                    if !isnothing(p)
                        Plots.plot!(plt, p[i],q[i], legend=true, labels=L"\mathrm{Simulation}", aspect_ratio=:equal, dpi=:400, lw=1, color=:darkorange)
                    end
                    if !isnothing(pObs)
                        Plots.plot!(plt, pObs[i],qObs[i], labels=L"\mathcal{B}_{k,j}", aspect_ratio=:equal, dpi=:400, lw=2, color=:cyan4)
                    end
                    if !isnothing(sim_pts_2d_a)
                        Plots.scatter!(plt, sim_pts_2d_a[i][1,:], sim_pts_2d_a[i][2,:], ms=:2, mc=:royalblue, ma=:0.7, legend=true, labels=L"\mathcal{U}_{k,j}", aspect_ratio=:equal, dpi=:400)
                    end
                    if !isnothing(borders_a)
                        Plots.scatter!(plt, borders_a[i][1,:], borders_a[i][2,:], ms=:3, mc=:indianred2, legend=true, labels=L"\widehat{\mathcal{B}}_{k,j}", aspect_ratio=:equal, dpi=:400)
                    end
                end

                Plots.xlabel!(plt, L"x\;[\mathrm{px}]")
                Plots.ylabel!(plt, L"y\;[\mathrm{px}]")
                Plots.xlims!(plt, 0,2048)
                Plots.ylims!(plt, 0,1536)
                Plots.yflip!(plt, true)
                next!(pr)
                try
                    yield()
                catch e
                    @debug "Animation yield interrupted (safe to continue): $(e.msg)"
                end
            end
            gif(animation, string(filepath, "/2D_grid_angle_$(a).gif"), fps=5)
        end
    end
end


"""
    animate3D(fields)

Function to animate the 3D fields as a gif

# Arguments:
- `fields::Vector{Vector{Float64}}`: solution vector
"""
function animate3D(;fields=nothing, surface_pts=nothing, IEN=nothing, filepath="images/3D_grid.gif", cam_pose::Union{Nothing,Vector{Float64}}=nothing, height::Union{Nothing,Float64}=nothing)
    
    if !isnothing(fields)
        sz = length(fields)
        xmax = maximum(fields[1][1,:])
        xmin = minimum(fields[1][1,:])
        ymax = maximum(fields[1][2,:])
        ymin = minimum(fields[1][2,:])
        zmax = maximum(fields[1][3,:])
        zmin = minimum(fields[1][3,:])
    elseif !isnothing(surface_pts)
        sz = length(surface_pts)            
        xmax = maximum(surface_pts[1][1,:])
        xmin = minimum(surface_pts[1][1,:])
        ymax = maximum(surface_pts[1][2,:])
        ymin = minimum(surface_pts[1][2,:])
        zmax = maximum(surface_pts[1][3,:])
        zmin = minimum(surface_pts[1][3,:])
    end

    if cam_pose !== nothing
        # Adjust limits to include camera position
        cam_x, cam_y, cam_z = cam_pose
        xmax = max(xmax, cam_x*1.1) # extend limits by 50% to ensure camera is visible
        xmin = min(xmin, cam_x*1.1)
        ymax = max(ymax, cam_y*1.1)
        ymin = min(ymin, cam_y*1.1)
        zmax = max(zmax, cam_z*1.1)
        zmin = min(zmin, cam_z*1.1)
    end

    fz = 12
    iter = 1:sz
    fac = 0.2 # factor to extend the limits of the plot
    pr = Progress(sz; desc="Animating 3D fields...",showspeed=true)
    animation = @animate for i in iter

        plt = Plots.plot(1, 
                            xlims=(xmin-fac*(xmax-xmin), xmax+fac*(xmax-xmin)),
                            ylims=(ymin-fac*(ymax-ymin), ymax+fac*(ymax-ymin)),
                            zlims=(zmin, zmax+fac*(zmax-zmin)*0.1),

                            xlabel=L"x\;[\mathrm{mm}]",ylabel=L"y\;[\mathrm{mm}]",zlabel=L"z\;[\mathrm{mm}]",

                            tickfontsize = fz-2,
                            guidefontsize = fz,
                            legendfontsize = fz-2, 

                            aspect_ratio = :equal, 

                            size = (500,500), 

                            fontfamily = "computer modern",
                            framestyle = :semi, 
                            margin = 0mm, 

                            grid = :true, 
                            minorgrid = :false, 
                            camera = (45, 20),
                            label="")
        if !isnothing(fields)
            if isnothing(IEN)
                Plots.scatter3d!(plt,fields[i][1,:], fields[i][2,:], fields[i][3,:], mc=:indianred2, markersize=3, label=:"", dpi=:400)
            else
                ien_iter = 1:size(IEN,2)
                Plots.scatter3d!(plt,fields[i][1,:], fields[i][2,:], fields[i][3,:], mc=:royalblue, markersize=3, label=:"", dpi=:400)
                seq = []
                if size(IEN,1) == 27
                    seq = [1,9,2,10,3,11,4,12,1,17,5,13,6,14,7,15,8,16,5,13,6,18,2,10,3,19,7,15,8,20,4]
                elseif size(IEN,1) == 8
                    seq = [1,2,3,4,1,5,6,8,5,6,2,3,7,8,4]
                end
                for e in ien_iter
                    IEN_e = IEN[:,e]
                    x = fields[i][1,IEN_e[seq]]
                    y = fields[i][2,IEN_e[seq]]
                    z = fields[i][3,IEN_e[seq]]
                    Plots.plot!(x,y,z,color=:black,label="",lw=1.5)
                end
            end
        end

        if !isnothing(cam_pose)
            cam_frame = _set_cam_frame(cam_pose, height/2)

            frame_x = cam_frame[1:3,1] + cam_pose
            frame_y = cam_frame[1:3,2] + cam_pose
            frame_z = cam_frame[1:3,3] + cam_pose

            cam_x, cam_y, cam_z = cam_pose
            len = 1.0 # length of the camera frame axes
            Plots.plot!(plt, [cam_pose[1], frame_x[1]*len], [cam_pose[2], frame_x[2]*len], [cam_pose[3], frame_x[3]*len], color=:red, label="", lw=2)
            Plots.plot!(plt, [cam_pose[1], frame_y[1]*len], [cam_pose[2], frame_y[2]*len], [cam_pose[3], frame_y[3]*len], color=:green, label="", lw=2)
            Plots.plot!(plt, [cam_pose[1], frame_z[1]*len], [cam_pose[2], frame_z[2]*len], [cam_pose[3], frame_z[3]*len], color=:blue, label="", lw=2)

        end
        
        if !isnothing(surface_pts)
            Plots.scatter3d!(plt, surface_pts[i][1,:], surface_pts[i][2,:], surface_pts[i][3,:], mc=:indianred2, markersize=3, label=:"", dpi=:400)
        end
        next!(pr)
        try
            yield()
        catch e
            @debug "Animation yield interrupted (safe to continue): $(e.msg)"
        end
    end

    gif(animation, string(filepath,"3D_grid.gif"), fps=10)
end

function plot_matches(simborderfields, p, q, pObs, qObs, pairsList, filepath::String="None")
    set_file(filepath)
    sz = length(simborderfields)
    @argcheck length(simborderfields) == length(pairsList)
    pr = Progress(sz; desc="Plotting matches ...",showspeed=true)
    iter = 1:sz
    animation = @animate for i in iter
        borderp = simborderfields[i][1,:]
        borderq = simborderfields[i][2,:]
        plt = Plots.plot(1,xlims=(0,2048), ylims=(0,1536), xlabel="x",ylabel="y",title="Prospective Projection of the 3D Grid", label="", dpi=400)
        Plots.plot!(p[i],q[i], legend=true, labels="Simulation",  dpi=:400, lw=2)
        Plots.plot!(pObs[i],qObs[i], labels="Observation",  dpi=:400, lw=2)
        pairs = pairsList[i]
        pairtIter = 1:size(pairs,1)
        for j in pairtIter
            Plots.plot!([borderp[pairs[j,1]], pObs[i][pairs[j,2]]], [borderq[pairs[j,1]], qObs[i][pairs[j,2]]], marker=1.5, lw=0.5, label="", dpi=:400)
        end
        Plots.xlims!(0,2048)
        Plots.ylims!(0,1536) 
        Plots.xlabel!(L"x")
        Plots.ylabel!(L"y")
        Plots.title!("Point Correspondence")
        next!(pr)
        try
            yield()
        catch e
            @debug "Animation yield interrupted (safe to continue): $(e.msg)"
        end
    end
    gif(animation, string(filepath,"/matches.gif"), fps=10)
end

function set_plot(fs::Int; sz::Tuple{Int,Int}=(477,350), legend_column::Int=1, right_margin=0pt, left_margin=0pt, top_margin=0pt, bottom_margin=-20mm, frame=false, legend=:outerbottom)

    frame_brdr = nothing
    if frame
        frame_brdr = :match
    end
    plt = Plots.plot(1, 
                            xtickfontsize = round(Int,fs*0.85), ytickfontsize = round(Int,fs*0.85),
                            titlefontsize = fs,
                            xguidefontsize = fs,
                            yguidefontsize = fs,
                            legendfontsize = round(Int,fs*0.8), 

                            size = sz, 
                            fontfamily = "computer modern",
                            framestyle = :box, 
                            top_margin=top_margin,
                            left_margin=left_margin,
                            bottom_margin=bottom_margin, # negative margin to move legend down
                            right_margin=right_margin,

                            grid = :true, 
                            minorgrid = :true, 
                            lw = 1,
                            legend=legend, 
                            legend_column=legend_column, 
                            foreground_color_legend = frame_brdr,

                            label=false)
    return plt
end


function set_subplot(fs::Int; sz::Tuple{Int,Int}=(1000,750), layout=(1,1), legend=nothing, legend_column=nothing, bottom_margin=nothing)

    # Build optional kwargs for legend placement if provided
    kw = ()
    if legend !== nothing
        kw = (kw..., legend=legend)
    end
    if legend_column !== nothing
        kw = (kw..., legend_column=legend_column)
    end
    if bottom_margin !== nothing
        kw = (kw..., bottom_margin=bottom_margin)
    end

    plt = Plots.plot(layout=layout,
                            xtickfontsize = fs, ytickfontsize = fs,
                            titlefontsize = fs,
                            xguidefontsize = fs,
                            yguidefontsize = fs,
                            legendfontsize = fs, 

                            size = sz, 
                            fontfamily = "computer modern",
                            framestyle = :box, 
                            margin = 6mm, 

                            grid = :true, 
                            minorgrid = :true, 
                            lw = 4,
                            
                            label=""; kw...)
    return plt
end
"""
    normalize(q, IEN)

Function normalize the solution vector for plotting
    
# Arguments:
- `q`: solution vector
- `IEN::Matrix{Float64}{nElem, nNodes}`: IEN array

# Returns:
- `qList`: normalized list of solutions 
"""
function normalize(q, IEN)

    qList = zeros(size(IEN))
    max_val = maximum(q)
    min_val = minimum(q)
    iter = 1:size(IEN,1)
    for e in iter
        for n in 1:4
            qList[e,n] = (q[IEN[e,n]] - min_val) / (max_val - min_val)
        end
    end
    return qList
end

"""
    truncate_colormap(minval=0.0, maxval=1.0, n=100)

Function to truncate a colormap
    
# Arguments:
- `minval::Integer`: minimum value of the colormap
- `maxval::Integer`: maximum value of the colormap
- `n::Integer`: number of colors

# Returns:
- `new_cmap`: truncated colormap
"""
function truncate_colormap(minval=0.0, maxval=1.0, n=100)
    new_cmap = matplotlib.colors.LinearSegmentedColormap.from_list("mycmap", get_cmap("jet")(collect(range(maxval, minval, n))))
    return new_cmap
end

function plot_covariance!(plt::Plots.Plot, η_list::Vector{Float64}, β_list::Vector{Float64}, label::String)

    mean_η = mean(η_list)
    mean_β = mean(β_list)

    cov_η = cov(η_list)
    cov_β = cov(β_list)
    cov_ηβ = cov(η_list, β_list)

    cov_mat = [cov_η cov_ηβ; cov_ηβ cov_β]
    mean_vec = [mean_η; mean_β]

    # Plot the covariance matrix
    StatsPlots.covellipse!(plt, mean_vec, cov_mat, label=label, color=:red, alpha=0.5, bottom_margin = -15mm, legend=:outerbottom, legend_column=legend_column)
    StatsPlots.scatter!(plt, η_list, β_list, label="Data points", dpi=:400, ms=:2, markerstrokewidth=0.1)
    xlabel!(plt, L"\eta")
    ylabel!(plt, L"\beta")    

end

function plot_covariance!(plt::Plots.Plot, η_list::Vector{Float64}, β_list::Vector{Float64}; label::String="Covariance", legend_column::Int=1, color_ellipse=nothing, color_scatter=nothing)

    mean_η = mean(η_list)
    mean_β = mean(β_list)

    cov_η = cov(η_list)
    cov_β = cov(β_list)
    cov_ηβ = cov(η_list, β_list)

    cov_mat = [cov_η cov_ηβ; cov_ηβ cov_β]
    mean_vec = [mean_η; mean_β]

    # Build kwargs for covellipse
    covellipse_kwargs = (label=label, alpha=0.5, legend_column=legend_column, bottom_margin=-30mm)
    if !isnothing(color_ellipse)
        covellipse_kwargs = (covellipse_kwargs..., color=color_ellipse, lw=1)
    end

    # Build kwargs for scatter
    scatter_kwargs = (label="Data points", ms=:10, markerstrokewidth=0.1)
    if !isnothing(color_scatter)
        scatter_kwargs = (scatter_kwargs..., color=color_scatter)
    end

    # Plot the covariance matrix
    StatsPlots.covellipse!(plt, mean_vec, cov_mat; covellipse_kwargs...)

    # Plots.scatter!(plt, η_list, β_list; scatter_kwargs...)
end
