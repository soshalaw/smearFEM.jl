using Plots, Plots.PlotMeasures
using LaTeXStrings
using StatsPlots
using ProgressMeter
using ArgCheck

"""
    PlotGrid(IEN, NodeList)

Function to plot the grid

# Arguments:
- `IEN::Matrix{Float64}{nElem, nNodes}`: IEN array.
- `NodeList::Matrix{Float64}{nNodes, ndim}`: array of nodes.
"""
function PlotGrid(IEN, NodeList)

    fig1 = plt.figure()
    if size(IEN,2) == 4 # 2D element with 4 nodes 
        ax = fig1.add_subplot(111)
        iter = 1:size(IEN,1)
        for i in iter
            x = NodeList[1,IEN[i,:]]
            y = NodeList[2,IEN[i,:]]
            ax.plot(x, y, "-k", linewidth=0.5)
        end
    
        ax.scatter(NodeList[1,:],NodeList[2,:],s=10,c=red)
        ax.axis("equal")
        ax.grid("on")
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_title("2D Grid")
        
    elseif size(IEN,2) == 8 # 3D element with 8 nodes
        ax = fig1.add_subplot(111, projection="3d")
        iter = 1:size(IEN,1)
        for i in iter
            x = NodeList[1,IEN[i,:]]
            y = NodeList[2,IEN[i,:]]
            z = NodeList[3,IEN[i,:]]
            ax.plot(x, y, z,"-k", linewidth=0.5)
        end

        ax.scatter(NodeList[1,:],NodeList[2,:],NodeList[3,:],s=10,c=red)
        ax.axis("equal")
        ax.grid("on")
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_zlabel("z")
        ax.set_title("3D Grid")
    end
    gcf()
end

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
        Plots.scatter3d(NodeList[1,:], NodeList[2,:], NodeList[3,:], markersize=2, label="", dpi=:400, lw=3)
        iter = 1:size(IEN,1)
        for i in iter
            x = NodeList[1,IEN[i,:]]
            y = NodeList[2,IEN[i,:]]
            z = NodeList[3,IEN[i,:]]
            Plots.plot3d!(x, y, z,marker=1.5, lw=0.5, label="", dpi=:400, lw=3)
        end

        Plots.xlabel!(L"x")
        Plots.ylabel!(L"y")
        Plots.zlabel!(L"z")
        Plots.title!("3D Grid")
    end
end

"""
    animate_fields(;filepath=nothing, fields=nothing , IEN=nothing, BorderNodes2D=nothing, fields2D=nothing, p=nothing, q=nothing)

Function to animate the fields as a gif

# Arguments:
- `fields::Vector{Vector{Float64}}`: solution vector
- `fields2D::Vector{Vector{Float64}}`: 2D projection of the solution vector
- `BorderNodes2D::Vector{Vector{Float64}}`: 2D coordinates of the border nodes of the mesh
- `IEN::Matrix{Float64}{nElem, nNodes}`: IEN array
- `p::Vector{Float64}`: x coordinates of the extracted convex hull
- `q::Vector{Float64}`: y coordinates of the extracted convex hull
"""
function animate_fields(; filepath::String="None", Nodes=nothing , SurfaceNodes3D = nothing, IEN=nothing, BorderNodes2D=nothing, fields2D=nothing, p=nothing, q=nothing, pObs=nothing, qObs=nothing)
    set_file(filepath) # create the directory to store the VTK files
    if isnothing(Nodes) && isnothing(fields2D) && isnothing(BorderNodes2D) && isnothing(IEN) && isnothing(p) && isnothing(q)
        AssertionError("No fields provided")
        return
    elseif isnothing(Nodes) && isnothing(SurfaceNodes3D)
        animate2D(BorderNodes2D=BorderNodes2D, fields2D=fields2D, p=p, q=q, pObs=pObs, qObs=qObs, filepath=filepath)
        return
    elseif isnothing(fields2D) && isnothing(BorderNodes2D) && isnothing(p) && isnothing(q)
        animate3D(fields=Nodes, IEN=IEN, filepath=filepath)
        return
    else 
        if !isnothing(SurfaceNodes3D)
            animate3D(surface_pts=SurfaceNodes3D, filepath=filepath)
        else
            animate3D(fields=Nodes, IEN=IEN, filepath=filepath)
        end
        animate2D(BorderNodes2D=BorderNodes2D, fields2D=fields2D, p=p, q=q, pObs=pObs, qObs=qObs, filepath=filepath)
    end
end 

"""
    animate2D(;BorderNodes2D=nothing, fields2D=nothing, p=nothing, q=nothing)

Function to animate the 2D fields as a gif

# Arguments:
- `BorderNodes2D::Vector{Vector{Float64}}`: 2D coordinates of the border nodes of the mesh
- `fields2D::Vector{Vector{Float64}}`: 2D projection of the solution vector
- `p::Vector{Float64}`: x coordinates of the extracted convex hull
- `q::Vector{Float64}`: y coordinates of the extracted convex hull
"""
function animate2D(;BorderNodes2D=nothing, fields2D=nothing, p=nothing, q=nothing, pObs=nothing, qObs=nothing, filepath="images/2D_grid.gif")
    
    if isnothing(BorderNodes2D) && isnothing(fields2D) && isnothing(p) && isnothing(pObs)
        AssertionError("No fields provided")
        return
    else
        if isnothing(fields2D) && isnothing(p)
            sz = length(BorderNodes2D)
        elseif isnothing(BorderNodes2D) && isnothing(fields2D)
            sz = length(p)
        elseif isnothing(BorderNodes2D) && isnothing(p)
            sz = length(fields2D)
        elseif isnothing(fields2D) && isnothing(p)    
            sz = length(BorderNodes2D)
        elseif isnothing(BorderNodes2D) && isnothing(fields2D) && isnothing(p)
            sz = length(pObs)
        elseif !isnothing(pObs)
            sz = length(pObs)
        else
            sz = length(fields2D)
        end
        pr = Progress(sz; desc="Animating 2D fields...",showspeed=true)
        iter = 1:sz
        animation2 = @animate for i in iter
            
            plt = set_plot(22)

            if !isnothing(p)   
                Plots.plot!(p[i],q[i], legend=true, labels="Simulation", aspect_ratio = :equal, dpi=:400, lw=3) #, marker = :circle, markersize = 2)
            end
            if !isnothing(pObs)
                Plots.plot!(pObs[i],qObs[i], labels="Observation", aspect_ratio = :equal, dpi=:400, lw=2)
            end
            if !isnothing(fields2D)
                Plots.scatter!(fields2D[i][1,:], fields2D[i][2,:], ms=:4, mc=:royalblue, ma=:0.7, legend=true, labels="Surface Nodes", aspect_ratio = :equal, 
                                dpi=:400)
            end
            if !isnothing(BorderNodes2D)
                Plots.scatter!(BorderNodes2D[i][1,:], BorderNodes2D[i][2,:], ms=:6, mc=:indianred2, legend=true, labels="Border Nodes", aspect_ratio = :equal, 
                                dpi=:400)
            end

            Plots.xlabel!(L"x")
            Plots.ylabel!(L"y")
            Plots.xlims!(0,2048)
            Plots.ylims!(0,1536)
            next!(pr)
        end
        gif(animation2, string(filepath,"/2D_grid.gif"), fps=10)
    end
end

function animate_2D_comp(;borders=nothing, filepath="images/2D_grid.gif")
    
    if isnothing(BorderNodes2D) && isnothing(fields2D) && isnothing(p) && isnothing(pObs)
        AssertionError("No fields provided")
        return
    else
        if isnothing(fields2D) && isnothing(p)
            sz = length(BorderNodes2D)
        elseif isnothing(BorderNodes2D) && isnothing(fields2D)
            sz = length(p)
        elseif isnothing(BorderNodes2D) && isnothing(p)
            sz = length(fields2D)
        elseif isnothing(fields2D) && isnothing(p)    
            sz = length(BorderNodes2D)
        elseif isnothing(BorderNodes2D) && isnothing(fields2D) && isnothing(p)
            sz = length(pObs)
        elseif !isnothing(pObs)
            sz = length(pObs)
        else
            sz = length(fields2D)
        end
        pr = Progress(sz; desc="Animating 2D fields...",showspeed=true)
        iter = 1:sz
        animation2 = @animate for i in iter
            
            plt = set_plot(22)

            if !isnothing(p)   
                Plots.plot!(p[i],q[i], legend=true, labels="Simulation", aspect_ratio = :equal, dpi=:400, lw=2)
            end
            if !isnothing(pObs)
                Plots.plot!(pObs[i],qObs[i], labels="Observation", aspect_ratio = :equal, dpi=:400, lw=2)
            end
            if !isnothing(fields2D)
                Plots.scatter!(fields2D[i][1,:], fields2D[i][2,:], ms=:4, mc=:royalblue, ma=:0.7, legend=true, labels="Surface Nodes", aspect_ratio = :equal, 
                                dpi=:400)
            end
            if !isnothing(BorderNodes2D)
                Plots.scatter!(BorderNodes2D[i][1,:], BorderNodes2D[i][2,:], ms=:6, mc=:indianred2, legend=true, labels="Border Nodes", aspect_ratio = :equal, 
                                dpi=:400)
            end

            Plots.xlabel!(L"x")
            Plots.ylabel!(L"y")
            Plots.xlims!(0,2048)
            Plots.ylims!(0,1536)
            next!(pr)
        end
        gif(animation2, string(filepath,"2D_grid.gif"), fps=10)
    end
end

"""
    animate3D(fields)

Function to animate the 3D fields as a gif

# Arguments:
- `fields::Vector{Vector{Float64}}`: solution vector
"""
function animate3D(;fields=nothing, surface_pts=nothing, IEN=nothing, filepath="images/3D_grid.gif")
    
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
    fz = 22
    iter = 1:sz
    fac = 0.2 # factor to extend the limits of the plot
    pr = Progress(sz; desc="Animating 3D fields...",showspeed=true)
    animation = @animate for i in iter
        if !isnothing(fields)
            zmin = minimum(fields[i][3,:])
        elseif !isnothing(surface_pts)
            zmin = minimum(surface_pts[i][3,:])
        end
        plt = Plots.plot(1, 
                            xlims=(xmin-fac*(xmax-xmin), xmax+fac*(xmax-xmin)),
                            ylims=(ymin-fac*(ymax-ymin), ymax+fac*(ymax-ymin)),
                            zlims=(zmin, zmax+fac*(zmax-zmin)),

                            xlabel=L"x",ylabel=L"y",zlabel=L"z",

                            tickfontsize = fz-4,
                            guidefontsize = fz,
                            legendfontsize = fz-2, 

                            aspect_ratio = :equal, 

                            size = (1000,750), 

                            fontfamily = "computer modern",
                            framestyle = :semi, 
                            margin = 0mm, 

                            grid = :true, 
                            minorgrid = :false, 
                            camera = (45, 20),
                            label="")
        if !isnothing(fields)
            if isnothing(IEN)
                Plots.scatter3d!(fields[i][1,:], fields[i][2,:], fields[i][3,:], mc=:indianred2, markersize=4, label=:"", dpi=:400)
            else
                ien_iter = 1:size(IEN,2)
                Plots.scatter3d!(fields[i][1,:], fields[i][2,:], fields[i][3,:], mc=:royalblue, markersize=4, label=:"", dpi=:400)
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
        if !isnothing(surface_pts)
            Plots.scatter3d!(surface_pts[i][1,:], surface_pts[i][2,:], surface_pts[i][3,:], mc=:indianred2, markersize=4, label=:"", dpi=:400)
        end
        next!(pr)
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
    end
    gif(animation, string(filepath,"/matches.gif"), fps=10)
end

function plot_matches_h(Exptx, Expty, Obsptx, p, q, pObs, qObs, filepath::String="None")
    set_file(filepath)
    @argcheck length(Exptx) == length(Obsptx)
    sz = length(p)
    sz == length(Obsptx) || AssertionError("The number of simulated and observed data should be the same")
    pr = Progress(sz; desc="Plotting matches...",showspeed=true)
    iter = 1:sz
    animation = @animate for i in iter
        plt = Plots.plot(1,xlims=(0,2048), ylims=(0,1536), xlabel="x",ylabel="y",title="Prospective Projection of the 3D Grid", dpi=400,label="")
        Plots.plot!(p[i],q[i], legend=true, labels="Simulation",  dpi=:400)
        Plots.plot!(pObs[i],qObs[i], labels="Observation",  dpi=:400)
        iterj = 1:length(Exptx[i])
        for j in iterj
            Plots.plot!([Exptx[i][j], Obsptx[i][j]], [Expty[i][j], Expty[i][j]], marker=1.5, lw=0.5, dpi=:400, label="")
        end
        Plots.xlims!(0,2048)
        Plots.ylims!(0,1536) 
        Plots.xlabel!(L"x")
        Plots.ylabel!(L"y")
        Plots.title!("Point Correspondence")
        next!(pr)
    end
    gif(animation, string(filepath,"/matches_h.gif"), fps=10)
end

function set_plot(fs::Int)

    plt = Plots.plot(1, 
                            xtickfontsize = fs, ytickfontsize = fs,
                            titlefontsize = fs,
                            xguidefontsize = fs,
                            yguidefontsize = fs,
                            legendfontsize = fs-2, 

                            size = (1000,750), 
                            fontfamily = "computer modern",
                            framestyle = :box, 
                            margin = 6mm, 

                            grid = :true, 
                            minorgrid = :true, 
                            lw = 3,
                            
                            label="")
    return plt
end
