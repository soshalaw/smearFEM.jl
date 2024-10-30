using Plots
using ProgressMeter

"""
    PlotGrid(IEN, NodeList)

Function to plot the grid

# Arguments:
- `IEN::Matrix{Float64}{nElem, nNodes}`: IEN array.
- `NodeList::Matrix{Float64}{nNodes, ndim}`: array of nodes.
"""
function PlotGrid(IEN, NodeList)

    fig1 = plt.figure()
    println(size(IEN))
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
        Plots.scatter3d(NodeList[1,:], NodeList[2,:], NodeList[3,:], markersize=2, label="", dpi=:400)
        iter = 1:size(IEN,1)
        for i in iter
            x = NodeList[1,IEN[i,:]]
            y = NodeList[2,IEN[i,:]]
            z = NodeList[3,IEN[i,:]]
            Plots.plot3d!(x, y, z,marker=1.5, lw=0.5, label="", dpi=:400)
        end

        Plots.xlabel!("x")
        Plots.ylabel!("y")
        Plots.zlabel!("z")
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
function animate_fields(; filepath=nothing, fields=nothing , IEN=nothing, BorderNodes2D=nothing, fields2D=nothing, p=nothing, q=nothing, pObs=nothing, qObs=nothing)

    if isnothing(fields) && isnothing(fields2D) && isnothing(BorderNodes2D) && isnothing(IEN) && isnothing(p) && isnothing(q)
        AssertionError("No fields provided")
        return
    elseif isnothing(fields)
        animate2D(BorderNodes2D=BorderNodes2D, fields2D=fields2D, p=p, q=q, pObs=pObs, qObs=qObs, filepath=filepath)
        return
    elseif isnothing(fields2D) && isnothing(BorderNodes2D) && isnothing(p) && isnothing(q)
        animate3D(fields, filepath=filepath)
        return
    else
        animate3D(fields, filepath=filepath)
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
    
    if isnothing(BorderNodes2D) && isnothing(fields2D) && isnothing(p) && isnothing(q)
        AssertionError("No fields provided")
        return
    else
        if isnothing(fields2D) || isnothing(BorderNodes2D)
            sz = length(p)
        elseif isnothing(BorderNodes2D) || isnothing(p)
            sz = length(fields2D)
        elseif isnothing(fields2D) || isnothing(p)    
            sz = length(BorderNodes2D)
        else
            sz = length(fields2D)
        end

        pr = Progress(sz; desc="Animating 2D fields...",showspeed=true)
        iter = 1:sz
        animation2 = @animate for i in iter
            plt = Plots.plot(1,xlims=(0,2048), ylims=(0,1536), xlabel="x",ylabel="y",title="Prospective Projection of the 3D Grid", label="", dpi=400)
            if !isnothing(p)   
                Plots.plot!(p[i],q[i], legend=true, labels="Simulation",  dpi=:400)
            end
            if !isnothing(pObs)
                Plots.plot!(pObs[i],qObs[i], labels="Observation",  dpi=:400)
            end
            if !isnothing(BorderNodes2D)
                Plots.scatter!(BorderNodes2D[i][1,:], BorderNodes2D[i][2,:], ms=:2, mc=:red, legend=true, labels="Border Nodes", dpi=:400)
            end
            if !isnothing(fields2D)
                Plots.scatter!(fields2D[i][1,:], fields2D[i][2,:], ms=:1, mc=:blue, ma=:0.5, legend=true, labels="Surface Nodes", dpi=:400)
            end
            Plots.xlims!(0,2048)
            Plots.ylims!(0,1536) 
            Plots.xlabel!("x")
            Plots.ylabel!("y")
            Plots.title!("Prospective Projection of the 3D Grid")
            next!(pr)
        end

        gif(animation2, string(filepath,"/2D_grid.gif"), fps=10)
    end
end

"""
    animate3D(fields)

Function to animate the 3D fields as a gif

# Arguments:
- `fields::Vector{Vector{Float64}}`: solution vector
"""
function animate3D(fields; filepath="images/3D_grid.gif")
    sz = length(fields)
    pr = Progress(sz; desc="Animating 3D fields...",showspeed=true)
    iter = 1:sz
    animation = @animate for i in iter
        Plots.scatter3d(fields[i][1,:], fields[i][2,:], fields[i][3,:], markersize=2, label=:"", dpi=:400)
        Plots.xlims!(-1,1)
        Plots.ylims!(-1,1)
        Plots.zlims!(0,1)
        Plots.xlabel!("x")
        Plots.ylabel!("y")
        Plots.zlabel!("z")
        Plots.title!("3D Grid")
        next!(pr)
    end

    gif(animation, string(filepath,"/3D_grid.gif"), fps=10)
end

function plot_matches(simborderfields, p, q, pObs, qObs, pairs, filepath="images")
    sz = length(simborderfields)
    pr = Progress(sz; desc="Plotting matches...",showspeed=true)
    iter = 1:sz
    animation = @animate for i in iter
        borderp = simborderfields[i][1,:]
        borderq = simborderfields[i][2,:]
        plt = Plots.plot(1,xlims=(0,2048), ylims=(0,1536), xlabel="x",ylabel="y",title="Prospective Projection of the 3D Grid", label="", dpi=400)
        Plots.plot!(p[i],q[i], legend=true, labels="Simulation",  dpi=:400)
        Plots.plot!(pObs[i],qObs[i], labels="Observation",  dpi=:400)
        for pair in pairs
            Plots.plot!([borderp[pair[1]], pObs[i][pair[2]]], [borderq[pair[1]], qObs[i][pair[2]]], marker=1.5, lw=0.5, label="", dpi=:400)
        end
        Plots.xlims!(0,2048)
        Plots.ylims!(0,1536) 
        Plots.xlabel!("x")
        Plots.ylabel!("y")
        Plots.title!("Point Correspondence")
        next!(pr)
    end
    gif(animation, string(filepath,"/matches.gif"), fps=10)
end

function plot_matches_h(Exptx, Expty, Obsptx, p, q, pObs, qObs, filepath="images")
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
        Plots.xlabel!("x")
        Plots.ylabel!("y")
        Plots.title!("Point Correspondence")
        next!(pr)
    end
    gif(animation, string(filepath,"/matches_h.gif"), fps=10)
end
