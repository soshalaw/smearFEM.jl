module PostProcess
    using LinearAlgebra
    using PyPlot
    using WriteVTK
    using ProgressMeter
    using LazySets
    using Polyhedra
    using DataInterpolations
    using Plots

    function greet()
        println("Hello, I am the PostProcess module")
    end

    function inflate_sphere(NodeList, x0, x1, y0, y1)
        """ Inflate the sphere to a unit sphere 
        
        Parameters:
            NodeList : {[ndims,nNodes], Matrix{Float64}} : The coordinates of the nodes
            x0 : Float64 : The lower bound of the x direction
            x1 : Float64 : The upper bound of the x direction
            y0 : Float64 : The lower bound of the y direction
            y1 : Float64 : The upper bound of the y direction
            
        Returns:
            NodeList : {[ndims,nNodes], Matrix{Float64}} :  The coordinates of the nodes after inflation
        """
        x_center = [0.5*(x0 + x1), 0.5*(y0 + y1)]
    
        for i in 1:size(NodeList,2)
            scale = maximum(abs.(NodeList[1:2,i] - x_center))
            if scale ≈ 0.
                NodeList[1:2,i] = [0 , 0]
            else
                r = sqrt((NodeList[1,i] - x_center[1])^2 + (NodeList[2,i] - x_center[2])^2)
                NodeList[1:2,i] = scale*(NodeList[1:2,i] - x_center)/r
            end
        end
        return NodeList
    end

    function extract_borders(NodeList, CameraMatrix, BorderNodesList, state, ne = nothing)
        """ Project the 3D mesh to 2D image plane and extract the border nodes (left and right)
        
        Args:
            NodeList: {[ndim,nNodes] ::Matrix{Float64}} : coordinates of the nodes
            CameraMatrix: {[3x3] ::Matrix{Float64}} : Camera matrix
            BorderNodesList: {[4xN] ::Vector{Vector{Any}}} : List of border nodes
            state: {::String} : State of the function (init:During the initialization of the mesh or update: when the mesh is updated)
            ne: {::Int} : Number of elements in each direction
        Returns:
            NodeList: {[ndim,nbNodes] ::Matrix{Float64}} : 2D coordinates of the border nodes
        """
        SideNodes = NodeList[:,BorderNodesList[1]]  # extract the border nodes 
    
        SideNodes2D = back_project(SideNodes, CameraMatrix) 

        # project the nodes to the image plane and extract the border nodes as an ordered list
        if state == "init"
            @assert !isnothing(ne) "Number of elements must be provided"

            LeftborderNodes = zeros(2,(ne+1))                  # vector to store indexes of the border nodes
            RightborderNodes = zeros(2,(ne+1))                 # vector to store indexes of the border nodes
            TopLayerList = []                             # vector to store indexes of the border nodes
            BottomLayerList = []                            # vector to store indexes of the border nodes
            szSide = size(SideNodes2D,2)÷(ne+1)                            # size of each layer
            
            for Layers in 1:ne+1                                        # loop through each layer
                nodes = SideNodes2D[:,(Layers-1)*szSide+1:Layers*szSide]
                minNode = (Layers-1)*szSide + argmin(nodes[1,:])
                maxNode = (Layers-1)*szSide + argmax(nodes[1,:])
                LeftborderNodes[:,Layers] = SideNodes2D[:,minNode]         # left border nodes
                RightborderNodes[:,Layers] = SideNodes2D[:,maxNode]        # right border nodes
                if Layers == ne + 1
                    for nodeId in 1:size(nodes,2)
                        if nodes[2,nodeId] > SideNodes2D[2,minNode]    
                            push!(TopLayerList, (Layers-1)*szSide+nodeId)
                        end
                    end 
                elseif Layers == 1
                    for nodeId in 1:size(nodes,2)
                        if nodes[2,nodeId] < SideNodes2D[2,minNode] 
                            push!(BottomLayerList, nodeId)
                        end
                    end 
                end
            end  
            TopLayer = sortslices(SideNodes2D[:,TopLayerList],dims=2)                     # top layer nodes
            BottomLayer = sortslices(SideNodes2D[:,BottomLayerList],dims=2)               # bottom layer nodes

            BorderPoints = hcat(LeftborderNodes, TopLayer, reverse(RightborderNodes,dims=2), reverse(BottomLayer,dims=2))         # concatenate the left and right border nodes

        elseif state == "update"
            p = Array{Vector{Float64}}(undef,0)
    
            for i in 1:size(SideNodes2D,2)
                push!(p, SideNodes2D[:,i])
            end
    
            hull = convex_hull(p)
            BorderPoints = zeros(2,size(hull,1))
            
            for points in 1:size(hull,1)
                BorderPoints[:,points] = hull[points]
            end
        end

        return BorderPoints, SideNodes2D # return the border nodes
    end
    
    function back_project(NodeList, CameraMatrix)
        """ Project the 3D mesh to 2D image plane
        
        Args:
            NodeList: 3D mesh grid
            CameraMatrix: Camera matrix
        
        Returns:
            NodeList: 2D coordinates of the nodes
        """
                    
        # transform point cloud wrt to camera frame 
        R = [1 0 0; 0 0 1; 0 -1 0]     # rotation matrix
        t = [0; -0.5; 3]               # translation vector

        NodeListTrans = R*NodeList .+ t
        
        NodeListNorm = zeros(3,size(NodeListTrans,2))
    
        for i in 1:size(NodeListNorm,2)
            NodeListNorm[1,i] = NodeListTrans[1,i]/NodeListTrans[3,i]
            NodeListNorm[2,i] = NodeListTrans[2,i]/NodeListTrans[3,i]
            NodeListNorm[3,i] = NodeListTrans[3,i]/NodeListTrans[3,i]
        end
        
        NodeListProj = CameraMatrix'*NodeListNorm   # project to image plane
    
        NodeList2D = NodeListProj[1:2,:]            # extract x and y coordinates
    
        return NodeList2D
    end 

    function fit_curve(border)
        """
            Fit a curve to the border nodes of the 2D mesh
        """

        len = size(border,2)
        seq = 1:(len+1)
        x = push!(border[1,:],border[1,1])
        y = push!(border[2,:],border[2,1])
        p = CubicSpline(x,seq)
        q = CubicSpline(y,seq)
        pi = [p(i) for i in 1:0.01:len+1]
        qi = [q(i) for i in 1:0.01:len+1]

    
        return pi, qi
    end

    function noramlize(q, IEN)
        """Function normalize the solution vector for plotting
        
        Parameters:
        q: solution vector
        IEN: IEN array type [nElem, nNodes]
        
        Returns:
        qList: normalized list of solutions 
        """

        qList = zeros(size(IEN))
        max = maximum(q)
        min = minimum(q)
        for e in 1:size(IEN,1)
            for n in 1:4
                qList[e,n] = (q[IEN[e,n]] - min) / (max - min)
            end
        end
        return qList
    end

    function animate(fields, fields2D, BorderNodes2D, IEN, p, q, op, oq)

        animation = @animate for i in 1:length(fields)
            Plots.scatter3d(fields[i][1,:], fields[i][2,:], fields[i][3,:], markersize=2, legend=:false, dpi=:400)
            Plots.xlims!(-1,1)
            Plots.ylims!(-1,1)
            Plots.zlims!(0,1)
            Plots.xlabel!("x")
            Plots.ylabel!("y")
            Plots.zlabel!("z")
            Plots.title!("3D Grid")
        end
    
        animation2 = @animate for i in 1:length(fields2D)
            Plots.plot(p[i],q[i], dpi=:400, legend=:false, color=:green)
            Plots.plot!(op[i],oq[i], dpi=:400, legend=:false, color=:pink)
            Plots.scatter!(fields2D[i][1,:], fields2D[i][2,:], ms=:1, mc=:blue, ma=:0.5, labels="Surface Nodes", dpi=:400)
            Plots.scatter!(BorderNodes2D[i][1,:], BorderNodes2D[i][2,:], ms=:2, mc=:red, labels="Border Nodes", dpi=:400)
            Plots.xlims!(0,2048)
            Plots.ylims!(0,1536) 
            Plots.xlabel!("x")
            Plots.ylabel!("y")
            Plots.title!("Prospective Projection of the 3D Grid")
        end
        gif(animation, "3D_grid.gif", fps=10)
        gif(animation2, "2D_grid.gif", fps=10)
    end 

    function truncate_colormap(minval=0.0, maxval=1.0, n=100)
        """Function to truncate a colormap
        
        Parameters:
        minval: minimum value of the colormap
        maxval: maximum value of the colormap
        n: number of colors
        
        Returns:
        new_cmap: truncated colormap
        
        """
        new_cmap = matplotlib.colors.LinearSegmentedColormap.from_list("mycmap", get_cmap("jet")(collect(range(maxval, minval, n))))
        return new_cmap
    end

    function write_vtk(fileName, fieldName, NodeList, IEN, ne, ndim, q)
        """Function to write the solution to a VTK file
        
        Parameters:
        fileName: name of the VTK file {string}
        NodeList: array of nodes {[nNodes, ndim]}
        IEN: IEN array {[nElem, nNodes]}
        ne: number of elements in each direction {int}
        ndim: number of dimensions {int}
        q: solution vector {Vector{Float64}}
        
        """

        if ndim == 1
            cellType = VTKCellTypes.VTK_LINE
        elseif ndim == 2
            cellType = VTKCellTypes.VTK_QUAD
        elseif ndim == 3
            cellType = VTKCellTypes.VTK_HEXAHEDRON
        end
        
        cells = [MeshCell(cellType,IEN[e,:]) for e in 1:ne^ndim]

        vtk_grid(fileName, NodeList, cells) do vtk
            vtk[fieldName] = q

        end
    end 

    function write_scene(fileName, NodeList, IEN, ne, ndim, fields)
        """Function to write the solution to a VTK file

        Parameters:
        fileName: name of the VTK file {string}
        NodeList: array of nodes {[nNodes, ndim]}
        IEN: IEN array {[nElem, nNodes]}
        ne: number of elements in each direction {int}
        ndim: number of dimensions {int}
        q: solution vector {Vector{Float64}}
        q2: solution vector {Vector{Float64}}

        """
        cellType = VTKCellTypes.VTK_HEXAHEDRON

        cells = [MeshCell(cellType,IEN[e,:]) for e in 1:ne^ndim]

        paraview_collection(string(fileName,"displacement")) do pvd # create a paraview collection
            @showprogress "Writing out to VTK..." for i in 1:length(fields)
                vtk_grid(string(fileName,"timestep_$i"), NodeList, cells) do vtk # write out the fields to VTK
                    vtk["u"] = fields[i]
                    time = (i - 1)
                    pvd[time] = vtk
                end
            end
        end
    end 

    function PlotGrid(IEN, NodeList)

        fig1 = plt.figure()
        println(size(IEN))
        if size(IEN,2) == 4 # 2D element with 4 nodes 
            ax = fig1.add_subplot(111)

            for i in 1:size(IEN,1)
                x = NodeList[1,IEN[i,:]]
                y = NodeList[2,IEN[i,:]]
                ax.plot(x, y, "-k", linewidth=0.5)
            end
        
            ax.scatter(NodeList[1,:],NodeList[2,:],s=10,c="r")
            ax.axis("equal")
            ax.grid("on")
            ax.set_xlabel("x")
            ax.set_ylabel("y")
            ax.set_title("2D Grid")
            
        elseif size(IEN,2) == 8 # 3D element with 8 nodes
            ax = fig1.add_subplot(111, projection="3d")

            for i in 1:size(IEN,1)
                x = NodeList[1,IEN[i,:]]
                y = NodeList[2,IEN[i,:]]
                z = NodeList[3,IEN[i,:]]
                ax.plot(x, y, z,"-k", linewidth=0.5)
            end

            ax.scatter(NodeList[1,:],NodeList[2,:],NodeList[3,:],s=10,c="r")
            ax.axis("equal")
            ax.grid("on")
            ax.set_xlabel("x")
            ax.set_ylabel("y")
            ax.set_zlabel("z")
            ax.set_title("3D Grid")
        end
        gcf()
    end

    function PlotMesh(NodeList, IEN)
        """Function to plot the mesh

        Parameters:
        NodeList: array of nodes {[nNodes, ndim]}
        IEN: IEN array {[nElem, nNodes]}
        
        """
        sz = size(NodeList,1)
        if sz == 2
            fig1 = plt.figure()
            ax = fig1.add_subplot(111)
    
            for i in 1:size(IEN,1)
                x = NodeList[1,IEN[i,:]]
                y = NodeList[2,IEN[i,:]]
                ax.plot(x, y, "-k", linewidth=0.5)
            end
            
            ax.scatter(NodeList[1,:],NodeList[2,:],s=10,c="r")
            ax.axis("equal")
            ax.grid("on")
            ax.set_xlabel("x")
            ax.set_ylabel("y")
        elseif sz == 3
            fig1 = plt.figure()
            ax = fig1.add_subplot(111, projection="3d")
    
            for i in 1:size(IEN,1)
                x = NodeList[1,IEN[i,:]]
                y = NodeList[2,IEN[i,:]]
                z = NodeList[3,IEN[i,:]]
                ax.plot(x, y, z,"-k", linewidth=0.5)
            end
    
            ax.scatter(NodeList[1,:],NodeList[2,:],NodeList[3,:],s=10,c="r")
            ax.axis("equal")
            ax.grid("on")
            ax.set_xlabel("x")
            ax.set_ylabel("y")
            ax.set_zlabel("z")
        end
        gcf()
    end
end # module PostProcess