module PostProcess
    using LinearAlgebra
    using PyPlot
    using WriteVTK

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
end # module PostProcess