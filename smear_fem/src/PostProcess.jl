module PostProcess
    using LinearAlgebra
    using PyPlot
    using WriteVTK

    function greet()
        println("Hello, I am the PostProcess module")
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

end # module PostProcess