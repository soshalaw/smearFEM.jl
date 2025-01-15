using WriteVTK
using ProgressMeter
using DelimitedFiles
using JSON3
using HDF5

"""
    write_vtk(fileName, fieldName, NodeList, IEN, ne, ndim, q)

Function to write the solution to a VTK file
    
# Arguments:
- `fileName::String`: name of the VTK file.
- `fieldName::String`:name of the field
- `NodeList::Matrix[nNodes, ndim]`: Spacial coordinates of the nodes of the mesh.
- `IEN::Matrix{Float64}{nNodes, nElem}`: Connectivity array.
- `ne::Integer`: Number of elements in each direction.
- `ndim::Integer`: Number of dimensions.
- `q::Vector{Float64}`: Solution field.
"""
function write_vtk(filePath::String, fieldName::String, NodeList, IEN, ne::Int64, ndim::Int64, q; ID=nothing, FunctionClass::String="Q1")
    if FunctionClass == "Q1"
        if ndim == 1
            cellType = VTKCellTypes.VTK_LINE
        elseif ndim == 2
            cellType = VTKCellTypes.VTK_QUAD
        elseif ndim == 3
            cellType = VTKCellTypes.VTK_HEXAHEDRON
        end
    elseif FunctionClass == "Q2"
        if ndim == 1
            cellType = VTKCellTypes.VTK_LAGRANGE_CURVE
        elseif ndim == 2
            cellType = VTKCellTypes.VTK_LAGRANGE_QUADRILATERAL
        elseif ndim == 3
            cellType = VTKCellTypes.VTK_LAGRANGE_HEXAHEDRON
        end

        IEN = rearrange(ndim, IEN)  # rearrange the solution
    end

    cells = [MeshCell(cellType,IEN[:,e]) for e in 1:ne^ndim]

    vtk_grid(string(filePath,"/vtkFiles/",fieldName), NodeList, cells) do vtk
        vtk[fieldName] = q
    end

end 

"""
    write_scene(fileName, NodeList, IEN, ne, ndim, fields)

Function to write the solution to a VTK file

# Arguments:
- `fileName::String`: Name and path to the VTK file.
- `NodeList::Matrix{Float64}{nNodes, ndim}`: Spacial coordinates of the nodes of the mesh.
- `IEN::Matrix{nElem, nNodes}`: Connectivity array.
- `ne::Integer`: Number of elements in x, y, z direction.
- `ndim::Integer`: Number of dimensions
- `fields::Vector{Vector{Float64}}`: Solution fields.
"""
function write_scene(fileName::String, NodeList, IEN, ne::Int64, ndim::Int64, fields; ID=nothing, FunctionClass::String="Q1")

    if FunctionClass == "Q1"
        if ndim == 1
            cellType = VTKCellTypes.VTK_LINE
        elseif ndim == 2
            cellType = VTKCellTypes.VTK_QUAD
        elseif ndim == 3
            cellType = VTKCellTypes.VTK_HEXAHEDRON
        end
    elseif FunctionClass == "Q2"
        if ndim == 1
            cellType = VTKCellTypes.VTK_LAGRANGE_CURVE
        elseif ndim == 2
            cellType = VTKCellTypes.VTK_LAGRANGE_QUADRILATERAL
        elseif ndim == 3
            cellType = VTKCellTypes.VTK_LAGRANGE_HEXAHEDRON
        end

        IEN = rearrange(ndim, IEN)  # rearrange the solution
    end

    cells = [MeshCell(cellType,IEN[:,e]) for e in 1:ne^ndim]

    fieldIter = 1:length(fields)
    paraview_collection(string(fileName,"/vtkFiles/displacement")) do pvd # create a paraview collection
        @showprogress "Writing out to VTK..." for i in fieldIter
            vtk_grid(string(fileName,"/vtkFiles/timestep_$i"), NodeList, cells) do vtk # write out the fields to VTK
                vtk["u"] = fields[i]
                time = (i - 1)
                pvd[time] = vtk
            end
        end
    end
end 

"""
    read_csv(csv_path)

Function to read the CSV files in the directory and fit a curve to the border observations in the CSV files.

# Arguments:
- `csv_path::String`: path to the directory containing the CSV files.

# Returns:
- `splinep::Vector{Float64}`: x coordinates samples of the spline parameters of the border nodes.
- `splineq::Vector{Float64}`: y coordinates samples of the spline parameters of the border nodes.
"""
function read_csv(csv_path::String; NOISE::Bool=false, nProfile::Int64=1, nFactor=0)   

    csv_files = readdir(csv_path, join=true)        # get the list of the csv files in the directory
    ObsDataList = AbstractArray[]                                  # store the observation data
    splinex = AbstractArray[]
    spliney = AbstractArray[]

    obsData_ = readdlm(csv_files[1], ',', Float64, '\n', header=false)  # read the observation data
    w = zeros(Float64,size(obsData_))
    if NOISE && nProfile == 1
        w = nFactor*randn(size(obsData_))
    end
    for file in csv_files
        obsData = readdlm(file, ',', Float64, '\n', header=false)  # read the observation data
        if NOISE && nProfile == 2
            w = nFactor*randn(size(obsData))
        end
        obsData = obsData + w
        push!(ObsDataList, obsData') # store the transpose of the observation data to fit the comparison function
        push!(splinex, obsData[:,1])
        push!(spliney, obsData[:,2])
    end

    return ObsDataList, splinex, spliney
end

"""
    write_csv(fileName, borders)

Function to write the border data to a CSV file

# Arguments:
- `fileName::String`: name of the CSV file.
- `Data::Vector{Vector{Float64}}`: vector of data.
"""
function write_csv(fileName::String, Data)
    open(string(fileName,".csv"), "w") do io
        writedlm(io, Data, ',')
    end
end

"""
    write_contour_data(fileName, borders)

Function to write the contour data to a CSV file

# Arguments:
- `fileName::String`: name of the file.
- `borders::Vector{Matrix{Float64}}`: vector of border data.
"""
function write_contour_data(fileName::String, borders)
    println("Writing contour files...")
    counter = 1
    for border in borders
        cStr = string(counter,pad=3)
        write_csv(string(fileName,"/contour_data/",cStr),[border[1,:] border[2,:]])
        counter += 1
    end
end

"""
    set_file(filepath)

Create the directories to store the data

# Arguments:
- `filepath::String` : path to the file
"""
function set_file(filepath)
    if !isdir(filepath)
        mkpath(filepath)
        mkdir(string(filepath,"/Results"))
        mkdir(string(filepath,"/Results/contour_data"))
        mkdir(string(filepath,"/Results/vtkFiles"))
        mkdir(string(filepath,"/Results/images"))
        mkdir(string(filepath,"/Results/cost"))
    end
end

"""
    read_sparse_mat(filename)

Function to read the sparse matrix from a JSON file

# Arguments:
- `filename::String`: name of the JSON file.

# Returns:
- `sparse_mat::SparseArrays`: sparse matrix.
"""
function read_sparse_mat(filename)
    sparse_mat = JSON3.read(filename)
    row_ptr = [i + 1 for i in sparse_mat[:row]] # -1 for Python numbering
    col_ptr = [i + 1 for i in sparse_mat[:col]]
    values = sparse_mat[:data]
    return sparse(col_ptr, row_ptr, values)
end

"""
    write_json(filename, sparse_mat)

Function to write data to a JSON file

# Arguments:
- `filename::String`: name of the JSON file.
- `data::Dict`: dictionary to write to the JSON file.
"""
function write_json(filename, data)
    open(string(filename,".json"), "w") do io
        JSON3.pretty(io, data)
    end
end

"""
    read_h5(filename)

Function to read mesh data from a .h5 file for IGA.

# Arguments:
- `filename::String`:Path to the the .h5 file.

# Returns:
-`CPointList::Vector{Float64}`:Control point list.
-`W`:
-`C`:
-`IEN::Matrix{Float64}{nNodes, nElem}``:Connectivity array.
-`IEN_top::Matrix{Float64}{nNodes, nElem}``:Connectivity array for the top surface mesh.
-`IEN_btm::Matrix{Float64}{nNodes, nElem}``:Connectivity array for the bottom surface mesh.
"""
function read_h5(filename::String)
    # Open the HDF5 file
    h5file = h5open(filename, "r")

    # Read mesh
    CPointList = read(h5file, "X")  # Control points
    W = read(h5file, "W")  # Control point weights
    IEN = read(h5file, "IEN").+1  # Element connectivity

    # IEN_top = read(h5file, "IEN_top").+1  # Element connectivity
    # IEN_btm = read(h5file, "IEN_btm").+1  # Element connectivity

    C = read(h5file, "C")  # Extraction operators
    vol_BSpline = read(h5file, "BSpline_vol")
    vol_NURBS = read(h5file, "NURBS_vol")
    C_new = permutedims(C,[2,1,3])

    # Close the HDF5 file after reading
    close(h5file)

    return CPointList, W, C_new, IEN, vol_BSpline, vol_NURBS #, IEN_top, IEN_btm
end