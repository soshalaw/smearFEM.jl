using WriteVTK
using ProgressMeter
using DelimitedFiles
using JSON3
using Distributions
using HDF5

"""
    write_vtk(filepath, fieldName, NodeList, IEN, ne, ndim, q)

Function to write the solution to a VTK file
    
# Arguments:
- `filepath::String`: name of the VTK file.
- `fieldName::String`:name of the field
- `NodeList::Matrix[nNodes, ndim]`: Spacial coordinates of the nodes of the mesh.
- `IEN::Matrix{Float64}{nNodes, nElem}`: Connectivity array.
- `ne::Integer`: Number of elements in each direction.
- `ndim::Integer`: Number of dimensions.
- `q::Vector{Float64}`: Solution field.
"""
function write_vtk(filePath::String, fieldName::String, NodeList, IEN, ne::Int64, ndim::Int64, q; ID=nothing, FunctionClass::String="Q1")

    set_file(string(filePath,"/vtkFiles")) # create the directory to store the VTK files
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
    write_scene(filepath, NodeList, IEN, ne, ndim, fields)

Function to write the solution to a VTK file

# Arguments:
- `filepath::String`: Name and path to the VTK file.
- `NodeList::Matrix{Float64}{nNodes, ndim}`: Spacial coordinates of the nodes of the mesh.
- `IEN::Matrix{nElem, nNodes}`: Connectivity array.
- `ne::Integer`: Number of elements in x, y, z direction.
- `ndim::Integer`: Number of dimensions
- `fields::Vector{Vector{Float64}}`: Solution fields.
"""
function write_scene(filepath::String, NodeList, IEN, ne::Int64, ndim::Int64, fields; ID=nothing, FunctionClass::String="Q1")

    set_file(string(filepath,"/vtkFiles")) # create the directory to store the VTK files
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

    # Create a VTK collection
    cells = [MeshCell(cellType,IEN[:,e]) for e in 1:ne^ndim]

    fieldIter = 1:length(fields)
    paraview_collection(string(filepath,"/vtkFiles/displacement")) do pvd # create a paraview collection
        @showprogress "Writing out to VTK..." for i in fieldIter
            vtk_grid(string(filepath,"/vtkFiles/timestep_$i"), NodeList[i], cells) do vtk # write out the fields to VTK
                vtk["u"] = fields[i]
                time = (i - 1)
                pvd[time] = vtk
            end
        end
    end
end 

"""
    read_csv(filepath)

Function to read the CSV files in the directory and fit a curve to the border observations in the CSV files.

# Arguments:
- `filepath::String`: path to the directory containing the CSV files.

# Returns:
- `ObsDataList::Vector{Matrix{Float64}}`: list of observation data.
- `splinep::Vector{Vector{Float64}}`: x coordinates samples of the spline parameters of the border nodes.
- `splineq::Vector{Vector{Float64}}`: y coordinates samples of the spline parameters of the border nodes.
"""
function read_csv(filepath::String)   
    if !isdir(filepath)
        throw(SystemError("Trying to read from $filepath, the directory does not exist."))
    end

    csv_files = readdir(filepath, join=true)        # get the list of the csv files in the directory
    ObsDataList = AbstractArray[]                                  # store the observation data
    splinex = AbstractArray[]
    spliney = AbstractArray[]

    for file in csv_files
        if !endswith(file, ".csv")  # check if the file is a CSV file
            continue
        end
        obsData = readdlm(file, ',', Float64, '\n', header=false)  # read the observation data
        push!(ObsDataList, obsData') # store the transpose of the observation data to fit the comparison function
        push!(splinex, obsData[:,1])
        push!(spliney, obsData[:,2])
    end

    return ObsDataList, splinex, spliney
end

"""
    write_csv(filepath, borders)

Function to write the border data to a CSV file

# Arguments:
- `filepath::String`: name of the CSV file.
- `Data::Vector{Vector{Float64}}`: vector of data.
"""
function write_csv(filepath::String, Data)
    set_file(dirname(filepath)) # create the directory to store the CSV files
    open(string(filepath,".csv"), "w") do io
        writedlm(io, Data,',')
    end
end

"""
    write_contour_data(filepath, borders)

Function to write the contour data to a CSV file

# Arguments:
- `filepath::String`: name of the file.
- `borders::Vector{Matrix{Float64}}`: vector of border data.
"""
function write_data(filepath::String, borders)
    @info "Writing contour files..."
    counter = 1
    set_file(string(filepath))
    for border in borders
        cStr = string(counter,pad=3)
        if size(border, 1) == 2
            border = hcat(border[1,:], border[2,:])
        elseif size(border, 1) == 3
            border = hcat(border[1,:], border[2,:], border[3,:])
        end
        write_csv(string(filepath,"/",cStr),border)
        counter += 1
    end
end

"""
    set_file(filepath)

Create the directories to store the data

# Arguments:
- `filepath::String` : path to the file
"""
function set_file(filepath::String)
    if filepath == "None"
        throw(ArgumentError("File path not defined."))
    elseif !isdir(filepath)
        @info "Specified file path $(filepath) not found, creating directories ..."
        mkpath(filepath)
    end
end

"""
    read_sparse_mat(filepath)

Function to read the sparse matrix from a JSON file

# Arguments:
- `filepath::String`: name of the JSON file.

# Returns:
- `sparse_mat::SparseArrays`: sparse matrix.
"""
function read_sparse_mat(filepath::String)
    if !isfile(filepath)
        throw(SystemError("Trying to read from $filepath, the directory does not exist."))
    end
    # Read the JSON file
    sparse_mat = JSON3.read(filepath)
    row_ptr = [i + 1 for i in sparse_mat[:row]] # -1 for Python numbering
    col_ptr = [i + 1 for i in sparse_mat[:col]]
    values = sparse_mat[:data]
    return sparse(col_ptr, row_ptr, values)
end

"""
    write_json(filepath, sparse_mat)

Function to write data to a JSON file

# Arguments:
- `filepath::String`: name of the JSON file.
- `data::Dict`: dictionary to write to the JSON file.
"""
function write_json(filepath::String, data::Dict)
    set_file(dirname(filepath))
    open(string(filepath,".json"), "w") do io
        JSON3.pretty(io, data)
    end
    @info "Data written to $(filepath)"
end

function read_json(filepath::String)
    if !isfile(filepath)
        throw(SystemError("Trying to read from $filepath, the directory does not exist."))
    end
    # Read the JSON file
    data = JSON3.read(filepath)
    return data
end
"""
    read_h5(filepath, mode)

Function to read mesh data from a .h5 file for IGA.

# Arguments:
- `filepath::String`:Path to the the .h5 file.
- `mode::String`:Simulation (sim) or test (test) mode. With the test modes the volume of the mesh is returned

# Returns:
-`CPointList::Vector{Float64}`:Control point list.
-`W`:
-`C`:
-`IEN::Matrix{Float64}{nNodes, nElem}``:Connectivity array.
-`IEN_top::Matrix{Float64}{nNodes, nElem}``:Connectivity array for the top surface mesh.
-`IEN_btm::Matrix{Float64}{nNodes, nElem}``:Connectivity array for the bottom surface mesh.
"""
function read_h5(filepath::String, mode::String="sim")
    # Open the HDF5 file
    if !isfile(filepath)
        throw(SystemError("Trying to read from $filepath, the directory does not exist."))
    elseif mode != "sim" && mode != "test"
        throw(ArgumentError("Mode not defined."))
    end
    
    h5file = h5open(filepath, "r")

    # Read mesh
    CPointList = read(h5file, "X")  # Control points
    W = read(h5file, "W")  # Control point weights
    IEN = read(h5file, "IEN").+1  # Element connectivity
    C = read(h5file, "C")  # Extraction operators
    C_new = permutedims(C,[2,1,3])

    if mode == "sim"
        IEN_top = read(h5file, "IEN_top").+1  # Element connectivity
        IEN_btm = read(h5file, "IEN_bottom").+1  # Element connectivity
        C_top = read(h5file, "C_top")
        C_top_new = permutedims(C_top,[2,1,3])
        C_btm = read(h5file, "C_bottom")
        C_btm_new = permutedims(C_btm,[2,1,3])

        # Close the HDF5 file after reading
        close(h5file)

        return CPointList, W, C_new, IEN, IEN_top, C_top_new, IEN_btm, C_btm_new
    elseif mode == "test"
        vol_BSpline = read(h5file, "BSpline_vol")
        vol_NURBS = read(h5file, "NURBS_vol")

        # Close the HDF5 file after reading
        close(h5file)   

        return CPointList, W, C_new, IEN, vol_BSpline, vol_NURBS
    end
end