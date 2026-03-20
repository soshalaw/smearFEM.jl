using WriteVTK
using ProgressMeter
using DelimitedFiles
using JSON
using Distributions
using HDF5
using CSV
using DataFrames

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

    cells = [MeshCell(cellType,IEN[:,e]) for e in 1:size(IEN,2)]

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
                vtk["p"] = fields[i]
                time = (i - 1)
                pvd[time] = vtk
            end
        end
    end
end 

"""
    write_stokes_scene(filepath, NodeList_u, IEN_u, NodeList_p, IEN_p, ne, ndim, velocities, pressures)

Write Stokes results when velocity uses a Q2 mesh and pressure uses a Q1 mesh.

# Arguments:
- `filepath::String`: Name and path to the VTK output directory.
- `NodeList_u`: Q2 velocity node coordinates (single array or vector of arrays per timestep).
- `IEN_u::Matrix{nNodes, nElem}`: Q2 velocity connectivity array.
- `NodeList_p`: Q1 pressure node coordinates (single array or vector of arrays per timestep).
- `IEN_p::Matrix{nNodes, nElem}`: Q1 pressure connectivity array.
- `ne::Integer`: Number of elements in x, y, z direction.
- `ndim::Integer`: Number of dimensions.
- `velocities::Vector`: Velocity field per timestep.
- `pressures::Vector`: Pressure field per timestep.
"""
function write_stokes_scene(
    filepath::String,
    NodeList_u,
    IEN_u,
    NodeList_p,
    IEN_p,
    ne::Int64,
    ndim::Int64,
    velocities,
    pressures;
    ID=nothing,
    velocity_name::String="v",
    pressure_name::String="p",
    collection_name_velocity::String="velocity",
    collection_name_pressure::String="pressure"
)

    set_file(string(filepath, "/vtkFiles")) # create the directory to store the VTK files

    if ndim == 1
        cellType_p = VTKCellTypes.VTK_LINE
        cellType_u = VTKCellTypes.VTK_LAGRANGE_CURVE
    elseif ndim == 2
        cellType_p = VTKCellTypes.VTK_QUAD
        cellType_u = VTKCellTypes.VTK_LAGRANGE_QUADRILATERAL
    elseif ndim == 3
        cellType_p = VTKCellTypes.VTK_HEXAHEDRON
        cellType_u = VTKCellTypes.VTK_LAGRANGE_HEXAHEDRON
    else
        throw(ArgumentError("ndim must be 1, 2, or 3."))
    end

    IEN_u = rearrange(ndim, IEN_u)  # rearrange the Q2 solution

    if length(velocities) != length(pressures)
        throw(ArgumentError("velocities and pressures must have the same length."))
    end

    cells_u = [MeshCell(cellType_u, IEN_u[:, e]) for e in 1:size(IEN_u, 2)]
    cells_p = [MeshCell(cellType_p, IEN_p[:, e]) for e in 1:size(IEN_p, 2)]
    fieldIter = 1:length(pressures)

    paraview_collection(string(filepath, "/vtkFiles/", collection_name_velocity)) do pvd
        @showprogress "Writing out velocity to VTK..." for i in fieldIter
            node_list_u = isa(NodeList_u, AbstractVector) ? NodeList_u[i] : NodeList_u
            vtk_grid(string(filepath, "/vtkFiles/velocity_$i"), node_list_u, cells_u) do vtk
                vtk[velocity_name] = velocities[i]
                time = (i - 1)
                pvd[time] = vtk
            end
        end
    end

    paraview_collection(string(filepath, "/vtkFiles/", collection_name_pressure)) do pvd
        @showprogress "Writing out pressure to VTK..." for i in fieldIter
            node_list_p = isa(NodeList_p, AbstractVector) ? NodeList_p[i] : NodeList_p
            vtk_grid(string(filepath, "/vtkFiles/pressure_$i"), node_list_p, cells_p) do vtk
                vtk[pressure_name] = pressures[i]
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
function write_data(filepath::String, data_array::AbstractArray)
    @info "Writing contour files..."
    counter = 0
    set_file(string(filepath))
    for data in data_array
        cStr = string(counter,pad=3)
        if size(data, 1) == 2
            data = hcat(data[1,:], data[2,:])
        elseif size(data, 1) == 3
            data = hcat(data[1,:], data[2,:], data[3,:])
        end
        write_csv(string(filepath,"/",cStr),data)
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
    sparse_mat = JSON.parsefile(filepath)
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
        JSON.print(io, data, 4)
    end
    @info "Data written to $(filepath)"
end

function read_json(filepath::String)
    if !isfile(filepath)
        throw(SystemError("Trying to read from $filepath, the directory does not exist."))
    end
    # Read the JSON file
    data = JSON.parsefile(filepath)
    return data
end

"""
    read_h5(filename, mode)
    
    Function to read mesh data from a .h5 file for IGA.
    
    # Arguments:
    - `filename::String`:Path to the the .h5 file.
    - `mode::String`:Simulation (sim) or test (test) mode. With the test modes the volume of the mesh is returned
    
    # Returns:
    -`CPointList::Vector{Float64}`:Control point list.
    -`W`:
    -`C`:
    -`IEN::Matrix{Float64}{nNodes, nElem}``:Connectivity array.
    -`IEN_top::Matrix{Float64}{nNodes, nElem}``:Connectivity array for the top surface mesh.
    -`IEN_btm::Matrix{Float64}{nNodes, nElem}``:Connectivity array for the bottom surface mesh.
    """
function read_h5(filename::String, mode::String="sim")
    # Open the HDF5 file
    @info "reading from $filename"
    h5file = h5open(filename, "r")

    # Read mesh
    CPointList = read(h5file, "X")  # Control points
    W = read(h5file, "W")  # Control point weights
    IEN = read(h5file, "IEN").+1  # Element connectivity
    C = read(h5file, "C")  # Extraction operators
    C_new = permutedims(C,[2,1,3])
    
    IEN_top = read(h5file, "IEN_back").+1  # Element connectivity
    IEN_btm = read(h5file, "IEN_front").+1  # Element connectivity
    IEN_vis = read(h5file, "IEN_vis").+1  # Element connectivity
    IEN_cp = read(h5file, "IEN_cp").+1  # Element connectivity
    
    C_top = read(h5file, "C_back")
    C_top_new = permutedims(C_top,[2,1,3])
    C_btm = read(h5file, "C_front")
    C_btm_new = permutedims(C_btm,[2,1,3])
    C_vis = read(h5file, "C_vis")
    C_vis_new = permutedims(C_vis,[2,1,3])
    
    if mode == "test"
        vol_BSpline = read(h5file, "BSpline_vol")
        vol_NURBS = read(h5file, "NURBS_vol")
        area_BSpline = read(h5file, "BSpline_area")
        area_NURBS = read(h5file, "NURBS_area")
        
        # Close the HDF5 file after reading
        close(h5file)   
        
        return CPointList, W, C_new, IEN, IEN_top, C_top_new, IEN_btm, C_btm_new, vol_BSpline, vol_NURBS, area_BSpline, area_NURBS
    end
    
    # Close the HDF5 file after reading
    close(h5file)
    
    return CPointList, W, C_new, IEN, IEN_cp, IEN_top, C_top_new, IEN_btm, C_btm_new, IEN_vis, C_vis_new
end

function read_perception_data(filepath::String)
    if !isfile(filepath)
        throw(SystemError("Trying to read from $filepath, the directory does not exist."))
    end
    # Read the JSON file
    h5file = h5open(filepath, "r")
    
    pose = read(h5file, "poses")
    pose_plt_top = read(h5file, "pose_top_plate")
    pose_plt_btm = read(h5file, "pose_btm_plate")

    pose_new = permutedims(pose,[2,1,3])
    pose_top_plt_new = permutedims(pose_plt_top,[2,1,3])
    pose_btm_plt_new = permutedims(pose_plt_btm,[2,1,3])

    return pose_new, pose_top_plt_new, pose_btm_plt_new
end

function get_time_windows(file_path::String)
    df = CSV.read(file_path,DataFrame,header=false)
    windows = dataframe_2_vec(df)
    return windows
end

