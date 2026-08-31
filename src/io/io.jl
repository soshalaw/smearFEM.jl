using WriteVTK
using ProgressMeter
using DelimitedFiles
using JSON
using Distributions
using HDF5
using CSV
using DataFrames
using JLD2
using FileIO

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
function write_vtk(filePath::String, fieldName::String, NodeList, IEN, ne::Int64, ndim::Int64, q; ID=nothing, element_shape::Symbol=:Hex, basis_order::Int=1)

    set_file(string(filePath,"/vtkFiles")) # create the directory to store the VTK files
    if element_shape == :Tet
        if basis_order == 1
            cellType = ndim == 2 ? VTKCellTypes.VTK_TRIANGLE : VTKCellTypes.VTK_TETRA
        else
            cellType = ndim == 2 ? VTKCellTypes.VTK_LAGRANGE_TRIANGLE : VTKCellTypes.VTK_LAGRANGE_TETRAHEDRON
        end
    elseif basis_order == 1
        if ndim == 1
            cellType = VTKCellTypes.VTK_LINE
        elseif ndim == 2
            cellType = VTKCellTypes.VTK_QUAD
        elseif ndim == 3
            cellType = VTKCellTypes.VTK_HEXAHEDRON
        end
    elseif basis_order == 2
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
function write_scene(filepath::String, node_list_list, IEN, ndim::Int64, fields; element_shape::Symbol=:Hex, basis_order::Int=1)

    set_file(string(filepath,"/vtkFiles")) # create the directory to store the VTK files
    if element_shape == :Tet
        if basis_order == 1
            if ndim == 1
                cellType = VTKCellTypes.VTK_LINE
            elseif ndim == 2
                cellType = VTKCellTypes.VTK_TRIANGLE
            elseif ndim == 3
                cellType = VTKCellTypes.VTK_TETRA
            end
        else
            if ndim == 1
                cellType = VTKCellTypes.VTK_LAGRANGE_CURVE
            elseif ndim == 2
                cellType = VTKCellTypes.VTK_LAGRANGE_TRIANGLE
            elseif ndim == 3
                cellType = VTKCellTypes.VTK_LAGRANGE_TETRAHEDRON
            end
        end
    elseif basis_order == 1
        if ndim == 1
            cellType = VTKCellTypes.VTK_LINE
        elseif ndim == 2
            cellType = VTKCellTypes.VTK_QUAD
        elseif ndim == 3
            cellType = VTKCellTypes.VTK_HEXAHEDRON
        end
    elseif basis_order == 2
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
    e_iter = 1:size(IEN, 2)
    cells = [MeshCell(cellType,IEN[:,e]) for e in e_iter]
    fieldIter = 1:length(fields)
    paraview_collection(string(filepath,"/vtkFiles/displacement")) do pvd # create a paraview collection
        @showprogress "Writing out to VTK..." for i in fieldIter
            vtk_grid(string(filepath,"/vtkFiles/timestep_$i"), node_list_list[i], cells) do vtk # write out the fields to VTK
                vtk["p"] = fields[i]
                time = (i - 1)
                pvd[time] = vtk
            end
        end
    end
end 

function extract_p_from_u_nodes(NodeList_u, NodeList_p, IEN_p)
    # Map each pressure node to its closest velocity node by Euclidean distance.
    _ = IEN_p  # kept for API compatibility

    node_axis_u = (size(NodeList_u, 1) <= 3 && size(NodeList_u, 2) > size(NodeList_u, 1)) ? :cols : :rows
    node_axis_p = (size(NodeList_p, 1) <= 3 && size(NodeList_p, 2) > size(NodeList_p, 1)) ? :cols : :rows

    nNodes_u = node_axis_u == :cols ? size(NodeList_u, 2) : size(NodeList_u, 1)
    nNodes_p = node_axis_p == :cols ? size(NodeList_p, 2) : size(NodeList_p, 1)

    p_to_u_idx = Vector{Int}(undef, nNodes_p)
    min_dist2 = Vector{Float64}(undef, nNodes_p)

    @views for p_idx in 1:nNodes_p
        p_node = node_axis_p == :cols ? NodeList_p[:, p_idx] : NodeList_p[p_idx, :]
        best_u_idx = 1
        best_dist2 = Inf

        for u_idx in 1:nNodes_u
            u_node = node_axis_u == :cols ? NodeList_u[:, u_idx] : NodeList_u[u_idx, :]
            dist2 = sum((p_node .- u_node) .^ 2)
            if dist2 < best_dist2
                best_dist2 = dist2
                best_u_idx = u_idx
            end
        end

        p_to_u_idx[p_idx] = best_u_idx
        min_dist2[p_idx] = best_dist2
    end

    node_list_p = node_axis_u == :cols ? NodeList_u[:, p_to_u_idx] : NodeList_u[p_to_u_idx, :]

    # Sanity check: alignment quality between original pressure nodes and mapped velocity nodes.
    min_dists = sqrt.(min_dist2)
    max_dist = maximum(min_dists)
    mean_dist = sum(min_dists) / length(min_dists)
    tol = 1e-8
    if max_dist > tol
        @warn "Nearest-neighbor mapping sanity check: max distance = $(max_dist), mean distance = $(mean_dist)."
    else
        @info "Nearest-neighbor mapping sanity check passed: max distance = $(max_dist), mean distance = $(mean_dist)."
    end

    println("Extracted $(size(node_list_p)) pressure nodes from $(size(NodeList_u)) velocity nodes using nearest-neighbor mapping.")

    return node_list_p, p_to_u_idx
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
    collection_name_pressure::String="pressure",
    pos3D::Vector{AbstractArray}=zeros(Float64, 0, 0),
    element_shape_u::Symbol=:Hex,
    basis_order_u::Int=2,
    element_shape_p::Symbol=:Hex,
    basis_order_p::Int=1
)

    set_file(string(filepath, "/vtkFiles")) # create the directory to store the VTK files

    if ndim == 1
        cellType_p = VTKCellTypes.VTK_LINE
        cellType_u = VTKCellTypes.VTK_LAGRANGE_CURVE
    elseif ndim == 2
        if element_shape_u == :Tet
            cellType_u = basis_order_u == 1 ? VTKCellTypes.VTK_TRIANGLE : VTKCellTypes.VTK_LAGRANGE_TRIANGLE
        else
            cellType_u = basis_order_u == 1 ? VTKCellTypes.VTK_QUAD : VTKCellTypes.VTK_LAGRANGE_QUADRILATERAL
        end
        cellType_p = element_shape_p == :Tet ? VTKCellTypes.VTK_TRIANGLE : VTKCellTypes.VTK_QUAD
    elseif ndim == 3
        if element_shape_u == :Tet
            cellType_u = basis_order_u == 1 ? VTKCellTypes.VTK_TETRA : VTKCellTypes.VTK_LAGRANGE_TETRAHEDRON
        else
            cellType_u = basis_order_u == 1 ? VTKCellTypes.VTK_HEXAHEDRON : VTKCellTypes.VTK_LAGRANGE_HEXAHEDRON
        end
        cellType_p = element_shape_p == :Tet ? VTKCellTypes.VTK_TETRA : VTKCellTypes.VTK_HEXAHEDRON
    else
        throw(ArgumentError("ndim must be 1, 2, or 3."))
    end

    if element_shape_u == :Hex && basis_order_u == 2
        IEN_u = rearrange(ndim, IEN_u)
    end

    if length(velocities) != length(pressures)
        throw(ArgumentError("velocities and pressures must have the same length."))
    end

    cells_u = [MeshCell(cellType_u, IEN_u[:, e]) for e in 1:size(IEN_u, 2)]
    cells_p = [MeshCell(cellType_p, IEN_p[:, e]) for e in 1:size(IEN_p, 2)]
    fieldIter = 1:length(pressures)

    
    node_list = length(pos3D) != 0 ? pos3D : NodeList_u
    paraview_collection(string(filepath, "/vtkFiles/", collection_name_velocity)) do pvd
        @showprogress "Writing out velocity to VTK..." for i in fieldIter
            node_list = isa(node_list, AbstractVector) ? node_list[i] : node_list
            vtk_grid(string(filepath, "/vtkFiles/velocity_$i"), node_list, cells_u) do vtk
                vtk[velocity_name] = velocities[i]
                time = (i - 1)
                pvd[time] = vtk
            end
        end
    end
    
    _node_list_p, p_to_u_idx = extract_p_from_u_nodes(NodeList_u, NodeList_p, IEN_p) 
    paraview_collection(string(filepath, "/vtkFiles/", collection_name_pressure)) do pvd
        @showprogress "Writing out pressure to VTK..." for i in fieldIter
            node_list_p = isa(pos3D, AbstractVector) ? pos3D[i] : NodeList_p
            node_list_p = length(pos3D) != 0 ? node_list_p[:, p_to_u_idx] : node_list_p
            vtk_grid(string(filepath, "/vtkFiles/pressure_$i"), node_list_p, cells_p) do vtk
                vtk[pressure_name] = pressures[i]
                time = (i - 1)
                pvd[time] = vtk
            end
        end
    end

    if length(pos3D) != 0
        paraview_collection(string(filepath, "/vtkFiles/geometry")) do pvd
            @showprogress "Writing out geometry to VTK..." for i in fieldIter
                vtk_grid(string(filepath, "/vtkFiles/geometry_$i"), pos3D[i], cells_u) do vtk
                    time = (i - 1)
                    pvd[time] = vtk
                end
            end
        end
    end

    @info "Finished writing VTK files to $(filepath)/vtkFiles"
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
    set_file(filepath)
    for (t, t_data) in enumerate(data_array)
        cStr = string(t - 1, pad=3)
        write_csv(joinpath(filepath, cStr), t_data)
    end
end

function write_2d_data(filepath::String, data_array::AbstractArray)
    @info "Writing contour files..."
    root_folder = dirname(filepath)
    folder_name = basename(filepath)
    for (t, t_data) in enumerate(data_array)
        cStr = string(t - 1, pad=3)
        for (a, angle_data) in enumerate(t_data)
            angle_dir = joinpath(root_folder, "view_$a", folder_name)
            set_file(angle_dir)
            write_csv(joinpath(angle_dir, cStr), angle_data)
        end
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

Function to read the sparse matrix from a JLD2 file

# Arguments:
- `filepath::String`: path to the file (with or without .jld2 extension).

# Returns:
- `sparse_mat::SparseMatrixCSC`: sparse matrix.
"""
function read_sparse_mat(filepath::String)
    filepath_with_ext = endswith(filepath, ".jld2") ? filepath : string(filepath, ".jld2")
    if !isfile(filepath_with_ext)
        throw(SystemError("File not found: $filepath_with_ext"))
    end
    @load filepath_with_ext sparse_mat
    return sparse_mat
end

"""                                                                                                                                                                    
    write_json(filepath, data)

Function to write data to a JLD2 file (preserves all Julia types and structures)

# Arguments:
- `filepath::String`: path to the file (with .json or .jld2 extension, or without extension).
- `data::Dict`: dictionary with mixed types (matrices, scalars, strings, etc.).

# Example:
```julia
data = Dict("matrix" => rand(3,3), "scalar" => 42, "string" => "hello")
write_json("results/mydata", data)  # saves as mydata.jld2
write_json("results/mydata.json", data)  # saves as mydata.json
```
"""
function write_json(filepath::String, data::Dict)
    set_file(dirname(filepath))
    # Add extension if not present
    if !endswith(filepath, ".jld2") && !endswith(filepath, ".json")
        filepath_with_ext = string(filepath, ".jld2")
    else
        filepath_with_ext = filepath
    end
    @save filepath_with_ext data
    @info "Data written to $(filepath_with_ext)"
end

"""
    read_json(filepath)

Function to read data from a JLD2 or JSON file (auto-detects format)

# Arguments:
- `filepath::String`: path to the file (with .json or .jld2 extension, or without extension).

# Returns:
- `data::Dict`: dictionary with original Julia types preserved.

# Example:
```julia
data = read_json("results/mydata")  # reads mydata.jld2 or mydata.json
data = read_json("results/mydata.json")  # reads mydata.json (JSON format)
```
"""
function read_json(filepath::String)
    # Check for file with various extensions
    filepath_with_ext = filepath
    if !isfile(filepath_with_ext)
        if !endswith(filepath, ".jld2") && !endswith(filepath, ".json")
            # Try .jld2 first
            filepath_with_ext = string(filepath, ".jld2")
            if !isfile(filepath_with_ext)
                # Try .json
                filepath_with_ext = string(filepath, ".json")
            end
        end
    end
    
    if !isfile(filepath_with_ext)
        throw(SystemError("File not found: $filepath"))
    end
    
    # Detect file format by extension or content
    if endswith(filepath_with_ext, ".json")
        # Read as JSON
        data = JSON.parsefile(filepath_with_ext)
        return data
    else
        # Read as JLD2 - use load() which returns a dictionary
        try
            file_data = load(filepath_with_ext)
            # If the file was saved with @save data, extract the data key
            if haskey(file_data, "data")
                return file_data["data"]
            else
                # If no "data" key, return the entire dictionary
                return file_data
            end
        catch e
            throw(SystemError("Failed to read JLD2 file $filepath_with_ext: $(string(e))"))
        end
    end
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
        throw(SystemError("File not found: $filepath"))
    end
    # Read the HDF5 file
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

