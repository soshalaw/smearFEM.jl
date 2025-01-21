using WriteVTK
using ProgressMeter
using DelimitedFiles
using JSON3
using Distributions

"""
    write_vtk(fileName, fieldName, NodeList, IEN, ne, ndim, q)

Function to write the solution to a VTK file
    
# Arguments:
- `fileName::String`: name of the VTK file.
- `NodeList::Matrix[nNodes, ndim]`: array of nodes.
- `IEN::Matrix{Float64}{nElem, nNodes}`: IEN array.
- `ne::Integer`: number of elements in each direction.
- `ndim::Integer`: number of dimensions.
- `q::Vector{Float64}`: vector of solution fields.
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

        IEN = rearrange(ne, ndim, IEN, ID)  # rearrange the solution
    end

    cells = [MeshCell(cellType,IEN[e,:]) for e in 1:ne^ndim]

    vtk_grid(string(filePath,"/vtkFiles",fieldName), NodeList, cells) do vtk
        vtk[fieldName] = q
    end

end 

"""
    write_scene(fileName, NodeList, IEN, ne, ndim, fields)

Function to write the solution to a VTK file

# Arguments:
- `fileName::String`: name of the VTK file.
- `NodeList::Matrix{Float64}{nNodes, ndim}`: array of nodes.
- `IEN::Matrix{nElem, nNodes}`: IEN array
- `ne::Integer`: number of elements in each direction.
- `ndim::Integer`: number of dimensions
- `fields::Vector{Vector{Float64}}`: vector of solution fields.
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

        IEN = rearrange(ne, ndim, IEN, ID)  # rearrange the solution
    end

    cells = [MeshCell(cellType,IEN[:,e]) for e in 1:ne^ndim]

    paraview_collection(string(fileName,"/vtkFiles/displacement")) do pvd # create a paraview collection
        @showprogress "Writing out to VTK..." for i in 1:length(fields)
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
- `splinep`: x coordinates samples of the spline parameters of the border nodes.
- `splineq`: y coordinates samples of the spline parameters of the border nodes.
"""
function read_csv(csv_path::String; nFactor=0)   

    csv_files = readdir(csv_path, join=true)        # get the list of the csv files in the directory
    ObsDataList = AbstractArray[]                                  # store the observation data
    splinex = AbstractArray[]
    spliney = AbstractArray[]
    pdf = Normal(0,nFactor)

    for file in csv_files
        obsData = readdlm(file, ',', Float64, '\n', header=false)  # read the observation data
        w = rand(pdf, size(obsData))
        obsData = obsData + w
        push!(ObsDataList, obsData') # store the transpose of the observation data to fit the comparison function
        push!(splinex, obsData[:,1])
        push!(spliney, obsData[:,2])
    end

    return ObsDataList, splinex, spliney, pdf
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