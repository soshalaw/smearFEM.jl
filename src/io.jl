using WriteVTK
using ProgressMeter
using DelimitedFiles

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
function write_vtk(filePath, fieldName, NodeList, IEN, ne, ndim, q)

    if ndim == 1
        cellType = VTKCellTypes.VTK_LINE
    elseif ndim == 2
        cellType = VTKCellTypes.VTK_QUAD
    elseif ndim == 3
        cellType = VTKCellTypes.VTK_HEXAHEDRON
    end
    
    println(size(q))
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
function write_scene(fileName, NodeList, IEN, ne, ndim, fields)

    cellType = VTKCellTypes.VTK_HEXAHEDRON

    cells = [MeshCell(cellType,IEN[e,:]) for e in 1:ne^ndim]

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
    readCSV(csv_path)

Function to read the CSV files in the directory and fit a curve to the border observations in the CSV files.

# Arguments:
- `csv_path::String`: path to the directory containing the CSV files.

# Returns:
- `splinep`: x coordinates samples of the spline parameters of the border nodes.
- `splineq`: y coordinates samples of the spline parameters of the border nodes.

"""

function readCSV(csv_path)

    csv_files = readdir(csv_path, join=true)        # get the list of the csv files in the directory
    ObsDataList = []                                  # store the observation data
    splinex = []
    spliney = []

    for file in csv_files
        obsData = readdlm(file, ',', Float64, '\n', header=false)  # read the observation data
        push!(ObsDataList, obsData') # store th*e transpose of the observation data to fit the comparison function
        push!(splinex, obsData[:,1])
        push!(spliney, obsData[:,2])
    end

    return ObsDataList, splinex, spliney
end

"""
    writeCSV(fileName, borders)

Function to write the border data to a CSV file.

# Arguments:
- `fileName::String`: name of the file.
- `borders::Vector{Matrix{Float64}}`: vector of border data.
"""
function writeCSV(fileName, borders)
    println("Writing CSV files...")
    counter = 1
    for border in borders
        cStr = string(counter,pad=3)
        open(string(fileName,"/contour_data/",cStr,".csv"), "w") do io
            writedlm(io, [border[1,:] border[2,:]], ',')
        end
        counter += 1
    end
end