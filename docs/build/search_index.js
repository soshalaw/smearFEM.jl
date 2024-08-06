var documenterSearchIndex = {"docs":
[{"location":"api/#API-Reference","page":"API Reference","title":"API Reference","text":"","category":"section"},{"location":"api/#Core-functions","page":"API Reference","title":"Core functions","text":"","category":"section"},{"location":"api/","page":"API Reference","title":"API Reference","text":"Modules = [smearFEM]\nPages = [\"fem.jl\"]","category":"page"},{"location":"api/#smearFEM.assemble_system","page":"API Reference","title":"smearFEM.assemble_system","text":"assemble_system(ne, NodeList, IEN, ndim, FunctionClass=\"Q1\", nDof=1, ID=nothing, Young=1, ν=0.3)\n\nAssembles the finite element system. # Returns the global stiffness matrix\n\nArguments:\n\nne::Interger: number of elements in each direction\nNodeList::Matrix{Float64}{ndim,nNodes} : coordinates of the nodes\nIEN::Matrix{Int}{nElements,nLocalNodes} : connectivity matrix\nndim::Interger: number of dimensions\nnDof::Interger: number of degree of freedom per node\nFunctionClass::String: type of basis functions to be considered (Q1:quadratic or Q2:Lagrange)\nID::Matrix{Int}{nNodes,nDof} : matrix that maps the global degrees of freedom to the local degrees of freedom\nYoung::Float64: Young's modulus\nν::Float64: Poisson's ratio\n\nReturns:\n\nK::SparseMatrixCSC{Float64,Int64}{ndof,ndof} : sparse stiffness matrix \n\n\n\n\n\n","category":"function"},{"location":"api/#smearFEM.basis_function","page":"API Reference","title":"smearFEM.basis_function","text":"basis_function(ξ,η=nothing,ζ=nothing,FunctionClass = \"Q1\")\n\nDefine the basis functions and the gradients for a master element\n\nArguments:\n\nξ::Float64: ξ coordinate of the point where the basis function is evaluated\nη::Float64: η coordinate of the point where the basis function is evaluated\nζ::Float64: ζ coordinate of the point where the basis function is evaluated\nFunctionClass::String: type of basis functions to be considered (Q1:quadratic or Q2:Lagrange)\n\nReturns:\n\nN::Vector{Float64}{,ndof}: basis functions\nDelta_N::Matrix{Float64}{ndof,ndim}: gradient of the basis functions \n\n\n\n\n\n","category":"function"},{"location":"api/#smearFEM.gaussian_quadrature","page":"API Reference","title":"smearFEM.gaussian_quadrature","text":"smearFEM.gaussian_quadrature(a,b,nGaussPoints)\n\nCompute the nodes and weights for the Gaussian quadrature of order 2\n\nArguments:\n\na,b::Integer : the limits of the integration interval\nnGaussPoints::Integer : number of Gauss points to be considered (2 or 3)\n\nReturns:\n\nξ::Vector{Float64}{,nGaussPoints}: nodes.\nw::Vector{Float64}{,nGaussPoints}: weights \n\n\n\n\n\n","category":"function"},{"location":"api/#Post-processing-functions","page":"API Reference","title":"Post processing functions","text":"","category":"section"},{"location":"api/","page":"API Reference","title":"API Reference","text":"Modules = [smearFEM]\nPages = [\"PostProcess.jl\"]","category":"page"},{"location":"api/#smearFEM.PlotMesh-Tuple{Any, Any}","page":"API Reference","title":"smearFEM.PlotMesh","text":"PlotMesh(NodeList, IEN)\n\nFunction to plot the mesh\n\nArguments:\n\nNodeList::Matrix{Float64}{nNodes,ndim}: array of nodes.\nIEN::Matrix{nElem,nNodes}: IEN array.\n\n\n\n\n\n","category":"method"},{"location":"api/#smearFEM.back_project-Tuple{Any, Any}","page":"API Reference","title":"smearFEM.back_project","text":"back_project(NodeList, CameraMatrix)\n\nProject the 3D mesh to 2D image plane\n\nArguments:\n\nNodeList::Matrix{Float64}{3,nbNodes}: 3D mesh grid\nCameraMatrix::Matrix{Float64}{3,3}: Camera matrix\n\nReturns:\n\nNodeList2D::Matrix{Float64}{2,nbNodes}: 2D coordinates of the nodes\n\n\n\n\n\n","category":"method"},{"location":"api/#smearFEM.extract_borders","page":"API Reference","title":"smearFEM.extract_borders","text":"Extract_borders(NodeList, CameraMatrix, BorderNodesList, state, ne = nothing)\n\nProject the 3D mesh to 2D image plane and extract the border nodes (left and right)\n\nArguments:\n\nNodeList::Matrix{Float64}{ndim,nNodes} : coordinates of the nodes\nCameraMatrix::Matrix{Float64}{3,3} : Camera matrix\nBorderNodesList::Vector{Vector{Any}{4,N}:  : List of border nodes\nstate::String : State of the function (init:During the initialization of the mesh or update: when the mesh is updated)\nne::Integer: Number of elements in each direction\n\nReturns:\n\nNodeList::Matrix{Float64}{ndim,nbNodes}: 2D coordinates of the border nodes\n\n\n\n\n\n","category":"function"},{"location":"api/#smearFEM.fit_curve-Tuple{Any}","page":"API Reference","title":"smearFEM.fit_curve","text":"fit_curve(border)\n\nFit a curve to the border nodes of the 2D mesh\n\n\n\n\n\n","category":"method"},{"location":"api/#smearFEM.inflate_sphere-NTuple{5, Any}","page":"API Reference","title":"smearFEM.inflate_sphere","text":"inflate_sphere(NodeList, x0, x1, y0, y1)\n\nInflate the sphere to a unit sphere \n\nArguments:\n\nNodeList::Matrix : {[ndims,nNodes], Matrix{Float64}} : The coordinates of the nodes\nx0 : Float64 : The lower bound of the x direction\nx1 : Float64 : The upper bound of the x direction\ny0 : Float64 : The lower bound of the y direction\ny1 : Float64 : The upper bound of the y direction\n\nReturns:\n\nNodeList::Matrix{Float64}{ndims,nNodes}:  The coordinates of the nodes after inflation\n\n\n\n\n\n","category":"method"},{"location":"api/#smearFEM.noramlize-Tuple{Any, Any}","page":"API Reference","title":"smearFEM.noramlize","text":"noramlize(q, IEN)\n\nFunction normalize the solution vector for plotting\n\nArguments:\n\nq: solution vector\nIEN::Matrix{Float64}{nElem, nNodes}: IEN array\n\nReturns:\n\nqList: normalized list of solutions \n\n\n\n\n\n","category":"method"},{"location":"api/#smearFEM.truncate_colormap","page":"API Reference","title":"smearFEM.truncate_colormap","text":"truncate_colormap(minval=0.0, maxval=1.0, n=100)\n\nFunction to truncate a colormap\n\nArguments:\n\nminval::Integer: minimum value of the colormap\nmaxval::Integer: maximum value of the colormap\nn::Integer: number of colors\n\nReturns:\n\nnew_cmap: truncated colormap\n\n\n\n\n\n","category":"function"},{"location":"api/#smearFEM.write_scene-NTuple{6, Any}","page":"API Reference","title":"smearFEM.write_scene","text":"write_scene(fileName, NodeList, IEN, ne, ndim, fields)\n\nFunction to write the solution to a VTK file\n\nArguments:\n\nfileName::String: name of the VTK file.\nNodeList::Matrix{Float64}{nNodes, ndim}: array of nodes.\nIEN::Matrix{nElem, nNodes}: IEN array {[nElem, nNodes]}\nne::Integer: number of elements in each direction.\nndim::Integer: number of dimensions\nfields::Vector{Vector{Float64}}: solution vector.\n\n\n\n\n\n","category":"method"},{"location":"api/#smearFEM.write_vtk-NTuple{7, Any}","page":"API Reference","title":"smearFEM.write_vtk","text":"write_vtk(fileName, fieldName, NodeList, IEN, ne, ndim, q)\n\nFunction to write the solution to a VTK file\n\nArguments:\n\nfileName::String: name of the VTK file.\nNodeList::Matrix[nNodes, ndim]: array of nodes.\nIEN::Matrix{Float64}{nElem, nNodes}: IEN array.\nne::Integer: number of elements in each direction.\nndim::Integer: number of dimensions.\nq::Vector{Float64}: solution vector.\n\n\n\n\n\n","category":"method"},{"location":"examples/#Examples","page":"Examples","title":"Examples","text":"","category":"section"},{"location":"#smearFEM","page":"Home","title":"smearFEM","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for smearFEM.","category":"page"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"using Pkg\npkg\"add smearFEM\"","category":"page"},{"location":"#Testing","page":"Home","title":"Testing","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"using Pkg\npkg\"test smearFEM\"","category":"page"}]
}
