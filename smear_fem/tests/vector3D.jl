using LinearAlgebra
using ProgressMeter
using WriteVTK
using SparseArrays

include("../src/fem.jl")
include("../src/PostProcess.jl")

# set up mesh grid
function meshgrid(x0,x1,y0,y1,z0,z1,ne,ndim)

    NodeList = zeros(ndim,(ne+1)^ndim)
    IEN = zeros(Int64,ne^ndim,2^ndim) # IEN for the 3D mesh
    IEN_top = zeros(Int64,ne^(ndim-1),2^(ndim-1)) # IEN for the top surface
    IEN_btm = zeros(Int64,ne^(ndim-1),2^(ndim-1)) # IEN for the bottom surface
    ID = zeros(Int64,(ne+1)^ndim,ndim)

    if ndim == 2
        x = collect(range(x0, x1, length=ne+1))
        y = collect(range(y0, y1, length=ne+1))
    
        m = 1
        for j in 1:ne+1 # y direction
            for i in 1:ne+1 # x direction
                NodeList[1,m] = x[i]
                NodeList[2,m] = y[j]
                for l in 1:ndim
                    ID[m,l] = ndim*(m-1) + l
                end
                m = m + 1
            end
        end 
        
        n = 1
        for j in 1:ne # y direction
            for i in 1:ne # x direction
                IEN[n,1] = (j-1)*(ne+1) + i
                IEN[n,2] = (j-1)*(ne+1) + i + 1
                IEN[n,3] = j*(ne+1) + i + 1
                IEN[n,4] = j*(ne+1) + i
                if j == 1 # populate the IEN for the bottom surface
                    IEN_1[i,1] = IEN[n,1]
                    IEN_1[i,2] = IEN[n,2]
                elseif j == ne # populate the IEN for the top surface
                    IEN_2[i,1] = IEN[n,4]
                    IEN_2[i,2] = IEN[n,3]
                end
                n = n + 1
            end
        end

    elseif ndim == 3

        x = collect(range(x0, x1, length=ne+1))
        y = collect(range(y0, y1, length=ne+1))
        z = collect(range(z0, z1, length=ne+1))
        
        m = 1
        for k in 1:ne+1 # z direction
            for j in 1:ne+1 # y direction
                for i in 1:ne+1 # x direction
                    NodeList[1,m] = x[i]
                    NodeList[2,m] = y[j]
                    NodeList[3,m] = z[k]
                    for l in 1:ndim
                        ID[m,l] = ndim*(m-1) + l
                    end
                    m = m + 1
                end
            end
        end
        
        n = 1 # element number on 3D mesh
        nt = 1 # element number on 2D mesh of top surface
        nb = 1 # element number on 2D mesh of bottom surface
        for k in 1:ne # z direction
            for j in 1:ne # y direction
                for i in 1:ne # x direction
                    IEN[n,1] = (k-1)*(ne+1)^2 + (j-1)*(ne+1) + i
                    IEN[n,2] = (k-1)*(ne+1)^2 + (j-1)*(ne+1) + i + 1
                    IEN[n,3] = (k-1)*(ne+1)^2 + j*(ne+1) + i + 1
                    IEN[n,4] = (k-1)*(ne+1)^2 + j*(ne+1) + i
                    IEN[n,5] = k*(ne+1)^2 + (j-1)*(ne+1) + i
                    IEN[n,6] = k*(ne+1)^2 + (j-1)*(ne+1) + i + 1
                    IEN[n,7] = k*(ne+1)^2 + j*(ne+1) + i + 1
                    IEN[n,8] = k*(ne+1)^2 + j*(ne+1) + i
                    if k == 1 # populate the IEN for the bottom surface
                        IEN_btm[nb,1] = IEN[n,1]
                        IEN_btm[nb,2] = IEN[n,2]
                        IEN_btm[nb,3] = IEN[n,3]
                        IEN_btm[nb,4] = IEN[n,4]
                        nb = nb + 1
                    elseif k == ne # populate the IEN for the top surface
                        IEN_top[nt,1] = IEN[n,5]
                        IEN_top[nt,2] = IEN[n,6]
                        IEN_top[nt,3] = IEN[n,7]
                        IEN_top[nt,4] = IEN[n,8]
                        nt = nt + 1
                    end
                    n = n + 1
                end
            end
        end
    end
    return NodeList, IEN, ID, IEN_top, IEN_btm
end


function setboundaryCond(NodeList, ne, ndim, FunctionClass, d, nDof=1)
    """
        Set the Dirichlet boundary conditions for the problem

        Parameters:
        NodeList: {[ndim,nNodes] Matrix{Float64}} : matrix of the coordinates of the nodes
        ne: {Int} : number of elements
        ndim: {Int} : number of dimensions
        FunctionClass: {String} : type of basis function
        d: {Float64} : displacement increment
        nDOf: {Int} : number of degree of freedom per node

        Returns:
        q_d: {[ndof] Vector{Float64}} : Dirichlet boundary conditions
    """

    # if FunctionClass == "Q1"
    #     q_d = zeros(nDof*(ne+1)^ndim,1)                       # initialize the vector of the Dirichlet boundary conditions
    #     C = Matrix{Int}(I,ndim*(ne+1)^ndim,ndim*(ne+1)^ndim)  # definition of the constraint matrix
    # end

    # z0Bound = 0
    # z1Bound = 1

    # rCol = Array{Int}(undef,0)

    # for nNode in 1:size(NodeList,2)
    #     coord = NodeList[:,nNode]    # get the coordinates of the node
    #     if coord[3] == z0Bound       # bottom boundary
    #         q_d[3*nNode] = 0         # constraint the z displacement to be zero at the bottom boundary
    #         push!(rCol,3*nNode)
    #     elseif coord[3] == z1Bound   # top boundary
    #         q_d[3*nNode] = -d        # constraint the z displacement to be -d
    #         push!(rCol,3*nNode)
    #     end
    # end

    # C = C[:,setdiff(1:size(C,2),rCol)]
        
    # return q_d, C

    q_d = zeros(ndim*(ne+1)^ndim,1)
    q_n = zeros(ndim*(ne+1)^ndim,1)

    C = Matrix{Int}(I,ndim*(ne+1)^ndim,ndim*(ne+1)^ndim) # definition of the constraint matrix

    z0Bound = 0
    z1Bound = 1

    yBBound = 0
    yTBound = 1

    xLBound = 0
    xRBound = 1

    rCol = Array{Int}(undef,0)

    for nNode in 1:size(NodeList,2)
        coord = NodeList[:,nNode] # get the coordinates of the node
        if ndim == 2
            if coord[2] == yBBound # bottom boundary
                if coord[1] == xLBound
                    q_d[2*nNode-1] = 0 # constraint the x displacement and y displacement to be zero at the bottom left corner
                    q_d[2*nNode] = 0

                    push!(rCol,2*nNode-1)
                    push!(rCol,2*nNode)
                end

                q_d[2*nNode] = 0 # constraint the y displacement to be zero at the bottom boundary
                push!(rCol,2*nNode)

            elseif coord[2] == yTBound # top boundary
                q_d[2*nNode] = -d # constraint the y displacement to be -delta
                push!(rCol,2*nNode)
            end
        elseif ndim == 3
            if coord[3] == z0Bound # bottom boundary
                # if nNode == 1
                if nNode == (ne+1)^2÷2 +1
                    q_d[3*(nNode-1)+1] = 0.
                    q_d[3*(nNode-1)+2] = 0.
                    push!(rCol, 3*(nNode-1)+1)
                    push!(rCol, 3*(nNode-1)+2)
                elseif nNode==ne+1 
                    q_d[3*(nNode-1)+2] = 0.
                    push!(rCol, 3*(nNode-1)+2)
                end

                q_d[3*nNode] = 0 # constraint the z displacement to be zero at the bottom boundary
                push!(rCol,3*nNode)
            end
            if coord[3] == z1Bound # top boundary
                q_d[3*nNode] = -d # constraint the z displacement to be -delta
                push!(rCol,3*nNode)
            end
        end

    end

    C = C[:,setdiff(1:size(C,2),rCol)]
        
    return q_d, C
end

function apply_boundary_conditions(ne, NodeList, IEN, IEN_top, IEN_btm, ndim, FunctionClass, ID=None, nDof=3) # , nDof=1, FunctionClass="Q1", ID=None, Young=1, ν=0.3, q_d=zeros(0), q_n=zeros(0))
    """ Apply the Neumann slip boundary conditions to the global stiffness matrix

        Parameters:
        K: {[ndof,ndof] SparseMatrixCSC{Float64,Int64}} : sparse stiffness matrix 
        ID: {[nNodes,nDof] Matrix{Int}} : matrix that maps the global degrees of freedom to the local degrees of freedom
        q_d: {[ndof] Vector{Float64}} : Dirichlet boundary conditions
        q_n: {[ndof] Vector{Float64}} : Neumann boundary conditions

        Returns:
        K: {[ndof,ndof] SparseMatrixCSC{Float64,Int64}} : sparse stiffness matrix with the boundary conditions applied
        F: {[ndof] Vector{Float64}} : force vector
    """

    E = zeros(  Int64, ne^(ndim-1)*(size(ID,2)*size(IEN_btm,2))^2*2) # *2 because we have two surfaces
    J = zeros(  Int64, ne^(ndim-1)*(size(ID,2)*size(IEN_btm,2))^2*2) # *2 because we have two surfaces
    V = zeros(Float64, ne^(ndim-1)*(size(ID,2)*size(IEN_btm,2))^2*2) # *2 because we have two surfaces

    for e in 1:ne^(ndim-1)
    
        coords_top = NodeList[:,IEN_top[e,:]] # get the coordinates of the nodes of the element
        coords_btm = NodeList[:,IEN_btm[e,:]] # get the coordinates of the nodes of the element

        if ndim == 2
            # gaussian quadrature points for the element [-1,1] 
            ξ, w_ξ = fem.gaussian_quadrature(-1,1)

            wpoints = [w_ξ[1], w_ξ[2]]
            
            x = [ξ[1], ξ[2]]
        elseif ndim == 3
            # gaussian quadrature points for the element [-1,1]x[-1,1] 
            ξ, w_ξ = fem.gaussian_quadrature(-1,1)
            η, w_η = fem.gaussian_quadrature(-1,1)
            
            wpoints = [w_ξ[1]*w_η[1], w_ξ[2]*w_η[1], w_ξ[2]*w_η[2], w_ξ[1]*w_η[2]]
            
            x = [ξ[1], ξ[2], ξ[2], ξ[1]]
            y = [η[1], η[1], η[2], η[2]]
        end 

        for gp in 1:2^(ndim-1)

            if ndim == 2
                N, ΔN = fem.basis_function(x[gp], nothing, nothing, FunctionClass)
            elseif ndim == 3
                N, ΔN = fem.basis_function(x[gp], y[gp], nothing, FunctionClass) 
            end

            dxdξ_top = coords_top*ΔN         # Jacobian matrix [dx/dxi dx/deta; dy/dxi dy/deta; dz/dxi dz/deta]
            dxdξ_btm = coords_btm*ΔN         # Jacobian matrix [dx/dxi dx/deta; dy/dxi dy/deta; dz/dxi dz/deta]

            w_top = wpoints[gp]*norm(cross(dxdξ_top[:,1],dxdξ_top[:,2]))     # weight of the quadrature point top surface
            w_btm = wpoints[gp]*norm(cross(dxdξ_btm[:,1],dxdξ_btm[:,2]))     # weight of the quadrature point bottom surface
            
            M = zeros(3, ndim*length(N))
            M[1,1:nDof:end] = N
            M[2,2:nDof:end] = N
            M[3,3:nDof:end] = N

            be = M'*M
            display(be)
            # loop between basis functions of the element
            for iNode in 1:size(be,1)÷nDof
                for jNode in 1:size(be,2)÷nDof
                    for iDof in 1:size(ID,2)
                        for jDof in 1:size(ID,2)
                            i = (iNode-1)*nDof + iDof
                            j = (jNode-1)*nDof + jDof
                            
                            inz_btm = length(be)*(e-1) + (iNode-1)*nDof*size(be,2) + (jNode-1)*nDof^2 + (iDof-1)*nDof + jDof # index for the COO sparse matrix
                            inz_top = length(be)*(ne)^(ndim-1)+ length(be)*((e-1)) + (iNode-1)*nDof*size(be,2) + (jNode-1)*nDof^2 + (iDof-1)*nDof + jDof # index for the COO sparse matrix
                            
                            E[inz_top] = ID[IEN_top[e,iNode],iDof]    # row index 
                            J[inz_top] = ID[IEN_top[e,jNode], jDof]   # column index
                            V[inz_top] += w_top*be[i,j] 

                            E[inz_btm] = ID[IEN_btm[e,iNode],iDof]    # row index 
                            J[inz_btm] = ID[IEN_btm[e,jNode], jDof]   # column index
                            V[inz_btm] += w_btm*be[i,j] 
                        end
                    end
                end
            end
            # TODO include in tha assembly function
        end
    end
    
    K = sparse(E,J,V)
    display(Matrix(K))

    return  K
end

# function main()

    # test case 
    x0 = 0
    x1 = 1
    y0 = 0
    y1 = 1
    z0 = 0
    z1 = 1
    ne = 2
    Young = 40
    ν = 0.4
    ndim = 3
    FunctionClass = "Q1"
    nDof = ndim  # number of degree of freedom per node
    β = 0.1

    NodeList, IEN, ID, IEN_top,  IEN_btm = meshgrid(x0,x1,y0,y1,z0,z1,ne,ndim) # generate the mesh grid
 
    NodeListCylinder = PostProcess.inflate_sphere(NodeList, x0, x1, y0, y1) # inflate the sphere to a unit sphere

    K = fem.assemble_system(ne, NodeListCylinder, IEN, ndim, nDof, FunctionClass, ID, ν, Young) # assemble the stiffness matrix

    b = apply_boundary_conditions(ne, NodeListCylinder, IEN, IEN_top, IEN_btm, ndim, FunctionClass, ID) # apply the neumann boundary conditions

    K_bar = K + β*b

    # delta = 0.1:0.01:0.3 # set the z displacement increment
    fields = [] # store the solution fields

    d = 0.1
    # for d in delta
    q_d, C = setboundaryCond(NodeList, ne, ndim, FunctionClass, d, nDof)

    # transpose the constraint matrix
    C_t = transpose(C)

    # extract the free part of the stiffness matrix
    K_free = C_t*K_bar*C

    invK = inv(K_free)

    # solve the system
    q_f = invK*C_t*(-K_bar*q_d)

    # assemble the solution 
    q = q_d + C*q_f;

    # post process the solution
    u = q[ID[:,1]]
    v = q[ID[:,2]] 
    w = q[ID[:,3]]

    # println("ν = ", -v[end]/w[end])
    push!(fields, [u v w]')

    println("d = ", d)
    # end

    cellType = VTKCellTypes.VTK_HEXAHEDRON

    cells = [MeshCell(cellType,IEN[e,:]) for e in 1:ne^ndim]

    paraview_collection("vtkFiles/displacement") do pvd # create a paraview collection
        @showprogress "Writing out to VTK..." for i in 1:length(fields)
            vtk_grid("timestep_$i", NodeList, cells) do vtk # write out the fields to VTK
                vtk["u"] = fields[i]
                time = (i - 1)*0.5
                pvd[time] = vtk
            end
        end
    end
    # return IEN, NodeList, ID, IEN_top, IEN_btm
# end
# IEN, NodeList, ID, IEN_top, IEN_btm = main()


