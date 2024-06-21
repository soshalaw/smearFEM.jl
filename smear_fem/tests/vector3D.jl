using LinearAlgebra
using ProgressMeter
using WriteVTK

include("../src/fem.jl")
include("../src/PostProcess.jl")

# set up mesh grid
function meshgrid(x0,x1,y0,y1,z0,z1,ne,ndim)

    NodeList = zeros(ndim,(ne+1)^ndim)
    IEN = zeros(Int64,ne^ndim,2^ndim) # IEN for the 3D mesh
    IEN_1 = zeros(Int64,ne^(ndim-1),2^(ndim-1)) # IEN for the top surface
    IEN_2 = zeros(Int64,ne^(ndim-1),2^(ndim-1)) # IEN for the bottom surface
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
        o = 1 # element number on 2D mesh of bottom surface
        p = 1 # element number on 2D mesh of top surface
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
                        IEN_1[p,1] = IEN[n,1]
                        IEN_1[p,2] = IEN[n,2]
                        IEN_1[p,3] = IEN[n,3]
                        IEN_1[p,4] = IEN[n,4]
                        p = p + 1
                    elseif k == ne # populate the IEN for the top surface
                        IEN_2[o,1] = IEN[n,5]
                        IEN_2[o,2] = IEN[n,6]
                        IEN_2[o,3] = IEN[n,7]
                        IEN_2[o,4] = IEN[n,8]
                        o = o + 1
                    end
                    n = n + 1
                end
            end
        end
    end
    return NodeList, IEN, ID, IEN_1, IEN_2
end


function setboundaryCond(NodeList, ne, ndim, FunctionClass, d)
    # set dirichlet boundary conditions

    if FunctionClass == "Q1"
        q_d = zeros(ndim*(ne+1)^ndim,1)
        q_n = zeros(ndim*(ne+1)^ndim,1)

        C = Matrix{Int}(I,ndim*(ne+1)^ndim,ndim*(ne+1)^ndim) # definition of the constraint matrix
    end

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
                # elseif nNode==ne+1 
                #     q_d[3*(nNode-1)+2] = 0.
                #     push!(rCol, 3*(nNode-1)+2)
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
        
    return q_d, q_n, C
end

function main()

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

    NodeList, IEN, ID, IEN_1,  IEN_2 = meshgrid(x0,x1,y0,y1,z0,z1,ne,ndim) # generate the mesh grid

    NodeListCylinder = PostProcess.inflate_sphere(NodeList, x0, x1, y0, y1) # inflate the sphere to a unit sphere

    K = fem.assemble_system(ne, NodeListCylinder, IEN, ndim, nDof, FunctionClass, ID, ν, Young) # assemble the stiffness matrix

    delta = 0.1:0.01:0.3 # set the displacement increment
    fields = [] # store the solution fields

    for d in delta
        q_d, q_n, C = setboundaryCond(NodeListCylinder, ne, ndim, FunctionClass, d);

        # transpose the constraint matrix
        C_t = transpose(C)

        # extract the free part of the stiffness matrix
        K_free = C_t*K*C

        b = q_n - K*q_d

        invK = inv(K_free)

        # solve the system
        q_f = invK*C_t*b

        # assemble the solution 
        q = q_d + C*q_f;

        # post process the solution
        u = q[ID[:,1]]
        v = q[ID[:,2]] 
        w = q[ID[:,3]]

        # println("ν = ", -u[end]/w[end])
        push!(fields, [u v w]')

        println("d = ", d)
    end

    cellType = VTKCellTypes.VTK_HEXAHEDRON

    cells = [MeshCell(cellType,IEN[e,:]) for e in 1:ne^ndim]

    paraview_collection("displacement") do pvd # create a paraview collection
        @showprogress "Writing out to VTK..." for i in 1:length(fields)
            vtk_grid("timestep_$i", NodeList, cells) do vtk # write out the fields to VTK
                vtk["u"] = fields[i]
                time = (i - 1)*0.5
                pvd[time] = vtk
            end
        end
    end

end
main()


