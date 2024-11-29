using LinearAlgebra
using ProgressMeter
using SparseArrays
using LinearAlgebra
using ProgressMeter
using SparseArrays

"""
Set the Dirichlet boundary conditions for the problem

# Arguments:
- `NodeList::Matrix{Float64}{nNodes,ndim}` : array of nodes
- `ne::Int` : number of elements
- `ndim::Int` : number of dimensions
- `FunctionClass::String` : type of basis function

# Returns:
- `q_upper::Vector{Float64}` : vector of the Dirichlet boundary conditions (for ndof = 1) / Dirichlet boundary conditions upper surface (for ndof > 1)
- `q_lower::Vector{Float64}` : vector of the Neumann boundary conditions (for ndof = 1) / Dirichlet boundary conditions lower surface (for ndof > 1)
- `C_uc::SparseMatrixCSC{Float64,Int64}` : onstraint matrix
"""
function set_boundary_cond_flow_cube(NodeList, ne, ndim, FunctionClass, nDof=1)
    if FunctionClass == "Q1"
        q_upper = zeros(nDof*(ne+1)^ndim,1)                  # initialize the vector of the Dirichlet boundary conditions (for ndof = 1) / Dirichlet boundary conditions upper surface (for ndof > 1)
        q_lower = zeros(nDof*(ne+1)^ndim,1)                  # initialize the vector of the Neumann boundary conditions (for ndof = 1) / Dirichlet boundary conditions lower surface (for ndof > 1)
        q_side = zeros(nDof*(ne+1)^ndim,1)                   # initialize the vector of the Neumann boundary conditions (for ndof = 1) / Dirichlet boundary conditions lower surface (for ndof > 1)
        C = sparse(I,nDof*(ne+1)^ndim,nDof*(ne+1)^ndim)      # definition of the constraint matrix
    elseif FunctionClass == "Q2"
        q_upper = zeros(nDof*(2*ne+1)^ndim,1)                # initialize the vector of the Dirichlet boundary conditions (for ndof = 1) / Dirichlet boundary conditions upper surface (for ndof > 1)
        q_lower = zeros(nDof*(2*ne+1)^ndim,1)                # initialize the vector of the Neumann boundary conditions (for ndof = 1) / Dirichlet boundary conditions lower surface (for ndof > 1)
        q_side = zeros(nDof*(2*ne+1)^ndim,1)                  # initialize the vector of the Neumann boundary conditions (for ndof = 1) / Dirichlet boundary conditions lower surface (for ndof > 1)
        C = sparse(I,nDof*(2*ne+1)^ndim,nDof*(2*ne+1)^ndim)  # definition of the constraint matrix
    end
    if nDof == 1
        if ndim == 3
            Dbound1 = 0
            Dbound2 = 1
            iter = 1:size(NodeList,2)
            for n in iter
                coord = NodeList[:,n] # get the coordinates of the node
                if coord[3] == Dbound1 # bottom boundary
                    q_upper[n] = 0
                elseif coord[3] == Dbound2 # top boundary
                    q_upper[n] = -1
                end
            end
        elseif ndim == 2
            Dbound1 = 0
            Dbound2 = 1
            iter = 1:size(NodeList,2)
            for n in iter
                coord = NodeList[:,n] # get the coordinates of the node
                if coord[2] == Dbound1 # bottom boundary
                    q_upper[n] = 0
                elseif coord[2] == Dbound2 # top boundary
                    q_upper[n] = -1
                end
            end
        end

        if FunctionClass == "Q1"
            C_uc = C[:,((ne+1)^(ndim-1)+1):((ne+1)^ndim-(ne+1)^(ndim-1))]
        elseif FunctionClass == "Q2"
            C_uc = C[:,((2*ne+1)^(ndim-1)+1):((2*ne+1)^ndim-(2*ne+1)^(ndim-1))]
        end
    else
        z0Bound = -0.5
        z1Bound = 0.5
        x0Bound = -0.5
        x1Bound = 0.5
        y0Bound = -0.5
        y1Bound = 0.5

        rCol = Array{Int}(undef,0)
        iter = 1:size(NodeList,2)
        for nNode in iter
            coord = NodeList[:,nNode]     # get the coordinates of the node
            if coord[1] == x0Bound || coord[1] == x1Bound || coord[2] == y0Bound || coord[2] == y1Bound
                q_side[3*nNode-2] = 1     # constraint the x displacement
                q_side[3*nNode-1] = 1     # constraint the y displacement
                q_side[3*nNode] = 1       # constraint the z displacement
                
                push!(rCol,3*nNode-2)
                push!(rCol,3*nNode-1)
                push!(rCol,3*nNode)
            elseif coord[3] == z1Bound         # top boundary
                q_upper[3*nNode-2] = 0     # constraint the x displacement
                q_upper[3*nNode-1] = 0     # constraint the y displacement
                q_upper[3*nNode] = (1-(2*coord[1])^2)*(1-(2*coord[2])^2)     # constraint the z displacement to be -d

                push!(rCol,3*nNode-2)
                push!(rCol,3*nNode-1)
                push!(rCol,3*nNode)
            elseif coord[3] == z0Bound     # bottom boundary
                q_lower[3*nNode-2] = 0     # constraint the x displacement
                q_lower[3*nNode-1] = 0     # constraint the y displacement
            #     q_lower[3*nNode] = (1-(2*coord[1])^2)*(1-(2*coord[2])^2)     # constraint the z displacement to be -d

                push!(rCol,3*nNode-2)
                push!(rCol,3*nNode-1)
                # push!(rCol,3*nNode)
            end
        end
        C_uc = C[:,setdiff(1:size(C,2),rCol)]
    end
    return q_upper, q_side, q_lower, C_uc
end

function set_boundary_cond_flow_cyl(NodeList, ne, ndim, FunctionClass, nDof=1)
    if FunctionClass == "Q1"
        q_upper = zeros(nDof*(ne+1)^ndim,1)                  # initialize the vector of the Dirichlet boundary conditions (for ndof = 1) / Dirichlet boundary conditions upper surface (for ndof > 1)
        q_lower = zeros(nDof*(ne+1)^ndim,1)                  # initialize the vector of the Neumann boundary conditions (for ndof = 1) / Dirichlet boundary conditions lower surface (for ndof > 1)
        q_side = zeros(nDof*(ne+1)^ndim,1)                   # initialize the vector of the Neumann boundary conditions (for ndof = 1) / Dirichlet boundary conditions lower surface (for ndof > 1)
        C = sparse(I,nDof*(ne+1)^ndim,nDof*(ne+1)^ndim)      # definition of the constraint matrix
    elseif FunctionClass == "Q2"
        q_upper = zeros(nDof*(2*ne+1)^ndim,1)                # initialize the vector of the Dirichlet boundary conditions (for ndof = 1) / Dirichlet boundary conditions upper surface (for ndof > 1)
        q_lower = zeros(nDof*(2*ne+1)^ndim,1)                # initialize the vector of the Neumann boundary conditions (for ndof = 1) / Dirichlet boundary conditions lower surface (for ndof > 1)
        q_side = zeros(nDof*(2*ne+1)^ndim,1)                  # initialize the vector of the Neumann boundary conditions (for ndof = 1) / Dirichlet boundary conditions lower surface (for ndof > 1)
        C = sparse(I,nDof*(2*ne+1)^ndim,nDof*(2*ne+1)^ndim)  # definition of the constraint matrix
    end
    if nDof == 1
        if ndim == 3
            Dbound1 = 0
            Dbound2 = 1
            iter = 1:size(NodeList,2)
            for n in iter
                coord = NodeList[:,n] # get the coordinates of the node
                if coord[3] == Dbound1 # bottom boundary
                    q_upper[n] = 0
                elseif coord[3] == Dbound2 # top boundary
                    q_upper[n] = -1
                end
            end
        elseif ndim == 2
            Dbound1 = 0
            Dbound2 = 1
            iter = 1:size(NodeList,2)
            for n in iter
                coord = NodeList[:,n] # get the coordinates of the node
                if coord[2] == Dbound1 # bottom boundary
                    q_upper[n] = 0
                elseif coord[2] == Dbound2 # top boundary
                    q_upper[n] = -1
                end
            end
        end

        if FunctionClass == "Q1"
            C_uc = C[:,((ne+1)^(ndim-1)+1):((ne+1)^ndim-(ne+1)^(ndim-1))]
        elseif FunctionClass == "Q2"
            C_uc = C[:,((2*ne+1)^(ndim-1)+1):((2*ne+1)^ndim-(2*ne+1)^(ndim-1))]
        end
    else
        z0Bound = -0.5
        z1Bound = 0.5
        x0Bound = -0.5
        x1Bound = 0.5
        y0Bound = -0.5
        y1Bound = 0.5

        rCol = Array{Int}(undef,0)
        iter = 1:size(NodeList,2)
        for nNode in iter
            coord = NodeList[:,nNode]     # get the coordinates of the node
            if coord[1] == x0Bound || coord[1] == x1Bound || coord[2] == y0Bound || coord[2] == y1Bound
                q_side[3*nNode-2] = 1     # constraint the x displacement
                q_side[3*nNode-1] = 1     # constraint the y displacement
                q_side[3*nNode] = 1       # constraint the z displacement
                
                push!(rCol,3*nNode-2)
                push!(rCol,3*nNode-1)
                push!(rCol,3*nNode)
            elseif coord[3] == z1Bound         # top boundary
                q_upper[3*nNode-2] = 0     # constraint the x displacement
                q_upper[3*nNode-1] = 0     # constraint the y displacement
                q_upper[3*nNode] = (1 - (coord[1]^2+coord[2]^2))     # constraint the z displacement to be -d

                push!(rCol,3*nNode-2)
                push!(rCol,3*nNode-1)
                push!(rCol,3*nNode)
            elseif coord[3] == z0Bound     # bottom boundary
                q_lower[3*nNode-2] = 0     # constraint the x displacement
                q_lower[3*nNode-1] = 0     # constraint the y displacement
            #     q_lower[3*nNode] = (1-(2*coord[1])^2)*(1-(2*coord[2])^2)     # constraint the z displacement to be -d

                push!(rCol,3*nNode-2)
                push!(rCol,3*nNode-1)
                # push!(rCol,3*nNode)
            end
        end
        C_uc = C[:,setdiff(1:size(C,2),rCol)]
    end
    return q_upper, q_side, q_lower, C_uc
end
