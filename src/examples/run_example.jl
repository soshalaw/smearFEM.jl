using LinearSolve
using SparseArrays
using IterativeSolvers
using Dates

"""
    simulate_single_tstep(r, h, ne, c1, c2, ndim, FunctionClass, nDof, β, μ_tp, μ_btm; mode="lame", GRAD=false, DENSE=false, CG=false)

Simulates the deformation of the mesh for a single time step using the linear elasticity model.

# Arguments
- `r::Number`: Radius of the cylinder.
- `h::Number`: Height of the cylinder.
- `ne::Int64`: Number of elements.
- `c1::Number`: First material parameter (e.g., Young's modulus or Lame's first parameter).
- `c2::Number`: Second material parameter (e.g., Poisson's ratio or Lame's second parameter).
- `ndim::Int64`: Number of dimensions.
- `FunctionClass::String`: Type of basis function.
- `nDof::Int64`: Number of degrees of freedom per node.
- `β::Number`: Friction parameter.
- `μ_tp::Number`: Top boundary condition.
- `μ_btm::Number`: Bottom boundary condition.

# Keyword Arguments
- `mode::String`: Type of constitutive matrix ("lame" by default).
- `GRAD::Bool`: Whether to compute gradients (default: `false`).
- `DENSE::Bool`: Whether to use dense matrices (default: `false`).
- `CG::Bool`: Whether to use the conjugate gradient method (default: `false`).

# Returns
- `q_out::Matrix{Float64}`: Displacement fields.
- `dqdθ_out::Matrix{Float64}` (if `GRAD=true`): Gradients wrt to model parameters at the nodal points.
- `mdl::model`: Model object containing mesh and material properties.
"""
function simulate_single_tstep(r::Number, h::Number, ne::Int64, c1::Number, c2::Number, ndim::Int64, FunctionClass::String, nDof::Int64, β::Number, μ_tp::Number, 
                            μ_btm::Number; mode::String="lame", GRAD::Bool=false, DENSE::Bool=false, CG::Bool=false)

    cMat = get_cMat(mode, c1, c2)
    dcdλ = get_cMat(mode, 1.0 , 0.0)

    NodeList, IEN, ID, IEN_top, IEN_btm, BorderNodesList = meshgrid_cylinder(r, h, ne, FunctionClass=FunctionClass)  # generate the mesh grid
    
    mdl = def_model("linear_elasticity", ne=ne, NodeList=NodeList, IEN=IEN, IEN_top=IEN_top, IEN_btm=IEN_btm, ndim=ndim, nDof=nDof, ID = ID,
                        FunctionClass=FunctionClass, θ1=c1, θ2=c2, cMat=cMat, dcMatdθ1=dcdλ)

    if DENSE == true
        q_tp, q_btm, C_uc = set_boundary_conditions_dense(mdl)
        K = assemble_system_dense(mdl)                   # assemble the stiffness matrix
        b = set_slip_conditions_dense(mdl)         # apply the neumann boundary conditions
    else
        q_tp, q_btm, C_uc = set_boundary_conditions(mdl)
        b = set_slip_conditions(mdl)         # apply the neumann boundary conditions
        if GRAD
            K, dKdλ = assemble_system(mdl, GRAD=true) 
            dKdβ = b             
        else
            K = assemble_system(mdl)                   # assemble the stiffness matrix
        end
    end
    
    q_d = (μ_btm*q_btm + μ_tp*q_tp)            # apply the Dirichlet boundary conditions
    K_bar = K + β*b
    C_T = transpose(C_uc)                      # transpose the constraint matrix
    K_free = C_T*K_bar*C_uc                    # extract the free part of the stiffness matrix

    if CG
        q_f = cg(K_free, C_T*(-K_bar*q_d))
    else
        q_f = K_free\(C_T*(-K_bar*q_d))         # solve the system of equations
    end

    q = q_d + C_uc*q_f                 # assemble the solution 
    q_out = hcat([q[ID[1,:]] q[ID[2,:]] q[ID[3,:]]])'    # update the nodal positions

    if GRAD && CG==false

        dqfdλ = -K_free\(C_T*dKdλ*C_uc*q_f + C_T*dKdλ*q_d)
        dqfdβ = -K_free\(C_T*dKdβ*C_uc*q_f + C_T*dKdβ*q_d)

        dqdλ = zeros(size(q_d)) + C_uc*dqfdλ
        dqdβ = zeros(size(q_d)) + C_uc*dqfdβ

        dqdλ_out = hcat(dqdλ[ID[1,:]], dqdλ[ID[2,:]], dqdλ[ID[3,:]])'
        dqdβ_out = hcat(dqdβ[ID[1,:]], dqdβ[ID[2,:]], dqdβ[ID[3,:]])'

        dqdθ_out = cat(dqdλ_out,dqdβ_out,dims=(3,3)) # concatenate the gradients in to a tensor
        return q_out, dqdθ_out, mdl
    else
        return q_out, mdl
    end
end

"""
    simulate_single_tstep_stokes(r, h, ne, η, ndim, FunctionClass_u, FunctionClass_p, nDof_u, nDof_p, β, μu_tp, μu_btm, μu_side; GRAD=false, DENSE=false)

Simulates the deformation of the mesh for a single time step using the Stokes model.

# Arguments
- `r::Number`: Radius of the cylinder.
- `h::Number`: Height of the cylinder.
- `ne::Int64`: Number of elements.
- `η::Number`: Viscosity parameter.
- `ndim::Int64`: Number of dimensions.
- `FunctionClass_u::String`: Basis function for the velocity field.
- `FunctionClass_p::String`: Basis function for the pressure field.
- `nDof_u::Int64`: Degrees of freedom per node for the velocity field.
- `nDof_p::Int64`: Degrees of freedom per node for the pressure field.
- `β::Number`: Friction parameter.
- `μu_tp::Number`: Top boundary condition.
- `μu_btm::Number`: Bottom boundary condition.
- `μu_side::Number`: Side boundary condition.

# Keyword Arguments
- `GRAD::Bool`: Whether to compute gradients (default: `false`).
- `DENSE::Bool`: Whether to use dense matrices (default: `false`).

# Returns
- `q_out::Matrix{Float64}`: Displacement fields.
- `dqdθ_out::Matrix{Float64}` (if `GRAD=true`): Gradients wrt to model parameters at the nodal points.
- `mdl::model`: Model object containing mesh and material properties.
"""
function simulate_single_tstep_stokes(r::Number, h::Number, ne::Int64, η::Number, ndim::Int64, FunctionClass_u::String, FunctionClass_p::String, nDof_u::Int64,
                                    nDof_p::Int64, β::Number, μu_tp::Number, μu_btm::Number, μu_side::Number; FunctionClass_x::String=FunctionClass_u, GRAD::Bool=false, DENSE::Bool=false)
    
    filePath = "/home/soshala/SMEAR-PhD/smear-modules/smearFEM.jl/cylindergen"

    mesh_x = meshgrid_cylinder(r, h, ne, FunctionClass=FunctionClass_x, filePath=filePath)  # generate the mesh grid for geometry
    mesh_u = meshgrid_cylinder(r, h, ne, FunctionClass=FunctionClass_u)  # generate the mesh grid
    mesh_p = meshgrid_cylinder(r, h, ne, FunctionClass=FunctionClass_p)  # generate the mesh grid

    mdl = Stokes(ndim=ndim, mesh_x=mesh_x, mesh_u=mesh_u, nDof_u=nDof_u, mesh_p=mesh_p, nDof_p=nDof_p, η=[η])

    ID_u = mdl.mesh_u.ID
    
    q_tp, q_side, q_btm, C_uc = set_boundary_cond(mdl)

    if DENSE == true
        A_bar = assemble_system_A(mdl)               # assemble the stiffness matrix
        B = assemble_system_B(mdl)                   # assemble the stiffness matrix
        b = apply_boundary_conditions(mdl)           # apply the neumann boundary conditions
    else
        A_bar = assemble_system_A(mdl)               # assemble the stiffness matrix
        B = assemble_system_B(mdl)                   # assemble the stiffness matrix
        b = apply_boundary_conditions(mdl)           # apply the neumann boundary conditions
    end

    q_d = (μu_btm*q_btm + μu_tp*q_tp + μu_side*q_side)      # apply the Dirichlet boundary conditions

    if GRAD
        dAdη = A_bar
        dAdβ = b        
    end

    A = η*A_bar + β*b

    C_Tu = transpose(C_uc)           # transpose the constraint matrix

    A_free = C_Tu*A*C_uc        # extract the free part of the stiffness matrix
    B_free = C_Tu*B             # extract the free part of the stiffness matrix

    K_free = [A_free B_free; B_free' zeros(size(B_free,2),size(B_free,2))]     # assemble the system of equations

    r = [C_Tu*A*q_d; B'*q_d]    # assemble the system of equations
    sol = -K_free\Matrix(r)                 # solve the system of equations

    q_f = sol[1:size(A_free,1)]     # extract the free part of the solution
    p_f = sol[size(A_free,1)+1:end] # extract the free part of the solution 

    q = q_d + C_uc*q_f;                 # assemble the solution 
    p = p_f;

    q_out = [q[ID_u[1,:]] q[ID_u[2,:]] q[ID_u[3,:]]]'

    if GRAD
        dKdη = [C_Tu*dAdη*C_uc zeros(size(B_free)); zeros(size(B_free')) zeros(size(B,2),size(B,2))] # assemble the system of equations
        dKdβ = [C_Tu*dAdβ*C_uc zeros(size(B_free)); zeros(size(B_free')) zeros(size(B,2),size(B,2))] # assemble the system of equations
    
        drdη = [C_Tu*dAdη*q_d; zeros(size(B,2),size(q_d,2))] # solve the system of equations
        drdβ = [C_Tu*dAdβ*q_d; zeros(size(B,2),size(q_d,2))] # solve the system of equations

        dsoldη = -K_free\(drdη + dKdη*sol) # solve the system of equations
        dsoldβ = -K_free\(drdβ + dKdβ*sol) # solve the system of equations

        dqfdη = dsoldη[1:size(A_free,1)] # extract the free part of the solution
        dqfdβ = dsoldβ[1:size(A_free,1)] # extract the free part of the solution

        dpfdη = dsoldη[size(A_free,1)+1:end] # extract the free part of the solution 
        dpfdβ = dsoldβ[size(A_free,1)+1:end] # extract the free part of the solution

        dqdη = C_uc*dqfdη;              # assemble the solution
        dqdβ = C_uc*dqfdβ;              # assemble the solution

        dpdη = dpfdη;                  # assemble the solution
        dpdβ = dpfdβ;                  # assemble the solution

        dqdη_out = hcat(dqdη[ID_u[1,:]], dqdη[ID_u[2,:]], dqdη[ID_u[3,:]])'
        dqdβ_out = hcat(dqdβ[ID_u[1,:]], dqdβ[ID_u[2,:]], dqdβ[ID_u[3,:]])'

        dqdθ_out = cat(dqdη_out,dqdβ_out,dims=(3,3)) # concatenate the gradients in to a tensor
        return q_out, dqdθ_out, mdl
    else
        return q_out, mdl
    end
end

function stokes_single_step_force(mdl::Stokes, scene::SqueezeFlow, conditions::Conditions)

    start_time = Dates.now()
    reset_model!(mdl)
    
    @unpack FunctionClass, IEN, IEN_cp, ID, NodeList, C_vol, W = mdl.mesh_x
    FunctionClass_x_cached::String = FunctionClass
    NodeList_x_cached::Matrix{Float64} = NodeList
    IEN_x_cached::Matrix{Int} = IEN
    IEN_x_cp_cached::Matrix{Int} = IEN_cp
    ID_x_cached::Matrix{Int} = ID
    C_vol_x_cached = C_vol
    W_x_cached = W

    @unpack IEN, ID, FunctionClass, top_nodes, bottom_nodes, side_nodes, nNodes = mdl.mesh_u
    FunctionClass_u_cached::String = FunctionClass
    IEN_u_cached::Matrix{Int} = IEN
    ID_u_cached::Matrix{Int} = ID
    nNodes_u_cached::Int = nNodes
    NodeList_u_cached::Matrix{Float64} = NodeList
    top_node_list_cached::Vector{Int} = top_nodes # top nodes
    bottom_node_list_cached::Vector{Int} = bottom_nodes # bottom nodes
    side_node_list_cached::Vector{Int} = side_nodes

    @unpack η, nDof_u, nDof_p, ndim = mdl
    η_cached::Any = η
    nDof_u_cached::Int = nDof_u
    nDof_p_cached::Int = nDof_p
    ndim_cached::Int = ndim

    @unpack β, viscosity_type, q_d, C_uc, t_steps, sim_time, control, cParam = scene
    β_cached::Any = β
    viscosity_type_cached::String = viscosity_type
    q_d_cached::Dict{Symbol, Matrix{Float64}} = q_d
    q_d_cached_top::SparseMatrixCSC{Bool, Int64} = q_d_cached[:top]
    q_d_cached_btm::SparseMatrixCSC{Bool, Int64} = q_d_cached[:bottom]
    q_d_cached_brdr::SparseMatrixCSC{Bool, Int64} = q_d_cached[:border]
    C_uc_cached::SparseMatrixCSC{Bool, Int64} = C_uc
    t_steps_cached::Float64 = t_steps
    sim_time_cached::Float64 = sim_time
    control_cached::String = control
    cParam_cached::Vector{Float64} = scene.cParam

    @unpack nNodes = mdl.mesh_p
    nNodes_p_cached::Int = nNodes

    @unpack camera_matrix, obj_pose, SIDES = conditions    
    camera_matrix_cached::Matrix{Float64} = camera_matrix
    obj_pose_cached::Matrix{Float64} = obj_pose
    SIDES_cached::Bool = SIDES
    
    NodeList_cached::Matrix{Float64} = NodeList_u_cached
    ID_cached::Matrix{Int} = ID_u_cached
    T = Matrix{Float64}(I, size(NodeList_x_cached,2), size(NodeList_u_cached,2)) # projection matrix from geometry to field mesh

    time = collect(Float64, range(start=t_steps_cached, stop=sim_time_cached, step=t_steps_cached))
    len_t = length(time)

    if FunctionClass_x_cached == "S2" && FunctionClass_u_cached != "S2"
        T = get_nurbs_2_lagrange_proj(IEN_x_cached, IEN_u_cached, C_vol_x_cached, NodeList_x_cached, W_x_cached)
        NodeList_cached = NodeList_x_cached
    end

    T_ = T'*inv(T*T')
    C_Tu = transpose(C_uc_cached) # transpose the constraint matrix

    if conditions.filepath != ""
        isnothing(conditions.filepath) || AssertionError("Please provide a filepath to write the data")
        set_file(conditions.filepath)
    end

    μu_btm = 0  
    μu_side = 0
    
    NodeList_proj = NodeList_cached*T # project the motion on the geometry mesh grid
            
    BorderPts2D, SurfacePts2D = extract_borders(NodeList_proj, camera_matrix_cached, obj_pose_cached, nNodes_u_cached, BorderNodesList=side_node_list_cached)
    pi, qi = fit_curve(border=BorderPts2D)
    
    dqdη = zeros(Float64, size(q_d_cached_top))
    dqdβ = zeros(Float64, size(q_d_cached_top))

    displacement = AbstractArray[zeros(Float64,size(NodeList_proj,1),size(NodeList_proj,2))] # store the displacement of the mesh in 3D
    surface_fields = AbstractArray[]
    surface_pts_3D = AbstractArray[vcat(NodeList_proj[:,top_node_list_cached]', 
                                        NodeList_proj[:,bottom_node_list_cached]', 
                                        NodeList_proj[:,side_node_list_cached]')']      # store the solution fields of the mesh in 3D
    gradList = AbstractArray[zeros(Float64, size(BorderPts2D,1),size(BorderPts2D,2),2)] # store the solution fields of the border nodes in 2D 
    pos3D = AbstractArray[NodeList_proj]         # store the solution fields of the mesh in 3D
    pos3D_cp = AbstractArray[NodeList_cached]  
    pos2D = AbstractArray[SurfacePts2D]          # store the solution fields of the mesh in 2D
    borderPts2DList = AbstractArray[BorderPts2D] # store the solution fields of the surfaces in 2D
    splinep = AbstractArray[BorderPts2D[1,:]]    # store the x coordinates samples of the spline parameters of the border nodes
    splineq = AbstractArray[BorderPts2D[2,:]]    # store the y coordinates samples of the spline parameters of the border nodes
    output = Float64[] 
    writeborderList = [vcat(pi', qi')]
    iter = 1
    
    if control_cached == "force"

        A_bar = SparseMatrixCSC{Float64,Int}(I, nDof_u_cached*(nNodes_u_cached)^ndim_cached, nDof_u_cached*(nNodes_u_cached)^ndim_cached)  # initialize the stiffness matrix
        B = SparseMatrixCSC{Float64,Int}(I, nDof_u_cached*(nNodes_u_cached)^ndim_cached, nDof_p_cached*(nNodes_p_cached)^ndim_cached)      # initialize the stiffness matrix
        b = SparseMatrixCSC{Float64,Int}(I, nDof_u_cached*(nNodes_u_cached)^ndim_cached, nDof_u_cached*(nNodes_u_cached)^ndim_cached)      # initialize the stiffness matrix
        q_d = spzeros(nDof_u_cached*(nNodes_u_cached)^ndim_cached,1)                                                                       # initialize the vector of the Dirichlet boundary conditions (for ndof = 1) / Dirichlet boundary conditions upper surface (for ndof > 1)
        A = similar(A_bar)

        A_free = SparseMatrixCSC{Float64, Int64}(I, size(C_Tu,1),size(C_uc_cached,2))   # convert to sparse matrix
        B_free = SparseMatrixCSC{Float64, Int64}(I, size(C_Tu,1),size(B,2))             # convert to sparse matrix

        dA_freedη = similar(A_free)                         
        dA_freedβ = similar(A_free)                         
        dB_free = spzeros(size(B_free))                     
        zero = spzeros(size(B_free,2),size(B_free,2))

        dAdη = similar(A_bar)
        dAdβ = similar(A_bar)
        dB = spzeros(size(B))

        q = similar(q_d)
        dqfdη = similar(q)
        dqfdβ = similar(q)

        M = spzeros((size(A_free,1)+size(B_free,2)+1),(size(A_free,2)+size(B_free,2)+1))
        dMdη = spzeros(size(M))
        dMdβ = spzeros(size(M))

        A_bar .= assemble_system_A(mdl)     # assemble the stiffness matrix
        B .= assemble_system_B(mdl)         # assemble the stiffness matrix
        b .= apply_boundary_conditions(mdl) # apply the neumann boundary conditions
    
        q_d .= (μu_btm*q_d_cached_btm + μu_side*q_d_cached_brdr)      # apply the Dirichlet boundary conditions
        
        if viscosity_type_cached == "bulk_viscosity"
            if length(β_cached) == 1
                A = η_cached[iter]*A_bar + β_cached[1]*b
            else
                A = η_cached[iter]*A_bar + β_cached[iter]*b
            end
        else
            A .= η_cached[1]*A_bar + β_cached[1]*b
        end
    
        dAdη .= A_bar
        dAdβ .= b

        A_free .= C_Tu*A*C_uc_cached # extract the free part of the stiffness matrix
        B_free .= C_Tu*B             # extract the free part of the stiffness matrix

        dA_freedη .= C_Tu*dAdη*C_uc_cached # extract the free part of the stiffness matrix
        dA_freedβ .= C_Tu*dAdβ*C_uc_cached # extract the free part of the stiffness matrix

        M[1:size(A_free,1),1:size(A_free,2)] = A_free
        M[(size(A_free,1)+1):(size(A_free,1)+size(B_free,2)),1:size(A_free,2)] = B_free'
        M[end,1:size(A_free,2)] = q_d_cached_top'*A*C_uc_cached

        M[1:size(A_free,1),(size(A_free,2)+1):(size(A_free,2)+size(B_free,2))] = B_free
        M[(size(A_free,1)+1):(size(A_free,1)+size(B_free,2)),(size(A_free,2)+1):(size(A_free,2)+size(B_free,2))] = zero
        M[end,(size(A_free,2)+1):(size(A_free,2)+size(B_free,2))] = q_d_cached_top'*B

        M[1:size(A_free,1),end] = C_Tu*A*q_d_cached_top
        M[(size(A_free,1)+1):(size(A_free,1)+size(B_free,2)),end] = B'*q_d_cached_top
        M[end,end] = (q_d_cached_top'A*q_d_cached_top)[end]
        
        dMdη[1:size(A_free,1),1:size(A_free,2)] = dA_freedη
        dMdη[(size(A_free,1)+1):(size(A_free,1)+size(B_free,2)),1:size(A_free,2)] = dB_free'
        dMdη[end,1:size(A_free,2)] = q_d_cached_top'*dAdη*C_uc_cached

        dMdη[1:size(A_free,1),(size(A_free,2)+1):(size(A_free,2)+size(B_free,2))] = dB_free
        dMdη[(size(A_free,1)+1):(size(A_free,1)+size(B_free,2)),(size(A_free,2)+1):(size(A_free,2)+size(B_free,2))] = zero
        dMdη[end,(size(A_free,2)+1):(size(A_free,2)+size(B_free,2))] = q_d_cached_top'*dB

        dMdη[1:size(A_free,1),end] = C_Tu*dAdη*q_d_cached_top
        dMdη[(size(A_free,1)+1):(size(A_free,1)+size(B_free,2)),end] = dB'*q_d_cached_top
        dMdη[end,end] = (q_d_cached_top'dAdη*q_d_cached_top)[end]
        
        dMdβ[1:size(A_free,1),1:size(A_free,2)] = dA_freedβ
        dMdβ[(size(A_free,1)+1):(size(A_free,1)+size(B_free,2)),1:size(A_free,2)] = dB_free'
        dMdβ[end,1:size(A_free,2)] = q_d_cached_top'*dAdβ*C_uc_cached

        dMdβ[1:size(A_free,1),(size(A_free,2)+1):(size(A_free,2)+size(B_free,2))] = dB_free
        dMdβ[(size(A_free,1)+1):(size(A_free,1)+size(B_free,2)),(size(A_free,2)+1):(size(A_free,2)+size(B_free,2))] = zero
        dMdβ[end,(size(A_free,2)+1):(size(A_free,2)+size(B_free,2))] = q_d_cached_top'*dB

        dMdβ[1:size(A_free,1),end] = C_Tu*dAdβ*q_d_cached_top
        dMdβ[(size(A_free,1)+1):(size(A_free,1)+size(B_free,2)),end] = dB'*q_d_cached_top
        dMdβ[end,end] = (q_d_cached_top'dAdβ*q_d_cached_top)[end]

        r = [-C_Tu*A*q_d; -B'*q_d; cParam_cached[iter].-q_d_cached_top'A*q_d]    # assemble the system of equations
        drdη = -[C_Tu*dAdη*q_d; zeros(Float64, size(B,2),size(q_d,2)); q_d_cached_top'dAdη*q_d] # solve the system of equations
        drdβ = -[C_Tu*dAdβ*q_d; zeros(Float64, size(B,2),size(q_d,2)); q_d_cached_top'dAdβ*q_d] # solve the system of equations

        lum = lu(M) # LU decomposition of the system of equations

        sol = lum\Matrix(r)             # solve the system of equations
        dsoldη = lum\(drdη - dMdη*sol)  # solve the system of equations
        dsoldβ = lum\(drdβ - dMdβ*sol)  # solve the system of equations

        q_f = view(sol, 1:size(A_free, 1))
        dqfdη = view(dsoldη, 1:size(A_free, 1))
        dqfdβ = view(dsoldβ, 1:size(A_free, 1))

        p_f = view(sol, (size(A_free, 1) + 1):(size(A_free, 1) + size(B_free, 2)))
        dpfdη = view(dsoldη, (size(A_free, 1) + 1):(size(A_free, 1) + size(B_free, 2)))
        dpfdβ = view(dsoldβ, (size(A_free, 1) + 1):(size(A_free, 1) + size(B_free, 2)))

        μ_tp = sol[end]
        dμdη = dsoldη[end]
        dμdβ = dsoldβ[end]
        
        q .= q_d + C_uc_cached*q_f + μ_tp*q_d_cached_top;       # assemble the solution 
        dqdη .= dqdη + C_uc_cached*dqfdη + dμdη*q_d_cached_top; # assemble the solution
        dqdβ .= dqdβ + C_uc_cached*dqfdβ + dμdβ*q_d_cached_top; # assemble the solution

        p = p_f;
        dpdη = dpfdη; # assemble the solution
        dpdβ = dpfdβ; # assemble the solution

        motion_y = @views hcat(q[ID_cached[1,:]], q[ID_cached[2,:]], q[ID_cached[3,:]])'*t_steps_cached # extract the motion of the mesh grid
        dmdη_out_y = @views hcat(dqdη[ID_cached[1,:]], dqdη[ID_cached[2,:]], dqdη[ID_cached[3,:]])'*t_steps_cached
        dmdβ_out_y = @views hcat(dqdβ[ID_cached[1,:]], dqdβ[ID_cached[2,:]], dqdβ[ID_cached[3,:]])'*t_steps_cached

        motion =  motion_y # extract the motion of the mesh grid

        NodeList_cached = NodeList_cached + motion # update the mesh grid
        mdl.mesh_x.NodeList = NodeList_cached      # update the mesh grid

        NodeList_proj = NodeList_cached # project the motion on the geometry mesh grid
        dmdη_out_proj = dmdη_out_y
        dmdβ_out_proj = dmdβ_out_y
        motion_proj = motion
        
        mat_nan_inf_check(dmdη_out_y)
        mat_nan_inf_check(dmdβ_out_y)
        
        dmdθ_out = @views cat(dmdη_out_proj,dmdβ_out_proj,dims=3) # concatenate the gradients in to a tensor

        BorderPts2D, dudθ, SurfacePts2D, ∇SurfacePts2D = extract_borders(NodeList_proj, camera_matrix_cached, obj_pose_cached, BorderNodesList=side_node_list_cached, GRAD=true, dqdθ=dmdθ_out, SIDES=SIDES_cached)
        pi, qi = fit_curve(border=BorderPts2D)

        push!(output, μ_tp*t_steps_cached) # store displacement at the top surface
        push!(displacement, motion_proj)
        push!(surface_fields, motion_proj[:,side_node_list_cached])
        push!(surface_pts_3D, vcat(NodeList_proj[:,top_node_list_cached]', NodeList_proj[:,bottom_node_list_cached]', NodeList_proj[:,side_node_list_cached]')')
        push!(gradList,dudθ)
        push!(pos2D, SurfacePts2D)
        push!(pos3D, NodeList_proj)
        push!(pos3D_cp, NodeList_cached)
        push!(borderPts2DList, BorderPts2D)
        push!(splinep, BorderPts2D[1,:])
        push!(splineq, BorderPts2D[2,:])
        push!(writeborderList, vcat(pi', qi'))

    elseif control_cached == "velocity"

        A_bar = assemble_system_A(mdl)     # assemble the stiffness matrix
        B = assemble_system_B(mdl)         # assemble the stiffness matrix
        b = apply_boundary_conditions(mdl) # apply the neumann boundary conditions
    
        q_d = (μu_btm*q_d_cached_btm + cParam_cached[iter]*q_d_cached_top + μu_side*q_d_cached_brdr)      # apply the Dirichlet boundary conditions
    
        if viscosity_type_cached == "bulk_viscosity"
            if length(β_cached) == 1
                A = η_cached[iter]*A_bar + β_cached[1]*b
            else
                A = η_cached[iter]*A_bar + β_cached[iter]*b
            end
        else
            A = η_cached[1]*A_bar + β_cached[1]*b
        end

        dAdη = A_bar
        dAdβ = b
        dB = zeros(Float64, size(B))

        C_Tu = transpose(C_uc_cached)      # transpose the constraint matrix
    
        A_free = C_Tu*A*C_uc_cached        # extract the free part of the stiffness matrix
        B_free = C_Tu*B                    # extract the free part of the stiffness matrix

        dA_freedη = C_Tu*dAdη*C_uc_cached        # extract the free part of the stiffness matrix
        dA_freedβ = C_Tu*dAdβ*C_uc_cached        # extract the free part of the stiffness matrix
        dB_free = zeros(Float64, size(B_free))
    
        K_free = [A_free B_free; B_free' zeros(Float64, size(B_free,2),size(B_free,2))]      # assemble the system of equations
        dKdη = [C_Tu*dAdη*C_uc_cached dB_free; dB_free' zeros(Float64, size(B,2),size(B,2))] # assemble the system of equations
        dKdβ = [C_Tu*dAdβ*C_uc_cached dB_free; dB_free' zeros(Float64, size(B,2),size(B,2))] # assemble the system of equations
        
        invK = inv(Matrix(K_free))
    
        r = [C_Tu*A*q_d; B'*q_d]    # assemble the system of equations
        drdη = [C_Tu*dAdη*q_d; zeros(Float64, size(B,2),size(q_d,2))] # solve the system of equations
        drdβ = [C_Tu*dAdβ*q_d; zeros(Float64, size(B,2),size(q_d,2))] # solve the system of equations

        sol = -invK*r                    # solve the system of equations
        dsoldη = -invK*(drdη + dKdη*sol) # solve the system of equations
        dsoldβ = -invK*(drdβ + dKdβ*sol) # solve the system of equations
    
        q_f = sol[1:size(A_free,1)]      # extract the free part of the solution
        dqfdη = dsoldη[1:size(A_free,1)] # extract the free part of the solution
        dqfdβ = dsoldβ[1:size(A_free,1)] # extract the free part of the solution

        p_f = sol[size(A_free,1)+1:end]      # extract the free part of the solution
        dpfdη = dsoldη[size(A_free,1)+1:end] # extract the free part of the solution 
        dpfdβ = dsoldβ[size(A_free,1)+1:end] # extract the free part of the solution
    
        q = q_d + C_uc_cached*q_f;         # assemble the solution 
        dqdη = dqdη + C_uc_cached*dqfdη;   # assemble the solution
        dqdβ = dqdβ + C_uc_cached*dqfdβ;   # assemble the solution

        p = p_f;
        dpdη = dpfdη;  # assemble the solution
        dpdβ = dpfdβ;  # assemble the solution

        motion = hcat(q[ID_cached[1,:]], q[ID_cached[2,:]], q[ID_cached[3,:]])'*t_steps_cached # get the motion of the mesh
        dmdη_out = hcat(dqdη[ID_cached[1,:]], dqdη[ID_cached[2,:]], dqdη[ID_cached[3,:]])'
        dmdβ_out = hcat(dqdβ[ID_cached[1,:]], dqdβ[ID_cached[2,:]], dqdβ[ID_cached[3,:]])'
        
        NodeList_cached = NodeList_cached + motion # update the mesh grid
        mdl.mesh_u.NodeList = NodeList_cached # update the mesh grid
    
        NodeList_proj = NodeList_cached*T # project the motion on the geometry mesh grid
        dmdη_out_proj = dmdη_out*T
        dmdβ_out_proj = dmdβ_out*T
        
        dmdθ_out = @views cat(dmdη_out_proj,dmdβ_out_proj,dims=3) # concatenate the gradients in to a tensor
    
        BorderPts2D, dudθ, SurfacePts2D, ∇SurfacePts2D = extract_borders(NodeList_proj, camera_matrix_cached, obj_pose_cached, side_node_list_cached, GRAD=true, dqdθ=dmdθ_out, SIDES=SIDES_cached)
        pi, qi = fit_curve(border=BorderPts2D)

        mat_nan_inf_check(dudθ[:,:,1])
        mat_nan_inf_check(dudθ[:,:,2])

        push!(output, μ_tp*t_steps_cached) # store displacement at the top surface
        push!(displacement, motion)
        push!(surface_fields, motion[:,side_node_list_cached])
        push!(surface_pts_3D, NodeList_proj[:,side_node_list_cached]')
        push!(gradList,dudθ)
        push!(pos2D, SurfacePts2D)
        push!(pos3D, NodeList_proj)
        push!(pos3D_cp, NodeList_cached)
        push!(borderPts2DList, BorderPts2D)
        push!(splinep, BorderPts2D[1,:])
        push!(splineq, BorderPts2D[2,:]) 
        push!(writeborderList, vcat(pi', qi'))
    else
            throw(ArgumentError("Control type not unknown"))
    end
    end_time = Dates.now()
    elapsed_time = end_time - start_time
    return output, gradList, borderPts2DList, displacement, surface_pts_3D, pos2D, splinep, splineq, elapsed_time
end
                                    """
write_data(exp_params::Dict)

    Writes simulation data to files based on the provided experiment parameters.

# Arguments
- `exp_params::Dict`: A dictionary containing experiment parameters such as mesh size, material properties, file paths, and simulation settings.

# Returns
None.
"""

function write_gt_data(exp_params::Dict)

    ndim::Int = 3
    FunctionClass_u::String = exp_params["FunctionClass_u"]
    FunctionClass_p::String = exp_params["FunctionClass_p"]    
    FunctionClass_x::String = exp_params["FunctionClass_x"]
    nDof_u::Int = ndim  # number of degree of freedom per node
    nDof_p::Int = 1  # number of degree of freedom per node

    # simulation parameters for the ground truth
    r::Float64 = exp_params["r"]  # radius of the cylinder in mm
    h::Float64 = exp_params["h"]  # height of the cylinder in mm

    control::String = exp_params["control"]
    viscosity_type::String = exp_params["viscosity_type"] # "constant" or "bulk_viscosity"
    filepath_gt::String = exp_params["filepath_gt"]
    obj_pose::AbstractArray = exp_params["obj_pose_gt"]
    camera_matrix::AbstractArray = exp_params["camera_matrix"]

    sim_time_gt::Float64 = exp_params["sim_time_gt"]
    steps_gt::Float64 = exp_params["steps_gt"]
    t_steps_gt::Float64 = sim_time_gt/steps_gt

    β_gt::Float64 = exp_params["β_gt"]
    η_gt::Float64 = exp_params["η_gt"]
    ne_gt::Int = exp_params["ne_gt"] # number of elements in the mesh for the ground truth

    F_ext::Float64 = exp_params["F_ext"]

    F = -F_ext*ones(Float64, round(Int, (sim_time_gt/t_steps_gt))) # force applied to the cylinder in N

    # Write the ground truth
    printstyled("Ground truth η: $(η_gt), ground truth β: $(β_gt)\n"; color = :green)
    model_gt, scene_gt = def_problem(r, h, ne_gt, η_gt, ndim, FunctionClass_u, nDof_u, FunctionClass_p, nDof_p, FunctionClass_x, β_gt, F, control, viscosity_type, 
                    sim_time_gt, t_steps_gt)

    @info "Writing ground truth gt data to with $ne_gt elements to $filepath_gt"
    write_sim_data(model_gt, scene_gt, camera_matrix, obj_pose, filepath_gt)

  end

"""
    write_sim_data(r, h, ne, Young, ν, ndim, FunctionClass, nDof, β, CameraMatrix, endTime, tSteps, Control, filename; mode="lame")

Initializes the simulation and writes the data to a file.

# Arguments
- `r::Number`: Radius of the cylinder.
- `h::Number`: Height of the cylinder.
- `ne::Int64`: Number of elements.
- `Young::Number`: Young's modulus.
- `ν::Number`: Poisson's ratio.
- `ndim::Int64`: Number of dimensions.
- `FunctionClass::String`: Type of basis function.
- `nDof::Int64`: Number of degrees of freedom per node.
- `β::Number`: Friction parameter.
- `CameraMatrix::AbstractMatrix{Float64}`: Camera matrix.
- `endTime::Number`: End time of the simulation.
- `tSteps::Number`: Number of time steps.
- `Control::String`: Type of control ("force" or "displacement").
- `filename::String`: Path to the output file.

# Keyword Arguments
- `mode::String`: Type of constitutive matrix ("lame" by default).

# Returns
None.
"""
function write_sim_data(_model::AbstractModel, _scene::AbstractScenario, camera_matrix::AbstractMatrix{Float64}, obj_pose::AbstractMatrix{Float64}, 
                        filepath::String="nothing"; ANIMATE=true, WRITECONTOUR=true, RENDER=true, WRITEVTK=true)

    # copy the model and scene to avoid modifying the original objects
    model::AbstractModel = deepcopy(_model)
    scene::AbstractScenario = deepcopy(_scene)
    
    # set the model parameters
    conditions = Conditions(ANIMATE=ANIMATE, WRITECONTOUR=WRITECONTOUR, RENDER=RENDER, WRITEVTK=WRITEVTK, camera_matrix=camera_matrix, obj_pose=obj_pose, filepath=filepath)
    
    # remove the previous results folder if it exists
    if isdir(string(filepath))
        @info "Removing previous results folder: $filepath"
        rm(string(filepath), recursive=true, force=true) # remove the previous results folder if it exists
    end

    h_, gradList, borderPts2DList, displacement, surface_pts_3D, pos2D, pos3D, splinep, splineq, velocity, pressure = simulate(model, scene, conditions) # run the simulation
    
    h = get_height(h_, model.mesh_u.h) # get the mesh height with time
    display(scene.cParam)
    params = Dict("r"=>model.mesh_u.r, "h"=>model.mesh_u.h, "η" => model.η, "β" => scene.β, "camera_matrix" => conditions.camera_matrix, "obj_pose" => conditions.obj_pose, 
                    "control_type"=>scene.control, "cParam"=>scene.cParam, "simulation_time" => scene.sim_time, "time_steps" => scene.t_steps, 
                    "viscosity_type"=>scene.viscosity_type)

    # write results to files
    write_json(string(filepath,"/data/sim_params"), params)
    write_csv(string(filepath,"/data/h"), h)
    write_data(string(conditions.filepath,"/data/sim_data/2D_surface_points"), pos2D)
    write_data(string(filepath,"/data/sim_data/node_points"), pos3D)
    write_data(string(filepath,"/data/sim_data/displacement_fields"), displacement)
    write_data(string(filepath,"/data/sim_data/2D_border_points"), borderPts2DList)
    write_data(string(filepath,"/data/sim_data/velocity_fields"), velocity)
    write_data(string(filepath,"/data/sim_data/pressure_fields"), pressure)

    if _scene.viscosity_type == "bulk_viscosity"
        time = collect(Float64, range(scene.t_steps, stop=scene.sim_time, step=scene.t_steps))
        # println(size(model.η))
        # println(size(time))
        plt = set_plot(22)
        Plots.plot!(plt, time, model.η, lw=3, label="Shear Viscosity", color=:red)
        Plots.xlabel!(plt, L"\mathrm{Time\;(s)}")
        Plots.ylabel!(plt, L"\eta(\mathrm{t})\;\mathrm{(KPa\cdot s)}")
        savefig(plt, string(filepath,"/Results/images/shear_viscosity_time.pdf"))
    end

    @info "Data written to $filepath"

end

"""
    compare(r, h, ne, Young, ν, ndim, FunctionClass, nDof, β, CameraMatrix, endTime, tSteps, Control, mode, ObsData; SIDES=false, PLOT=false, filepath="nothing")

Compares the simulation results with observation data.

# Arguments
- `r::Number`: Radius of the cylinder.
- `h::Number`: Height of the cylinder.
- `ne::Int64`: Number of elements.
- `Young::Number`: Young's modulus.
- `ν::Number`: Poisson's ratio.
- `ndim::Int64`: Number of dimensions.
- `FunctionClass::String`: Type of basis function.
- `nDof::Int64`: Number of degrees of freedom per node.
- `β::Number`: Friction parameter.
- `CameraMatrix::AbstractMatrix{Float64}`: Camera matrix.
- `endTime::Number`: End time of the simulation.
- `tSteps::Number`: Number of time steps.
- `Control::String`: Type of control ("force" or "displacement").
- `mode::String`: Type of constitutive matrix.
- `ObsData::AbstractArray`: Observation data.

# Keyword Arguments
- `SIDES::Bool`: Whether to include side boundaries (default: `false`).
- `PLOT::Bool`: Whether to generate plots (default: `false`).
- `filepath::String`: Path to the output file (default: `"nothing"`).

# Returns
- `d_cp::Vector{Float64}`: Closest point distances between borders at each timestep.
"""
function compare(r::Number, h::Number, ne::Int64, Young::Number, ν::Number, ndim::Int64, FunctionClass::String, nDof::Int64, β::Number, CameraMatrix::AbstractMatrix{Float64}, 
                endTime::Number, tSteps::Number, Control::String, mode::String, ObsData::AbstractArray, SIDES::Bool=false, PLOT::Bool=false, 
                filepath::String="nothing")

    μ_list, gradList, simBorderPts, splinex, spliney, mdl, SurfacePt2D = test(r, h, ne, Young, ν, ndim, FunctionClass, nDof, β, CameraMatrix, 
                                                                            endTime, tSteps, Control, mode=mode, SIDES=SIDES)    

    obsBorderPts = ObsData[1]
    splinexObs = ObsData[2]
    splineyObs = ObsData[3]
    
    # test the closest point function
    
    d_cp, pairs = closest_point(simBorderPts, obsBorderPts)
    
    if PLOT
        @assert !isnothing(filepath) "File path not provided"
        plot_matches(simBorderPts, splinex, spliney, splinexObs, splineyObs, pairs, string(filepath,"/Results/images/"))
    end
    d_h = 0
    return d_h, d_cp
end

"""
    initialize_mesh_test(r, h, ne, FunctionClass, CameraMatrix; filepath="nothing", SIDES=false)

Initializes the mesh and writes the data to a file.

# Arguments
- `r::Number`: Radius of the cylinder.
- `h::Number`: Height of the cylinder.
- `ne::Int64`: Number of elements.
- `FunctionClass::String`: Type of basis function.
- `CameraMatrix::AbstractMatrix{Float64}`: Camera matrix.

# Keyword Arguments
- `filepath::String`: Path to the output file (default: `"nothing"`).
- `SIDES::Bool`: Whether to include side boundaries (default: `false`).

# Returns
- `borderPts2DList::Vector{Matrix{Float64}}`: List of 2D border points at each timestep.
- `surfaceNodesList::Vector{Matrix{Float64}}`: List of surface nodes of the mesh at each timestep.
- `splinep::Vector{Vector{Float64}}`: x coordinates of the border observation at each timestep interpolated.
- `splineq::Vector{Vector{Float64}}`: y coordinates of the border observation at each timestep interpolated.
"""
function initialize_mesh(r::Number, h::Number, ne::Int64, FunctionClass::String, camera_matrix::AbstractMatrix{Float64}, obj_pose::AbstractMatrix{Float64},filepath::String="nothing", 
                            SIDES::Bool=false)

    isnothing(filepath) || AssertionError("Please provide a filepath to write the data")
    set_file(filepath)
    
    ndim::Int = 3
    # mesh = meshgrid_cylinder(r, h, ne, FunctionClass=FunctionClass)
    filepath_mesh = joinpath("/home", "soshala", "SMEAR-PhD", "smear-modules", "smear-meshes")
    mesh = get_gmsh_cylinder(joinpath(filepath_mesh, "cylinder_x_$ne.msh"), 1, r, h, FunctionClass)

    BorderPts2D, SurfacePts2D = extract_borders(mesh.NodeList, camera_matrix, obj_pose, 13, BorderNodesList= mesh.side_nodes)
    pi, qi = fit_curve(border=BorderPts2D)
        
                                               # store the solution fields of the border nodes in 2D 
    pos3D = AbstractArray[mesh.NodeList]                                                             # store the solution fields of the mesh in 3D
    surface_pts_3D = [vcat(mesh.NodeList[:,mesh.top_nodes]', mesh.NodeList[:,mesh.bottom_nodes]', mesh.NodeList[:,mesh.side_nodes]')'] # store the solution fields of the mesh in 3D
    pos2D = AbstractArray[SurfacePts2D]                                                                   # store the solution fields of the mesh in 2D
    borderPts2DList = AbstractArray[BorderPts2D]                                                          # store the solution fields of the surfaces in 2D
    splinep = AbstractArray[pi]                                                                           # store the x coordinates samples of the spline parameters of the border nodes
    splineq = AbstractArray[qi]                                                                           # store the y coordinates samples of the spline parameters of the border nodes 
    writeborderList = [vcat(pi', qi')]

    write_scene(joinpath(filepath,"data"), pos3D, mesh.IEN, ndim, pos3D, function_class = mesh.FunctionClass)
    animate_fields(filepath = joinpath(filepath,"Results","images"), Nodes=pos3D , IEN=mesh.IEN, BorderNodes2D=borderPts2DList, fields2D=pos2D)
    write_data(joinpath(filepath,"Results"), writeborderList)

    return borderPts2DList, pos2D, splinep, splineq
end 

"""
    readData(filepath)

Reads the data from a file.

# Arguments
- `filepath::String`: Path to the file.

# Returns
- `obsBorderPts::Vector{Matrix{Float64}}`: list of observation data.
- `splinexObs::Vector{Vector{Float64}}`: x-coordinates of the border observation at each time step.
- `splineyObs::Vector{Vector{Float64}}`: y-coordinates of the border observation at each time step.
"""
function readData(filepath::String)
    obsBorderPts, splinexObs, splineyObs, pd = read_csv(string(filepath,"/Results/contour_data"))   
    animate_fields(filepath = string(filepath,"/Results/cost"),pObs=splinexObs, qObs=splineyObs) # animate the fields
    plot(x->pdf(pd, x))
    savefig(string(filepath,"/Results"))
    return obsBorderPts, xObs, yObs
end

"""
    plot_rad_vel(file_path)

Plots the mean radial velocity against the slip parameter.

# Arguments
- `file_path::String`: Path to the directory containing the simulation data.

# Returns
- A plot of mean radial velocity vs slip parameter.
"""
function plot_rad_norm_vel_vs_slip(file_path::String)
    dirs = readdir(file_path)
    v_rad_slip_list = zeros(Float64, (length(dirs)-1))
    v_norm_slip_list = zeros(Float64, (length(dirs)-1))
    β_list = zeros(Float64, (length(dirs)-1))
    for (i, dir) in enumerate(dirs)
        if dir == "analysis" || dir == "Results" || !isdir(joinpath(file_path, dir))
            continue
        end
        file = joinpath(file_path, dir, "data", "sim_data", "velocity_fields")
        exp_params = read_json(joinpath(file_path, dir, "data", "sim_params.json"))
        β = exp_params["β"]
        β_list[i] = β[1]
        csv_file = readdir(file)
        v_rad_mean_list = zeros(Float64, length(csv_file))
        v_norm_mean_list = zeros(Float64, length(csv_file))
        for (j, csv) in enumerate(csv_file)
            if csv == "000.csv"
                continue
            end
            data = readdlm(joinpath(file, csv), ',')
            v1 = data[2028:end,1]
            v2 = data[2028:end,2]
            v3 = data[2028:end,3]
            v_rad = sqrt.(v1.^2 + v2.^2)
            v_norm_mean = mean(abs.(v3))
            vrad_mean = mean(v_rad)
            v_rad_mean_list[j] = vrad_mean
            v_norm_mean_list[j] = v_norm_mean
        end
        v_rad_time_mean = mean(v_rad_mean_list)
        v_norm_time_mean = mean(v_norm_mean_list)
        v_rad_slip_list[i] = v_rad_time_mean
        v_norm_slip_list[i] = v_norm_time_mean
    end
    sort_idx = sortperm(β_list)
    β_list = β_list[sort_idx]
    v_rad_slip_list = v_rad_slip_list[sort_idx]
    v_norm_slip_list = v_norm_slip_list[sort_idx]
    
    plot_path = joinpath(file_path, "analysis")
    if !isdir(plot_path)
        mkdir(plot_path)
    end

    # println("Slip parameters: ", β_list)

    plt_rad = set_plot(10, sz=(400,360), bottom_margin=0.0mm)
    Plots.plot!(plt_rad, β_list, v_rad_slip_list, marker=:o, ms=2, xlabel=latexstring("\$\\beta\$ [kPa s mm\$^{-1}\$]"), ylabel="Mean Radial Velocity [mm/s]", xscale=:log10, yscale=:log10, legend=false)
    Plots.xticks!(plt_rad, [1, 1e2, 1e4, 1e6, 1e8, 1e10])
    Plots.savefig(string(plot_path,"/mean_radial_velocity_vs_slip_parameter.pdf"))

    plt_norm = set_plot(10, sz=(400,360), bottom_margin=0.0mm)
    Plots.plot!(plt_norm, β_list, v_norm_slip_list, marker=:o, ms=2, xlabel=latexstring("\$\\beta\$ [kPa s mm\$^{-1}\$]"), ylabel="Mean Normal Velocity [mm/s]", xscale=:log10, yscale=:log10, legend=false)
    Plots.xticks!(plt_norm, [1, 1e2, 1e4, 1e6, 1e8, 1e10])
    Plots.savefig(string(plot_path,"/mean_normal_velocity_vs_slip_parameter.pdf"))

end

function plot_rad_norm_vel_vs_visc(file_path::String)
    dirs = readdir(file_path)
    v_rad_slip_list = zeros(Float64, (length(dirs)-1))
    v_norm_slip_list = zeros(Float64, (length(dirs)-1))
    η_list = zeros(Float64, (length(dirs)-1))
    for (i, dir) in enumerate(dirs)
        if dir == "analysis" || dir == "Results" || !isdir(joinpath(file_path, dir))
            continue
        end
        file = joinpath(file_path, dir, "data", "sim_data", "velocity_fields")
        exp_params = read_json(joinpath(file_path, dir, "data", "sim_params.json"))
        η = exp_params["η"]
        η_list[i] = η[1]
        csv_file = readdir(file)
        v_rad_mean_list = zeros(Float64, length(csv_file)-1)
        v_norm_mean_list = zeros(Float64, length(csv_file)-1) 
        for (j, csv) in enumerate(csv_file)
            if csv == "000.csv"
                continue
            end
            j = j - 1 # adjust index to account for skipping the first file
            data = readdlm(joinpath(file, csv), ',')
            v1 = data[2028:end,1]
            v2 = data[2028:end,2]
            v3 = data[2028:end,3]
            v_rad = sqrt.(v1.^2 + v2.^2)
            v_norm_mean = mean(abs.(v3))
            vrad_mean = mean(v_rad)
            v_rad_mean_list[j] = vrad_mean
            v_norm_mean_list[j] = v_norm_mean
        end
        v_rad_time_mean = mean(v_rad_mean_list)
        v_norm_time_mean = mean(v_norm_mean_list)
        v_rad_slip_list[i] = v_rad_time_mean
        v_norm_slip_list[i] = v_norm_time_mean
    end
    sort_idx = sortperm(η_list)
    η_list = η_list[sort_idx]
    v_rad_slip_list = v_rad_slip_list[sort_idx]
    v_norm_slip_list = v_norm_slip_list[sort_idx]
    inv_eta_list = 1 ./ η_list
    plot_path = joinpath(file_path, "analysis")
    if !isdir(plot_path)
        mkdir(plot_path)
    end

    # println("Viscosity parameters: ", η_list)

    plt_rad = set_plot(10, sz=(400,360), bottom_margin=0.0mm)
    Plots.plot!(plt_rad, inv_eta_list, v_rad_slip_list, marker=:o, ms=2, xlabel=latexstring("\$\\eta\$ [kPa s]"), ylabel="Mean Radial Velocity [mm/s]", legend=false)
    # Plots.xticks!(plt_rad, [1, 1e2, 1e4, 1e6, 1e8, 1e10])
    Plots.savefig(string(plot_path,"/mean_radial_velocity_vs_viscosity_parameter.pdf"))

    plt_norm = set_plot(10, sz=(400,360), bottom_margin=0.0mm)
    Plots.plot!(plt_norm, inv_eta_list, v_norm_slip_list, marker=:o, ms=2, xlabel=latexstring("\$1 / \\eta\$ [kPa s]"), ylabel="Mean Normal Velocity [mm/s]", legend=false)
    # Plots.xticks!(plt_norm, [1, 1e2, 1e4, 1e6, 1e8, 1e10])
    Plots.savefig(string(plot_path,"/mean_normal_velocity_vs_viscosity_parameter.pdf"))

end