using smearFEM
using LinearAlgebra
using SparseArrays
using WriteVTK

using ProgressMeter
using Dates

function main()
  dateTime = Dates.now()

  x0 = 0
  x1 = 1
  y0 = 0
  y1 = 1
  z0 = 0
  z1 = 1
  ne = 4
  ndim = 3
  FunctionClass_u = "Q2"
  FunctionClass_p = "Q1"
  nDof_u = ndim  # number of degree of freedom per node
  nDof_p = 1
  β = 100
  ν = 1
  Control = "displacement"

  filepath = string("/home/soshala/SMEAR-PhD/SMEAR/Data/sim_experiments/single_simulation/Stokes/fem_runs/",Control,"/",dateTime,"/")
  writeData = true

  if writeData
    isnothing(filepath) || AssertionError("Please provide a filepath to write the data")
    set_file(filepath)
  end

  μu_btm = 0  
  μu_tp = -0.02
  μu_side = 0

  μp_btm = 0
  μp_tp = 0
  μp_side = 0

  NodeList_u, IEN_u, ID_u, IEN_u_top, IEN_u_btm, BorderNodesList_u = meshgrid_cube(x0,x1,y0,y1,z0,z1,ne,ndim,FunctionClass=FunctionClass_u)  # generate the mesh grid
  NodeListCylinder = inflate_cylinder(NodeList_u, x0, x1, y0, y1)
  q_tp, q_side, q_btm, C_uc = set_boundary_cond_stokes(NodeList_u, ne, ndim, FunctionClass_u, nDof_u)

  NodeList_p, IEN_p, ID_p, IEN_p_top, IEN_p_btm, BorderNodesList_p = meshgrid_cube(x0,x1,y0,y1,z0,z1,ne,ndim,FunctionClass=FunctionClass_p)  # generate the mesh grid
  NodeListCylinderp = inflate_cylinder(NodeList_p, x0, x1, y0, y1)

  mdl = def_model("stokes", ne=ne, NodeList=NodeListCylinder, IEN=IEN_u, IEN_top=IEN_u_top, IEN_btm=IEN_u_btm, ndim=ndim, nDof=nDof_u, 
                      FunctionClass=FunctionClass_u, ID=ID_u, IEN_2=IEN_p, IEN_2_top=IEN_p_top, IEN_2_btm=IEN_p_btm, 
                          ndim_2=ndim, nDof_2=nDof_u, FunctionClass_2=FunctionClass_p ) # define the model

  A_bar = assemble_system_A(mdl)                   # assemble the stiffness matrix
  B = assemble_system_B(mdl)                   # assemble the stiffness matrix
  b = apply_boundary_conditions_stokes(mdl)           # apply the neumann boundary conditions

  q_d = (μu_btm*q_btm + μu_tp*q_tp + μu_side*q_side)      # apply the Dirichlet boundary conditions

  A = A_bar + β*b

  C_Tu = transpose(C_uc)           # transpose the constraint matrix

  A_free = C_Tu*A*C_uc        # extract the free part of the stiffness matrix
  B_free = C_Tu*B             # extract the free part of the stiffness matrix

  K_free = [A_free B_free; B_free' zeros(size(B_free,2),size(B_free,2))]     # assemble the system of equations

  invK = inv(Matrix(K_free))

  r = -[C_Tu*A*q_d; B'*q_d]    # assemble the system of equations
  sol = invK*r                 # solve the system of equations

  q_f = sol[1:size(A_free,1)]     # extract the free part of the solution
  p_f = sol[size(A_free,1)+1:end] # extract the free part of the solution 

  q = q_d + C_uc*q_f;                 # assemble the solution 
  p = p_f;

  motion = [q[ID_u[:,1]] q[ID_u[:,2]] q[ID_u[:,3]]]'

  write_vtk(string(filepath,"/Results"), "u", NodeListCylinder, IEN_u, ne, ndim, motion, ID = ID_u, FunctionClass=FunctionClass_u)
  write_vtk(string(filepath,"/Results"), "p", NodeListCylinderp, IEN_p, ne, ndim, p)
end

main()