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
  element_shape_u = :Hex
  basis_order_u = 2
  element_shape_p = :Hex
  basis_order_p = 1
  nDof_u = ndim  # number of degree of freedom per node
  nDof_p = 1
  β = 100
  η = 1
  Control = "displacement"

    filepath = string(resolve_data_path("sim_experiments/single_simulation/Stokes/fem_runs/"), Control, "/", Date(dateTime), "/", Time(dateTime), "/")
  writeData = true

  if writeData
    isnothing(filepath) && throw(AssertionError("Please provide a filepath to write the data"))
    set_file(filepath)
  end

  μu_btm = 0  
  μu_tp = -0.02
  μu_side = 0

  mesh_u = meshgrid_cuboid(x1-x0, y1-y0, z1-z0; ne=ne, element_shape=element_shape_u, basis_order=basis_order_u, ndof=nDof_u)
  IEN_u, ID_u, IEN_u_top, IEN_u_btm = mesh_u.IEN, mesh_u.ID, mesh_u.IEN_top, mesh_u.IEN_bottom
  NodeListCylinder = _inflate_cylinder(mesh_u.NodeList, x0, x1, y0, y1)
  q_tp, q_side, q_btm, C_uc = set_boundary_cond_stokes(mesh_u.NodeList, ne, ndim, element_shape_u, basis_order_u, nDof_u)

  mesh_p = meshgrid_cuboid(x1-x0, y1-y0, z1-z0; ne=ne, element_shape=element_shape_p, basis_order=basis_order_p, ndof=nDof_p)
  IEN_p, ID_p, IEN_p_top, IEN_p_btm = mesh_p.IEN, mesh_p.ID, mesh_p.IEN_top, mesh_p.IEN_bottom
  NodeListCylinderp = _inflate_cylinder(mesh_p.NodeList, x0, x1, y0, y1)

  mdl = def_model("stokes", ne=ne, NodeList=NodeListCylinder, IEN=IEN_u, IEN_top=IEN_u_top, IEN_btm=IEN_u_btm, ndim=ndim, nDof=nDof_u,
                      element_shape=element_shape_u, basis_order=basis_order_u, ID=ID_u, IEN_2=IEN_p, IEN_2_top=IEN_p_top, IEN_2_btm=IEN_p_btm,
                          ndim_2=ndim, nDof_2=nDof_p, element_shape_2=element_shape_p, basis_order_2=basis_order_p, ν=η) # define the model

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

  write_vtk(string(filepath,"/Results"), "u", NodeListCylinder, IEN_u, ne, ndim, motion, ID = ID_u, element_shape=element_shape_u, basis_order=basis_order_u)
  write_vtk(string(filepath,"/Results"), "p", NodeListCylinderp, IEN_p, ne, ndim, p)
end

main()