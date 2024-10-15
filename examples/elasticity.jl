using LinearAlgebra
using ProgressMeter
using SparseArrays
using Plots

using smearFEM

# function test_meshgrid()
    x0 = 0
    x1 = 1
    y0 = 0
    y1 = 1
    z0 = 0
    z1 = 1
    ne = 8
    ndim = 3
    nDof = 1 # number of degree of freedom per node
    FunctionClass = "Q2"
    mode = "standard" # "standard" or "lame"
    Young = 30
    ν = 0.4

    μ_btm = 0 
    μ_tp = -1
    β = 100

    cMat = get_cMat(mode, Young, ν)

    NodeList, IEN, ID, IEN_top, IEN_btm, BorderNodesList = meshgrid(x0,x1,y0,y1,z0,z1,ne,ndim,FunctionClass=FunctionClass)  # generate the mesh grid
    NodeListCylinder = inflate_sphere(NodeList, x0, x1, y0, y1)    
    mdl = def_model(ne=ne, NodeList=NodeListCylinder, IEN=IEN, IEN_top=IEN_top, IEN_btm=IEN_btm, ndim=ndim, nDof=nDof, FunctionClass=FunctionClass, ID=ID, Young=Float64(Young), ν=ν, cMat=cMat)

    K = assemble_system(mdl)                   # assemble the stiffness matrix

    q_d, q_n, C = setboundaryCond(NodeList, ne, ndim, FunctionClass, nDof)

    # transpose the constraint matrix
    C_t = transpose(C)

    # extract the free part of the stiffness matrix
    K_free = C_t*K*C

    b = q_n - K*q_d

    invK = inv(Matrix(K_free))

    # solve the system
    q_f = invK*C_t*b

    # assemble the solution 
    q = q_d + C*q_f

    q_new, IEN_new = rearrange(q, ne, ndim, IEN, FunctionClass)

    filePath = "/home/soshala/SMEAR-PhD/SMEAR/Data/sim_experiments/Lame/experiment_single_run/Results/vtkFiles/test"
    write_vtk(filePath ,"u",NodeList,IEN_new, ne, ndim, q)
# end
    #######################################
    # test basis functions
    #######################################
    # x0 = -1
    # x1 = 1
    # y0 = -1
    # y1 = 1
    # z0 = -1
    # z1 = 1
    # ne = 1

    # NodeList, IEN, ID, IEN_top, IEN_btm, BorderNodesList = meshgrid(x0,x1,y0,y1,z0,z1,ne,ndim,FunctionClass=FunctionClass)

    # iter = 1:size(IEN,2)
    # for i in iter
    #     coord = NodeList[:,IEN[i]]
    #     N, dN = basis_function(coord[1],coord[2],coord[3], FunctionClass)
    #     println(findall(x->x==1,N)==[i])
    # end
    #######################################

