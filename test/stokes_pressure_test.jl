using smearFEM
using LinearAlgebra

r = 0.5
h = 1
ne = 5
ndim = 3
FunctionClass_u = "Q2"
FunctionClass_p = "Q1"
nDof_u = ndim  # number of degree of freedom per node
nDof_p = 1
β = 1e-5

h = z1 - z0
μu_btm = 0  
μu_tp = -0.1
μu_side = 0

μp_btm = 0
μp_tp = 0
μp_side = 0

NodeList_p, IEN_p, ID_p, IEN_p_top, IEN_p_btm, BorderNodesList_p = meshgrid_cylinder(r, h, ne, FunctionClass=FunctionClass)  # generate the mesh grid

p = zeros(size(NodeListCylinderp, 2))

iter = 1:size(NodeListCylinderp, 2)
for i in iter
    R = 0.5
    r_ = NodeListCylinderp[1,i]^2 + NodeListCylinderp[2,i]^2
    # h_ = NodeListCylinderp[3,i] - z0
    h_ = 0.5
    p[i] = 3*μu_tp*(r_ - R^2)
end

filePath = "/home/soshala/SMEAR-PhD/SMEAR/Data/sim_experiments/stokes_test/vtkFiles"

write_vtk(filePath, "p_an", NodeListCylinderp, IEN_p, ne, ndim, p)