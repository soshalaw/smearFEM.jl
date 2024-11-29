using smearFEM
using LinearAlgebra

x0 = -0.5
x1 = 0.5
y0 = -0.5
y1 = 0.5
z0 = -0.5
z1 = 0.5
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

NodeList_p, IEN_p, ID_p, IEN_p_top, IEN_p_btm, BorderNodesList_p = meshgrid_cube(x0,x1,y0,y1,z0,z1,ne,ndim,FunctionClass=FunctionClass_p)  # generate the mesh grid
NodeListCylinderp = inflate_cylinder(NodeList_p, x0, x1, y0, y1)

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