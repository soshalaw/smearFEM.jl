# using WriteVTK

# R = 3.0

# points = permutedims(
#     [
#       1.     1.     1.   
#      -1.     1.     1.   
#      -1.     1.    -1.   
#       1.     1.    -1.   
#       1.    -1.     1.   
#      -1.    -1.     1.   
#      -1.    -1.    -1.   
#       1.    -1.    -1.   
#       0.     1.     1.   
#      -1.     1.     0.   
#       0.     1.    -1.   
#       1.     1.     0.   
#       0.    -1.     1.   
#      -1.    -1.     0.   
#       0.    -1.    -1.   
#       1.    -1.     0.   
#       1.     0.     1.   
#      -1.     0.     1.   
#      -1.     0.    -1.   
#       1.     0.    -1.   
#       0.     0.     0.   
#       0.     1.     0.   
#       0.    -1.     0.   
#       0.     0.     1.   
#      -1.     0.     0.   
#       0.     0.    -1.   
#       1.     0.     0.   
#     ]
# )

# data = zeros(size(points,2))

# for pts = 1:size(points,2)
#   if sum(points[:,pts].^2) > 1e-10
#     points[:,pts] = R*points[:,pts]/sqrt(sum(points[:,pts].^2))
#     data[pts] = sqrt(sum(points[:,pts].^2))
#   end
# end
# # My Element Node Numbering
# elem = [ 1, 2, 3, 4,        # Vertex   Y =  1
#          5, 6, 7, 8,        # Vertex   Y = -1
#          9,10,11,12,        # MidEdge  Y =  1
#         13,14,15,16,        # MidEdge  Y = -1
#         17,18,19,20,        # MidEdge  Y =  0
#         21,                 # Body     X=Y=Z=0
#         22,23,24,25,26,27]  # Face Centroids

# # The following shows the mapping of ParaView's node numbering
# # generated via Sources to my node numbering:
# #   VTK NID          My NID
# #  0  1  3  2  =>  7  8  4  3    # Bottom   Z = -1
# #  4  5  7  6  =>  6  5  1  2    # Top      Z =  1
# #  8 11 12  9  => 15 20 11 19    # MidEdge  Z = -1
# # 22 25 26 23  => 13 17  9 18    # MidEdge  Z =  1
# # 13 15 21 19  => 14 16 12 10    # MidEdge  Z =  0
# # 16           => 25             # Face     X = -1
# # 18           => 27             # Face     X =  1
# # 14           => 23             # Face     Y = -1
# # 20           => 22             # Face     Y =  1
# # 10           => 26             # Face     Z = -1
# # 24           => 24             # Face     Z =  1
# # 17           => 21             # Body    X=Y=Z = 0

# # From what I've output from ParaView, this is what 
# # I THINK the permutation for VTK should be:
# perm27 = [ 7,  8,  4,  3,
#            6,  5,  1,  2,
#           15, 20, 11, 19,  
#           13, 17,  9, 18,
#           14, 16, 10, 12,
#           25, 27, 23, 22, 26, 24, 
#           21]

# # Define cells for WriteVTK (imitating my larger program's logic)
# cells = Vector{MeshCell}(undef,1)
# for eid = 1:1
#     cells[eid] = MeshCell(VTKCellTypes.VTK_LAGRANGE_HEXAHEDRON, elem[perm27,eid])
# end

# # vtk_model = vtk_grid("HEX27_ex", points, cells)
#  vtk_grid("HEX27_ex", points, cells) do vtk
#     vtk["data"] = data
#  end

using smearFEM
using LinearAlgebra
using SparseArrays
using WriteVTK

x0 = 0
x1 = 1
y0 = 0
y1 = 1
z0 = 0
z1 = 1
ne = 5
ndim = 3
FunctionClass_u = "Q2"
FunctionClass_p = "Q1"
nDof_u = ndim  # number of degree of freedom per node
nDof_p = 1
β = 0.0001
Young = 40
ν = 0.35

μu_btm = 0  
μu_tp = -1

μp_btm = 0
μp_tp = 0

NodeList_u, IEN_u, ID_u, IEN_u_top, IEN_u_btm, BorderNodesList = meshgrid_cube(x0,x1,y0,y1,z0,z1,ne,ndim,FunctionClass=FunctionClass_u)  # generate the mesh grid
NodeListCylinder = inflate_cylinder(NodeList_u, x0, x1, y0, y1)
q_tp, q_btm, C_uc = set_boundary_cond_stokes(NodeListCylinder, ne, ndim, FunctionClass_u, nDof_u)

q = (μu_btm*q_btm + μu_tp*q_tp)      # apply the Dirichlet boundary conditions

q_ = [q[ID_u[:,1]] q[ID_u[:,2]] q[ID_u[:,3]]]'

q_new, IEN_new = rearrange(q, ne, ndim, IEN_u, FunctionClass_u, ID_u) # rearrange the solution

filePath = "/home/soshala/SMEAR-PhD/SMEAR/Data/sim_experiments/stokes_test/vtkFiles"

write_vtk(filePath, "u", NodeListCylinder, IEN_new, ne, ndim, q_, FunctionClass=FunctionClass_u)

