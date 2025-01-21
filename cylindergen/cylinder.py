#! /usr/bin/env python3
#
# This is a Python 3 script based on the open source
# spline library Splipy (github.com/SINTEF/Splipy) and the open source
# Nutils (v7.x, http://www.nutils.org) finite element library.
# This script generates a Bezier extraction file for a cylinder.

import splipy, splipy.volume_factory, splipy.surface_factory, splipy.utils.nutils
from nutils import mesh, function, export
import numpy, h5py
import tests

_ = numpy.newaxis

# problem parameters
r = 1 # Radius
h = 1 # Height

# discretization parameters
p = 2 # Basis order
nref = 1 # Number of uniform refinements

# creation of the Splipy NURBS object
disk = splipy.surface_factory.disc(r, type='square')
cylinder = splipy.volume_factory.extrude(disk, [0.,0.,h])
cylinder.raise_order(p-2,p-2,p-1)
cylinder.refine(2**nref-1)

# extaction of the control points and weights
Xh = numpy.reshape(cylinder.controlpoints[:,:,:,[0,1,2]].swapaxes(0,2), (-1,3), order='F')
W  = numpy.reshape(cylinder.controlpoints[:,:,:,3].swapaxes(0,2), (-1,), order='F')
X  = Xh / W[:,numpy.newaxis]

knotvalues         = cylinder.knots()
knotmultiplicities = splipy.utils.nutils.multiplicities(cylinder)
degree             = splipy.utils.nutils.degree(cylinder)

# define the domain with coordinate ξ
domain, ξ = mesh.rectilinear(knotvalues)
print(domain.shape)
# B-spline basis
bspline_basis = domain.basis('spline', degree=degree, knotvalues=knotvalues, knotmultiplicities=knotmultiplicities)
w             = bspline_basis.dot(W)
nurbs_basis   = (W*bspline_basis)/w

# B-spline geometry interpolation
geom = function.matmat(nurbs_basis, X)

vol_BSpline = domain.integrate(function.J(function.matmat(bspline_basis, X)), degree=2*p)
vol_NURBS = domain.integrate(function.J(geom), degree=2*p)
print(f'BSpline volume = {vol_BSpline}')
print(f'NURBS volume = {vol_NURBS}')

cp_domain, cp_geom = mesh.rectilinear([numpy.linspace(0,1,n) for n in cylinder.controlpoints.shape[:-1]])
cp_basis = cp_domain.basis("std",1)
cp_geom = function.matmat(cp_basis, X)

# plot the CAD object
bezier = domain.sample('bezier', 3)
x = bezier.eval(geom)
export.vtk('/home/soshala/SMEAR-PhD/smear-modules/smearFEM.jl/cylindergen/vtkFiles/cylinder', bezier.tri, x)

# plot the CAD object
cp_bezier = cp_domain.sample('bezier', 3)
x_cp = cp_bezier.eval(cp_geom)
export.vtk('/home/soshala/SMEAR-PhD/smear-modules/smearFEM.jl/cylindergen/vtkFiles/cylinder_cp', cp_bezier.tri, x_cp)

map_cp = [0, 4, 6, 2, 1, 5, 7, 3]
IEN_cp = numpy.empty(shape=(numpy.prod(cp_domain.shape),(2)**3),dtype=int)
for e, ref  in enumerate(cp_domain.references):
    IEN_e = cp_basis.get_dofs(e)
    IEN_cp[e,:] = IEN_e[map_cp]

# get the Lagrange extraction operators
lagrange_basis = domain.basis('lagrange', degree=p)
Afun = lagrange_basis[:,_]*lagrange_basis[_,:]
bfun = lagrange_basis[:,_]*bspline_basis[_,:]
A, b = domain.integrate_elementwise([Afun*function.J(geom), bfun*function.J(geom)], degree=2*p)

ne = A.shape[0]
assert ne==2**(3*nref), 'Incorrect number of elements'

# map node order of nutils to fit with smearFEM
map = [0, 18, 24, 6, 2, 20, 26, 8, 9, 21, 15, 3, 11, 23, 17, 5, 1, 19, 25, 7, 10, 22, 16, 4, 12, 14, 13]

IEN = numpy.empty(shape=(ne,(p+1)**3),dtype=int)
C = numpy.empty(shape=(ne,(p+1)**3,(p+1)**3))
for e, (Ae,be) in enumerate(zip(A,b)):
    IEN[e,:] = bspline_basis.get_dofs(e)
    Ce = numpy.transpose(numpy.linalg.inv(Ae[numpy.ix_(lagrange_basis.get_dofs(e),lagrange_basis.get_dofs(e))]).dot(be[numpy.ix_(lagrange_basis.get_dofs(e),IEN[e,:])]))
    C[e,:,:] = Ce[:,map]

tests.test_extraction_operators(C)

boundaries = ['top', 'bottom', 'front', 'back', 'left' , 'right']
IEN_bound_lst = {}
C_bound_lst = {}

for bound in boundaries:
    # Boundary mesh (1 boundary at a time)
    boundary = domain.boundary[bound]
    lagrange_basis = boundary.basis('lagrange', degree=p)
    Afun = lagrange_basis[:,_]*lagrange_basis[_,:]
    bfun = lagrange_basis[:,_]*bspline_basis[_,:]
    A, b = boundary.integrate_elementwise([Afun*function.J(geom), bfun*function.J(geom)], degree=2*p)

    ne = A.shape[0]
    assert ne==2**(2*nref), 'Incorrect number of elements'

    map_boundary = [0, 6, 8, 2, 3, 7, 5, 1, 4]

    IEN_boundary = numpy.empty(shape=(ne,(p+1)**2),dtype=int)
    C_boundary = numpy.empty(shape=(ne,(p+1)**2,(p+1)**2))
    for e, (Ae,be) in enumerate(zip(A,b)):

        # Only consider basis functions that are supported on the boundary
        supp = (numpy.sum(be,axis=0)>1e-10)
        assert sum(supp)==(p+1)**2

        i, tail = domain.transforms.index_with_tail(boundary.transforms[e])
        print(f'boundary element {e} corresponds to volume element {i}')

        IEN_boundary[e,:] = [dof for dof in bspline_basis.get_dofs(i) if supp[dof]]
        Ce_boundary = numpy.transpose(numpy.linalg.inv(Ae[numpy.ix_(lagrange_basis.get_dofs(e),lagrange_basis.get_dofs(e))]).dot(be[numpy.ix_(lagrange_basis.get_dofs(e),IEN_boundary[e,:])]))
        C_boundary[e,:,:] = Ce_boundary[:,map_boundary]
         
    tests.test_extraction_operators(C_boundary)

    IEN_bound_lst[bound] = IEN_boundary.copy()
    C_bound_lst[bound] = C_boundary.copy()

# Save to an HDF5 file
with h5py.File('/home/soshala/SMEAR-PhD/smear-modules/smearFEM.jl/cylindergen/cylinder.h5', 'w') as f:
    f.create_dataset('X', data=X)
    f.create_dataset('W', data=W)
    f.create_dataset('C', data=C)
    f.create_dataset('IEN', data=IEN)
    f.create_dataset('IEN_cp', data=IEN_cp)
    for key in IEN_bound_lst:
        f.create_dataset(('IEN_'+str(key)), data=IEN_bound_lst[key])
    for key in C_bound_lst:
        f.create_dataset(('C_'+str(key)), data=C_bound_lst[key])
    f.create_dataset('BSpline_vol', data=vol_BSpline)
    f.create_dataset('NURBS_vol', data=vol_NURBS)