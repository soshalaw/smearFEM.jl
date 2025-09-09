#! /usr/bin/env python3
#
# This is a Python 3 script based on the open source
# spline library Splipy (github.com/SINTEF/Splipy) and the open source
# Nutils (v7.x, http://www.nutils.org) finite element library.
# This script generates a Bezier extraction file for a cylinder.

import splipy, splipy.volume_factory, splipy.surface_factory, splipy.utils.nutils
from nutils import mesh, function, cli, export
import treelog
import numpy, h5py
import stokes as stokes, extraction
from pathlib import Path

_ = numpy.newaxis

def main(case='nurbs geometry',
         r=25,
         h=50,
         p=2,
         nref=3):

    # Creation of the Splipy NURBS object
    disk = splipy.surface_factory.disc(r, type='square')
    cylinder = splipy.volume_factory.extrude(disk, [0.,0.,h])
    cylinder.raise_order(p-2,p-2,p-1)
    cylinder.refine(2**nref-1)

    # Extraction of the control points and weights
    Xh = numpy.reshape(cylinder.controlpoints[:,:,:,[0,1,2]].swapaxes(0,2), (-1,3), order='F')
    W  = numpy.reshape(cylinder.controlpoints[:,:,:,3].swapaxes(0,2), (-1,), order='F')
    X  = Xh / W[:,numpy.newaxis]

    knotvalues         = cylinder.knots()
    knotmultiplicities = splipy.utils.nutils.multiplicities(cylinder)
    degree             = splipy.utils.nutils.degree(cylinder)

    # Define the domain with coordinate ξ
    domain, ξ = mesh.rectilinear(knotvalues)

    # B-spline basis
    bspline_basis = domain.basis('spline', degree=degree, knotvalues=knotvalues, knotmultiplicities=knotmultiplicities)
    w             = bspline_basis.dot(W)
    nurbs_basis   = (W*bspline_basis)/w

    # NURBS geometry interpolation
    geom = function.matmat(nurbs_basis, X)

    vol_BSpline = domain.integrate(function.J(function.matmat(bspline_basis, X)), degree=2*p)
    vol_NURBS = domain.integrate(function.J(geom), degree=2*p)

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

    treelog.user(f'BSpline volume = {vol_BSpline}')
    treelog.user(f'NURBS volume = {vol_NURBS}')

    # Stokes test case
    if case == 'nurbs geometry':
        x, lhsu = stokes.test_stokes(domain, geom, u_basis=domain.basis('lagrange', degree=2), p_basis=domain.basis('lagrange', degree=1))
    elif case == 'iga':
        x, lhsu = stokes.test_stokes(domain, geom, u_basis=domain.basis('spline', degree=3, continuity=-2), p_basis=domain.basis('spline', degree=2, continuity=-1))

    # Interior lagrange extraction
    IEN, C = extraction.get_lagrange_extraction(extraction_topo=domain, topo=domain, geom=ξ, basis=bspline_basis, degree=p)
    
    # Boundary lagrange extraction
    boundaries = {'front': 'front', 'back': 'back', 'sides': 'left, right, top, bottom'}
    boundaries_IEN = {}
    boundaries_C   = {}

    for boundary_name, boundary_tags in boundaries.items():
        boundary = domain.boundary[boundary_tags]
        bIEN, bC = extraction.get_lagrange_extraction(extraction_topo=boundary, topo=domain, geom=ξ, basis=bspline_basis, degree=p)
        bnd_renumbering = numpy.swapaxes(numpy.arange(bIEN.shape[0]).reshape(boundary.shape),0,1).ravel()
        boundaries_IEN[boundary_name] = bIEN[bnd_renumbering,:] 
        boundaries_C[boundary_name] = bC[bnd_renumbering,:,:]

    print(numpy.shape(x))
    print(numpy.shape(lhsu))
    # Save to an HDF5 file
    elem_renumbering = numpy.swapaxes(numpy.arange(IEN.shape[0]).reshape(domain.shape),0,2).ravel()
    print("IEN shape :",numpy.shape(IEN))
    script_dir = Path( __file__ ).parent.absolute()

    save_dir = Path(script_dir,'slip_1')
    save_dir.mkdir(parents=True, exist_ok=True)
    
    numpy.savetxt(str(Path(save_dir,'node_list.csv')),x.astype(numpy.float64),delimiter=",")
    numpy.savetxt(str(Path(save_dir,'sol_u.csv')),lhsu.astype(numpy.float64),delimiter=",")

    with h5py.File(Path(save_dir.parent, ('cylinder_'+str(nref))).with_suffix('.h5'), 'w') as f:
        f.create_dataset('X', data=X)
        f.create_dataset('W', data=W)
        f.create_dataset('C', data=C[elem_renumbering,:,:])
        f.create_dataset('IEN', data=IEN[elem_renumbering,:])
        f.create_dataset('IEN_cp', data=IEN_cp)
        for boundary_name, boundary_IEN in boundaries_IEN.items():
            f.create_dataset(('IEN_'+str(boundary_name)), data=boundary_IEN)
        for boundary_name, boundary_C in boundaries_C.items():
            f.create_dataset(('C_'+str(boundary_name)), data=boundary_C)
        f.create_dataset('BSpline_vol', data=vol_BSpline)
        f.create_dataset('NURBS_vol', data=vol_NURBS)

if __name__ == '__main__':
    # cli.run( main )
    main()