#! /usr/bin/env python3
#
# This is a Python 3 script based on the open source
# spline library Splipy (github.com/SINTEF/Splipy) and the open source
# Nutils (v7.x, http://www.nutils.org) finite element library.
# This script generates a Bezier extraction file for a cylinder.

import splipy, splipy.volume_factory, splipy.surface_factory, splipy.utils.nutils
from nutils import mesh, function, cli
import treelog
import numpy, h5py
import cylindergen.stokes as stokes, extraction
from pathlib import Path

_ = numpy.newaxis

def main(case = 'nurbs geometry', # alternative: 'nurbs geometry'
         r = 1,        # Radius
         h = 1,        # Height
         p = 2,        # Basis order
         nref = 1):    # Number of uniform refinements

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

    treelog.user(f'BSpline volume = {vol_BSpline}')
    treelog.user(f'NURBS volume = {vol_NURBS}')

    # Stokes test case
    if case == 'nurbs geometry':
        stokes.test_stokes(domain, geom, u_basis=domain.basis('lagrange', degree=2), p_basis=domain.basis('lagrange', degree=1))
    elif case == 'iga':
        stokes.test_stokes(domain, geom, u_basis=domain.basis('spline', degree=3, continuity=-2), p_basis=domain.basis('spline', degree=2, continuity=-1))

    # Interior lagrange extraction
    IEN, C = extraction.get_lagrange_extraction(extraction_topo=domain, topo=domain, geom=ξ, basis=bspline_basis, degree=p)

    # Boundary lagrange extraction
    boundaries = {'front': 'front', 'back': 'back', 'sides': 'left, right, top, bottom'}
    boundaries_IEN = {}
    boundaries_C   = {}

    for boundary_name, boundary_tags in boundaries.items():
        boundary = domain.boundary[boundary_tags]
        boundaries_IEN[boundary_name], boundaries_C[boundary_name] = extraction.get_lagrange_extraction(extraction_topo=boundary, topo=domain, geom=ξ, basis=bspline_basis, degree=p)
        
    # Save to an HDF5 file

    elem_renumbering = numpy.swapaxes(numpy.arange(IEN.shape[0]).reshape(domain.shape),0,2).ravel()
    #TODO: RENUMBERING FOR BOUNDARIES?
    script_dir = Path( __file__ ).parent.absolute()
    with h5py.File((script_dir / 'cylinder').with_suffix('.h5'), 'w') as f:
        f.create_dataset('X', data=X)
        f.create_dataset('W', data=W)
        f.create_dataset('C', data=C[elem_renumbering,:,:])
        f.create_dataset('IEN', data=IEN[elem_renumbering,:])
        for boundary_name, boundary_IEN in boundaries_IEN.items():
            f.create_dataset(('IEN_'+str(boundary_name)), data=boundary_IEN)
        for boundary_name, boundary_C in boundaries_C.items():
            f.create_dataset(('C_'+str(boundary_name)), data=boundary_C)
        f.create_dataset('BSpline_vol', data=vol_BSpline)
        f.create_dataset('NURBS_vol', data=vol_NURBS)

if __name__ == '__main__':
    cli.run(main)