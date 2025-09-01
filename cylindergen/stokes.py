import numpy
from nutils import function, solver, export
from nutils.expression_v2 import Namespace

def test_stokes(domain, geom, u_basis, p_basis):

    ns = Namespace()
    
    ns.ν = 1.0 # viscosity
    ns.β = 1e4 # slip

    ns.x = geom
    ns.define_for('x', gradient='∇', normal='n', jacobians=('dV', 'dS'))
    
    ns.Nu = u_basis.vector(domain.ndims)
    ns.Np = p_basis
    ns.u = function.dotarg('lhs_u', ns.Nu)
    ns.p = function.dotarg('lhs_p', ns.Np)

    # Weak formulation for Stokes flow
    A = domain.integrate('2 ν (∇_j(Nu_mi) + ∇_i(Nu_mj)) (∇_j(Nu_ni) + ∇_i(Nu_nj)) dV' @ ns, degree=4).export('dense')
    B = domain.integrate('-Np_m ∇_i(Nu_ni) dV' @ ns, degree=4).export('dense')

    top_boundary = domain.boundary['back']
    bot_boundary = domain.boundary['front']
    slp_boundary = domain.boundary['front,back']

    Stop, Qtop = top_boundary.integrate(['dS','x_2 dS'] @ ns, degree=4)
    Sbot, Qbot = bot_boundary.integrate(['dS','x_2 dS'] @ ns, degree=4)
    
    print(f'Top boundary area = {Stop}, height = {Qtop/Stop}')
    print(f'Bottom boundary area = {Sbot}, height = {Qbot/Sbot}')

    # Boundary condition matrix
    A_bcs = slp_boundary.integrate('β Nu_mi Nu_ni dS' @ ns, degree=4).export('dense')

    top_sqr = top_boundary.integral('( (u_2 + 0.1)^2 ) dS' @ ns, degree=4)
    top_cons = solver.optimize('lhs_u', top_sqr, droptol=1e-9)

    bot_sqr = bot_boundary.integral('( u_2^2 ) dS' @ ns, degree=4)
    bot_cons = solver.optimize('lhs_u', bot_sqr, droptol=1e-9)

    # Constraint matrix
    cons = numpy.full_like(bot_cons, numpy.nan)
    cons[~numpy.isnan(bot_cons)] = bot_cons[~numpy.isnan(bot_cons)]
    cons[~numpy.isnan(top_cons)] = top_cons[~numpy.isnan(top_cons)]

    C = numpy.eye(len(cons))[:, numpy.isnan(cons)]
    up = cons.copy()
    up[numpy.isnan(up)] = 0.
    
    # Solving
    print(numpy.min(A), numpy.max(A))

    A_total = A  + A_bcs
    Ac = C.T @ A_total @ C
    fc = -C.T @ A_total @ up

    Bc = B @ C
    gc = -B @ up
    
    M = numpy.block([[Ac, Bc.T],
                     [Bc, numpy.zeros((Bc.shape[0], Bc.shape[0]))]])
    b = numpy.concatenate([fc, gc])

    # Solve the augmented system
    lhs = numpy.linalg.solve(M, b)
    
    # Extract solutions
    lhsuc = lhs[:Ac.shape[0]]
    lhsp = lhs[Ac.shape[0]:]
    
    # Reconstruct full velocity solution
    lhsu = C @ lhsuc + up

    # Plot
    bezier = domain.sample('bezier', 3)
    x, pvals, uvals = bezier.eval(['x_i', 'p', 'u_i'] @ ns, arguments={'lhs_u': lhsu, 'lhs_p': lhsp})
    export.vtk('stokes', bezier.tri, x, u=uvals, p=pvals)