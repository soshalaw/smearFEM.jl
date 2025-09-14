import unittest, numpy, itertools
from nutils import mesh, function

_ = numpy.newaxis

maps = {
    (2,2): [0, 6, 8, 2, 3, 7, 5, 1, 4],
    (2,3): [0, 18, 24, 6, 2, 20, 26, 8, 9, 21, 15, 3, 11, 23, 17, 5, 1, 19, 25, 7, 10, 22, 16, 4, 12, 14, 13],
    (3,3): [0, 48, 51, 3, 12, 60, 63, 15, 1, 49, 50, 2, 13, 61, 62, 14, 16, 32, 35, 19, 28, 44, 47, 31, 4, 52, 55, 7, 8, 56, 59, 11, 
            5, 53, 54, 6, 9, 57, 58, 10, 20, 36, 39, 23, 24, 40, 43, 27, 17, 33, 34, 18, 29, 45, 46, 30, 21, 37, 38, 22, 25, 41, 42, 25]
}

invert_map = lambda map: [map.index(i) for i in map]

def get_lagrange_extraction(extraction_topo, topo, geom, basis, degree, increase_degree_by_vis=-1):

    # get the Lagrange extraction operators
    new_degree = (degree)
    ndims = extraction_topo.ndims
    n_elems = len(extraction_topo)
    n_lagrange = (new_degree + 1) ** ndims

    lagrange_basis = extraction_topo.basis('lagrange', degree=new_degree)
    
    A = numpy.empty(shape=(n_elems,n_lagrange,n_lagrange))
    b = numpy.empty(shape=(n_elems,n_lagrange,basis.ndofs))
    for e, idx in enumerate(numpy.ndindex(extraction_topo.shape)):

        edofs = lagrange_basis.get_dofs(e)
        Afun = lagrange_basis[edofs,_]*lagrange_basis[_,edofs]
        bfun = lagrange_basis[edofs,_]*basis[_,:]

        if ndims==3:
            elem = extraction_topo[idx[0]:idx[0]+1, idx[1]:idx[1]+1, idx[2]:idx[2]+1]
        elif ndims==2:
            elem = extraction_topo[idx[0]:idx[0]+1, idx[1]:idx[1]+1]
        else:
            raise NotImplementedError('Extraction only implemented for 2D and 3D volumes')
        
        quad_degree = 2 * new_degree

        Ae, be = elem.integrate([Afun*function.J(geom), bfun*function.J(geom)], degree=quad_degree)

        A[e,:,:] = Ae.export('dense')
        b[e,:,:] = be.export('dense')

    # map node order of nutils to fit with smearFEM+
    if not (new_degree, ndims) in maps:
        raise NotImplementedError(f'Lagrange extraction not implemented for ndims={ndims} and p={new_degree}')
    else:
        map = maps[(new_degree, ndims)]

    ne  = A.shape[0]
    IEN = numpy.empty(shape=(ne,n_lagrange),dtype=int)
    C   = numpy.empty(shape=(ne,n_lagrange,n_lagrange))

    for e_extract, (Ae,be) in enumerate(zip(A,b)):

        supp = (numpy.sum(be,axis=0)>1e-10)
        assert sum(supp)==n_lagrange, 'Incorrect number of supported basis functions'

        e_bspline, tail = topo.transforms.index_with_tail(extraction_topo.transforms[e_extract])

        assert len(tail) == (topo.ndims - ndims), f'Tail length mismatch: {len(tail)} != {topo.ndims - ndims}'
        
        IEN[e_extract,:] = [dof for dof in basis.get_dofs(e_bspline) if supp[dof]]
        Ce = numpy.linalg.solve(Ae, be[:,IEN[e_extract,:]]).T
        C[e_extract,:,:] = Ce[:,map]

 # --- VISUALIZATION BASIS ---
    if increase_degree_by_vis >= 0:
            # get the Lagrange extraction operators
        new_degree = degree + increase_degree_by_vis
        n_lagrange_vis = (new_degree + 1) ** ndims

        lagrange_basis_vis = extraction_topo.basis('lagrange', degree=new_degree)
        
        A = numpy.empty(shape=(n_elems,n_lagrange_vis,n_lagrange_vis))
        b = numpy.empty(shape=(n_elems,n_lagrange_vis,basis.ndofs))
        for e, idx in enumerate(numpy.ndindex(extraction_topo.shape)):

            edofs = lagrange_basis_vis.get_dofs(e)
            Afun = lagrange_basis_vis[edofs,_]*lagrange_basis_vis[_,edofs]
            bfun = lagrange_basis_vis[edofs,_]*basis[_,:]

            if ndims==3:
                elem = extraction_topo[idx[0]:idx[0]+1, idx[1]:idx[1]+1, idx[2]:idx[2]+1]
            elif ndims==2:
                elem = extraction_topo[idx[0]:idx[0]+1, idx[1]:idx[1]+1]
            else:
                raise NotImplementedError('Extraction only implemented for 2D and 3D volumes')
            
            quad_degree = 2 * new_degree

            Ae, be = elem.integrate([Afun*function.J(geom), bfun*function.J(geom)], degree=quad_degree)

            A[e,:,:] = Ae.export('dense')
            b[e,:,:] = be.export('dense')

        # map node order of nutils to fit with smearFEM+
        if not (new_degree, ndims) in maps:
            raise NotImplementedError(f'Lagrange extraction not implemented for ndims={ndims} and p={new_degree}')
        else:
            vis_map = maps[(new_degree, ndims)]

        ne  = A.shape[0]
        n_nurbs = len(basis.get_dofs(0))  # assume regular mesh
        IEN_vis = numpy.empty(shape=(ne,n_nurbs),dtype=int)
        C_vis   = numpy.empty(shape=(ne,n_nurbs,n_lagrange_vis))

        print(numpy.shape(C_vis))

        for e_extract, (Ae,be) in enumerate(zip(A,b)):

            supp_vis = (numpy.sum(be,axis=0)>1e-10)
            assert sum(supp)==n_nurbs, f'Incorrect number of supported basis functions sum(supp_vis)={sum(supp_vis)}, n_nurbs={n_nurbs}'

            e_bspline, tail = topo.transforms.index_with_tail(extraction_topo.transforms[e_extract])
            assert len(tail) == (topo.ndims - ndims), f'Tail length mismatch: {len(tail)} != {topo.ndims - ndims}'
            
            IEN_vis[e_extract,:] = [dof for dof in basis.get_dofs(e_bspline) if supp_vis[dof]]
            Ce = numpy.linalg.solve(Ae, be[:,IEN[e_extract,:]]).T

            print(numpy.shape(Ce))
            C_vis[e_extract,:,:] = Ce[:,vis_map]
        
        return IEN, C, IEN_vis, C_vis
        # If no visualization degree increase requested, return only analysis extraction
    return IEN, C


# def project_control_to_vis(extraction_topo, geom, nurbs_basis, nurbs_degree, vis_degree, int_degree=None):
#     """
#     Project from a coarse NURBS control net (27 DOFs/element) to a visualization mesh (64 DOFs/element).

#     Parameters:
#         extraction_topo: Nutils topology for mesh (same as NURBS control mesh)
#         geom: geometry function
#         nurbs_basis: Nutils global basis (degree 2 for NURBS/B-spline example)
#         vis_degree: visualization Lagrange degree (e.g. 3 for 64/element)
#         int_degree: (optional) integration degree, otherwise 2 * max(degrees)

#     Returns:
#         IEN_vis: Per-element (27,) global control net indices (for each visualization element)
#         C_vis:   Per-element (64, 27) projection matrices; multiply C_vis[e] @ u[IEN_vis[e]]
#     """
#     ndims = extraction_topo.ndims
#     n_elems = len(extraction_topo)
#     n_lagrange = (vis_degree + 1) ** ndims

#     # Determine integration degree
#     if int_degree is None:
#         int_degree = 2 * max(nurbs_degree, vis_degree)

#     # Visualization Lagrange basis on each extraction element
#     vis_lagrange_basis = extraction_topo.basis('lagrange', degree=vis_degree)

#     if (vis_degree, ndims) not in maps:
#         raise NotImplementedError(f'Mapping for degree {vis_degree} and {ndims}D not implemented.')
#     node_map_vis = maps[(vis_degree, ndims)]

#     # Allocate outputs
#     n_nurbs = len(nurbs_basis.get_dofs(0))  # assume regular mesh
#     IEN_vis = numpy.empty((n_elems, n_nurbs), dtype=int)
#     C_vis = numpy.empty((n_elems, n_lagrange, n_nurbs))

#     # For each visualization element (same as control net element)
#     for e, idx in enumerate(numpy.ndindex(extraction_topo.shape)):
#         # Identify the local control net DOFs (27 for degree-2, 3D mesh)
#         dofs = nurbs_basis.get_dofs(e)
#         assert len(dofs) == n_nurbs, f"Expected {n_nurbs} DOFs but got {len(dofs)} in element {e}"

#         IEN_vis[e, :] = dofs

#         # Define symbolic integrands for Lagrange basis coupling with NURBS basis
#         edofs_vis = vis_lagrange_basis.get_dofs(e)
#         assert len(edofs_vis) == n_lagrange

#         element = extraction_topo[idx[0]:idx[0]+1, idx[1]:idx[1]+1, idx[2]:idx[2]+1] if ndims == 3 else \
#                   extraction_topo[idx[0]:idx[0]+1, idx[1]:idx[1]+1]
        
#         edofs = vis_lagrange_basis.get_dofs(e)
#         Afun = vis_lagrange_basis[edofs,_]*vis_lagrange_basis[_,edofs]
#         bfun = vis_lagrange_basis[edofs,_]*nurbs_basis[_,:]


#         # Integrate over just this element; for more accuracy, use a patch
#         Ae, be = element.integrate([Afun * function.J(geom), bfun * function.J(geom)], degree=int_degree)

#         print('Ae:', Ae.shape)  # should be (64,64)
#         print('be:', be.shape)  # should be (64,27)

#         A_vis = Ae.export('dense')        # shape (64, 64)
#         b_vis = be.export('dense')        # shape (64, 27)

#         # Solve for projection matrix: from NURBS to vis Lagrange (64, 27)
#         Ce_vis = numpy.linalg.solve(A_vis, b_vis)   # shape (64, 27)

#         C_vis[e, :, :] = Ce_vis[:,node_map_vis]

#     return IEN_vis, C_vis


###########
# Testing #
###########

julia_quad2 = numpy.array( [[0.0,0.0],
                           [1.0,0.0],
                           [1.0,1.0],
                           [0.0,1.0],
                           [0.5,0.0],
                           [1.0,0.5],
                           [0.5,1.0],
                           [0.0,0.5],
                           [0.5,0.5]])

julia_hex2 = numpy.array([[0.0,0.0,0.0],
                          [1.0,0.0,0.0],
                          [1.0,1.0,0.0],
                          [0.0,1.0,0.0],
                          [0.0,0.0,1.0],
                          [1.0,0.0,1.0],
                          [1.0,1.0,1.0],
                          [0.0,1.0,1.0],
                          [0.5,0.0,0.0],
                          [1.0,0.5,0.0],
                          [0.5,1.0,0.0],
                          [0.0,0.5,0.0],
                          [0.5,0.0,1.0],
                          [1.0,0.5,1.0],
                          [0.5,1.0,1.0],
                          [0.0,0.5,1.0],
                          [0.0,0.0,0.5],
                          [1.0,0.0,0.5],
                          [1.0,1.0,0.5],
                          [0.0,1.0,0.5],
                          [0.5,0.0,0.5],
                          [1.0,0.5,0.5],
                          [0.5,1.0,0.5],
                          [0.0,0.5,0.5],
                          [0.5,0.5,0.0],
                          [0.5,0.5,1.0],
                          [0.5,0.5,0.5]])

julia_top = [3, 2, 6, 7, 10, 18, 14, 19, 22]

class TestExtraction:

    def test_IEN_shape(self):
        self.assertEqual(self.IEN.shape, (self.nelems, self.ndofs))

    def test_extraction_shape(self):
        self.assertEqual(self.C.shape, (self.nelems, self.ndofs, self.ndofs))

    def test_column_sums(self):
        numpy.testing.assert_allclose(numpy.sum(self.C, axis=1), 1.0)

    def test_lagrange_points(self):
        for e in range(self.nelems):
            numpy.testing.assert_allclose((self.X[self.IEN[e,:],:].T @ self.C[e,:,:]).T, self.Y[e], atol=1e-10, rtol=1e-10)

class TestUnitSquareVolume(unittest.TestCase,TestExtraction):
    def setUp(self):
        self.ndims  = 2
        self.degree = 2
        self.topo, self.geom = mesh.rectilinear([2]*self.ndims)
        self.nelems = len(self.topo)
        self.basis = self.topo.basis('spline', degree=self.degree)
        self.ndofs = (self.degree+1)**self.ndims
        self.IEN, self.C = get_lagrange_extraction(self.topo, self.topo, self.geom, self.basis, degree=2)
        self.X = self.topo.project(self.geom, onto=self.basis.vector(self.ndims), geometry=self.geom, degree=2*self.degree).reshape(-1,self.ndims)
        self.Y = [julia_quad2 + numpy.array(offset) for offset in itertools.product(range(self.topo.shape[0]),repeat=self.ndims)]

class TestUnitCubeVolume(unittest.TestCase,TestExtraction):
    def setUp(self):
        self.ndims = 3
        self.degree = 2
        self.topo, self.geom = mesh.rectilinear([2]*self.ndims)
        self.nelems = len(self.topo)
        self.basis = self.topo.basis('spline', degree=self.degree)
        self.ndofs = (self.degree+1)**self.ndims        
        self.IEN, self.C = get_lagrange_extraction(self.topo, self.topo, self.geom, self.basis, degree=2)
        self.X = self.topo.project(self.geom, onto=self.basis.vector(self.ndims), geometry=self.geom, degree=2*self.degree).reshape(-1,self.ndims)
        self.Y = [julia_hex2 + numpy.array(offset) for offset in itertools.product(range(self.topo.shape[0]),repeat=self.ndims)]        

class TestUnitCubeTop(unittest.TestCase,TestExtraction):
    def setUp(self):
        self.ndims = 3
        self.degree = 2
        self.topo, self.geom = mesh.rectilinear([2]*self.ndims)
        self.nelems = len(self.topo.boundary['top'])
        self.basis = self.topo.basis('spline', degree=self.degree)
        self.ndofs = (self.degree+1)**self.topo.boundary['top'].ndims
        self.X = self.topo.boundary['top'].project(self.geom, onto=self.basis.vector(self.ndims), geometry=self.geom, degree=2*self.degree).reshape(-1,self.ndims)
        self.IEN, self.C = get_lagrange_extraction(self.topo.boundary['top'], self.topo, self.geom, self.basis, degree=2)
        self.Y =[julia_hex2[julia_top] + numpy.insert(numpy.array(offset), 1, 1) for offset in itertools.product(range(self.topo.shape[0]),repeat=self.ndims-1)]

if __name__ == '__main__':
    unittest.main(verbosity=2)


# Node order maps (unchanged from your original maps)
# maps = {
#     (2,2): [0, 6, 8, 2, 3, 7, 5, 1, 4],
#     (2,3): [0, 18, 24, 6, 2, 20, 26, 8, 9, 21, 15, 3, 11, 23, 17, 5, 1, 19, 25, 7, 10, 22, 16, 4, 12, 14, 13],
#     (3,3): [0, 48, 51, 3, 12, 60, 63, 15, 1, 49, 50, 2, 13, 61, 62, 14, 16, 32, 35, 19, 28, 44, 47, 31, 4, 52, 55, 7, 8, 56, 59, 11,
#             5, 53, 54, 6, 9, 57, 58, 10, 20, 36, 39, 23, 24, 40, 43, 27, 17, 33, 34, 18, 29, 45, 46, 30, 21, 37, 38, 22, 25, 41, 42, 25]
# }

# def get_lagrange_extraction(extraction_topo, topo, geom, basis, degree, oversample_factor=1):
#     """
#     Compute Lagrange extraction operator C and IEN array.
    
#     Parameters:
#     - extraction_topo : finer topology for extraction (over-sampled if oversample_factor > 1)
#     - topo : original topology (NURBS mesh)
#     - geom : geometry mapping
#     - basis : global NURBS/B-spline basis
#     - degree : polynomial degree of Lagrange basis
#     - oversample_factor : int > 0, factor to increase Lagrange mesh density
    
#     Returns:
#     - IEN : Element-to-global-DOF connectivity for Lagrange mesh
#     - C : Extraction operators per element mapping Lagrange -> NURBS
#     """

#     # Determine oversampled degree and topology
#     ndims = extraction_topo.ndims

#     # Adjust Lagrange mesh if oversampling is requested
#     if oversample_factor > 1:
#         # Refine extraction topology by subdividing uniformly oversample_factor times each dimension
#         refined_shape = tuple(oversample_factor * s for s in extraction_topo.shape)
#         refined_topo = extraction_topo.refined(refined_shape)
#         lagrange_basis = refined_topo.basis('lagrange', degree=degree)
#         use_topo = refined_topo
#     else:
#         # Use extraction_topo as is (no oversampling)
#         lagrange_basis = extraction_topo.basis('lagrange', degree=degree)
#         use_topo = extraction_topo

#     n_elems = len(use_topo)
#     n_lagrange = (degree + 1) ** ndims
#     n_nurbs = basis.ndofs

#     A = numpy.empty((n_elems, n_lagrange, n_lagrange))
#     b = numpy.empty((n_elems, n_lagrange, n_nurbs))

#     for e, idx in enumerate(numpy.ndindex(use_topo.shape)):
#         edofs = lagrange_basis.get_dofs(e)

#         Afun = Afun = lagrange_basis[edofs,_]*lagrange_basis[_,edofs]
#         bfun = lagrange_basis[edofs,_]*basis[_,:]

#         if ndims == 3:
#             elem = use_topo[idx[0]:idx[0]+1, idx[1]:idx[1]+1, idx[2]:idx[2]+1]
#         elif ndims == 2:
#             elem = use_topo[idx[0]:idx[0]+1, idx[1]:idx[1]+1]
#         else:
#             raise NotImplementedError('Extraction only implemented for 2D and 3D')

#         # Use higher integration degree for oversampled cases to ensure accuracy
#         integrate_degree = 2 * degree if oversample_factor == 1 else 3 * degree

#         Ae, be = elem.integrate([Afun * function.J(geom), bfun * function.J(geom)], degree=integrate_degree)

#         A[e, :, :] = Ae.export('dense')
#         b[e, :, :] = be.export('dense')

#     # Apply node reordering map if needed
#     if not (degree, ndims) in maps:
#         raise NotImplementedError(f'Lagrange extraction not implemented for ndims={ndims} and degree={degree}')
#     else:
#         node_map = maps[(degree, ndims)]

#     IEN = numpy.empty((n_elems, n_lagrange), dtype=int)
#     C = numpy.empty((n_elems, n_lagrange, n_nurbs))

#     for e_extract, (Ae, be) in enumerate(zip(A, b)):

#         supp = (numpy.sum(be, axis=0) > 1e-10)
#         assert sum(supp) == n_nurbs, 'Unexpected number of supported NURBS basis functions'

#         e_bspline, tail = topo.transforms.index_with_tail(use_topo.transforms[e_extract])
#         assert len(tail) == topo.ndims - ndims, 'Tail length mismatch'

#         IEN[e_extract, :] = [dof for dof in basis.get_dofs(e_bspline) if supp[dof]]

#         Ce = numpy.linalg.solve(Ae, be[:, IEN[e_extract, :]]).T
#         C[e_extract, :, :] = Ce[:, node_map]

#     return IEN, C