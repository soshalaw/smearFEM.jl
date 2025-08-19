import unittest, numpy, itertools
from nutils import mesh, function

_ = numpy.newaxis

maps = {
    (2,2): [0, 6, 8, 2, 3, 7, 5, 1, 4],
    (2,3): [0, 18, 24, 6, 2, 20, 26, 8, 9, 21, 15, 3, 11, 23, 17, 5, 1, 19, 25, 7, 10, 22, 16, 4, 12, 14, 13]
}

invert_map = lambda map: [map.index(i) for i in map]

def get_lagrange_extraction(extraction_topo, topo, geom, basis, degree):

    # get the Lagrange extraction operators
    lagrange_basis = extraction_topo.basis('lagrange', degree=degree)
    Afun = lagrange_basis[:,_]*lagrange_basis[_,:]
    bfun = lagrange_basis[:,_]*basis[_,:]
    A, b = extraction_topo.integrate_elementwise([Afun*function.J(geom), bfun*function.J(geom)], degree=2*degree)

    # map node order of nutils to fit with smearFEM
    if not (degree, extraction_topo.ndims) in maps:
        raise NotImplementedError(f'Lagrange extraction not implemented for ndims={extraction_topo.ndims} and p={degree}')
    else:
        map = maps[(degree, extraction_topo.ndims)]

    ne  = A.shape[0]
    IEN = numpy.empty(shape=(ne,(degree+1)**extraction_topo.ndims),dtype=int)
    C   = numpy.empty(shape=(ne,(degree+1)**extraction_topo.ndims,(degree+1)**extraction_topo.ndims))

    for e_extract, (Ae,be) in enumerate(zip(A,b)):

        supp = (numpy.sum(be,axis=0)>1e-10)
        assert sum(supp)==(degree+1)**extraction_topo.ndims, 'Incorrect number of supported basis functions'

        e_bspline, tail = topo.transforms.index_with_tail(extraction_topo.transforms[e_extract])
        assert len(tail)==topo.ndims-extraction_topo.ndims, 'Tail length mismatch'

        IEN[e_extract,:] = [dof for dof in basis.get_dofs(e_bspline) if supp[dof]]
        Ce = numpy.transpose(numpy.linalg.inv(Ae[numpy.ix_(lagrange_basis.get_dofs(e_extract),lagrange_basis.get_dofs(e_extract))]).dot(be[numpy.ix_(lagrange_basis.get_dofs(e_extract),IEN[e_extract,:])]))
        C[e_extract,:,:] = Ce[:,map]

    return IEN, C

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
