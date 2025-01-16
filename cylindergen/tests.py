import numpy

def test_extraction_operators(C):
    for Ce in C:
        numpy.testing.assert_allclose(numpy.sum(Ce, axis=0), 1.0)