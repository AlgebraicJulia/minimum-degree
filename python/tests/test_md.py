# Orders the graph of the README with each interface of `md`.
import numpy as np
import scipy.sparse

from minimum_degree import check, md

XADJ = [0, 1, 4, 9, 10, 12, 15, 18, 20]
ADJNCY = [1, 0, 2, 5, 1, 3, 5, 6, 7, 2, 5, 6, 1, 2, 4, 2, 4, 7, 2, 6]
INVP = [1, 2, 7, 0, 4, 3, 6, 5]


def test_arrays():
    assert md(XADJ, ADJNCY).tolist() == INVP


def test_matrix():
    n = len(XADJ) - 1
    a = scipy.sparse.csr_matrix((np.ones(len(ADJNCY)), ADJNCY, XADJ), shape=(n, n))
    # the pattern of A + A^T: a triangle, or the matrix with its diagonal, is the same graph
    assert md(a).tolist() == INVP
    assert md(scipy.sparse.triu(a)).tolist() == INVP
    assert md(a + scipy.sparse.eye(n)).tolist() == INVP


def test_options():
    invp = md(XADJ, ADJNCY, external=False, mass=False, twin=False)
    assert sorted(invp.tolist()) == list(range(8))
    assert md(XADJ, ADJNCY, vwght=np.ones(8, dtype=int)).tolist() == INVP


def test_empty():
    assert md([0], []).tolist() == []


def test_check():
    check(XADJ, ADJNCY)
    bad = list(ADJNCY)
    bad[0] = 2  # vertex 0 lists 2, but 2 does not list 0

    try:
        check(XADJ, bad)
        raised = False
    except ValueError:
        raised = True

    assert raised


if __name__ == "__main__":
    test_arrays()
    test_matrix()
    test_options()
    test_empty()
    test_check()
    print("ok")
