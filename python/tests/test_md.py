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
    invp = md(XADJ, ADJNCY, external=False, mass=False, twin=False, compress=False)
    assert sorted(invp.tolist()) == list(range(8))
    assert md(XADJ, ADJNCY, vwght=np.ones(8, dtype=int)).tolist() == INVP


def test_heavy():
    # a path on 6 vertices, every weight 100: the space for a reachable set
    # is bounded by its number of nodes, not by its weight
    xadj = [0, 1, 3, 5, 7, 9, 10]
    adjncy = [1, 0, 2, 1, 3, 2, 4, 3, 5, 4]
    invp = md(xadj, adjncy, vwght=np.full(6, 100))
    assert sorted(invp.tolist()) == list(range(6))


def test_dense():
    # two hubs (40 and 41) adjacent to every vertex of a path on 40 vertices
    # are dense, removed and ordered last
    adj = [[] for _ in range(42)]

    for u in range(39):
        adj[u].append(u + 1)
        adj[u + 1].append(u)

    for h in (40, 41):
        for v in range(40):
            adj[h].append(v)
            adj[v].append(h)

    xadj = np.cumsum([0] + [len(a) for a in adj])
    adjncy = [v for a in adj for v in a]

    for external in (False, True):
        invp = md(xadj, adjncy, external=external, dense=0.0).tolist()
        assert sorted(invp) == list(range(42)) and invp[40] == 40 and invp[41] == 41


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
    test_heavy()
    test_dense()
    test_empty()
    test_check()
    print("ok")
