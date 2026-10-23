"""The exact minimum degree algorithm.

The function `md` computes a minimum degree permutation of a graph with n
vertices and m edges in O(mn) time and O(m + n) space. It returns the inverse
permutation `invp`. It calls the C library in ../c (see md.h).
"""

import ctypes
import importlib.util
import os

import numpy as np

__all__ = ["md", "check"]
__version__ = "0.1.0"

# return values of the C library
_MD_OK = 0
_MD_INVALID = -1
_MD_TOO_LARGE = -2
_MD_OUT_OF_MEMORY = -3


class _Options(ctypes.Structure):
    _fields_ = [
        ("external", ctypes.c_int),
        ("dense", ctypes.c_double),
        ("mass", ctypes.c_int),
        ("twin", ctypes.c_int),
    ]


def _load():
    # the library built with the package, or the one named by MD_LIBRARY
    path = os.environ.get("MD_LIBRARY")

    if path is None:
        spec = importlib.util.find_spec("minimum_degree._libmd")

        if spec is not None:
            path = spec.origin

    if path is None:
        raise ImportError("minimum_degree: the C library was not built; "
                          "install the package with pip, or set MD_LIBRARY")

    lib = ctypes.CDLL(path)

    for name, ctype in (("md_order", ctypes.c_int32), ("md_l_order", ctypes.c_int64)):
        fun = getattr(lib, name)
        ptr = ctypes.POINTER(ctype)
        fun.argtypes = [ctype, ptr, ptr, ptr, ptr, ctypes.POINTER(_Options)]
        fun.restype = ctypes.c_int

    for name, ctype in (("md_check", ctypes.c_int32), ("md_l_check", ctypes.c_int64)):
        fun = getattr(lib, name)
        ptr = ctypes.POINTER(ctype)
        fun.argtypes = [ctype, ptr, ptr]
        fun.restype = ctypes.c_int

    return lib


_lib = _load()


def _pattern(matrix):
    # the graph of a square sparse matrix: the pattern of A + A^T
    import scipy.sparse

    a = scipy.sparse.csr_matrix(matrix)

    if a.shape[0] != a.shape[1]:
        raise ValueError("md: the matrix must be square")

    s = scipy.sparse.csr_matrix((np.ones(a.nnz), a.indices, a.indptr), shape=a.shape)
    s = (s + s.T).tocsr()
    s.sum_duplicates()
    return s.indptr, s.indices


def _call(fun, dtype, n, vwght, xadj, adjncy, options):
    ptr = ctypes.POINTER(np.ctypeslib.as_ctypes_type(dtype))
    x = np.ascontiguousarray(xadj, dtype=dtype)
    a = np.ascontiguousarray(adjncy, dtype=dtype)
    invp = np.empty(n, dtype=dtype)
    w = None

    if vwght is not None:
        w = np.ascontiguousarray(vwght, dtype=dtype).ctypes.data_as(ptr)

    status = fun(n, w, x.ctypes.data_as(ptr), a.ctypes.data_as(ptr),
                 invp.ctypes.data_as(ptr), ctypes.byref(options))
    return status, invp


def md(xadj, adjncy=None, vwght=None, *, external=True, dense=10.0, mass=True, twin=True):
    """Compute a minimum degree permutation of a graph with n vertices and m
    edges in O(mn) time and O(m + n) space. Returns the inverse permutation
    `invp`: `invp[v]` is the position of vertex `v` (so the permutation is
    `numpy.argsort(invp)`).

    - `xadj`, `adjncy`: the graph, indexed from 0; the neighbors of `v` are
      `adjncy[xadj[v]:xadj[v + 1]]`. Each edge must appear in both
      directions, and at most once in each (see `check`); self loops are
      ignored. Or, with `adjncy` omitted, `xadj` is a square SciPy sparse
      matrix (or a dense array), and the graph is that of its pattern plus
      its transpose
    - `vwght`: vertex weights, positive integers; a vertex of weight w stands
      for w indistinguishable vertices (default 1)
    - `external`: minimize the external degree instead of the degree
    - `dense`: vertices with more than `max(16, dense * sqrt(n))` neighbors
      are handled as dense (a negative value disables this)
    - `mass`: mass elimination
    - `twin`: false-twin detection
    """
    if adjncy is None:
        xadj, adjncy = _pattern(xadj)

    xadj = np.asarray(xadj)
    adjncy = np.asarray(adjncy)

    if xadj.ndim != 1 or adjncy.ndim != 1 or xadj.size < 1:
        raise ValueError("md: xadj and adjncy must be vectors, with len(xadj) >= 1")

    n = xadj.size - 1

    if vwght is not None:
        vwght = np.asarray(vwght)

        if vwght.shape != (n,) or not np.all(vwght == np.floor(vwght)) or np.any(vwght < 1):
            raise ValueError("md: vwght must hold n positive integers")

    options = _Options(int(bool(external)), float(dense), int(bool(mass)), int(bool(twin)))
    status = _MD_TOO_LARGE
    small = np.iinfo(np.int32).max

    # 32-bit integers when they suffice
    if n < small and xadj[-1] < small:
        status, invp = _call(_lib.md_order, np.int32, n, vwght, xadj, adjncy, options)

    if status == _MD_TOO_LARGE:
        status, invp = _call(_lib.md_l_order, np.int64, n, vwght, xadj, adjncy, options)

    if status == _MD_INVALID:
        raise ValueError("md: invalid graph or weights")
    elif status == _MD_OUT_OF_MEMORY:
        raise MemoryError("md: out of memory")
    elif status != _MD_OK:
        raise RuntimeError(f"md: error {status}")

    return invp.astype(np.int64)


def check(xadj, adjncy):
    """Check that the graph `(xadj, adjncy)` (see `md`) is valid: indices in
    range, each edge in both directions, and at most once in each. `md` does
    not check this, and its result is undefined if it does not hold. Raises
    ValueError if it does not."""
    xadj = np.asarray(xadj)
    adjncy = np.asarray(adjncy)

    if xadj.ndim != 1 or adjncy.ndim != 1 or xadj.size < 1:
        raise ValueError("check: xadj and adjncy must be vectors, with len(xadj) >= 1")

    n = xadj.size - 1
    small = np.iinfo(np.int32).max

    if n < small and xadj[-1] < small:
        dtype, fun = np.int32, _lib.md_check
    else:
        dtype, fun = np.int64, _lib.md_l_check

    ptr = ctypes.POINTER(np.ctypeslib.as_ctypes_type(dtype))
    x = np.ascontiguousarray(xadj, dtype=dtype)
    a = np.ascontiguousarray(adjncy, dtype=dtype)
    status = fun(n, x.ctypes.data_as(ptr), a.ctypes.data_as(ptr))

    if status == _MD_INVALID:
        raise ValueError("check: the graph is not valid")
    elif status == _MD_OUT_OF_MEMORY:
        raise MemoryError("check: out of memory")
