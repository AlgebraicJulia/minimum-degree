# minimum-degree

The exact minimum degree algorithm, in Julia, Python and C.

The function `md` computes a minimum degree permutation of a graph with n
vertices and m edges in O(mn) time and O(m + n) space. It returns the inverse
permutation `invp`: `invp[v]` is the position of vertex `v`.

The three versions compute the same orderings.

- `julia/`: the Julia package `MDLib`
- `c/`: the C library (`md.h`, `md.c`, `md_impl.h`)
- `python/`: the Python package `minimum_degree`, which calls the C library

## Julia

```julia-repl
julia> using Pkg; Pkg.develop(path = "julia")

julia> using MDLib, SparseArrays

julia> graph = sparse([
           0 1 0 0 0 0 0 0
           1 0 1 0 0 1 0 0
           0 1 0 1 0 1 1 1
           0 0 1 0 0 0 0 0
           0 0 0 0 0 1 1 0
           0 1 1 0 1 0 0 0
           0 0 1 0 1 0 0 1
           0 0 1 0 0 0 1 0
       ]);

julia> invp = md(8, graph.colptr, graph.rowval)
8-element Vector{Int64}:
 2
 3
 8
 1
 5
 4
 7
 6
```

Tests: `julia --project=julia -e 'using Pkg; Pkg.test()'`.

## Python

The package needs NumPy (and SciPy, for sparse matrices), and a C compiler
to build the library. From this directory:

```sh
pip install ./python
```

```python
>>> import numpy as np, scipy.sparse
>>> from minimum_degree import md, check
>>> xadj = [0, 1, 4, 9, 10, 12, 15, 18, 20]
>>> adjncy = [1, 0, 2, 5, 1, 3, 5, 6, 7, 2, 5, 6, 1, 2, 4, 2, 4, 7, 2, 6]
>>> md(xadj, adjncy)
array([1, 2, 7, 0, 4, 3, 6, 5])
>>> a = scipy.sparse.csr_matrix((np.ones(20), adjncy, xadj))
>>> md(a)                       # the graph of A + A^T
array([1, 2, 7, 0, 4, 3, 6, 5])
>>> perm = np.argsort(md(a))    # the permutation
```

Indices start at 0. Tests: `python python/tests/test_md.py`.

## C

```sh
cd c
make          # libmd.a, and libmd.so (libmd.dylib on macOS)
make test     # builds and runs example.c
```

Or compile `md.c` with the rest of a program (`cc -O3 -c md.c`, and link
with `-lm`). With MSVC, build a DLL with `cl /O2 /LD /DMD_BUILD_DLL md.c`,
and compile programs that use it with `/DMD_USE_DLL`.

```c
#include "md.h"

int32_t xadj[9] = {0, 1, 4, 9, 10, 12, 15, 18, 20};
int32_t adjncy[20] = {1, 0, 2, 5, 1, 3, 5, 6, 7, 2, 5, 6, 1, 2, 4, 2, 4, 7, 2, 6};
int32_t invp[8];
int status = md_order(8, NULL, xadj, adjncy, invp, NULL);   /* MD_OK */
/* invp = {1, 2, 7, 0, 4, 3, 6, 5} */
```

Indices start at 0. `md_order` takes 32-bit integers, and `md_l_order` 64-bit
integers (`md_order` returns `MD_TOO_LARGE` when 32 bits do not suffice). The
last argument holds the options (`NULL` for the defaults; see `md_options`
and `md_defaults` in `md.h`). The functions return `MD_OK`, `MD_INVALID`,
`MD_TOO_LARGE` or `MD_OUT_OF_MEMORY`.

`md_order` checks its vectors, but not that the graph is symmetric, with no
repeated entries; `md_check` (`md_l_check`) checks that too, in O(m + n) time.
Compiled with `-DMD_BENCH`, the library also provides `md_core_seconds()`, the
time of the last ordering without the checks, the allocation and the copy of
the graph.

## Parameters

    md(neqns[, vwght], xadj, adjncy; external=true, dense=10.0, mass=true, twin=true)   # Julia
    md(xadj, adjncy=None, vwght=None, *, external=True, dense=10.0, mass=True, twin=True)  # Python
    md_order(neqns, vwght, xadj, adjncy, invp, options)                                # C

  - `neqns`: the number of vertices
  - `vwght`: vertex weights, positive integers; a vertex of weight w stands
    for w indistinguishable vertices (default 1)
  - `xadj`, `adjncy`: the graph; the neighbors of `v` are
    `adjncy[xadj[v]:xadj[v + 1] - 1]`, and self loops are ignored. Each edge
    appears in both directions, and at most once in each
  - `external`: minimize the external degree instead of the degree
  - `dense`: vertices with more than `max(16, dense * sqrt(neqns))` neighbors
    are handled as dense (a negative value disables this)
  - `mass`: mass elimination
  - `twin`: false-twin detection
