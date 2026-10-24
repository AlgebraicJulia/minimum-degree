// md.c -- the exact minimum degree algorithm (see md.h).
//
// The degree of every node is kept as a lower bound. Exact degrees are
// computed only for the nodes of minimum bound, once no exact node of
// minimum bound is left, in batches that share their most recent element,
// so each element is scanned once per batch.
// With external degree, an eager pass at each elimination also recomputes
// the nodes whose bound lies within their own weight of the floor.
// For information on the minimum degree algorithm see the articles:
// The evolution of the minimum degree algorithm by Alan George and
// Joseph Liu, SIAM Rev. 31 pp. 1 - 19, 1989
// A fast minimum degree algorithm and matching lower bound by Robert
// Cummings, Matthew Fahrbach, and Animesh Fatehpuria, SODA 2021
//
// The algorithm is in md_impl.h, which is compiled here once for each choice
// of integer types.

#include "md.h"

#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#ifdef MD_BENCH
#include <time.h>
#endif

// ---------------------------------------------------------------------
// helpers
// ---------------------------------------------------------------------

#define MD_MAX(a, b) ((a) > (b) ? (a) : (b))
#define MD_MIN(a, b) ((a) < (b) ? (a) : (b))
#define MD_ABS(a) ((a) < 0 ? -(a) : (a))

// This routine asks the processor to bring `A[i]` into its caches. It is a
// hint: it never faults, so `i` may lie outside `A`.
#if defined(__GNUC__) || defined(__clang__)
#define md_prefetch(A, i) __builtin_prefetch((const void *) &(A)[i], 0, 3)
#else
#define md_prefetch(A, i) ((void) 0)
#endif

// number of ones in `x`
static inline int64_t md_popcount(uint64_t x) {
#if defined(__GNUC__) || defined(__clang__)
    return (int64_t) __builtin_popcountll(x);
#else
    int64_t cnt = 0;

    while (x != 0) {
        x &= x - 1; cnt += 1;
    }

    return cnt;
#endif
}

// number of trailing zeros of `x`, which is not 0
static inline int64_t md_ctz(uint64_t x) {
#if defined(__GNUC__) || defined(__clang__)
    return (int64_t) __builtin_ctzll(x);
#else
    int64_t cnt = 0;

    while ((x & 1) == 0) {
        x >>= 1; cnt += 1;
    }

    return cnt;
#endif
}

// the least power of 2 at least `x`, and the greatest at most `x` (x >= 1)
static int64_t md_nextpow2(int64_t x) {
    int64_t p = 1;

    while (p < x) {
        p *= 2;
    }

    return p;
}

static int64_t md_prevpow2(int64_t x) {
    int64_t p = 1;

    while (2 * p <= x) {
        p *= 2;
    }

    return p;
}

// a list is hashed as the sum of a mix of its entries, so that lists with
// the same sum of entries rarely collide
static inline uint64_t md_mix(uint64_t x) {
    uint64_t z = x * UINT64_C(0x9e3779b97f4a7c15);
    return z ^ (z >> 29);
}

// malloc(count * size), or NULL if that overflows
static void *md_malloc(size_t count, size_t size) {
    void *p = NULL;

    if (size == 0 || count <= SIZE_MAX / size) {
        p = malloc(count * size);
    }

    return p;
}

// the candidates for the dense nodes: sorted by key (minus the degree), then
// by node
typedef struct {
    int64_t key;
    int64_t node;
} md_pair;

static int md_pair_compare(const void *a, const void *b) {
    const md_pair *x = (const md_pair *) a;
    const md_pair *y = (const md_pair *) b;
    int cmp = 0;

    if (x->key < y->key) {
        cmp = -1;
    } else if (x->key > y->key) {
        cmp = 1;
    } else if (x->node < y->node) {
        cmp = -1;
    } else if (x->node > y->node) {
        cmp = 1;
    }

    return cmp;
}

static int md_pair_compare_node(const void *a, const void *b) {
    const md_pair *x = (const md_pair *) a;
    const md_pair *y = (const md_pair *) b;
    int cmp = 0;

    if (x->node < y->node) {
        cmp = -1;
    } else if (x->node > y->node) {
        cmp = 1;
    }

    return cmp;
}

#ifdef MD_BENCH
// the time of the last ordering, without the allocation and the copy of the
// graph (a benchmarking build only)
static double md_seconds = 0.0;

static double md_clock(void) {
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return (double) t.tv_sec + 1e-9 * (double) t.tv_nsec;
}

double md_core_seconds(void) {
    return md_seconds;
}
#endif

#ifdef MD_COUNT
// the work counters of the last ordering (a counting build only; see md.h).
// `md_count_rc` is the first counter of the recomputation being run
// (`MD_C_PASS_MARK` or `MD_C_WATER_MARK`), and `md_count_h` the position of
// the water walk (-1 before the first operation). The water of an operation
// is the key of its level minus 1; after an elimination, the walk is at the
// intermediate position (see `mdupdate`).
static int64_t md_count[MD_NCOUNT];
static int md_count_rc = MD_C_WATER_MARK;
static int64_t md_count_h = -1;

#define MD_ADD(k, x) (md_count[k] += (int64_t) (x))
#define MD_ADD_RC(i, x) (md_count[md_count_rc + (i)] += (int64_t) (x))
#define MD_SET_RC(k) (md_count_rc = (k))
#define MD_WALK(w) md_count_walk((int64_t) (w))

// an operation at water `w`: the walk rises (or falls) to `w`
static void md_count_walk(int64_t w) {
    if (md_count_h >= 0 && w > md_count_h) {
        md_count[MD_C_RISE] += w - md_count_h;
    }

    md_count_h = w;
}

void md_counts(int64_t out[MD_NCOUNT]) {
    for (int k = 0; k < MD_NCOUNT; k++) {
        out[k] = md_count[k];
    }
}
#else
#define MD_ADD(k, x) ((void) 0)
#define MD_ADD_RC(i, x) ((void) 0)
#define MD_SET_RC(k) ((void) 0)
#define MD_WALK(w) ((void) 0)
#endif

// ---------------------------------------------------------------------
// the algorithm, for each choice of integer types
// ---------------------------------------------------------------------

// 32-bit vertices, positions and tags
#define V int32_t
#define E int32_t
#define I int32_t
#define V_MAX INT32_MAX
#define I_MAX INT32_MAX
#define MDF(f) f##_32
#include "md_impl.h"
#undef V
#undef E
#undef I
#undef V_MAX
#undef I_MAX
#undef MDF

// 32-bit vertices and positions, 64-bit tags
#define V int32_t
#define E int32_t
#define I int64_t
#define V_MAX INT32_MAX
#define I_MAX INT64_MAX
#define MDF(f) f##_32l
#include "md_impl.h"
#undef V
#undef E
#undef I
#undef V_MAX
#undef I_MAX
#undef MDF

// 64-bit vertices, positions and tags
#define V int64_t
#define E int64_t
#define I int64_t
#define V_MAX INT64_MAX
#define I_MAX INT64_MAX
#define MDF(f) f##_64
#include "md_impl.h"
#undef V
#undef E
#undef I
#undef V_MAX
#undef I_MAX
#undef MDF

// ---------------------------------------------------------------------
// the entry points
// ---------------------------------------------------------------------

void md_defaults(md_options *options) {
    options->external = 1;
    options->dense = 10.0;
    options->mass = 1;
    options->twin = 1;
}

// This routine checks the vectors of the graph and the weights (not the
// symmetry of the graph: see `md_check`), and computes the total weight and
// the length of `iw`. Returns MD_OK or MD_INVALID.
#define MD_VALID(T)                                                                  \
    static int md_valid_##T(T neqns, const T *vwght, const T *xadj, const T *adjncy, \
            int64_t *total, int64_t *iwlen) {                                        \
        int status = MD_OK;                                                          \
        int64_t nnz = 0;                                                             \
        *total = 0; *iwlen = 0;                                                      \
                                                                                     \
        if (neqns < 0 || xadj == NULL || xadj[0] != 0) {                             \
            status = MD_INVALID;                                                     \
        }                                                                            \
                                                                                     \
        /* (the loops over xadj and adjncy do not stop early, so that they           \
           vectorize) */                                                             \
        if (status == MD_OK) {                                                       \
            int bad = 0;                                                             \
                                                                                     \
            for (int64_t v = 0; v < (int64_t) neqns; v++) {                          \
                bad |= xadj[v + 1] < xadj[v];                                        \
            }                                                                        \
                                                                                     \
            nnz = (int64_t) xadj[neqns];                                             \
            bad |= nnz > 0 && adjncy == NULL;                                        \
                                                                                     \
            if (bad) {                                                               \
                status = MD_INVALID;                                                 \
            }                                                                        \
        }                                                                            \
                                                                                     \
        if (status == MD_OK) {                                                       \
            uint64_t bad = 0;                                                        \
                                                                                     \
            for (int64_t i = 0; i < nnz; i++) {                                      \
                bad |= (uint64_t) (int64_t) adjncy[i] >= (uint64_t) neqns;           \
            }                                                                        \
                                                                                     \
            if (bad) {                                                               \
                status = MD_INVALID;                                                 \
            }                                                                        \
        }                                                                            \
                                                                                     \
        if (status == MD_OK && vwght == NULL) {                                      \
            *total = (int64_t) neqns;                                                \
        } else if (status == MD_OK) {                                                \
            for (int64_t v = 0; v < (int64_t) neqns && status == MD_OK; v++) {       \
                if (vwght[v] <= 0 || (int64_t) vwght[v] > INT64_MAX / 4 - *total) {  \
                    status = MD_INVALID;                                             \
                } else {                                                             \
                    *total += (int64_t) vwght[v];                                    \
                }                                                                    \
            }                                                                        \
        }                                                                            \
                                                                                     \
        /* the lists never take more than `nnz` entries of `iw`, except while a      \
           dense node rejoins, which adds fewer than `2 neqns`. The extra            \
           `max(nnz / 5, neqns) + neqns` covers that, and leaves room between        \
           compressions */                                                           \
        *iwlen = nnz + MD_MAX(nnz / 5, (int64_t) neqns) + (int64_t) neqns;           \
        return status;                                                               \
    }

MD_VALID(int32_t)
MD_VALID(int64_t)

int md_order(int32_t neqns, const int32_t *vwght, const int32_t *xadj,
        const int32_t *adjncy, int32_t *invp, const md_options *options) {
    md_options defaults;
    int64_t total, iwlen;
    md_defaults(&defaults);

    if (options == NULL) {
        options = &defaults;
    }

    int status = md_valid_int32_t(neqns, vwght, xadj, adjncy, &total, &iwlen);

    if (status == MD_OK && invp == NULL) {
        status = MD_INVALID;
    }

    // 32-bit integers when they suffice: positions in `iw` too, and 32-bit
    // tags when they leave room for many steps between resets
    if (status == MD_OK) {
        int64_t nnz = (int64_t) xadj[neqns];

        if (!(neqns < INT32_MAX / 8 && total < INT32_MAX / 4 && nnz < INT32_MAX && iwlen < INT32_MAX)) {
            status = MD_TOO_LARGE;
        } else if (neqns < INT32_MAX / 64) {
            status = md_32(neqns, (int32_t) total, iwlen, vwght, xadj, adjncy, options->external,
                options->dense, options->mass, options->twin, invp);
        } else {
            status = md_32l(neqns, (int32_t) total, iwlen, vwght, xadj, adjncy, options->external,
                options->dense, options->mass, options->twin, invp);
        }
    }

    return status;
}

int md_l_order(int64_t neqns, const int64_t *vwght, const int64_t *xadj,
        const int64_t *adjncy, int64_t *invp, const md_options *options) {
    md_options defaults;
    int64_t total, iwlen;
    md_defaults(&defaults);

    if (options == NULL) {
        options = &defaults;
    }

    int status = md_valid_int64_t(neqns, vwght, xadj, adjncy, &total, &iwlen);

    if (status == MD_OK && invp == NULL) {
        status = MD_INVALID;
    }

    if (status == MD_OK) {
        status = md_64(neqns, total, iwlen, vwght, xadj, adjncy, options->external,
            options->dense, options->mass, options->twin, invp);
    }

    return status;
}

// This routine checks that the graph is symmetric, with no repeated
// entries: it counts the entries of each column, lists them (the transpose),
// and compares each row with its column.
#define MD_CHECK(T, NAME)                                                        \
    int NAME(T neqns, const T *xadj, const T *adjncy) {                          \
        int64_t total, iwlen;                                                    \
        int status = md_valid_##T(neqns, NULL, xadj, adjncy, &total, &iwlen);    \
        int64_t n = (int64_t) neqns;                                             \
        int64_t *ptr = NULL;                                                     \
        int64_t *mark = NULL;                                                    \
        T *trans = NULL;                                                         \
                                                                                 \
        if (status == MD_OK) {                                                   \
            ptr = (int64_t *) md_malloc((size_t) n + 1, sizeof(int64_t));        \
            mark = (int64_t *) md_malloc((size_t) n + 1, sizeof(int64_t));       \
            trans = (T *) md_malloc((size_t) xadj[n] + 1, sizeof(T));            \
                                                                                 \
            if (ptr == NULL || mark == NULL || trans == NULL) {                  \
                status = MD_OUT_OF_MEMORY;                                       \
            }                                                                    \
        }                                                                        \
                                                                                 \
        if (status == MD_OK) {                                                   \
            for (int64_t v = 0; v <= n; v++) {                                   \
                ptr[v] = 0;                                                      \
            }                                                                    \
                                                                                 \
            for (int64_t i = 0; i < (int64_t) xadj[n]; i++) {                    \
                ptr[(int64_t) adjncy[i] + 1] += 1;                               \
            }                                                                    \
                                                                                 \
            for (int64_t v = 0; v < n && status == MD_OK; v++) {                 \
                if (ptr[v + 1] != (int64_t) (xadj[v + 1] - xadj[v])) {           \
                    status = MD_INVALID;                                         \
                }                                                                \
                                                                                 \
                ptr[v + 1] += ptr[v];                                            \
            }                                                                    \
        }                                                                        \
                                                                                 \
        if (status == MD_OK) {                                                   \
            for (int64_t v = 0; v < n; v++) {                                    \
                for (int64_t i = (int64_t) xadj[v]; i < (int64_t) xadj[v + 1]; i++) { \
                    trans[ptr[adjncy[i]]] = (T) v; ptr[adjncy[i]] += 1;          \
                }                                                                \
            }                                                                    \
                                                                                 \
            for (int64_t v = 0; v < n; v++) {                                    \
                mark[v] = -1;                                                    \
            }                                                                    \
                                                                                 \
            for (int64_t v = 0; v < n && status == MD_OK; v++) {                 \
                for (int64_t i = (int64_t) xadj[v]; i < (int64_t) xadj[v + 1]; i++) { \
                    if (mark[adjncy[i]] == v) {                                  \
                        status = MD_INVALID;                                     \
                    }                                                            \
                                                                                 \
                    mark[adjncy[i]] = v;                                         \
                }                                                                \
                                                                                 \
                for (int64_t i = (int64_t) xadj[v]; i < ptr[v] && status == MD_OK; i++) { \
                    if (mark[trans[i]] != v) {                                   \
                        status = MD_INVALID;                                     \
                    }                                                            \
                }                                                                \
            }                                                                    \
        }                                                                        \
                                                                                 \
        free(ptr); free(mark); free(trans);                                      \
        return status;                                                           \
    }

MD_CHECK(int32_t, md_check)
MD_CHECK(int64_t, md_l_check)
