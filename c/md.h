/* md.h -- the exact minimum degree algorithm.
 *
 * md_order computes a minimum degree ordering of a graph with n vertices and
 * m edges in O(mn) time and O(m + n) space, and returns the inverse
 * permutation: invp[v] is the position of vertex v. md_order takes 32-bit
 * integers and md_l_order 64-bit integers; otherwise they are the same.
 *
 * The graph is given in compressed form, indexed from 0: the neighbors of
 * vertex v are adjncy[xadj[v]], ..., adjncy[xadj[v + 1] - 1], with
 * xadj[0] = 0. Each edge must appear in both directions, and at most once in
 * each; self loops are ignored.
 */

#ifndef MD_H
#define MD_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#if defined(_WIN32) && defined(MD_BUILD_DLL)
#define MD_API __declspec(dllexport)
#elif defined(_WIN32) && defined(MD_USE_DLL)
#define MD_API __declspec(dllimport)
#else
#define MD_API
#endif

/* return values */
#define MD_OK 0             /* success */
#define MD_INVALID (-1)     /* invalid input */
#define MD_TOO_LARGE (-2)   /* too large for md_order: use md_l_order */
#define MD_OUT_OF_MEMORY (-3)

typedef struct {
    int external;   /* minimize the external degree instead of the degree
                       (default 1) */
    double dense;   /* vertices with more than max(16, dense * sqrt(n))
                       neighbors are handled as dense; a negative value
                       disables this (default 10.0) */
    int mass;       /* mass elimination (default 1) */
    int twin;       /* false-twin detection (default 1) */
    int compress;   /* merge the vertices with the same closed neighborhood
                       first (default 1) */
} md_options;

/* set the default options */
MD_API void md_defaults(md_options *options);

/* neqns: the number of vertices
 * vwght: vertex weights, positive integers; a vertex of weight w stands for w
 *        indistinguishable vertices (NULL: all 1)
 * xadj, adjncy: the graph
 * invp: on output, the inverse permutation (length neqns)
 * options: the options (NULL: the defaults)
 *
 * Returns MD_OK, or one of the error values above. */
MD_API int md_order(int32_t neqns, const int32_t *vwght, const int32_t *xadj,
    const int32_t *adjncy, int32_t *invp, const md_options *options);

MD_API int md_l_order(int64_t neqns, const int64_t *vwght, const int64_t *xadj,
    const int64_t *adjncy, int64_t *invp, const md_options *options);

/* md_order checks the vectors it is given, but not that the graph is
 * symmetric, with no repeated entries; if it is not, the result is undefined.
 * md_check checks that too, in O(m + n) time and space. Returns MD_OK,
 * MD_INVALID or MD_OUT_OF_MEMORY. */
MD_API int md_check(int32_t neqns, const int32_t *xadj, const int32_t *adjncy);

MD_API int md_l_check(int64_t neqns, const int64_t *xadj, const int64_t *adjncy);

#ifdef MD_BENCH
/* the time in seconds of the last ordering, without the checks, the
 * allocation and the copy of the graph (a benchmarking build, compiled with
 * -DMD_BENCH) */
MD_API double md_core_seconds(void);
#endif

#ifdef MD_COUNT
/* work counters of the last ordering (a counting build, compiled with
 * -DMD_COUNT; not for timing). The first eleven count entries of adjacency
 * lists read, each once by each loop that runs over its list:
 *
 *   MD_C_ELIM         the elimination (mdelim, mdupdate), outside the passes
 *   MD_C_COMPRESS     entries moved by the compression of iw
 *   MD_C_DENSE        dense nodes: rejoins, and the reachable sets read to
 *                     update the masks
 *   MD_C_PASS_*       the eager passes (recomputations from mdupdate):
 *                     marking the nodes of the element (MARK), partitioning
 *                     the other elements (SCAN), the unions (UNION), and the
 *                     lists of the recomputed nodes (VAR)
 *   MD_C_WATER_*      the same, for the recomputations at the water (mdvisit)
 *
 * and the others count events:
 *
 *   MD_C_STEPS        eliminations (nodes of degree 0 at the start, pivots,
 *                     and rejoining dense nodes)
 *   MD_C_TOUCHED      the principal nodes of each reachable set, after the
 *                     merges and before mass elimination
 *   MD_C_CERTIFIED    those of them exact at their initial bound
 *   MD_C_RECOMP_PASS  nodes recomputed by the eager passes
 *   MD_C_RECOMP_WATER nodes recomputed at the water
 *   MD_C_GROUPS       groups recomputed at the water
 *   MD_C_RISE         the total rise of the water walk */
enum { MD_C_ELIM, MD_C_COMPRESS, MD_C_DENSE,
       MD_C_PASS_MARK, MD_C_PASS_SCAN, MD_C_PASS_UNION, MD_C_PASS_VAR,
       MD_C_WATER_MARK, MD_C_WATER_SCAN, MD_C_WATER_UNION, MD_C_WATER_VAR,
       MD_C_STEPS, MD_C_TOUCHED, MD_C_CERTIFIED,
       MD_C_RECOMP_PASS, MD_C_RECOMP_WATER, MD_C_GROUPS,
       MD_C_RISE, MD_NCOUNT };

MD_API void md_counts(int64_t out[MD_NCOUNT]);
#endif

#ifdef __cplusplus
}
#endif

#endif
