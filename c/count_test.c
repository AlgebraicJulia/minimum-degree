// count_test.c -- orders random graphs with every option, and prints a hash
// of each ordering. Built twice, with and without -DMD_COUNT, it must print
// the same lines: the counters do not change the orderings (make count-test).
// The counting build also prints the counters, summed over the runs, to
// stderr.

#include <stdio.h>
#include <stdlib.h>

#include "md.h"

// a linear congruential generator
static uint64_t state = 2026;

static uint64_t next(void) {
    state = state * 6364136223846793005ULL + 1442695040888963407ULL;
    return state >> 33;
}

// a random symmetric graph on `n` vertices with about `n * deg / 2` edges,
// plus `nhubs` vertices adjacent to about half of the others; the neighbors
// of vertex v are adjncy[xadj[v]], ..., adjncy[xadj[v + 1] - 1]
static void graph(int32_t n, int32_t deg, int32_t nhubs, int32_t *xadj, int32_t **adjncy) {
    char *adj = calloc((size_t) n * (size_t) n, 1);

    for (int32_t k = 0; k < n * deg / 2; k++) {
        int32_t u = (int32_t) (next() % (uint64_t) n);
        int32_t v = (int32_t) (next() % (uint64_t) n);

        if (u != v) {
            adj[(size_t) u * n + v] = 1; adj[(size_t) v * n + u] = 1;
        }
    }

    for (int32_t h = 0; h < nhubs; h++) {
        int32_t u = (int32_t) (next() % (uint64_t) n);

        for (int32_t v = 0; v < n; v++) {
            if (v != u && next() % 2 == 0) {
                adj[(size_t) u * n + v] = 1; adj[(size_t) v * n + u] = 1;
            }
        }
    }

    xadj[0] = 0;

    for (int32_t v = 0; v < n; v++) {
        xadj[v + 1] = xadj[v];

        for (int32_t u = 0; u < n; u++) {
            xadj[v + 1] += adj[(size_t) v * n + u];
        }
    }

    *adjncy = malloc(sizeof(int32_t) * ((size_t) xadj[n] + 1));

    for (int32_t v = 0; v < n; v++) {
        int32_t p = xadj[v];

        for (int32_t u = 0; u < n; u++) {
            if (adj[(size_t) v * n + u]) {
                (*adjncy)[p] = u; p += 1;
            }
        }
    }

    free(adj);
}

int main(void) {
    double dense[3] = {10.0, 1.0, -1.0};
    int nruns = 0;

#ifdef MD_COUNT
    int64_t sum[MD_NCOUNT] = {0};
#endif

    for (int trial = 0; trial < 60; trial++) {
        int32_t n = 1 + (int32_t) (next() % 200);
        int32_t deg = 1 + (int32_t) (next() % 12);
        int32_t nhubs = (int32_t) (next() % 4);
        int32_t *xadj = malloc(sizeof(int32_t) * ((size_t) n + 1));
        int32_t *adjncy = NULL;
        int32_t *vwght = NULL;
        int32_t *invp = malloc(sizeof(int32_t) * (size_t) n);
        graph(n, deg, nhubs, xadj, &adjncy);

        // vertex weights in half of the trials
        if (trial % 2 == 1) {
            vwght = malloc(sizeof(int32_t) * (size_t) n);

            for (int32_t v = 0; v < n; v++) {
                vwght[v] = 1 + (int32_t) (next() % 3);
            }
        }

        for (int k = 0; k < 48; k++) {
            md_options options;
            options.external = k % 2;
            options.mass = (k / 2) % 2;
            options.twin = (k / 4) % 2;
            options.dense = dense[(k / 8) % 3];
            options.compress = (k / 24) % 2;

            if (md_order(n, vwght, xadj, adjncy, invp, &options) != MD_OK) {
                printf("error\n");
                return 1;
            }

            // a hash of the ordering
            uint64_t hash = 14695981039346656037ULL;

            for (int32_t v = 0; v < n; v++) {
                hash = (hash ^ (uint64_t) (uint32_t) invp[v]) * 1099511628211ULL;
            }

            printf("%d %d %016llx\n", trial, k, (unsigned long long) hash);
            nruns += 1;

#ifdef MD_COUNT
            int64_t count[MD_NCOUNT];
            md_counts(count);

            for (int c = 0; c < MD_NCOUNT; c++) {
                sum[c] += count[c];
            }
#endif
        }

        free(xadj); free(adjncy); free(vwght); free(invp);
    }

    printf("%d runs\n", nruns);

#ifdef MD_COUNT
    const char *names[MD_NCOUNT] = {"ELIM", "COMPRESS", "DENSE",
        "PASS_MARK", "PASS_SCAN", "PASS_UNION", "PASS_VAR",
        "WATER_MARK", "WATER_SCAN", "WATER_UNION", "WATER_VAR",
        "STEPS", "TOUCHED", "CERTIFIED", "RECOMP_PASS", "RECOMP_WATER", "GROUPS", "RISE"};

    fprintf(stderr, "counters, summed over %d runs:\n", nruns);

    for (int c = 0; c < MD_NCOUNT; c++) {
        fprintf(stderr, "  %-13s %lld\n", names[c], (long long) sum[c]);
    }
#endif

    return 0;
}
