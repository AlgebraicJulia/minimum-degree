// example.c -- orders the graph of the README, and checks the result.

#include <stdio.h>

#include "md.h"

int main(void) {
    // the neighbors of vertex v are adjncy[xadj[v]], ..., adjncy[xadj[v + 1] - 1]
    int32_t xadj[9] = {0, 1, 4, 9, 10, 12, 15, 18, 20};
    int32_t adjncy[20] = {1, 0, 2, 5, 1, 3, 5, 6, 7, 2, 5, 6, 1, 2, 4, 2, 4, 7, 2, 6};
    int32_t expected[8] = {1, 2, 7, 0, 4, 3, 6, 5};
    int32_t invp[8];
    int64_t xadj64[9], adjncy64[20], invp64[8];
    int ok = md_check(8, xadj, adjncy) == MD_OK;

    // the default options (NULL), with 32-bit and with 64-bit integers
    ok = ok && md_order(8, NULL, xadj, adjncy, invp, NULL) == MD_OK;

    for (int v = 0; v < 9; v++) {
        xadj64[v] = xadj[v];
    }

    for (int i = 0; i < 20; i++) {
        adjncy64[i] = adjncy[i];
    }

    ok = ok && md_l_order(8, NULL, xadj64, adjncy64, invp64, NULL) == MD_OK;
    printf("invp:");

    for (int v = 0; v < 8 && ok; v++) {
        printf(" %d", (int) invp[v]);
        ok = invp[v] == expected[v] && invp64[v] == expected[v];
    }

    printf("\n");
    int status = 0;

    if (ok) {
        printf("ok\n");
    } else {
        printf("unexpected ordering\n");
        status = 1;
    }

    return status;
}
