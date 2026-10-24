// md_impl.h -- the exact minimum degree algorithm, for one choice of integer
// types. md.c includes this file once for each choice, with these macros
// defined:
//
//   V       the integer type of the vertices, weights and degrees
//   E       the integer type of the positions in iw
//   I       the integer type of the tags in marker
//   V_MAX   the largest V
//   I_MAX   the largest I
//   MDF(f)  the name of the routine f for this choice
//
// The arrays are indexed from 1: entry 0 of each is not used. Only the input
// graph (xadj, adjncy) is indexed from 0.

static int MDF(md)(V n, V total, int64_t iwlen, const V *vwght, const V *xadj,
    const V *adjncy, int external, double dense, int mass, int twin, V *invp);
static V MDF(mdndense)(V neqns, const V *xadj, double dense);
static int MDF(mdselect)(V *dindex, V neqns, const V *xadj, const V *adjncy, double dense,
    V *ndense);
static void MDF(md_impl)(I *marker, V *superwght, V *degree, V *deghead, V *degnext,
    V *degprev, V *wext, V *dwght, V *extlen, E *pe, V *len, V *elen, V *iw, E iwlen,
    uint64_t *heads, V total, V neqns, int external, int mass, int twin,
    const V *vwght, const V *xadj, const V *adjncy, V ndense, int64_t nwords, V *dindex,
    V *dnode, V *ddeg, uint64_t *dmask);
static void MDF(mdmain)(I *marker, V *superwght, V *degree, V *deghead, V *degnext,
    V *degprev, V *wext, V *dwght, V *extlen, E *pe, V *len, V *elen, V *iw, E iwlen,
    uint64_t *heads, E pfree, V total, V neqns, int external, int mass,
    int twin, const V *vwght, V ndense, int64_t nwords, V *dindex, V *dnode, V *ddeg,
    uint64_t *dmask);
static void MDF(mdelim)(V mdnode, E *pe, V *len, V *elen, V *iw, E iwlen, E *pfreep,
    uint64_t *heads, V neqns, V *superwght, V *degree, V *deghead, V *degnext, V *degprev,
    I *marker, V *wext, V *extlen, V *wflgp, V *lemaxp, V dS, V offset, int external,
    int twin, V fnum, V nbypass, int64_t nwords, uint64_t *dmask, V *ddeg, V *rchszep,
    V *deltap, V *swghtp, V *fszep, V *fwghtp);
static void MDF(mdupdate)(V mdnode, V rchsze, V delta, V dS, V swght, V fwght, E *pe, V *len,
    V *elen, V *iw, V *superwght, V *degree, V *deghead, V *degnext, V *degprev, V *wext,
    I *marker, I *tagp, V neqns, V *extlen, V *dwght, int64_t nwords,
    uint64_t *dmask, V *ddeg, V nbypass, V *mindegp, int external, int mass, V offset,
    V *mwghtp);
static void MDF(mdvisit)(V key, E *pe, V *len, V *elen, V *iw, V *superwght, V *degree,
    V *deghead, V *degnext, V *degprev, I *marker, I *tagp, V neqns,
    V *extlen, V *dwght, V nbypass, int external, V offset, V *mindegp);
static void MDF(mdrecompute)(V zhead, V elmnt, E *pe, V *len, V *elen, V *iw, V *superwght,
    V *degree, V *deghead, V *degnext, V *degprev, I *marker, I *tagp,
    V neqns, V *extlen, V *dwght, V nbypass, V *mindegp, int external, V offset);
static void MDF(mdremove)(V node, V *deghead, V *degnext, V *degprev, V offset);
static V MDF(mdinsert)(V key, V mindeg, V node, int isexact, V *deghead, V *degnext,
    V *degprev, V offset);
static V MDF(mdkey)(V deg, V wght, int external);
static V MDF(mdvwght)(const V *vwght, V node);
static int64_t MDF(mdoffset)(int64_t nwords, V node);
static void MDF(mdsetbit)(uint64_t *dmask, int64_t nwords, V node, V j);
static int64_t MDF(mdcount)(const uint64_t *dmask, int64_t nwords, V node);
static V MDF(mdweight)(const uint64_t *dmask, int64_t nwords, const V *dnode, V node,
    const V *superwght);
static int MDF(mdsame)(const uint64_t *dmask, int64_t nwords, V inode, V jnode);
static void MDF(mddenseinit)(V *ddeg, const uint64_t *dmask, int64_t nwords, V ndense,
    const V *dnode, V neqns, const V *vwght);
static V MDF(mdlose)(V *ddeg, const uint64_t *dmask, int64_t nwords, V node, V wght);
static int64_t MDF(mddense)(V mdnode, E rstart, V rchsze, const V *iw, V *dwght,
    const V *superwght, int64_t nwords, const V *dnode, V *ddeg, uint64_t *dmask,
    int64_t nmissing);
static E MDF(mdreinsert)(V j, V neqns, E *pe, V *len, V *elen, V *iw, E iwlen, E pfree,
    uint64_t *heads, const V *superwght, int64_t nwords, const V *dindex, V *dnode,
    uint64_t *dmask, V *wext);
static E MDF(mdcompress)(V neqns, E *pe, const V *len, V *iw, E pfree, uint64_t *heads);
static void MDF(mdnumber)(V neqns, V *degnext, V *size, V *start);

// This routine allocates the working arrays, with vertices, weights and
// degrees of type `V`, positions in `iw` of type `E`, and tags of type `I`,
// and computes the ordering. Unless the allocation fails, `invp[v]` is then
// the position of vertex `v`, both indexed from 0.
static int MDF(md)(V n, V total, int64_t iwlen, const V *vwght, const V *xadj,
        const V *adjncy, int external, double dense, int mass, int twin, V *invp) {
    int status = MD_OK;
    I *marker = (I *) md_malloc((size_t) n + 1, sizeof(I));
    V *superwght = (V *) md_malloc((size_t) n + 1, sizeof(V));
    V *degree = (V *) md_malloc((size_t) n + 1, sizeof(V));
    V *deghead = (V *) md_malloc(2 * ((size_t) total + 1) + 1, sizeof(V));
    V *degnext = (V *) md_malloc((size_t) n + 1, sizeof(V));
    V *degprev = (V *) md_malloc((size_t) n + 1, sizeof(V));
    V *wext = (V *) md_malloc((size_t) n + 1, sizeof(V));
    V *extlen = (V *) md_malloc((size_t) n + 1, sizeof(V));
    E *pe = (E *) md_malloc((size_t) n + 1, sizeof(E));
    V *len = (V *) md_malloc((size_t) n + 1, sizeof(V));
    V *elen = (V *) md_malloc((size_t) n + 1, sizeof(V));
    V *iw = (V *) md_malloc((size_t) iwlen + 1, sizeof(V));
    uint64_t *heads = (uint64_t *) md_malloc((size_t) ((iwlen + 63) / 64) + 1, sizeof(uint64_t));
    V *weights = NULL;
    V *dindex = NULL;
    V *dwght = NULL;
    V *dnode = NULL;
    V *ddeg = NULL;
    uint64_t *dmask = NULL;
    V ndense = 0;
    int64_t nwords = 0;

    if (marker == NULL || superwght == NULL || degree == NULL || deghead == NULL
            || degnext == NULL || degprev == NULL || wext == NULL || extlen == NULL
            || pe == NULL || len == NULL || elen == NULL || iw == NULL || heads == NULL) {
        status = MD_OUT_OF_MEMORY;
    }

    // the vertex weights, indexed from 1 (none if they are all 1)
    if (status == MD_OK && vwght != NULL) {
        weights = (V *) md_malloc((size_t) n + 1, sizeof(V));

        if (weights == NULL) {
            status = MD_OUT_OF_MEMORY;
        } else {
            for (V node = 1; node <= n; node++) {
                weights[node] = vwght[node - 1];
            }
        }
    }

    if (status == MD_OK) {

        // select the dense nodes (`dindex` and `dwght` are needed only if
        // there are any)
        V dlen = 0;
        ndense = MDF(mdndense)(n, xadj, dense);

        if (ndense > 0) {
            dlen = n;
        }

        dindex = (V *) md_malloc((size_t) dlen + 1, sizeof(V));

        if (dindex == NULL) {
            status = MD_OUT_OF_MEMORY;
        } else {
            status = MDF(mdselect)(dindex, n, xadj, adjncy, dense, &ndense);
        }

        nwords = ((int64_t) ndense + 63) / 64;
        dwght = (V *) md_malloc((size_t) dlen + 1, sizeof(V));
        dnode = (V *) md_malloc((size_t) ndense + 1, sizeof(V));
        ddeg = (V *) md_malloc((size_t) ndense + 1, sizeof(V));
        dmask = (uint64_t *) md_malloc((size_t) (nwords * (int64_t) n) + 1, sizeof(uint64_t));

        if (dwght == NULL || dnode == NULL || ddeg == NULL || dmask == NULL) {
            status = MD_OUT_OF_MEMORY;
        }
    }

    if (status == MD_OK) {
        MDF(md_impl)(marker, superwght, degree, deghead, degnext, degprev,
            wext, dwght, extlen, pe, len,
            elen, iw, (E) iwlen, heads, total, n, external, mass,
            twin, weights, xadj,
            adjncy, ndense, nwords, dindex, dnode, ddeg, dmask);

        for (V node = 1; node <= n; node++) {
            invp[node - 1] = degnext[node] - 1;
        }
    }

    free(marker); free(superwght); free(degree); free(deghead); free(degnext); free(degprev);
    free(wext); free(extlen); free(pe); free(len); free(elen); free(iw); free(heads);
    free(weights); free(dindex); free(dwght); free(dnode); free(ddeg); free(dmask);
    return status;
}

// This routine counts the nodes with more than `max(16, dense * sqrt(neqns))`
// neighbors: the candidates of `mdselect` (0 if `dense` is negative).
static V MDF(mdndense)(V neqns, const V *xadj, double dense) {
    V ndense = 0;

    if (!(dense < 0)) {
        double thresh = MD_MAX(16.0, dense * sqrt((double) neqns));

        for (V node = 1; node <= neqns; node++) {
            if ((double) (xadj[node] - xadj[node - 1]) > thresh) {
                ndense += 1;
            }
        }
    }

    return ndense;
}

// This routine selects the dense nodes: the nodes with more than
// `max(16, dense * sqrt(neqns))` neighbors. To keep the space linear, the
// masks may use at most `nnz / 2` words, so if there are too many dense
// nodes, only those of highest degree are selected (the lowest-numbered
// among nodes of equal degree). `ndense` is the number of candidates (see
// `mdndense`). Returns the number of dense nodes in `ndense` and sets
// `dindex` (empty if there are no candidates).
static int MDF(mdselect)(V *dindex, V neqns, const V *xadj, const V *adjncy, double dense,
        V *ndensep) {
    int status = MD_OK;
    V ndense = *ndensep;
    V maxdense = ndense;
    (void) adjncy;

    if (ndense > 0) {
        double thresh = MD_MAX(16.0, dense * sqrt((double) neqns));

        for (V node = 1; node <= neqns; node++) {
            dindex[node] = 0;
        }

        // at most `64 * maxwords` dense nodes, where
        // `maxwords * neqns <= nnz / 2`
        int64_t nnz = (int64_t) xadj[neqns];
        int64_t maxwords = MD_MAX((int64_t) 1, nnz / (2 * (int64_t) neqns));
        maxdense = (V) MD_MIN((int64_t) ndense, 64 * maxwords);
        md_pair *cand = (md_pair *) md_malloc((size_t) ndense, sizeof(md_pair)); int64_t k = 0;

        if (cand == NULL) {
            status = MD_OUT_OF_MEMORY;
        } else {
            for (V node = 1; node <= neqns; node++) {
                if ((double) (xadj[node] - xadj[node - 1]) > thresh) {
                    cand[k].key = (int64_t) (xadj[node - 1] - xadj[node]); cand[k].node = (int64_t) node; k += 1;
                }
            }

            if (maxdense < ndense) {
                qsort(cand, (size_t) ndense, sizeof(md_pair), md_pair_compare);
                qsort(cand, (size_t) maxdense, sizeof(md_pair), md_pair_compare_node);
            }

            for (int64_t j = 1; j <= (int64_t) maxdense; j++) {
                dindex[cand[j - 1].node] = (V) j;
            }

            free(cand);
        }
    }

    *ndensep = maxdense;
    return status;
}

// This routine copies the adjacency structure into `iw`, and computes the
// ordering (`mdmain`).
static void MDF(md_impl)(
        I *marker,
        V *superwght,
        V *degree,
        V *deghead,
        V *degnext,
        V *degprev,
        V *wext,
        V *dwght,
        V *extlen,
        E *pe,
        V *len,
        V *elen,
        V *iw,
        E iwlen,
        uint64_t *heads,
        V total,
        V neqns,
        int external,
        int mass,
        int twin,
        const V *vwght,
        const V *xadj,
        const V *adjncy,
        V ndense,
        int64_t nwords,
        V *dindex,
        V *dnode,
        V *ddeg,
        uint64_t *dmask) {
#ifdef MD_COUNT
    for (int k = 0; k < MD_NCOUNT; k++) {
        md_count[k] = 0;
    }

    md_count_h = -1;
#endif

    // copy the adjacency structure into `iw`, dropping self loops (the graph
    // is indexed from 0)
    E pfree = 1;

    if (ndense == 0) {
        for (V node = 1; node <= neqns; node++) {
            pe[node] = pfree; V deg = 0;
            E istart = (E) xadj[node - 1];
            E istop = (E) xadj[node] - 1;

            for (E i = istart; i <= istop; i++) {
                V neighbor = adjncy[i] + 1;

                if (neighbor != node) {
                    iw[pfree] = neighbor; pfree += 1;
                    deg += MDF(mdvwght)(vwght, neighbor);
                }
            }

            len[node] = (V) (pfree - pe[node]); degree[node] = deg;
        }
    } else {
        // the dense nodes are left out of `iw`: the adjacency between any
        // node and a dense node is stored in the masks
        for (int64_t i = 1; i <= nwords * (int64_t) neqns; i++) {
            dmask[i] = 0;
        }

        for (V node = 1; node <= neqns; node++) {
            V j = dindex[node];

            if (j > 0) {
                dnode[j] = node;
            }

            pe[node] = pfree; V deg = 0;
            E istart = (E) xadj[node - 1];
            E istop = (E) xadj[node] - 1;

            for (E i = istart; i <= istop; i++) {
                V neighbor = adjncy[i] + 1;

                if (neighbor != node) {
                    V k = dindex[neighbor];

                    if (k > 0) {
                        MDF(mdsetbit)(dmask, nwords, node, k);
                    } else if (!(j > 0)) {
                        iw[pfree] = neighbor; pfree += 1;
                        deg += MDF(mdvwght)(vwght, neighbor);
                    }
                }
            }

            len[node] = (V) (pfree - pe[node]); degree[node] = deg;
        }
    }

#ifdef MD_BENCH
    double start = md_clock();
#endif

    MDF(mdmain)(marker, superwght, degree, deghead, degnext, degprev,
        wext, dwght, extlen, pe, len,
        elen, iw, iwlen, heads, pfree, total, neqns, external, mass, twin,
        vwght, ndense, nwords, dindex, dnode, ddeg, dmask);

#ifdef MD_BENCH
    md_seconds = md_clock() - start;
#endif

#ifdef MD_COUNT
    // the water rises fewer than n times (2n with external degree)
    if (external) {
        assert(total == 0 || md_count[MD_C_RISE] < 2 * (int64_t) total);
    } else {
        assert(total == 0 || md_count[MD_C_RISE] < (int64_t) total);
    }
#endif
}

// This routine implements the exact minimum degree algorithm. It makes use
// of the implicit representation of elimination graphs by quotient graphs,
// and the notion of indistinguishable nodes. It also implements the
// modification by minimum external degree.
//
// The degree of a node is kept as a lower bound, which is exact for some
// nodes. Let the key of a node be its bound plus 1 (with external degree)
// or plus its weight (with degree). The algorithm repeats the following.
// If an exact node has the minimum key, it is eliminated (`mdelim`,
// `mdupdate`), and the nodes of its reachable set get new bounds.
// Otherwise every node at the minimum key is recomputed: they are grouped
// by most recent element, and each group is recomputed at once
// (`mdrecompute`).
//
// A node of weight w stands for w indistinguishable nodes, so the degree of
// a node is the weight of its closed neighborhood minus one, and its external
// degree is the weight of its neighbors outside its supernode.
//
// input parameters:
//
//   - `neqns`: number of equations
//   - `total`: sum of the vertex weights
//   - `vwght`: vertex weights (none if they are all 1)
//   - `external`: minimize the external degree instead of the degree. Nodes
//     are indistinguishable if they lie in the reachable set of a pivot and
//     have the same elements and neighbors there; elements contained in a
//     reachable set are absorbed into the new element (aggressive
//     absorption). The ordering is exact for these supernodes.
//   - `mass`: mass elimination. Each node of the reachable set whose only
//     neighbors are those of the new element (it lies in no other element,
//     has no other neighbors, and has the same dense neighbors as the pivot)
//     is eliminated with the pivot. With degree such a node is
//     indistinguishable from the pivot, so the ordering remains exact. With
//     `external` it need not be, and the ordering is no longer exact.
//   - `twin`: false-twin elimination. Each node outside the reachable set
//     whose neighbors are exactly the reachable set of the pivot, which has
//     never been in a reachable set, and whose dense neighbors (and, with
//     degree, weight) are those of the pivot, is eliminated right after the
//     pivot's supernode. With degree it is then simplicial and of minimum
//     degree, so the ordering remains exact. With `external` it need not be,
//     and the ordering is no longer exact (see `mdelim`).
//   - `iwlen`: length of `iw`
//   - `(ndense, nwords, dindex, dnode, ddeg, dmask)`: the dense nodes (see
//     `mdselect`). Their lists are empty: their adjacency, and every
//     adjacency to them, is stored in the masks. They are kept out of the
//     quotient graph until one of them has minimum degree, so the ordering
//     remains exact; with `external`, indistinguishable nodes must also have
//     the same dense neighbors.
//
// updated parameters:
//
//   - `(pe, len, elen, iw)`: quotient graph. On input, `(pe, len, iw)` holds
//     the adjacency structure, and `degree` the weight of the neighbors in
//     each list. For an uneliminated node, the list
//     `iw[pe[node]:pe[node] + len[node] - 1]` holds `elen[node]` elements,
//     sorted by elimination order, followed by its uneliminated neighbors;
//     for an eliminated node (element), it holds the nodes in the element.
//     An element is absorbed, and the list of a merged node released, by
//     setting its length to 0. Destroyed on output.
//   - `pfree`: the quotient graph is stored in `iw[1:pfree - 1]`
//
// output parameters:
//
//   - `degnext`: the minimum degree ordering (repurposed from degree list
//     links)
//
// working arrays:
//
//   - `deghead`: `deghead[key]` points to the first exact node with a key,
//     and `deghead[key + total + 1]` to the first node that is not exact, or
//     0 if there are no such nodes
//   - `degnext`: during algorithm, points to the next node in the degree list
//     (or hash bucket) associated with node, or stores negative values for
//     eliminated/merged nodes (`-(neqns + rank)` for an eliminated principal
//     node, the `rank`-th supernode eliminated, and -parent for a merged
//     node). During a visit, the next node in a group. After algorithm,
//     stores the final ordering.
//   - `degprev`: points to the previous node in a degree list associated with
//     node, or the negative of the index in `deghead` of its list (if node
//     was the first in it), or 0 if the node is not in the degree lists.
//     While indistinguishable nodes are detected, the hash bucket of node.
//     During a visit, the first node of the group of an element.
//   - `degree`: the bound on the external degree of each uneliminated
//     supernode; for an element, the weight of the nodes in it
//   - `superwght`: the weight of the supernodes. It is positive exactly for
//     the uneliminated supernodes outside the reachable set being formed: it
//     is negated for the nodes of that set and for the pivot, and 0 for
//     merged and eliminated nodes.
//   - `wext`: for an element, 0 if it is absorbed; otherwise less than `wflg`
//     if it has not been met at the current elimination, and `wflg` plus the
//     weight of its nodes outside the reachable set if it has. For a node of
//     the reachable set, the first term of its bound (see `mdelim`).
//   - `dwght`: the weight of the dense neighbors of each node (empty if there
//     are no dense nodes)
//   - hash buckets, used to detect indistinguishable nodes, share `deghead`
//     with the degree lists. They use its last entries, whose lists (not
//     exact, with keys near `total`) are almost always empty: the first node
//     of the bucket at `deghead[b]` is `-deghead[b]` if no degree list starts
//     there, and otherwise `degprev[deghead[b]]`, which is restored to `-b`
//     once the bucket is processed. The hash of a node of the reachable set
//     is kept in `marker` (see `mdelim`).
//   - `extlen`: for an element, the number of its nodes outside the element
//     being recomputed, which come first in its list (see `mdrecompute`).
//     For a node of the reachable set, `extlen` is the second term of its
//     bound, and then its upper bound (see `mdupdate`).
//   - `marker`: temporary marker vector, used with increasing tags. A node
//     of the reachable set holds its hash there, as a negative number, until
//     the indistinguishable nodes are merged
//   - `heads`: the first entries of the lists, while `iw` is compressed (see
//     `mdcompress`)
static void MDF(mdmain)(
        I *marker,
        V *superwght,
        V *degree,
        V *deghead,
        V *degnext,
        V *degprev,
        V *wext,
        V *dwght,
        V *extlen,
        E *pe,
        V *len,
        V *elen,
        V *iw,
        E iwlen,
        uint64_t *heads,
        E pfree,
        V total,
        V neqns,
        int external,
        int mass,
        int twin,
        const V *vwght,
        V ndense,
        int64_t nwords,
        V *dindex,
        V *dnode,
        V *ddeg,
        uint64_t *dmask) {
    // initialization for the minimum degree algorithm.
    for (V node = 1; node <= neqns; node++) {
        marker[node] = 0;
        superwght[node] = MDF(mdvwght)(vwght, node);
        elen[node] = 0;
        degnext[node] = 0;
        degprev[node] = 0;
        wext[node] = 1;
    }

    for (V key = 1; key <= 2 * (total + 1); key++) {
        deghead[key] = 0;
    }

    // - `mindeg` is the current minimum key
    // - `rank` counts the eliminated supernodes plus 1, and `wleft` is the
    //   weight of the uneliminated nodes
    // - `tag` is used to facilitate marking nodes; the tags are reset
    //   before they overflow
    // - `wflg` is the flag of `wext`, and `lemax` the largest weight of an
    //   element
    // - `offset` separates the keys of the lists that are not exact
    V mindeg = total + 1; V rank = 1; V wleft = total; I tag = 0;
    V wflg = 1; V lemax = 0; V offset = total + 1;

    // - `nbypass` is the number of dense nodes outside the quotient graph
    // - `nmissing` is the number of ordered pairs of such nodes that are
    //   not adjacent
    V nbypass = ndense; int64_t nmissing = 0; V dmin = 0;

    if (ndense > 0) {
        for (V node = 1; node <= neqns; node++) {
            dwght[node] = 0;
        }

        MDF(mddenseinit)(ddeg, dmask, nwords, ndense, dnode, neqns, superwght);

        for (V j = 1; j <= ndense; j++) {
            nmissing += (int64_t) ndense - 1 - MDF(mdcount)(dmask, nwords, dnode[j]);
        }
    }

    // `degree` holds the weight of the neighbors in each list (summed while
    // the adjacency was copied)
    for (V node = 1; node <= neqns; node++) {
        // dense nodes stay out of the degree lists
        if (ndense == 0 || !(dindex[node] > 0)) {
            V deg = degree[node];

            if (ndense > 0) {
                dwght[node] = MDF(mdweight)(dmask, nwords, dnode, node, superwght);
                deg += dwght[node];
            }

            degree[node] = deg;
            V key = MDF(mdkey)(deg, superwght[node], external);

            if (key > 1) {
                mindeg = MDF(mdinsert)(key, mindeg, node, 1, deghead, degnext, degprev, offset);
            } else {
                // eliminate node of degree 0
                degnext[node] = -(neqns + rank); rank += 1; wleft -= superwght[node];
                MD_ADD(MD_C_STEPS, 1); MD_WALK(key - 1);
            }
        }
    }

    // each step uses fewer than `8 neqns + 64` tags: reset them before they
    // overflow
    I tagmax = I_MAX - (I) 8 * (I) neqns - (I) 64;

    while (wleft > 0) {
        // find the pivot `mdnode`. If no node at the minimum key is exact,
        // they are all recomputed, and the search starts over
        V mdnode = 0;

        while (mdnode == 0) {
            if (tag > tagmax) {
                for (V node = 1; node <= neqns; node++) {
                    marker[node] = 0;
                }

                tag = 0;
            }

            while (mindeg <= total && !(deghead[mindeg] > 0) && !(deghead[mindeg + offset] > 0)) {
                mindeg += 1;
            }

            // if a dense node has a smaller key than every bound, it rejoins
            // the quotient graph and is eliminated. `dmin` is a lower bound on
            // the keys of the dense nodes: they are scanned only if it is
            // below the minimum key
            V dj = 0;

            if (nbypass > 0 && !(mindeg <= total && mindeg <= dmin)) {
                V dkey = V_MAX;

                for (V j = 1; j <= ndense; j++) {
                    V dnd = dnode[j];

                    if (dnd > 0) {
                        V key = MDF(mdkey)(ddeg[j], superwght[dnd], external);

                        if (key < dkey) {
                            dkey = key; dj = j;
                        }
                    }
                }

                dmin = dkey;

                if (mindeg <= total && mindeg <= dkey) {
                    dj = 0;
                }
            }

            if (dj > 0) {
                mdnode = dnode[dj];
                nmissing -= 2 * ((int64_t) nbypass - 1 - MDF(mdcount)(dmask, nwords, mdnode));
                nbypass -= 1;
                pfree = MDF(mdreinsert)(dj, neqns, pe, len, elen, iw, iwlen, pfree, heads,
                    superwght, nwords, dindex, dnode, dmask, wext);
                degree[mdnode] = ddeg[dj];

                // its neighbors now reach it through their lists
                V dwgt = superwght[mdnode];
                E pstart = pe[mdnode];
                E pstop = pstart + (E) len[mdnode] - 1;
                MD_ADD(MD_C_DENSE, pstop - pstart + 1);

                for (E p = pstart; p <= pstop; p++) {
                    dwght[iw[p]] -= dwgt;
                }
            } else if (deghead[mindeg] > 0) {
                // remove `mdnode` from the degree structure. The exact nodes
                // at the minimum key go before the others are recomputed:
                // their eliminations may touch those nodes and give them new
                // bounds
                mdnode = deghead[mindeg];
                MDF(mdremove)(mdnode, deghead, degnext, degprev, offset);

                // the new head of the list is likely the next pivot: ask for
                // its fields and its list while `mdnode` is eliminated
                V nxnode = deghead[mindeg];

                if (nxnode > 0) {
                    md_prefetch(superwght, nxnode);
                    md_prefetch(degree, nxnode); md_prefetch(len, nxnode);
                    md_prefetch(elen, nxnode); md_prefetch(iw, pe[nxnode]);
                }
            } else {
                // recompute every node at the minimum key
                MD_WALK(mindeg - 1);
                MDF(mdvisit)(mindeg, pe, len, elen, iw, superwght,
                    degree, deghead, degnext, degprev, marker,
                    &tag, neqns, extlen, dwght, nbypass, external, offset, &mindeg);
            }
        }

        MD_ADD(MD_C_STEPS, 1); MD_WALK(MDF(mdkey)(degree[mdnode], superwght[mdnode], external) - 1);
        degprev[mdnode] = 0;
        degnext[mdnode] = -(neqns + rank);

        // `mdnode` is the last supernode
        if (superwght[mdnode] >= wleft) {
            break;
        }

        // the dense neighbors of `mdnode` lose it
        V dS;

        if (nbypass > 0) {
            // a key is at least the degree plus 1
            V dlow = MDF(mdlose)(ddeg, dmask, nwords, mdnode, superwght[mdnode]);
            dS = MDF(mdweight)(dmask, nwords, dnode, mdnode, superwght);

            if (dlow < V_MAX) {
                dmin = MD_MIN(dmin, dlow + 1);
            }
        } else {
            dS = 0;
        }

        // eliminate `mdnode` and perform quotient graph
        // transformation
        V rchsze, delta, swght, fsze, fwght;
        MDF(mdelim)(mdnode,
            pe, len, elen, iw, iwlen, &pfree, heads, neqns, superwght, degree,
            deghead, degnext, degprev, marker, wext, extlen, &wflg, &lemax,
            dS, offset, external, twin, rank + 1,
            nbypass, nwords, dmask, ddeg, &rchsze, &delta, &swght, &fsze, &fwght);

        // the nodes in the reachable set, and the dense neighbors of
        // `mdnode`, become adjacent to the dense neighbors of `mdnode`
        if (nbypass > 0) {
            nmissing = MDF(mddense)(mdnode, pe[mdnode], rchsze, iw, dwght, superwght,
                nwords, dnode, ddeg, dmask, nmissing);
        }

        // merge the indistinguishable nodes of the reachable set, eliminate
        // some of them with `mdnode` (mass elimination), assign the bounds
        // of the others, and return them to the degree structure. The
        // reachable set also loses the false twins of `mdnode`. Mass
        // elimination is skipped when there are any: a node of the reachable
        // set adjacent only to the new element and to them is not
        // indistinguishable from `mdnode`, and it must follow them (they
        // are numbered right after the supernode of `mdnode`)
        V mwght;
        MDF(mdupdate)(mdnode, rchsze, delta, dS, swght, fwght, pe,
            len, elen, iw, superwght, degree, deghead, degnext,
            degprev, wext, marker, &tag,
            neqns, extlen, dwght, nwords, dmask, ddeg, nbypass, &mindeg, external,
            mass && fsze == 0, offset, &mwght);

        rank += 1 + fsze; wleft -= swght + fwght + mwght;
        superwght[mdnode] = 0;

        if (dS > 0 && (fsze > 0 || mass)) {
            dmin = 0;
        }
    }

    MDF(mdnumber)(neqns, degnext, degprev, degree);
}

// This routine eliminates the node `mdnode` of minimum degree from the
// adjacency structure, which is stored in the quotient graph format. It
// forms the reachable set of `mdnode`, absorbs the elements adjacent to
// `mdnode`, and computes, for every other element adjacent to the reachable
// set, the weight of its nodes outside the reachable set. The elements with
// none are absorbed (aggressive absorption). The lists of the nodes in the
// reachable set are purged, and the terms of their bounds are gathered:
// `wext[u]` is the weight of their uneliminated neighbors plus their largest
// external part, and `extlen[u]` the sum of their other external parts (these
// slots are unused for uneliminated nodes). The nodes are placed in hash
// buckets, and the reachable set becomes the new element `mdnode`.
//
// With `twin`, the false twins of `mdnode` are found first: the nodes
// outside the reachable set whose neighbors are exactly the reachable set.
// Each lies in every list of the reachable set, so only the shortest is
// scanned. A node of it qualifies if it has never been in a reachable set
// (its adjacency is then all in its list, and its degree is exact), its
// degree and dense neighbors (and, with degree, its weight) are those of
// `mdnode`, and its list holds no node outside the reachable set. Once
// `mdnode` is eliminated, it is simplicial and has the key of `mdnode`. It
// leaves the quotient graph, numbered right after the supernode of
// `mdnode`: its weight is set to 0 before the lists are purged, so it drops
// out of them.
//
// The reachable set is written to `iw` as it is formed, and becomes the list of
// the new element. If `mdnode` is adjacent to no element, its list of
// uneliminated neighbors holds the reachable set, which is formed there, in
// place. Otherwise it is formed at `pfree`, after a compression if the tail of
// `iw` is too short: it holds at most `degree[mdnode] - dS` nodes (each node
// weighs at least 1). That fits after a compression: outside a rejoin, the
// lists take at most `nnz` entries, and `iw` has at least `nnz + 2 neqns`
// (see `md`); a rejoining dense node has no element.
//
// input parameters:
//
//   - `mdnode`: node of minimum degree
//   - `neqns`: number of equations
//   - `iwlen`: length of `iw`
//   - `dS`: the weight of the dense neighbors of `mdnode`
//   - `offset`: separates the keys of the lists that are not exact
//   - `external`: the keys are external degrees (the false twins need not
//     have the weight of `mdnode`)
//   - `twin`: find the false twins of `mdnode`
//   - `fnum`: the rank of the first false twin (the rank after that of
//     `mdnode`)
//   - `(nbypass, nwords, dmask)`: the dense nodes
//
// updated parameters:
//
//   - `(pe, len, elen, iw)`: quotient graph
//   - `pfree`: the quotient graph is stored in `iw[1:pfree - 1]`
//   - `heads`: working array of `mdcompress`
//   - `ddeg`: the degrees of the dense nodes, which lose the false twins
//   - `superwght`: weight of the supernodes; on output, negative for
//     `mdnode`
//   - `degree`: on output, `degree[mdnode]` is the weight of the nodes in the
//     new element
//   - `(deghead, degnext, degprev)`: degree lists; on output, the
//     nodes in the reachable set are removed from them, and the false twins
//     are numbered in `degnext`
//   - `marker`: on output, the nodes in the reachable set are placed in hash
//     buckets, which share `deghead` with the degree lists (`degprev` and
//     `degnext` hold their buckets and links), and each holds its hash in
//     `marker`, as a negative number, so that it is never taken for a tag
//   - `(wext, wflg, lemax)`: external parts of the elements
//
// output parameters:
//
//   - `rchsze`: size of the reachable set, which is the list of the new
//     element `mdnode`
//   - `(wext, extlen)`: for each node of the reachable set, the terms of its
//     bound
//   - `delta`: the weight of the reachable set, with the dense neighbors of
//     `mdnode`
//   - `swght`: the weight of `mdnode`
//   - `fsze`: the number of false twins of `mdnode` (supernodes)
//   - `fwght`: their weight
static void MDF(mdelim)(
        V mdnode,
        E *pe,
        V *len,
        V *elen,
        V *iw,
        E iwlen,
        E *pfreep,
        uint64_t *heads,
        V neqns,
        V *superwght,
        V *degree,
        V *deghead,
        V *degnext,
        V *degprev,
        I *marker,
        V *wext,
        V *extlen,
        V *wflgp,
        V *lemaxp,
        V dS,
        V offset,
        int external,
        int twin,
        V fnum,
        V nbypass,
        int64_t nwords,
        uint64_t *dmask,
        V *ddeg,
        V *rchszep,
        V *deltap,
        V *swghtp,
        V *fszep,
        V *fwghtp) {
    E pfree = *pfreep; V wflg = *wflgp; V lemax = *lemaxp;

    // the pivot leaves the quotient graph
    V swght = superwght[mdnode]; superwght[mdnode] = -swght;

    // the reachable set is written to `iw[rstart:rstop]` as it is formed:
    // in place if `mdnode` is adjacent to no element, and otherwise at the
    // end of `iw` (compress it first if the tail is too short)
    int inplace = elen[mdnode] == 0;
    E rstart;

    if (inplace) {
        rstart = pe[mdnode];
    } else {
        E need = (E) (degree[mdnode] - dS);

        if (pfree + need - 1 > iwlen) {
            pfree = MDF(mdcompress)(neqns, pe, len, iw, pfree, heads);
            assert(pfree + need - 1 <= iwlen);
        }

        rstart = pfree;
    }

    // - `iw[estart:estop]` is the list of elements adjacent to `mdnode`
    // - `iw[nstart:nstop]` is the list of its uneliminated neighbors
    E estart = pe[mdnode];
    E estop = estart + (E) elen[mdnode] - 1;
    E nstart = estop + 1;
    E nstop = estart + (E) len[mdnode] - 1;

    // form the reachable set: the weight of its nodes is negated, and they
    // are removed from the degree structure. With `twin`, `rsum` is the
    // weight of the set, and `rbest` its node with the fewest neighbors in
    // its list (`lbest`)
    V rchsze = 0; V rsum = 0; V rbest = 0; V lbest = V_MAX;
    MD_ADD(MD_C_ELIM, estop - estart + 1);

    for (E i = estart; i <= estop; i++) {
        V elmnt = iw[i];

        // ask for the length and the list of the next element
        V nxelmnt = iw[MD_MIN(i + 1, estop)];
        md_prefetch(len, nxelmnt); md_prefetch(iw, pe[nxelmnt]);

        E jstart = pe[elmnt];
        E jstop = jstart + (E) len[elmnt] - 1;
        MD_ADD(MD_C_ELIM, jstop - jstart + 1);

        for (E j = jstart; j <= jstop; j++) {
            V node = iw[j]; V wght = superwght[node];

            if (wght > 0) {
                superwght[node] = -wght;
                iw[rstart + (E) rchsze] = node; rchsze += 1;
                MDF(mdremove)(node, deghead, degnext, degprev, offset);

                // the fields and the list of `node` are read only after
                // the reachable set is formed: ask for them now
                md_prefetch(len, node); md_prefetch(pe, node); md_prefetch(elen, node);
                md_prefetch(degree, node); md_prefetch(iw, pe[node]);

                if (twin) {
                    rsum += wght; V nlen = len[node] - elen[node];

                    if (nlen < lbest) {
                        lbest = nlen; rbest = node;
                    }
                }
            }
        }

        // element absorption: `elmnt` is absorbed into the new element
        len[elmnt] = 0; wext[elmnt] = 0;
    }

    MD_ADD(MD_C_ELIM, nstop - nstart + 1);

    for (E i = nstart; i <= nstop; i++) {
        V node = iw[i]; V wght = superwght[node];

        if (wght > 0) {
            // (in place, the write trails the read)
            superwght[node] = -wght;
            iw[rstart + (E) rchsze] = node; rchsze += 1;
            MDF(mdremove)(node, deghead, degnext, degprev, offset);
            md_prefetch(len, node); md_prefetch(pe, node); md_prefetch(elen, node);
            md_prefetch(degree, node); md_prefetch(iw, pe[node]);

            if (twin) {
                rsum += wght; V nlen = len[node] - elen[node];

                if (nlen < lbest) {
                    lbest = nlen; rbest = node;
                }
            }
        }
    }

    // false twins of `mdnode`: they lie in the list of `rbest`. Each is
    // numbered (from `fnum`), leaves the degree lists and its dense
    // neighbors, and gets weight 0 so that the purge below drops it
    V fsze = 0; V fwght = 0;

    if (twin && lbest < V_MAX) {
        V tdeg = rsum + dS;
        E pstart = pe[rbest] + (E) elen[rbest];
        E pstop = pe[rbest] + (E) len[rbest] - 1;
        MD_ADD(MD_C_ELIM, pstop - pstart + 1);

        for (E p = pstart; p <= pstop; p++) {
            V node = iw[p]; V wght = superwght[node];

            // cheapest tests first: the weight is in hand
            if (wght > 0 && (external || wght == swght) &&
                    degree[node] == tdeg && elen[node] == 0 &&
                    (!(nbypass > 0) || MDF(mdsame)(dmask, nwords, node, mdnode))) {
                // every neighbor of `node` must lie in the reachable set
                int same = 1;
                E qstart = pe[node];
                E qstop = qstart + (E) len[node] - 1;

                for (E q = qstart; q <= qstop; q++) {
                    MD_ADD(MD_C_ELIM, 1);

                    if (superwght[iw[q]] > 0) {
                        same = 0;
                        break;
                    }
                }

                if (same) {
                    MDF(mdremove)(node, deghead, degnext, degprev, offset);

                    if (nbypass > 0) {
                        MDF(mdlose)(ddeg, dmask, nwords, node, wght);
                    }

                    degnext[node] = -(neqns + fnum + fsze); fsze += 1;
                    superwght[node] = 0; len[node] = 0; fwght += wght;
                }
            }
        }
    }

    // the list of `mdnode` is released (its elements were absorbed above):
    // the reachable set takes its place
    E rstop = rstart + (E) rchsze - 1;
    len[mdnode] = 0; elen[mdnode] = 0;

    // advance the flag of `wext` by more than the weight of any element
    // (reset it before it overflows)
    if ((int64_t) wflg + 2 * ((int64_t) lemax + 1) >= (int64_t) V_MAX) {
        for (V node = 1; node <= neqns; node++) {
            if (wext[node] != 0) {
                wext[node] = 1;
            }
        }

        wflg = 1;
    }

    wflg += lemax + 1;

    // for each element adjacent to a node in the reachable set, compute
    // `wext[elmnt] - wflg`, the weight of the nodes of `elmnt` outside the
    // reachable set
    MD_ADD(MD_C_ELIM, rchsze);

    for (E r = rstart; r <= rstop; r++) {
        V rnode = iw[r]; V rwght = -superwght[rnode];
        E jstart = pe[rnode];
        E jstop = jstart + (E) elen[rnode] - 1;
        MD_ADD(MD_C_ELIM, jstop - jstart + 1);

        for (E j = jstart; j <= jstop; j++) {
            V elmnt = iw[j]; V we = wext[elmnt];

            if (we != 0) {
                if (we < wflg) {
                    we = wflg + degree[elmnt];
                }

                wext[elmnt] = we - rwght;
            }
        }
    }

    // for each node in the reachable set, do the following... The hash
    // buckets are many (32 to 64 per node), so that a node is seldom in a
    // bucket with another
    int64_t hashmask = MD_MIN(md_nextpow2(MD_MAX(32 * (int64_t) rchsze, (int64_t) 1)), md_prevpow2(2 * ((int64_t) offset))) - 1;
    MD_ADD(MD_C_ELIM, rchsze);

    for (E r = rstart; r <= rstop; r++) {
        V rnode = iw[r];

        // purge inactive quotient neighbors of `rnode`: absorbed elements,
        // and eliminated nodes, merged nodes, and nodes in the reachable set
        // (they are now reachable through the new element). The elements of
        // `rnode` are `iw[p1:p2]` and its neighbors are `iw[p3:p4]`; `pn` is
        // the next free position, and `pv` the position of the first
        // neighbor after the purging. Along the way, `elmax` is the largest
        // external part, `elsum` their sum, and `nwght` the weight of the
        // neighbors.
        E p1 = pe[rnode];
        E p2 = p1 + (E) elen[rnode] - 1;
        E p3 = p2 + 1;
        E p4 = p1 + (E) len[rnode] - 1;
        MD_ADD(MD_C_ELIM, p4 - p1 + 1);
        E pn = p1; uint64_t hash = md_mix((uint64_t) mdnode);
        V elmax = 0; V elsum = 0; V nwght = 0;

        for (E p = p1; p <= p2; p++) {
            V elmnt = iw[p]; V we = wext[elmnt];

            if (we != 0) {
                V dext = we - wflg;

                if (dext == 0) {
                    // aggressive absorption: `elmnt` lies in the reachable set
                    len[elmnt] = 0; wext[elmnt] = 0;
                } else {
                    iw[pn] = elmnt; pn += 1; hash += md_mix((uint64_t) elmnt);
                    elmax = MD_MAX(elmax, dext); elsum += dext;
                }
            }
        }

        E pv = pn;

        for (E p = p3; p <= p4; p++) {
            V node = iw[p]; V wght = superwght[node];

            if (wght > 0) {
                iw[pn] = node; pn += 1; hash += md_mix((uint64_t) node);
                nwght += wght;
            }
        }

        wext[rnode] = nwght + elmax; extlen[rnode] = elsum - elmax;

        // append the new element `mdnode` to the element list, which stays
        // sorted by elimination order: move the first neighbor to the end
        // (`rnode` lost at least one entry, so there is room)
        iw[pn] = iw[pv]; iw[pv] = mdnode;
        elen[rnode] = (V) (pv - p1) + 1;
        len[rnode] = (V) (pn - p1) + 1;

        // and place `rnode` in a hash bucket to detect indistinguishable
        // nodes
        marker[rnode] = -(I) (hash & (uint64_t) I_MAX) - 1;
        V bucket = 2 * offset - (V) ((hash >> 32) & (uint64_t) hashmask);
        degprev[rnode] = bucket;
        V head = deghead[bucket];

        if (head > 0) {
            // a degree list starts at `deghead[bucket]`: its first node
            // holds the bucket
            V bnode = degprev[head];

            if (bnode > 0) {
                degnext[rnode] = bnode;
            } else {
                degnext[rnode] = 0;
            }

            degprev[head] = rnode;
        } else {
            degnext[rnode] = -head;
            deghead[bucket] = -rnode;
        }
    }

    // the reachable set is the list of the new element `mdnode`: make its
    // nodes active again
    V ewght = 0;
    MD_ADD(MD_C_ELIM, rchsze);

    for (E r = rstart; r <= rstop; r++) {
        V rnode = iw[r]; V rwght = -superwght[rnode]; superwght[rnode] = rwght;
        ewght += rwght;
    }

    if (!inplace) {
        pfree = rstop + 1;
    }

    pe[mdnode] = rstart; len[mdnode] = rchsze;
    degree[mdnode] = ewght; wext[mdnode] = 1;
    *rchszep = rchsze; *pfreep = pfree; *deltap = ewght + dS; *swghtp = swght; *wflgp = wflg; *lemaxp = MD_MAX(lemax, ewght); *fszep = fsze; *fwghtp = fwght;
}

// This routine updates the degree structure after the elimination of
// `mdnode`. Indistinguishable nodes in the reachable set are detected by
// hashing their quotient neighbors, and merged into supernodes. With mass
// elimination, each supernode whose list holds only `mdnode`, and whose
// dense neighbors are those of `mdnode`, is eliminated with `mdnode`; let
// `m` be their weight (0 otherwise). Each other supernode of the reachable
// set gets the bound
//
//     delta - m - w(u) + wext[u] + (weight of its dense neighbors not adjacent to `mdnode`),
//
// raised to its old bound minus `swght + fwght + m` if it was not formed by
// a merge. It is exact if it equals the upper bound, the same expression
// with `extlen[u]` added. With external degree, the supernodes whose bound
// is less than `max(0, delta - m - swght) + w(u) - 1` are recomputed at once
// (`mdrecompute`). The supernodes are returned to the degree lists.
//
// input parameters:
//
//   - `mdnode`: the eliminated node
//   - `rchsze`: the size of the reachable set of `mdnode`, which is its list
//   - `delta`: the weight of the reachable set, with the dense neighbors of
//     `mdnode`
//   - `dS`: the weight of the dense neighbors of `mdnode`
//   - `swght`: the weight of `mdnode`
//   - `fwght`: the weight of the false twins of `mdnode` (see `mdelim`)
//   - `(pe, len, elen, iw)`: quotient graph
//   - `(wext, extlen)`: the terms of the bounds (see `mdelim`); on output,
//     `extlen[u]` is the upper bound of each node `u` of the reachable set
//   - `external`: key the degree lists by external degree instead of degree
//   - `mass`: mass elimination
//   - `offset`: separates the keys of the lists that are not exact
//
// updated parameters:
//
//   - `mindeg`: new minimum key after degree update
//   - `ddeg`: the degrees of the dense nodes, which lose the nodes eliminated
//     with `mdnode`
//   - `(superwght, degree)`: supernodes and their bounds
//   - `(deghead, degnext, degprev)`: degree lists
//   - `marker`: the hashes of the nodes of the reachable set (see `mdelim`),
//     then tags; the hash buckets share `deghead` with the degree lists
//     (empty on output)
//   - `(marker, tag)`: marker vector
//
// output parameters:
//
//   - `mwght`: the weight of the nodes eliminated with `mdnode`
static void MDF(mdupdate)(
        V mdnode,
        V rchsze,
        V delta,
        V dS,
        V swght,
        V fwght,
        E *pe,
        V *len,
        V *elen,
        V *iw,
        V *superwght,
        V *degree,
        V *deghead,
        V *degnext,
        V *degprev,
        V *wext,
        I *marker,
        I *tagp,
        V neqns,
        V *extlen,
        V *dwght,
        int64_t nwords,
        uint64_t *dmask,
        V *ddeg,
        V nbypass,
        V *mindegp,
        int external,
        int mass,
        V offset,
        V *mwghtp) {
    I tag = *tagp; V mindeg = *mindegp;

    // the supernodes formed by a merge are marked with `mtag`
    tag += 1; I mtag = tag;

    // the reachable set is `iw[rstart:rstop]`
    E rstart = pe[mdnode]; E rstop = rstart + (E) rchsze - 1;

    // for each hash bucket of the reachable set, do the following...
    MD_ADD(MD_C_ELIM, rchsze);

    for (E r = rstart; r <= rstop; r++) {
        V hash = degprev[iw[r]];
        V head = deghead[hash];

        // take the first node of the bucket, and empty it (`inode` stays 0
        // if it is empty, or was taken before)
        V inode = 0;

        if (head < 0) {
            inode = -head; deghead[hash] = 0;
        } else if (head > 0 && degprev[head] > 0) {
            inode = degprev[head]; degprev[head] = -hash;
        }

        // for each node `inode` in the bucket, do the following...
        while (inode != 0 && degnext[inode] != 0) {
            // find the first node `jnode` after `inode` whose hash value and
            // list lengths match, with `jlast` before it: only then are the
            // lists compared
            // (the hash of `inode` is read once: a merge overwrites it)
            V jlast = inode; V jnode = degnext[inode]; I ihash = marker[inode];

            while (jnode != 0 && !(marker[jnode] == ihash && len[jnode] == len[inode] && elen[jnode] == elen[inode])) {
                jlast = jnode; jnode = degnext[jnode];
            }

            if (jnode != 0) {
                // mark the quotient neighbors of `inode`
                tag += 1;
                E pstart = pe[inode];
                E pstop = pstart + (E) len[inode] - 1;

                MD_ADD(MD_C_ELIM, pstop - pstart + 1);

                for (E p = pstart; p <= pstop; p++) {
                    marker[iw[p]] = tag;
                }

                // for each node `jnode` from there on in the bucket, do the
                // following...

                while (jnode != 0) {
                    V jnext = degnext[jnode];
                    int same = marker[jnode] == ihash && len[jnode] == len[inode] && elen[jnode] == elen[inode];

                    if (same) {
                        E qstart = pe[jnode];
                        E qstop = qstart + (E) len[jnode] - 1;

                        for (E q = qstart; q <= qstop; q++) {
                            MD_ADD(MD_C_ELIM, 1);

                            if (marker[iw[q]] != tag) {
                                same = 0;
                                break;
                            }
                        }
                    }

                    if (same && nbypass > 0) {
                        same = MDF(mdsame)(dmask, nwords, inode, jnode);
                    }

                    if (same) {
                        // `jnode` is indistinguishable from `inode`: merge
                        // them into a new supernode, and remove `jnode` from
                        // the bucket
                        superwght[inode] += superwght[jnode]; superwght[jnode] = 0;
                        degnext[jnode] = -inode;
                        len[jnode] = 0; elen[jnode] = 0;
                        degnext[jlast] = jnext;
                        marker[inode] = mtag;
                    } else {
                        jlast = jnode;
                    }

                    jnode = jnext;
                }
            }

            inode = degnext[inode];
        }
    }

#ifdef MD_COUNT
    // the principal nodes of the reachable set (these reads are not counted)
    for (E r = rstart; r <= rstop; r++) {
        if (superwght[iw[r]] > 0) {
            md_count[MD_C_TOUCHED] += 1;
        }
    }
#endif

    // keep the principal nodes of the reachable set. With mass elimination,
    // a node whose list holds only `mdnode`, and whose dense neighbors are
    // those of `mdnode`, is eliminated with `mdnode` instead: with degree, it
    // is indistinguishable from `mdnode`. Its dense neighbors lose it
    V mwght = 0;

    if (mass) {
        MD_ADD(MD_C_ELIM, rchsze);

        for (E r = rstart; r <= rstop; r++) {
            V rnode = iw[r]; V rwght = superwght[rnode];

            if (rwght > 0 && len[rnode] == 1 && (!(nbypass > 0) || dwght[rnode] == dS)) {
                if (nbypass > 0) {
                    MDF(mdlose)(ddeg, dmask, nwords, rnode, rwght);
                }

                superwght[rnode] = 0; degprev[rnode] = 0; degnext[rnode] = -mdnode;
                len[rnode] = elen[rnode] = 0; mwght += rwght;
            }
        }
    }

    // the other nodes lose the nodes eliminated with `mdnode`, and its false
    // twins, as they lose `mdnode`. The floor of the eager pass is that of
    // `mdnode`, without the nodes eliminated with it
    delta -= mwght; V floor = MD_MAX((V) 0, delta - swght);
    swght += mwght + fwght;
    degree[mdnode] -= mwght;

    // for each principal node in the reachable set, assign its bound. With
    // degree, return it to the degree structure at once; with external
    // degree, chain the nodes for the eager pass through `degnext`, in order,
    // and the others in a second chain
    V zhead = 0; V ztail = 0; V ohead = 0; V otail = 0;
    MD_ADD(MD_C_ELIM, rchsze);

    for (E r = rstart; r <= rstop; r++) {
        V rnode = iw[r]; V rwght = superwght[rnode];

        if (rwght > 0) {
            degprev[rnode] = degnext[rnode] = 0;
            V lo = delta - rwght + wext[rnode];

            if (nbypass > 0) {
                lo += dwght[rnode] - dS;
            }

            V hi = lo + extlen[rnode];

            // a supernode not formed by a merge lost at most `swght` neighbors
            if (marker[rnode] != mtag) {
                lo = MD_MAX(lo, degree[rnode] - swght);
            }

            degree[rnode] = lo; extlen[rnode] = hi;
            MD_ADD(MD_C_CERTIFIED, lo == hi);

            if (!external) {
                V key = MDF(mdkey)(lo, rwght, external);
                mindeg = MDF(mdinsert)(key, mindeg, rnode, lo == hi, deghead, degnext, degprev, offset);
            } else if (lo < hi && lo < floor + rwght - 1) {
                if (ztail == 0) {
                    zhead = rnode;
                } else {
                    degnext[ztail] = rnode;
                }

                ztail = rnode;
            } else {
                if (otail == 0) {
                    ohead = rnode;
                } else {
                    degnext[otail] = rnode;
                }

                otail = rnode;
            }
        }
    }

    // the eager pass
    if (zhead > 0) {
        MD_SET_RC(MD_C_PASS_MARK);
        MDF(mdrecompute)(zhead, mdnode, pe, len, elen, iw,
            superwght, degree, deghead, degnext, degprev, marker,
            &tag, neqns, extlen, dwght, nbypass, &mindeg, external, offset);
    }

    // with external degree, return the other nodes to the degree structure
    // (the eager pass leaves their links alone: they are in no list)
    V rnode = ohead;

    while (rnode > 0) {
        V rnext = degnext[rnode];
        V key = MDF(mdkey)(degree[rnode], superwght[rnode], external);
        mindeg = MDF(mdinsert)(key, mindeg, rnode, degree[rnode] == extlen[rnode], deghead, degnext, degprev, offset);
        rnode = rnext;
    }

#ifdef MD_COUNT
    // the end of the elimination: the walk is at `max(0, min(phi, l))`, where
    // `phi` is the floor of the eager pass (`max(0, delta - 1)` with degree)
    // and `l` the least key minus 1 in the reachable set, now (these reads are
    // not counted)
    int64_t phi = (int64_t) floor;

    if (!external) {
        phi = MD_MAX((int64_t) 0, (int64_t) delta - 1);
    }

    int64_t lmin = phi;

    for (E r = rstart; r <= rstop; r++) {
        V rnode = iw[r];

        if (superwght[rnode] > 0) {
            lmin = MD_MIN(lmin, (int64_t) MDF(mdkey)(degree[rnode], superwght[rnode], external) - 1);
        }
    }

    md_count_h = MD_MAX((int64_t) 0, lmin);
#endif

    *mindegp = mindeg; *tagp = tag; *mwghtp = mwght;
}

// This routine recomputes every node with key `key` that is not exact. The
// nodes are grouped by their most recent element, and each group is
// recomputed at once (`mdrecompute`). Returns the new minimum key and the
// tag.
//
// working arrays:
//
//   - `degprev`: for an element, the first node of its group; for the node
//     of a group taken first from the list, the corresponding node of the
//     next group (the groups are chained in order)
//   - `degnext`: for a node taken from the list, the next node in its group
static void MDF(mdvisit)(
        V key,
        E *pe,
        V *len,
        V *elen,
        V *iw,
        V *superwght,
        V *degree,
        V *deghead,
        V *degnext,
        V *degprev,
        I *marker,
        I *tagp,
        V neqns,
        V *extlen,
        V *dwght,
        V nbypass,
        int external,
        V offset,
        V *mindegp) {
    I tag = *tagp;

    // take the whole list, and group it by most recent element. The
    // groups are chained, in order, from `ghead` through `degprev` of the
    // node of each group taken first
    tag += 1; I gtag = tag; V mindeg = key; V ghead = 0; V gtail = 0;
    V node = deghead[key + offset]; deghead[key + offset] = 0;

    while (node > 0) {
        V nxnode = degnext[node];

        // ask for what the next node needs here, and for what this one
        // needs in its recomputation
        md_prefetch(pe, nxnode); md_prefetch(elen, nxnode);
        md_prefetch(len, node); md_prefetch(degree, node);

        degnext[node] = degprev[node] = 0;
        V elmnt = iw[pe[node] + (E) elen[node] - 1];
        MD_ADD(MD_C_WATER_VAR, 1);

        if (marker[elmnt] != gtag) {
            marker[elmnt] = gtag; degprev[elmnt] = 0;

            if (gtail == 0) {
                ghead = node;
            } else {
                degprev[gtail] = node;
            }

            gtail = node;
        }

        degnext[node] = degprev[elmnt]; degprev[elmnt] = node;
        node = nxnode;
    }

    // recompute each group (read the link first: the recomputation
    // inserts the node into a degree list)
    V gnode = ghead;

    while (gnode > 0) {
        V gnext = degprev[gnode];
        V elmnt = iw[pe[gnode] + (E) elen[gnode] - 1];
        MD_ADD(MD_C_WATER_VAR, 1); MD_ADD(MD_C_GROUPS, 1); MD_SET_RC(MD_C_WATER_MARK);

        MDF(mdrecompute)(degprev[elmnt], elmnt, pe, len, elen, iw,
            superwght, degree, deghead, degnext, degprev, marker,
            &tag, neqns, extlen, dwght, nbypass, &mindeg, external, offset);

        gnode = gnext;
    }

    *mindegp = mindeg; *tagp = tag;
}

// This routine recomputes the external degrees of the nodes of the chain that
// starts at `zhead` and is linked by `degnext`, whose most recent element is
// `elmnt`. The nodes of `elmnt` are marked, the list of each other element of
// the nodes is scanned once and partitioned in place, its nodes outside
// `elmnt` first, and the external degree of each node is the weight of
// `elmnt`, of the union of these nodes and of its uneliminated neighbors, and
// of its dense neighbors, minus its own weight. The nodes become exact and are
// inserted into the degree lists. Returns the new minimum key and the tag.
//
// The elements in the lists of all the nodes (the core) are united once, as
// they are scanned, so each node only adds its other elements. This keeps
// groups whose nodes share many overlapping elements cheap.
//
// working arrays:
//
//   - `extlen`: the list of each element `f` is partitioned: its `extlen[f]`
//     nodes outside `elmnt` come first, in order, then its nodes in `elmnt`.
//     Merged and eliminated nodes are dropped, so `len[f]` can shrink.
//   - `marker`: the nodes of `elmnt` are marked `etag`; while the elements
//     are counted, `marker[f] - base` is the number of nodes that list `f`;
//     the core and the nodes of its lists are marked `ctag`, the other
//     elements `stag`, and the union of each node `utag`.
static void MDF(mdrecompute)(
        V zhead,
        V elmnt,
        E *pe,
        V *len,
        V *elen,
        V *iw,
        V *superwght,
        V *degree,
        V *deghead,
        V *degnext,
        V *degprev,
        I *marker,
        I *tagp,
        V neqns,
        V *extlen,
        V *dwght,
        V nbypass,
        V *mindegp,
        int external,
        V offset) {
    I tag = *tagp; V mindeg = *mindegp;
    (void) neqns;

    // mark the nodes of `elmnt`
    tag += 1; I etag = tag;
    E estart = pe[elmnt];
    E estop = estart + (E) len[elmnt] - 1;
    MD_ADD_RC(0, estop - estart + 1);

    for (E j = estart; j <= estop; j++) {
        V node = iw[j];

        if (superwght[node] > 0) {
            marker[node] = etag;
        }
    }

    // count the nodes of the group
    int64_t gsize = 0; V znode = zhead;

    while (znode > 0) {
        gsize += 1; znode = degnext[znode];
    }

#ifdef MD_COUNT
    if (md_count_rc == MD_C_PASS_MARK) {
        md_count[MD_C_RECOMP_PASS] += gsize;
    } else {
        md_count[MD_C_RECOMP_WATER] += gsize;
    }
#endif

    // mark the other elements of the nodes with the number of nodes that list
    // them, `marker[f] - base` (elements and nodes are marked in the same
    // vector: an element is never an uneliminated node)
    tag += 1; I base = tag; tag += (I) gsize;

    znode = zhead;

    while (znode > 0) {
        V znext = degnext[znode];
        E p1 = pe[znode];
        E p2 = p1 + (E) elen[znode] - 1;
        MD_ADD_RC(3, p2 - p1 + 1);

        for (E p = p1; p <= p2; p++) {
            V f = iw[p];

            if (f != elmnt) {
                if (marker[f] <= base) {
                    marker[f] = base + 1;
                } else {
                    marker[f] += 1;
                }
            }
        }

        znode = znext;
    }

    // scan each other element of the nodes once, and partition its list. The
    // elements listed by every node of the group (the core) are marked
    // `ctag`, and the union of their nodes outside `elmnt`, also marked
    // `ctag`, is taken once here: its weight is `cwght`
    tag += 1; I stag = tag; tag += 1; I ctag = tag;
    I full = base + (I) gsize; V cwght = 0;

    znode = zhead;

    while (znode > 0) {
        V znext = degnext[znode];
        E p1 = pe[znode];
        E p2 = p1 + (E) elen[znode] - 1;
        MD_ADD_RC(3, p2 - p1 + 1);

        for (E p = p1; p <= p2; p++) {
            V f = iw[p];

            if (f != elmnt && base < marker[f] && marker[f] <= full) {
                int core = marker[f] == full;

                if (core) {
                    marker[f] = ctag;
                } else {
                    marker[f] = stag;
                }

                // the nodes outside `elmnt` go to the front of the list
                // (`pout` is the next position for them, and `pin` for the
                // nodes in `elmnt`), and merged and eliminated nodes are
                // dropped
                E jstart = pe[f];
                E jstop = jstart + (E) len[f] - 1;
                E pout = jstart; E pin = jstart;
                MD_ADD_RC(1, jstop - jstart + 1);

                for (E j = jstart; j <= jstop; j++) {
                    V node = iw[j]; V wght = superwght[node];

                    if (wght > 0) {
                        if (marker[node] != etag) {
                            iw[pin] = iw[pout]; iw[pout] = node;
                            pout += 1; pin += 1;

                            if (core && marker[node] != ctag) {
                                marker[node] = ctag; cwght += wght;
                            }
                        } else {
                            iw[pin] = node; pin += 1;
                        }
                    }
                }

                len[f] = (V) (pin - jstart); extlen[f] = (V) (pout - jstart);
            }
        }

        znode = znext;
    }

    // for each node, take the union of the nodes outside `elmnt` of its
    // elements, and of its uneliminated neighbors
    V ewght = degree[elmnt];

    znode = zhead;

    while (znode > 0) {
        V znext = degnext[znode];
        tag += 1; I utag = tag; V deg = cwght;
        E p1 = pe[znode];
        E p2 = p1 + (E) elen[znode] - 1;
        E p4 = p1 + (E) len[znode] - 1;
        MD_ADD_RC(3, p4 - p1 + 1);

        for (E p = p1; p <= p2; p++) {
            V f = iw[p];

            if (f != elmnt && marker[f] != ctag) {
                E qstart = pe[f];
                E qstop = qstart + (E) extlen[f] - 1;
                MD_ADD_RC(2, qstop - qstart + 1);

                for (E q = qstart; q <= qstop; q++) {
                    V node = iw[q];

                    if (marker[node] != utag && marker[node] != ctag) {
                        marker[node] = utag; deg += superwght[node];
                    }
                }
            }
        }

        for (E p = p2 + 1; p <= p4; p++) {
            V node = iw[p]; V wght = superwght[node];

            if (wght > 0 && marker[node] != etag && marker[node] != utag && marker[node] != ctag) {
                marker[node] = utag; deg += wght;
            }
        }

        deg += ewght - superwght[znode];

        if (nbypass > 0) {
            deg += dwght[znode];
        }

        degree[znode] = deg;
        V key = MDF(mdkey)(deg, superwght[znode], external);
        mindeg = MDF(mdinsert)(key, mindeg, znode, 1, deghead, degnext, degprev, offset);
        znode = znext;
    }

    *mindegp = mindeg; *tagp = tag;
}

// This routine removes `node` from the degree lists.
static void MDF(mdremove)(
        V node,
        V *deghead,
        V *degnext,
        V *degprev,
        V offset) {
    V pvnode = degprev[node];
    (void) offset;

    if (pvnode != 0) {
        V nxnode = degnext[node];

        if (nxnode > 0) {
            degprev[nxnode] = pvnode;
        }

        if (pvnode > 0) {
            degnext[pvnode] = nxnode;
        } else {
            deghead[-pvnode] = nxnode;
        }

        degprev[node] = 0;
    }
}

// This routine inserts `node` into the exact (`isexact`) or other degree
// list with key `key`, and returns the new minimum key.
static V MDF(mdinsert)(
        V key,
        V mindeg,
        V node,
        int isexact,
        V *deghead,
        V *degnext,
        V *degprev,
        V offset) {
    V index = key;

    if (!isexact) {
        index += offset;
    }

    V firstnode = deghead[index];
    deghead[index] = node;
    degnext[node] = firstnode;
    degprev[node] = -index;

    if (firstnode > 0) {
        degprev[firstnode] = node;
    }

    return MD_MIN(key, mindeg);
}

// the key of a node: its external degree plus 1, or its degree plus 1
static V MDF(mdkey)(V deg, V wght, int external) {
    V key;

    if (external) {
        key = deg + 1;
    } else {
        key = deg + wght;
    }

    return key;
}

// the weight of `node`: `vwght[node]`, or 1 if there are no weights
static V MDF(mdvwght)(const V *vwght, V node) {
    V wght = 1;

    if (vwght != NULL) {
        wght = vwght[node];
    }

    return wght;
}

// ---------------------------------------------------------------------
// dense nodes
// ---------------------------------------------------------------------
//
// The dense nodes are kept out of the quotient graph. Fill never passes
// through an uneliminated node, so the quotient graph of the remaining nodes
// represents their elimination graph exactly. The adjacency between any node
// and the dense nodes is kept in bit masks: when a node is eliminated, its
// neighbors become a clique, so every node in its reachable set (and every
// dense neighbor) acquires its dense neighbors. A dense node rejoins the
// quotient graph when its degree is minimum, and is then eliminated.
//
//   - `ndense`: number of dense nodes
//   - `nwords`: number of 64-bit words in a mask
//   - `dindex`: `dindex[node]` is the index of a dense node, or 0
//   - `dnode`: `dnode[j]` is the dense node with index `j`; it is negated
//     once the node rejoins the quotient graph
//   - `ddeg`: `ddeg[j]` is the weighted external degree of dense node `j`
//   - `dmask`: `dmask[(node - 1) * nwords + w]` is word `w` of the mask of
//     the dense neighbors of `node`

static int64_t MDF(mdoffset)(int64_t nwords, V node) {
    return ((int64_t) node - 1) * nwords;
}

static void MDF(mdsetbit)(uint64_t *dmask, int64_t nwords, V node, V j) {
    int64_t jj = (int64_t) j - 1;
    dmask[MDF(mdoffset)(nwords, node) + (jj >> 6) + 1] |= (uint64_t) 1 << (jj & 63);
}

// number of dense neighbors of `node`
static int64_t MDF(mdcount)(const uint64_t *dmask, int64_t nwords, V node) {
    int64_t off = MDF(mdoffset)(nwords, node); int64_t cnt = 0;

    for (int64_t w = 1; w <= nwords; w++) {
        cnt += md_popcount(dmask[off + w]);
    }

    return cnt;
}

// weight of the dense neighbors of `node`
static V MDF(mdweight)(const uint64_t *dmask, int64_t nwords, const V *dnode, V node, const V *superwght) {
    int64_t off = MDF(mdoffset)(nwords, node); V wght = 0;

    for (int64_t w = 1; w <= nwords; w++) {
        uint64_t word = dmask[off + w];

        while (word != 0) {
            int64_t j = ((w - 1) << 6) + md_ctz(word) + 1;
            wght += superwght[MD_ABS(dnode[j])];
            word &= word - 1;
        }
    }

    return wght;
}

// do `inode` and `jnode` have the same dense neighbors?
static int MDF(mdsame)(const uint64_t *dmask, int64_t nwords, V inode, V jnode) {
    int64_t ioff = MDF(mdoffset)(nwords, inode); int64_t joff = MDF(mdoffset)(nwords, jnode);
    int same = 1;

    for (int64_t w = 1; w <= nwords; w++) {
        if (dmask[ioff + w] != dmask[joff + w]) {
            same = 0;
            break;
        }
    }

    return same;
}

// This routine computes the degrees of the dense nodes. On input, the masks
// hold the original adjacency to the dense nodes.
static void MDF(mddenseinit)(V *ddeg, const uint64_t *dmask, int64_t nwords, V ndense, const V *dnode, V neqns, const V *vwght) {
    (void) dnode;

    for (V j = 1; j <= ndense; j++) {
        ddeg[j] = 0;
    }

    for (V node = 1; node <= neqns; node++) {
        int64_t off = MDF(mdoffset)(nwords, node); V wght = vwght[node];

        for (int64_t w = 1; w <= nwords; w++) {
            uint64_t word = dmask[off + w];

            while (word != 0) {
                int64_t j = ((w - 1) << 6) + md_ctz(word) + 1;
                ddeg[j] += wght;
                word &= word - 1;
            }
        }
    }
}

// This routine removes `node`, which is about to be eliminated, from the
// degrees of its dense neighbors. Returns the lowest of their new degrees
// (`V_MAX` if it has none).
static V MDF(mdlose)(V *ddeg, const uint64_t *dmask, int64_t nwords, V node, V wght) {
    int64_t off = MDF(mdoffset)(nwords, node); V dlow = V_MAX;

    for (int64_t w = 1; w <= nwords; w++) {
        uint64_t word = dmask[off + w];

        while (word != 0) {
            int64_t j = ((w - 1) << 6) + md_ctz(word) + 1;
            ddeg[j] -= wght; dlow = MD_MIN(dlow, ddeg[j]);
            word &= word - 1;
        }
    }

    return dlow;
}

// This routine updates the adjacency to the dense nodes after the elimination
// of `mdnode`: its neighbors become a clique, so every node in the reachable
// set, and every dense neighbor of `mdnode`, becomes adjacent to every dense
// neighbor of `mdnode`. The reachable set is `iw[rstart:rstart + rchsze - 1]`.
// Each new adjacency is added to the dense weight
// `dwght` of the node, or to the degree `ddeg` of the dense node. Returns the
// updated number `nmissing` of ordered pairs of nonadjacent dense nodes.
static int64_t MDF(mddense)(
        V mdnode,
        E rstart,
        V rchsze,
        const V *iw,
        V *dwght,
        const V *superwght,
        int64_t nwords,
        const V *dnode,
        V *ddeg,
        uint64_t *dmask,
        int64_t nmissing) {
    int64_t moff = MDF(mdoffset)(nwords, mdnode);

    // skip if `mdnode` has no dense neighbors
    int empty = 1;

    for (int64_t w = 1; w <= nwords; w++) {
        if (dmask[moff + w] != 0) {
            empty = 0;
            break;
        }
    }

    if (!empty) {
        MD_ADD(MD_C_DENSE, rchsze);

        for (E r = rstart; r <= rstart + (E) rchsze - 1; r++) {
            V rnode = iw[r];
            int64_t roff = MDF(mdoffset)(nwords, rnode); V rdeg = 0; V rwght = superwght[rnode];

            for (int64_t w = 1; w <= nwords; w++) {
                uint64_t mword = dmask[moff + w]; uint64_t rword = dmask[roff + w];
                uint64_t added = mword & ~rword;

                while (added != 0) {
                    int64_t j = ((w - 1) << 6) + md_ctz(added) + 1;
                    rdeg += superwght[dnode[j]]; ddeg[j] += rwght;
                    added &= added - 1;
                }

                dmask[roff + w] = mword | rword;
            }

            dwght[rnode] += rdeg;
        }

        // the dense neighbors of `mdnode` become pairwise adjacent
        if (nmissing > 0) {
            for (int64_t w = 1; w <= nwords; w++) {
                uint64_t word = dmask[moff + w];

                while (word != 0) {
                    int64_t j = ((w - 1) << 6) + md_ctz(word) + 1;
                    word &= word - 1;
                    V dj = dnode[j]; int64_t joff = MDF(mdoffset)(nwords, dj); V jdeg = 0;

                    for (int64_t v = 1; v <= nwords; v++) {
                        uint64_t mword = dmask[moff + v];

                        if (v == w) {
                            // a node is not its own neighbor
                            mword &= ~((uint64_t) 1 << ((j - 1) & 63));
                        }

                        uint64_t added = mword & ~dmask[joff + v];

                        while (added != 0) {
                            int64_t l = ((v - 1) << 6) + md_ctz(added) + 1;
                            jdeg += superwght[dnode[l]]; nmissing -= 1;
                            added &= added - 1;
                        }

                        dmask[joff + v] |= mword;
                    }

                    ddeg[j] += jdeg;
                }
            }
        }
    }

    return nmissing;
}

// This routine returns the dense node with index `j` to the quotient graph.
// Its neighbors are the uneliminated principal nodes, other than dense nodes,
// whose masks contain `j`; they become its list, and it is appended to each
// of their lists (which are moved to the end of `iw` to make room). The node
// is then removed from every mask. Returns the updated `pfree`.
//
// The lists fit after a compression. Outside a rejoin, each list is no
// longer than in the quotient graph that keeps the dense nodes, so the lists
// take at most `nnz` entries. Until the dense node is eliminated, the rejoin
// adds fewer than `2 neqns`: its list and one entry per neighbor, or, while a
// list is copied, the entries appended so far and that list.
//
// working arrays:
//
//   - `wext`: chains the neighbors of the dense node, in order (it is not in
//     use for uneliminated nodes between eliminations)
//   - `heads`: working array of `mdcompress`
static E MDF(mdreinsert)(
        V j,
        V neqns,
        E *pe,
        V *len,
        V *elen,
        V *iw,
        E iwlen,
        E pfree,
        uint64_t *heads,
        const V *superwght,
        int64_t nwords,
        const V *dindex,
        V *dnode,
        uint64_t *dmask,
        V *wext) {
    V dnd = dnode[j]; dnode[j] = -dnd;
    int64_t w = (((int64_t) j - 1) >> 6) + 1;
    uint64_t bit = (uint64_t) 1 << (((int64_t) j - 1) & 63);
    V cnt = 0; V head = 0; V tail = 0;

    for (V node = 1; node <= neqns; node++) {
        int64_t i = MDF(mdoffset)(nwords, node) + w;

        if ((dmask[i] & bit) != 0) {
            dmask[i] &= ~bit;

            if (!(dindex[node] > 0) && superwght[node] > 0) {
                // move the list of `node` to the end of `iw` and append
                // `dnd` to it
                V nlen = len[node];

                if (pfree + (E) nlen > iwlen) {
                    pfree = MDF(mdcompress)(neqns, pe, len, iw, pfree, heads);
                    assert(pfree + (E) nlen <= iwlen);
                }

                E pstart = pe[node]; pe[node] = pfree;
                MD_ADD(MD_C_DENSE, nlen);

                for (E p = pstart; p <= pstart + (E) nlen - 1; p++) {
                    iw[pfree] = iw[p]; pfree += 1;
                }

                iw[pfree] = dnd; pfree += 1;
                len[node] = nlen + 1;
                cnt += 1;

                if (tail == 0) {
                    head = node;
                } else {
                    wext[tail] = node;
                }

                tail = node;
            }
        }
    }

    if (pfree + (E) cnt - 1 > iwlen) {
        pfree = MDF(mdcompress)(neqns, pe, len, iw, pfree, heads);
        assert(pfree + (E) cnt - 1 <= iwlen);
    }

    pe[dnd] = pfree;

    V node = head;

    MD_ADD(MD_C_DENSE, cnt);

    for (V k = 1; k <= cnt; k++) {
        iw[pfree] = node; pfree += 1; node = wext[node];
    }

    len[dnd] = cnt; elen[dnd] = 0;
    return pfree;
}

// This routine compresses the quotient graph storage `iw`, releasing the
// space used by absorbed elements, merged nodes, and purged entries. The
// first entry of each list is marked in a bitmap, so the released space is
// skipped 64 entries at a time. Returns the new value of `pfree`.
//
// input parameters:
//
//   - `neqns`: number of equations
//   - `len`: lengths of the lists (0 for released lists)
//
// updated parameters:
//
//   - `(pe, iw)`: quotient graph storage
//   - `pfree`: the quotient graph is stored in `iw[1:pfree - 1]`
//
// working arrays:
//
//   - `heads`: the bitmap, of at least `cld(pfree - 1, 64)` words
static E MDF(mdcompress)(
        V neqns,
        E *pe,
        const V *len,
        V *iw,
        E pfree,
        uint64_t *heads) {
    int64_t nheads = ((int64_t) pfree - 1 + 63) / 64;

    for (int64_t w = 1; w <= nheads; w++) {
        heads[w] = 0;
    }

    // store `-node` in the first entry of the list of each node, saving the
    // entry in `pe[node]`, and mark the entry in `heads`
    for (V node = 1; node <= neqns; node++) {
        if (len[node] > 0) {
            E p = pe[node];
            pe[node] = (E) iw[p];
            iw[p] = -node;
            int64_t i = (int64_t) p - 1;
            heads[(i >> 6) + 1] |= (uint64_t) 1 << (i & 63);
        }
    }

    // move the lists to the front of `iw`, in order
    E pdst = 1;

    for (int64_t w = 1; w <= nheads; w++) {
        uint64_t bits = heads[w];

        while (bits != 0) {
            E psrc = (E) (((w - 1) << 6) + md_ctz(bits) + 1);
            bits &= bits - 1;
            V node = -iw[psrc]; psrc += 1;
            MD_ADD(MD_C_COMPRESS, len[node]);
            iw[pdst] = (V) pe[node]; pe[node] = pdst; pdst += 1;

            for (V k = 2; k <= len[node]; k++) {
                iw[pdst] = iw[psrc]; pdst += 1; psrc += 1;
            }
        }
    }

    return pdst;
}

// This routine computes the ordering. On input, `degnext[node]` is
// `-(neqns + rank)` for a principal node, the `rank`-th supernode eliminated,
// and `-parent` for a node merged into another (or eliminated with it). The
// nodes of each supernode are counted from these trees, and the supernodes laid
// out in the order of their ranks: the principal first, then its other nodes in
// increasing order. On output, `degnext[node]` is the position of `node`.
//
// working arrays:
//
//   - `size`: for a principal, the number of nodes in its supernode, then the
//     next position in it
//   - `start`: for a rank, the size of its supernode, then its first position
static void MDF(mdnumber)(V neqns, V *degnext, V *size, V *start) {
    for (V node = 1; node <= neqns; node++) {
        size[node] = 0;
    }

    // count the nodes of each supernode, pointing every merged node at its
    // principal (path compression)
    V nrank = 0;

    for (V node = 1; node <= neqns; node++) {
        V x = -degnext[node];

        if (x > neqns) {
            size[node] += 1; nrank = MD_MAX(nrank, x - neqns);
        } else {
            V root = x;

            while (-degnext[root] <= neqns) {
                root = -degnext[root];
            }

            V j = node;

            while (j != root) {
                V k = -degnext[j]; degnext[j] = -root; j = k;
            }

            size[root] += 1;
        }
    }

    // the first position of each rank
    for (V node = 1; node <= neqns; node++) {
        V x = -degnext[node];

        if (x > neqns) {
            start[x - neqns] = size[node];
        }
    }

    V pos = 1;

    for (V r = 1; r <= nrank; r++) {
        V sz = start[r]; start[r] = pos; pos += sz;
    }

    // the principals, then the other nodes in increasing order
    for (V node = 1; node <= neqns; node++) {
        V x = -degnext[node];

        if (x > neqns) {
            size[node] = degnext[node] = start[x - neqns];
        }
    }

    for (V node = 1; node <= neqns; node++) {
        V x = -degnext[node];

        if (x > 0) {
            degnext[node] = size[x] += 1;
        }
    }
}
