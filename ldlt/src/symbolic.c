/* Symbolic phase: ordering -> permute -> etree -> postorder -> repermute ->
   L-pattern -> fundamental supernodes -> per-supernode row lists, update DAG,
   panel offsets, flop estimate. */

#include "internal.h"
#include <limits.h>
#include <stdio.h>

static int cmp_i32(const void *a, const void *b) {
    int32_t x=*(const int32_t*)a, y=*(const int32_t*)b;
    return (x>y)-(x<y);
}

static int64_t super_panel_rows(int32_t fc, int32_t end,
                                const int32_t *Lp, const int32_t *Li,
                                int32_t *mark, int32_t tag,
                                int32_t *rows_out)
{
    int32_t nrows = 0;
    for (int32_t c = fc; c < end; ++c) {
        for (int32_t p = Lp[c]; p < Lp[c+1]; ++p) {
            int32_t r = Li[p];
            if (r < end) continue;
            if (mark[r] == tag) continue;
            mark[r] = tag;
            if (rows_out) rows_out[nrows] = r;
            nrows++;
        }
    }
    if (rows_out && nrows > 1) qsort(rows_out, (size_t)nrows, sizeof(int32_t), cmp_i32);
    return nrows;
}

static int64_t super_panel_size(int32_t fc, int32_t end,
                                const int32_t *Lp, const int32_t *Li,
                                int32_t *mark, int32_t tag)
{
    int64_t k = (int64_t)(end - fc);
    int64_t m = super_panel_rows(fc, end, Lp, Li, mark, tag, NULL);
    return (k + m) * k;
}

static int64_t super_panel_overlap(int32_t a0, int32_t a1, int32_t b0, int32_t b1,
                                   const int32_t *Lp, const int32_t *Li,
                                   int32_t *mark, int32_t tag)
{
    int64_t overlap = 0;
    for (int32_t c = a0; c < a1; ++c) {
        for (int32_t p = Lp[c]; p < Lp[c+1]; ++p) {
            int32_t r = Li[p];
            if (r < a1) continue;
            mark[r] = tag;
        }
    }
    for (int32_t c = b0; c < b1; ++c) {
        for (int32_t p = Lp[c]; p < Lp[c+1]; ++p) {
            int32_t r = Li[p];
            if (r < b1) continue;
            if (mark[r] == tag) {
                mark[r] = tag + 1;
                overlap++;
            }
        }
    }
    return overlap;
}

static int adaptive_relax_accept(int relax,
                                 int32_t cur_start, int32_t cur_end,
                                 int32_t next_start, int32_t next_end,
                                 int64_t cur_size, int64_t next_size,
                                 int64_t merged_size, int64_t cur_fund_size,
                                 int64_t next_fund_size,
                                 int64_t cur_rows, int64_t next_rows,
                                 int64_t overlap,
                                 const int32_t *parent, const int32_t *nchild)
{
    (void)cur_rows;
    (void)next_rows;
    int64_t added = merged_size - cur_size - next_size;
    int64_t width = (int64_t)(next_end - cur_start);
    int64_t fund_size = cur_fund_size + next_fund_size;
    int64_t merged_rows = merged_size / (width > 0 ? width : 1);
    int32_t boundary = cur_end - 1;

    int path_like = (width <= 16 && merged_rows <= width + 4 &&
                     parent[boundary] == next_start && nchild[next_start] <= 1);

    int64_t fill_budget = (int64_t)relax * width;
    int64_t growth_num = (int64_t)(relax + 1);
    int64_t growth_den = 1;

    if (path_like) {
        fill_budget = width / 2;
        growth_num = 2;
        growth_den = 1;
    } else {
        (void)overlap;
        if (merged_rows > 4096) {
            fill_budget /= 2;
            growth_num = (growth_num + 1) / 2;
        }
    }

    if (added > fill_budget) return 0;
    if (merged_size * growth_den > growth_num * fund_size) return 0;
    return 1;
}

/* Symbolic L pattern. For each row i, ereach(i) gives the cols k <= i with L[i,k]!=0.
   We collect those into per-column row lists: for each k in ereach(i), append i to col[k]. */
static ldlt_status build_L_pattern(int32_t n,
    const int32_t *Bp, const int32_t *Bi, const int32_t *parent,
    int32_t **Lp_out, int32_t **Li_out)
{
    /* Build upper-incidence: per-row i, the cols k<i with A[i,k]!=0 (= lower entries (i,k)). */
    int32_t *Up = (int32_t*)ldlt_xcalloc((size_t)n+1, sizeof(int32_t));
    for (int32_t k = 0; k < n; ++k)
        for (int32_t p = Bp[k]; p < Bp[k+1]; ++p) {
            int32_t i = Bi[p];
            if (i > k) {
                if (Up[i+1] == INT32_MAX) { free(Up); return LDLT_ERR_NOMEM; }
                Up[i+1]++;
            }
        }
    int64_t up_prefix = 0;
    for (int32_t i = 0; i < n; ++i) {
        up_prefix += Up[i+1];
        if (up_prefix > INT32_MAX) { free(Up); return LDLT_ERR_NOMEM; }
        Up[i+1] = (int32_t)up_prefix;
    }
    int32_t *Ui = (int32_t*)ldlt_xmalloc((size_t)(Up[n]>0?Up[n]:1) * sizeof(int32_t));
    int32_t *cur = (int32_t*)ldlt_xmalloc((size_t)n * sizeof(int32_t));
    memcpy(cur, Up, (size_t)n * sizeof(int32_t));
    for (int32_t k = 0; k < n; ++k)
        for (int32_t p = Bp[k]; p < Bp[k+1]; ++p) {
            int32_t i = Bi[p];
            if (i > k) Ui[cur[i]++] = k;
        }
    free(cur);

    int32_t *mark = (int32_t*)ldlt_xmalloc((size_t)n * sizeof(int32_t));
    int32_t *stk  = (int32_t*)ldlt_xmalloc((size_t)n * sizeof(int32_t));
    int32_t *colcount = (int32_t*)ldlt_xcalloc((size_t)n, sizeof(int32_t));
    for (int32_t i = 0; i < n; ++i) mark[i] = -1;

    /* Pass 1: count column lengths.
       For each row i, run ereach to get cols k < i with L[i,k] != 0. Each such k contributes
       "row i in column k's pattern". Plus diagonal at column i itself. */
    for (int32_t i = 0; i < n; ++i) colcount[i] += 1; /* diagonal */
    for (int32_t i = 0; i < n; ++i) {
        mark[i] = i;
        for (int32_t p = Up[i]; p < Up[i+1]; ++p) {
            int32_t s = Ui[p];
            while (s != -1 && mark[s] != i) {
                mark[s] = i;
                if (colcount[s] == INT32_MAX) {
                    free(colcount); free(mark); free(stk); free(Up); free(Ui);
                    return LDLT_ERR_NOMEM;
                }
                colcount[s] += 1; /* row i goes into col s */
                s = parent[s];
            }
        }
    }
    int32_t *Lp = (int32_t*)ldlt_xcalloc((size_t)n+1, sizeof(int32_t));
    int64_t lp_prefix = 0;
    for (int32_t i = 0; i < n; ++i) {
        lp_prefix += colcount[i];
        if (lp_prefix > INT32_MAX) {
            free(Lp); free(colcount); free(mark); free(stk); free(Up); free(Ui);
            return LDLT_ERR_NOMEM;
        }
        Lp[i+1] = (int32_t)lp_prefix;
    }
    int32_t *Li = (int32_t*)ldlt_xmalloc((size_t)(Lp[n]>0?Lp[n]:1) * sizeof(int32_t));

    /* Pass 2: fill. We can place rows in increasing-row order naturally by iterating
       i from 0 to n-1 and appending into each column. */
    int32_t *colcur = (int32_t*)ldlt_xmalloc((size_t)n * sizeof(int32_t));
    memcpy(colcur, Lp, (size_t)n * sizeof(int32_t));
    for (int32_t i = 0; i < n; ++i) {
        Li[colcur[i]++] = i; /* diagonal */
    }
    for (int32_t i = 0; i < n; ++i) mark[i] = -1;
    for (int32_t i = 0; i < n; ++i) {
        mark[i] = i;
        for (int32_t p = Up[i]; p < Up[i+1]; ++p) {
            int32_t s = Ui[p];
            while (s != -1 && mark[s] != i) {
                mark[s] = i;
                Li[colcur[s]++] = i;
                s = parent[s];
            }
        }
    }
    free(colcur); free(colcount); free(mark); free(stk); free(Up); free(Ui);
    *Lp_out = Lp; *Li_out = Li;
    return LDLT_OK;
}

ldlt_status ldlt_build_symbolic(int32_t n, const int32_t *Ap, const int32_t *Ai,
                                const ldlt_options *opt, ldlt_symbolic *S)
{
    if (n < 0 || !Ap || (n>0 && !Ai)) return LDLT_ERR_INPUT;
    S->n = n;

    int32_t *perm0 = (int32_t*)ldlt_xmalloc((size_t)n * sizeof(int32_t));
    int32_t *iperm0 = (int32_t*)ldlt_xmalloc((size_t)n * sizeof(int32_t));
    double t_ord0 = ldlt_wall_time_seconds();
    ldlt_status st = ldlt_compute_ordering(n, Ap, Ai, opt, perm0);
    if (st != LDLT_OK) { free(perm0); free(iperm0); return st; }
    printf("LDLT %-39s %8.3f s\n", "analyze: ordering", ldlt_wall_time_seconds() - t_ord0);
    for (int32_t i = 0; i < n; ++i) iperm0[perm0[i]] = i;

    int32_t *B0p=NULL, *B0i=NULL;
    st = ldlt_permute_lower_pattern(n, Ap, Ai, perm0, iperm0, &B0p, &B0i);
    if (st != LDLT_OK) { free(perm0); free(iperm0); return st; }

    int32_t *parent0 = (int32_t*)ldlt_xmalloc((size_t)n * sizeof(int32_t));
    st = ldlt_etree(n, B0p, B0i, parent0);
    if (st != LDLT_OK) { free(perm0); free(iperm0); free(B0p); free(B0i); free(parent0); return st; }
    int32_t *post = (int32_t*)ldlt_xmalloc((size_t)n * sizeof(int32_t));
    st = ldlt_postorder(n, parent0, post);
    if (st != LDLT_OK) {
        free(perm0); free(iperm0); free(B0p); free(B0i); free(parent0); free(post);
        return st;
    }
    free(B0p); free(B0i);

    /* Compose perm with postorder: final perm[k] = perm0[post[k]] */
    int32_t *perm = (int32_t*)ldlt_xmalloc((size_t)n * sizeof(int32_t));
    int32_t *iperm = (int32_t*)ldlt_xmalloc((size_t)n * sizeof(int32_t));
    for (int32_t k = 0; k < n; ++k) perm[k] = perm0[post[k]];
    for (int32_t k = 0; k < n; ++k) iperm[perm[k]] = k;
    free(perm0); free(iperm0);

    /* Permute A by composed perm to get final B (pattern only). */
    int32_t *Bp=NULL, *Bi=NULL;
    st = ldlt_permute_lower_pattern(n, Ap, Ai, perm, iperm, &Bp, &Bi);
    if (st != LDLT_OK) {
        free(perm); free(iperm); free(parent0); free(post); return st;
    }

    /* Relabel parent0 (in perm0 coords) into parent (in postorder coords) instead
       of recomputing the etree on B. parent_new[k] = inv_post[parent0[post[k]]]. */
    int32_t *inv_post = (int32_t*)ldlt_xmalloc((size_t)n * sizeof(int32_t));
    for (int32_t k = 0; k < n; ++k) inv_post[post[k]] = k;
    int32_t *parent = (int32_t*)ldlt_xmalloc((size_t)n * sizeof(int32_t));
    for (int32_t k = 0; k < n; ++k) {
        int32_t p = parent0[post[k]];
        parent[k] = (p < 0) ? -1 : inv_post[p];
    }
    free(inv_post); free(parent0); free(post);
    int32_t *postI = (int32_t*)ldlt_xmalloc((size_t)n * sizeof(int32_t));
    for (int32_t k = 0; k < n; ++k) postI[k] = k; /* identity after composition */

    int32_t *Lp=NULL, *Li=NULL;
    st = build_L_pattern(n, Bp, Bi, parent, &Lp, &Li);
    if (st != LDLT_OK) {
        free(perm); free(iperm); free(parent); free(postI); free(Bp); free(Bi);
        return st;
    }

    /* Detect fundamental supernodes */
    int32_t *nchild = (int32_t*)ldlt_xcalloc((size_t)n, sizeof(int32_t));
    for (int32_t i = 0; i < n; ++i) if (parent[i] != -1) nchild[parent[i]]++;

    int32_t *col_to_super = (int32_t*)ldlt_xmalloc((size_t)n * sizeof(int32_t));
    int32_t *fund_first = (int32_t*)ldlt_xmalloc((size_t)(n+1) * sizeof(int32_t));
    int32_t nfund = 0;
    if (n > 0) {
        fund_first[0] = 0;
        for (int32_t j = 1; j < n; ++j) {
            int32_t prev = j-1;
            int32_t plen = Lp[prev+1] - Lp[prev]; /* incl diagonal */
            int32_t clen = Lp[j+1] - Lp[j];
            int merge = (parent[prev] == j) && (nchild[j] == 1)
                        && (plen - 1 == clen);
            if (!merge) fund_first[++nfund] = j;
        }
        fund_first[++nfund] = n;
    }

    int32_t *sn_first = (int32_t*)ldlt_xmalloc((size_t)(n+1) * sizeof(int32_t));
    int32_t nsuper = 0;
    if (n > 0) {
        int relax = (opt && opt->relax_supernodes > 0) ? opt->relax_supernodes : 0;
        int32_t *mark = (int32_t*)ldlt_xmalloc((size_t)n * sizeof(int32_t));
        for (int32_t i = 0; i < n; ++i) mark[i] = -1;
        int32_t tag = 0;
        sn_first[0] = fund_first[0];
        int32_t cur_start = fund_first[0];
        int32_t cur_end = fund_first[1];
        int64_t cur_rows = super_panel_rows(cur_start, cur_end, Lp, Li, mark, tag++, NULL);
        int64_t cur_size = (int64_t)(cur_end - cur_start) * (cur_rows + cur_end - cur_start);
        int64_t cur_fund_size = cur_size;
        for (int32_t f = 1; f < nfund; ++f) {
            int32_t next_start = fund_first[f];
            int32_t next_end = fund_first[f+1];
            int64_t next_rows = super_panel_rows(next_start, next_end, Lp, Li, mark, tag++, NULL);
            int64_t next_size = (int64_t)(next_end - next_start) *
                                (next_rows + next_end - next_start);
            int64_t next_fund_size = next_size;
            int do_merge = 0;
            if (relax > 0 && parent[cur_end - 1] >= next_start && parent[cur_end - 1] < next_end) {
                int64_t merged_size = super_panel_size(cur_start, next_end, Lp, Li, mark, tag++);
                int64_t overlap = super_panel_overlap(cur_start, cur_end, next_start, next_end,
                                                      Lp, Li, mark, tag);
                tag += 2;
                int64_t fund_size = cur_fund_size + next_fund_size;
                if (adaptive_relax_accept(relax, cur_start, cur_end, next_start, next_end,
                                          cur_size, next_size, merged_size,
                                          cur_fund_size, next_fund_size,
                                          cur_rows, next_rows, overlap,
                                          parent, nchild)) {
                    do_merge = 1;
                    cur_end = next_end;
                    cur_size = merged_size;
                    cur_fund_size = fund_size;
                    cur_rows = merged_size / (cur_end - cur_start) - (cur_end - cur_start);
                }
            }
            if (!do_merge) {
                sn_first[++nsuper] = next_start;
                cur_start = next_start;
                cur_end = next_end;
                cur_size = next_size;
                cur_fund_size = next_fund_size;
                cur_rows = next_rows;
            }
        }
        sn_first[++nsuper] = n;
        free(mark);
    }
    free(fund_first);
    free(nchild);

    for (int32_t s = 0; s < nsuper; ++s)
        for (int32_t j = sn_first[s]; j < sn_first[s+1]; ++j)
            col_to_super[j] = s;

    S->nsuper = nsuper;
    S->perm = perm; S->iperm = iperm;
    S->perm_identity = 1;
    for (int32_t i = 0; i < n; ++i) {
        if (perm[i] != i) { S->perm_identity = 0; break; }
    }
    S->parent = parent; S->postorder = postI;
    S->col_to_super = col_to_super;
    S->super = (ldlt_super*)ldlt_xcalloc((size_t)nsuper, sizeof(ldlt_super));
    S->panel_off = (int64_t*)ldlt_xmalloc((size_t)(nsuper+1) * sizeof(int64_t));
    S->d_off     = (int64_t*)ldlt_xmalloc((size_t)(nsuper+1) * sizeof(int64_t));
    S->panel_off[0] = 0; S->d_off[0] = 0;

    /* Build per-supernode rows_below from the union of L patterns in the supernode.
       Relaxed supernodes can include structural zeros, so the first column pattern is
       not sufficient once amalgamation is enabled. */
    int32_t *rowmark = (int32_t*)ldlt_xmalloc((size_t)(n>0?n:1) * sizeof(int32_t));
    for (int32_t i = 0; i < n; ++i) rowmark[i] = -1;
    int32_t rowtag = 0;
    int64_t rows_pool_cap = 0;
    for (int32_t s = 0; s < nsuper; ++s) {
        int32_t fc = sn_first[s], k = sn_first[s+1] - sn_first[s];
        rows_pool_cap += super_panel_rows(fc, fc + k, Lp, Li, rowmark, rowtag++, NULL);
        (void)k;
    }
    if (rows_pool_cap < 1) rows_pool_cap = 1;
    S->rows_pool = (int32_t*)ldlt_xmalloc((size_t)rows_pool_cap * sizeof(int32_t));
    int64_t roff = 0;

    int64_t lnz = 0;
    double flops = 0.0;
    for (int32_t s = 0; s < nsuper; ++s) {
        int32_t fc = sn_first[s];
        int32_t k = sn_first[s+1] - fc;
        int32_t m = (int32_t)super_panel_rows(fc, fc + k, Lp, Li, rowmark, rowtag++,
                                              S->rows_pool + roff);
        S->super[s].first_col = fc;
        S->super[s].width = k;
        S->super[s].nrows_below = m;
        S->super[s].row_off = roff;
        S->super[s].sn_perm = NULL; /* local supernode permutation reserved */
        roff += m;
        S->panel_off[s+1] = S->panel_off[s] + (int64_t)(k + m) * k;
        S->d_off[s+1] = S->d_off[s] + k;
        lnz += (int64_t)(k * (k+1) / 2) + (int64_t)m * k; /* approx L nnz */
        /* flops: dense LDLT k^3/3 + 2*m*k^2 (TRSM) + sum over updates handled later;
           use simple lower bound */
        flops += (double)k*k*k/3.0 + 2.0*(double)m*k*k;
    }
    S->rows_pool_len = roff;
    S->lnz = lnz;
    S->flops = flops;

    S->super_parent = (int32_t*)ldlt_xmalloc((size_t)(nsuper>0?nsuper:1) * sizeof(int32_t));
    S->child_ptr = (int32_t*)ldlt_xcalloc((size_t)nsuper+1, sizeof(int32_t));
    for (int32_t s = 0; s < nsuper; ++s) {
        int32_t last_col = S->super[s].first_col + S->super[s].width - 1;
        int32_t p = parent[last_col];
        while (p != -1 && col_to_super[p] == s) p = parent[p];
        int32_t sp = (p == -1) ? -1 : col_to_super[p];
        S->super_parent[s] = sp;
        if (sp >= 0) S->child_ptr[sp+1]++;
    }
    for (int32_t s = 0; s < nsuper; ++s) S->child_ptr[s+1] += S->child_ptr[s];
    int32_t nchildren = S->child_ptr[nsuper];
    S->child_idx = (int32_t*)ldlt_xmalloc((size_t)(nchildren>0?nchildren:1) * sizeof(int32_t));
    int32_t *child_cur = (int32_t*)ldlt_xmalloc((size_t)(nsuper>0?nsuper:1) * sizeof(int32_t));
    memcpy(child_cur, S->child_ptr, (size_t)(nsuper>0?nsuper:1) * sizeof(int32_t));
    for (int32_t s = 0; s < nsuper; ++s) {
        int32_t p = S->super_parent[s];
        if (p >= 0) S->child_idx[child_cur[p]++] = s;
    }
    free(child_cur);

    int32_t *lvl = (int32_t*)ldlt_xcalloc((size_t)(nsuper>0?nsuper:1), sizeof(int32_t));
    int32_t maxlvl = 0;
    for (int32_t s = 0; s < nsuper; ++s) {
        int32_t p = S->super_parent[s];
        if (p >= 0 && lvl[p] < lvl[s] + 1) {
            lvl[p] = lvl[s] + 1;
            if (lvl[p] > maxlvl) maxlvl = lvl[p];
        }
    }
    S->nlevels = nsuper > 0 ? maxlvl + 1 : 0;
    int32_t *level_width = (int32_t*)ldlt_xcalloc((size_t)(S->nlevels>0?S->nlevels:1),
                                                 sizeof(int32_t));
    for (int32_t s = 0; s < nsuper; ++s) {
        level_width[lvl[s]]++;
        if (level_width[lvl[s]] > S->max_level_width)
            S->max_level_width = level_width[lvl[s]];
    }
    free(level_width);
    free(lvl);

    S->edge_ptr = (int32_t*)ldlt_xmalloc((size_t)(nsuper+1) * sizeof(int32_t));
    S->edge_ptr[0] = 0;
    for (int32_t s = 0; s < nsuper; ++s) {
        int64_t next_edge = (int64_t)S->edge_ptr[s] + S->super[s].nrows_below;
        if (next_edge > INT32_MAX) {
            free(rowmark); free(Lp); free(Li); free(Bp); free(Bi); free(sn_first);
            return LDLT_ERR_NOMEM;
        }
        S->edge_ptr[s+1] = (int32_t)next_edge;
    }
    int32_t edge_total = S->edge_ptr[nsuper];
    S->edge_pos = (int32_t*)ldlt_xmalloc((size_t)(edge_total>0?edge_total:1) * sizeof(int32_t));
    int32_t *front_map = (int32_t*)ldlt_xmalloc((size_t)(n>0?n:1) * sizeof(int32_t));
    for (int32_t i = 0; i < n; ++i) front_map[i] = -1;
    for (int32_t p = 0; p < nsuper; ++p) {
        const ldlt_super *Sp = &S->super[p];
        int32_t pfc = Sp->first_col, pk = Sp->width, pm = Sp->nrows_below;
        const int32_t *prb = S->rows_pool + Sp->row_off;
        for (int32_t i = 0; i < pk; ++i) front_map[pfc + i] = i;
        for (int32_t i = 0; i < pm; ++i) front_map[prb[i]] = pk + i;
        for (int32_t cp = S->child_ptr[p]; cp < S->child_ptr[p+1]; ++cp) {
            int32_t c = S->child_idx[cp];
            const ldlt_super *Sc = &S->super[c];
            const int32_t *crb = S->rows_pool + Sc->row_off;
            int32_t *emap = S->edge_pos + S->edge_ptr[c];
            for (int32_t i = 0; i < Sc->nrows_below; ++i) {
                int32_t pos = front_map[crb[i]];
                if (pos < 0) {
                    free(front_map);
                    free(Lp); free(Li); free(Bp); free(Bi); free(sn_first);
                    return LDLT_ERR_INPUT;
                }
                emap[i] = pos;
            }
        }
        for (int32_t i = 0; i < pk; ++i) front_map[pfc + i] = -1;
        for (int32_t i = 0; i < pm; ++i) front_map[prb[i]] = -1;
    }
    free(front_map);

    /* Update DAG: for each supernode d, for each row r in rows_below, target supernode = col_to_super[r] */
    S->upd_ptr = (int32_t*)ldlt_xcalloc((size_t)nsuper+1, sizeof(int32_t));
    int32_t *last = (int32_t*)ldlt_xmalloc((size_t)nsuper * sizeof(int32_t));
    for (int32_t i = 0; i < nsuper; ++i) last[i] = -1;
    /* Count */
    for (int32_t d = 0; d < nsuper; ++d) {
        int64_t off = S->super[d].row_off; int32_t m = S->super[d].nrows_below;
        for (int32_t t = 0; t < m; ++t) {
            int32_t r = S->rows_pool[off + t];
            int32_t s = col_to_super[r];
            if (last[s] != d) { last[s] = d; S->upd_ptr[s+1]++; }
        }
    }
    for (int32_t s = 0; s < nsuper; ++s) S->upd_ptr[s+1] += S->upd_ptr[s];
    int32_t total = S->upd_ptr[nsuper];
    S->upd_idx = (int32_t*)ldlt_xmalloc((size_t)(total>0?total:1) * sizeof(int32_t));
    int32_t *cur = (int32_t*)ldlt_xmalloc((size_t)nsuper * sizeof(int32_t));
    memcpy(cur, S->upd_ptr, (size_t)nsuper * sizeof(int32_t));
    for (int32_t i = 0; i < nsuper; ++i) last[i] = -1;
    for (int32_t d = 0; d < nsuper; ++d) {
        int64_t off = S->super[d].row_off; int32_t m = S->super[d].nrows_below;
        for (int32_t t = 0; t < m; ++t) {
            int32_t r = S->rows_pool[off + t];
            int32_t s = col_to_super[r];
            if (last[s] != d) { last[s] = d; S->upd_idx[cur[s]++] = d; }
        }
    }
    free(cur); free(last); free(rowmark);

    free(Lp); free(Li); free(Bp); free(Bi); free(sn_first);
    return LDLT_OK;
}

ldlt_status ldlt_analyze(int32_t n, const int32_t *Ap, const int32_t *Ai,
                         const ldlt_options *opt, ldlt_symbolic **out)
{
    if (!out) return LDLT_ERR_INPUT;
    *out = NULL;
    ldlt_options local; ldlt_options_default(&local);
    if (opt) local = *opt;
    if (local.pivot != LDLT_PIVOT_NONE &&
        local.pivot != LDLT_PIVOT_BUNCH_KAUFMAN)
        return LDLT_ERR_UNSUPPORTED;
    ldlt_symbolic * volatile S = NULL;
    ldlt_error_trap trap;
    ldlt_push_error_trap(&trap);
    if (setjmp(trap.env) != 0) {
        ldlt_pop_error_trap(&trap);
        ldlt_free_symbolic((ldlt_symbolic*)S);
        return trap.status == LDLT_OK ? LDLT_ERR_NOMEM : trap.status;
    }

    S = (ldlt_symbolic*)ldlt_xcalloc(1, sizeof(*S));
    ldlt_status st = ldlt_build_symbolic(n, Ap, Ai, &local, S);
    if (st != LDLT_OK) {
        ldlt_free_symbolic((ldlt_symbolic*)S);
        ldlt_pop_error_trap(&trap);
        return st;
    }
    *out = (ldlt_symbolic*)S;
    ldlt_pop_error_trap(&trap);
    return LDLT_OK;
}
