#include "internal.h"

static int cmp_int(const void *a, const void *b) {
    int32_t x = *(const int32_t*)a, y = *(const int32_t*)b;
    return (x>y) - (x<y);
}

/* Pattern-only variant: compute the structure of B = P A P^T (lower triangle)
   without touching values. Used during the symbolic phase, where build_L_pattern
   only needs Bp/Bi. Avoids ~nnz(A) extra allocations and value sorts. */
ldlt_status ldlt_permute_lower_pattern(int32_t n,
    const int32_t *Ap, const int32_t *Ai,
    const int32_t *perm, const int32_t *iperm,
    int32_t **Bp_out, int32_t **Bi_out)
{
    (void)perm;
    int32_t *Bp = (int32_t*)ldlt_xcalloc((size_t)n+1, sizeof(int32_t));
    for (int32_t j = 0; j < n; ++j) {
        for (int32_t p = Ap[j]; p < Ap[j+1]; ++p) {
            int32_t i = Ai[p];
            int32_t ni = iperm[i], nj = iperm[j];
            int32_t a = ni < nj ? ni : nj;
            Bp[a+1]++;
        }
    }
    for (int32_t i = 0; i < n; ++i) Bp[i+1] += Bp[i];
    int32_t nnz = Bp[n];
    int32_t *Bi = (int32_t*)ldlt_xmalloc((size_t)(nnz>0?nnz:1) * sizeof(int32_t));
    int32_t *cur = (int32_t*)ldlt_xmalloc((size_t)n * sizeof(int32_t));
    memcpy(cur, Bp, (size_t)n * sizeof(int32_t));
    for (int32_t j = 0; j < n; ++j) {
        for (int32_t p = Ap[j]; p < Ap[j+1]; ++p) {
            int32_t i = Ai[p];
            int32_t ni = iperm[i], nj = iperm[j];
            int32_t newcol = ni < nj ? ni : nj;
            int32_t newrow = ni < nj ? nj : ni;
            Bi[cur[newcol]++] = newrow;
        }
    }
    /* sort each column's row indices (insertion sort: typically short columns) */
    for (int32_t j = 0; j < n; ++j) {
        int32_t s = Bp[j], e = Bp[j+1];
        for (int32_t a = s+1; a < e; ++a) {
            int32_t ki = Bi[a];
            int32_t b = a-1;
            while (b >= s && Bi[b] > ki) { Bi[b+1]=Bi[b]; --b; }
            Bi[b+1] = ki;
        }
    }
    free(cur);
    *Bp_out = Bp; *Bi_out = Bi;
    return LDLT_OK;
}

/* Compute B = P A P^T but keep only the lower triangle (rows >= cols).
   Inputs: A is symmetric, lower-triangle stored in CSC.
   perm: perm[k] = original index placed at new position k.
   iperm[i] = new position of original index i. */
ldlt_status ldlt_permute_lower(int32_t n,
    const int32_t *Ap, const int32_t *Ai, const double *Ax,
    const int32_t *perm, const int32_t *iperm,
    int32_t **Bp_out, int32_t **Bi_out, double **Bx_out)
{
    (void)perm;
    int32_t *Bp = (int32_t*)ldlt_xcalloc((size_t)n+1, sizeof(int32_t));
    /* Count entries per new-column */
    for (int32_t j = 0; j < n; ++j) {
        for (int32_t p = Ap[j]; p < Ap[j+1]; ++p) {
            int32_t i = Ai[p];
            int32_t ni = iperm[i], nj = iperm[j];
            int32_t a = ni < nj ? ni : nj; /* new column */
            (void)nj;
            Bp[a+1]++;
        }
    }
    for (int32_t i = 0; i < n; ++i) Bp[i+1] += Bp[i];
    int32_t nnz = Bp[n];
    int32_t *Bi = (int32_t*)ldlt_xmalloc((size_t)nnz * sizeof(int32_t));
    double  *Bx = (double*)ldlt_xmalloc((size_t)nnz * sizeof(double));
    int32_t *cur = (int32_t*)ldlt_xmalloc((size_t)n * sizeof(int32_t));
    memcpy(cur, Bp, (size_t)n * sizeof(int32_t));
    for (int32_t j = 0; j < n; ++j) {
        for (int32_t p = Ap[j]; p < Ap[j+1]; ++p) {
            int32_t i = Ai[p];
            int32_t ni = iperm[i], nj = iperm[j];
            int32_t newcol = ni < nj ? ni : nj;
            int32_t newrow = ni < nj ? nj : ni;
            int32_t k = cur[newcol]++;
            Bi[k] = newrow;
            Bx[k] = Ax[p];
        }
    }
    /* sort each column's row indices */
    for (int32_t j = 0; j < n; ++j) {
        int32_t s = Bp[j], e = Bp[j+1];
        /* simple sort with parallel value array */
        for (int32_t a = s+1; a < e; ++a) {
            int32_t ki = Bi[a]; double kx = Bx[a];
            int32_t b = a-1;
            while (b >= s && Bi[b] > ki) { Bi[b+1]=Bi[b]; Bx[b+1]=Bx[b]; --b; }
            Bi[b+1] = ki; Bx[b+1] = kx;
        }
    }
    free(cur);
    *Bp_out = Bp; *Bi_out = Bi; *Bx_out = Bx;
    (void)cmp_int;
    return LDLT_OK;
}

void ldlt_apply_perm(int32_t n, const int32_t *perm, const double *x, double *y) {
    for (int32_t k = 0; k < n; ++k) y[k] = x[perm[k]];
}
void ldlt_apply_iperm(int32_t n, const int32_t *perm, const double *x, double *y) {
    for (int32_t k = 0; k < n; ++k) y[perm[k]] = x[k];
}
