/* Minimal Matrix Market reader for symmetric coordinate real/pattern matrices.
   Returns lower-triangle CSC (rows >= cols), 0-indexed.
   Pattern matrices are loaded with value 1.0 for each structural entry; callers
   should not treat symmetric pattern-only files as validated numeric matrices. */

#include "internal.h"
#include <stdio.h>
#include <ctype.h>
#include <string.h>

static int starts_with(const char *s, const char *pre) {
    size_t n = strlen(pre);
    return strncmp(s, pre, n) == 0;
}

static int cmp_pair(const void *a, const void *b) {
    const int32_t *p = (const int32_t*)a, *q = (const int32_t*)b;
    if (p[1] != q[1]) return p[1] - q[1];
    return p[0] - q[0];
}

ldlt_status ldlt_read_mm_symmetric(const char *path, int32_t *n_out,
                                   int32_t **Ap_out, int32_t **Ai_out, double **Ax_out)
{
    FILE *f = fopen(path, "r");
    if (!f) return LDLT_ERR_INPUT;
    char line[1024];
    int symmetric = 0;
    int pattern = 0;
    if (!fgets(line, sizeof(line), f)) { fclose(f); return LDLT_ERR_INPUT; }
    if (!starts_with(line, "%%MatrixMarket")) { fclose(f); return LDLT_ERR_INPUT; }
    if (strstr(line, "symmetric")) symmetric = 1;
    if (strstr(line, "pattern")) pattern = 1;
    if (symmetric && pattern) {
        fprintf(stderr,
                "ldlt warning: loading symmetric pattern Matrix Market file '%s' "
                "with implicit value 1.0 entries; this is a structural pattern, "
                "not a validated numeric matrix.\n",
                path);
    }
    do {
        if (!fgets(line, sizeof(line), f)) { fclose(f); return LDLT_ERR_INPUT; }
    } while (line[0] == '%');
    long M, N, nnz;
    if (sscanf(line, "%ld %ld %ld", &M, &N, &nnz) != 3) { fclose(f); return LDLT_ERR_INPUT; }
    if (M != N) { fclose(f); return LDLT_ERR_INPUT; }

    /* read triplets, keep lower-tri */
    int32_t *trip = (int32_t*)ldlt_xmalloc((size_t)nnz * 3 * sizeof(int32_t));
    double *vals  = (double*)ldlt_xmalloc((size_t)nnz * sizeof(double));
    long actual = 0;
    for (long t = 0; t < nnz; ++t) {
        long i, j; double v = 1.0;
        int nread = pattern ? fscanf(f, "%ld %ld", &i, &j)
                            : fscanf(f, "%ld %ld %lf", &i, &j, &v);
        if (nread < 2) { fclose(f); free(trip); free(vals); return LDLT_ERR_INPUT; }
        i--; j--;
        if (i < j) { long tmp=i; i=j; j=tmp; } /* keep lower tri (row >= col) */
        if (symmetric || i >= j) {
            trip[3*actual+0] = (int32_t)i;
            trip[3*actual+1] = (int32_t)j;
            trip[3*actual+2] = (int32_t)actual;
            vals[actual] = v;
            actual++;
        }
    }
    fclose(f);

    /* sort by (col, row) */
    /* Use a temporary index array to sort vals along with trip */
    int32_t *idx = (int32_t*)ldlt_xmalloc((size_t)actual * sizeof(int32_t));
    for (long t = 0; t < actual; ++t) idx[t] = (int32_t)t;
    /* simple in-place sort: pack into pair-array then sort, then permute vals */
    /* Build pair array (row,col,orig_idx) */
    int32_t *pairs = (int32_t*)ldlt_xmalloc((size_t)actual * 3 * sizeof(int32_t));
    for (long t = 0; t < actual; ++t) {
        pairs[3*t+0] = trip[3*t+0]; pairs[3*t+1] = trip[3*t+1]; pairs[3*t+2] = (int32_t)t;
    }
    qsort(pairs, (size_t)actual, 3*sizeof(int32_t), cmp_pair);
    int32_t *Ap = (int32_t*)ldlt_xcalloc((size_t)N+1, sizeof(int32_t));
    int32_t *Ai = (int32_t*)ldlt_xmalloc((size_t)actual * sizeof(int32_t));
    double  *Ax = (double*)ldlt_xmalloc((size_t)actual * sizeof(double));
    for (long t = 0; t < actual; ++t) Ap[pairs[3*t+1]+1]++;
    for (long j = 0; j < N; ++j) Ap[j+1] += Ap[j];
    for (long t = 0; t < actual; ++t) {
        Ai[t] = pairs[3*t+0];
        Ax[t] = vals[pairs[3*t+2]];
    }
    free(pairs); free(idx); free(trip); free(vals);
    *n_out = (int32_t)N; *Ap_out = Ap; *Ai_out = Ai; *Ax_out = Ax;
    return LDLT_OK;
}
