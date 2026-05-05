/* Elimination tree and postorder.
   Lower-triangle CSC input; etree built via Liu's path-compression algorithm. */

#include "internal.h"

/* Liu/Davis etree with path compression. Iterates rows in increasing order
   (outer = "k"), walking upward from each upper-triangle row index of column k.
   Input is lower-CSC, so we first build per-row lists of cols j<i (Up/Ui), then
   process each row's upper entries against the in-progress ancestor[] chain. */
ldlt_status ldlt_etree(int32_t n, const int32_t *Ap, const int32_t *Ai, int32_t *parent) {
    int32_t *Up = (int32_t*)ldlt_xcalloc((size_t)n+1, sizeof(int32_t));
    for (int32_t k = 0; k < n; ++k)
        for (int32_t p = Ap[k]; p < Ap[k+1]; ++p) {
            int32_t i = Ai[p];
            if (i > k) Up[i+1]++;
        }
    for (int32_t i = 0; i < n; ++i) Up[i+1] += Up[i];
    int32_t *Ui = (int32_t*)ldlt_xmalloc((size_t)(Up[n]>0?Up[n]:1) * sizeof(int32_t));
    int32_t *cur = (int32_t*)ldlt_xmalloc((size_t)n * sizeof(int32_t));
    memcpy(cur, Up, (size_t)n * sizeof(int32_t));
    for (int32_t k = 0; k < n; ++k)
        for (int32_t p = Ap[k]; p < Ap[k+1]; ++p) {
            int32_t i = Ai[p];
            if (i > k) Ui[cur[i]++] = k;
        }
    free(cur);

    int32_t *anc = (int32_t*)ldlt_xmalloc((size_t)n * sizeof(int32_t));
    for (int32_t k = 0; k < n; ++k) { parent[k] = -1; anc[k] = -1; }
    for (int32_t i = 0; i < n; ++i) {
        for (int32_t p = Up[i]; p < Up[i+1]; ++p) {
            int32_t s = Ui[p];                  /* col, s < i */
            while (s != -1 && s < i) {
                int32_t nxt = anc[s];
                anc[s] = i;
                if (nxt == -1) { parent[s] = i; break; }
                s = nxt;
            }
        }
    }
    free(anc); free(Up); free(Ui);
    return LDLT_OK;
}

static int32_t tdfs(int32_t root, int32_t k, int32_t *head, int32_t *next,
                    int32_t *post, int32_t *stack) {
    int32_t top = 0;
    stack[0] = root;
    while (top >= 0) {
        int32_t p = stack[top];
        int32_t i = head[p];
        if (i == -1) { post[k++] = p; --top; }
        else { head[p] = next[i]; stack[++top] = i; }
    }
    return k;
}

ldlt_status ldlt_postorder(int32_t n, const int32_t *parent, int32_t *post) {
    int32_t *head = (int32_t*)ldlt_xmalloc((size_t)n * sizeof(int32_t));
    int32_t *next = (int32_t*)ldlt_xmalloc((size_t)n * sizeof(int32_t));
    int32_t *stk  = (int32_t*)ldlt_xmalloc((size_t)n * sizeof(int32_t));
    for (int32_t i = 0; i < n; ++i) head[i] = -1;
    for (int32_t j = n-1; j >= 0; --j) {
        if (parent[j] == -1) continue;
        next[j] = head[parent[j]];
        head[parent[j]] = j;
    }
    int32_t k = 0;
    for (int32_t j = 0; j < n; ++j) {
        if (parent[j] != -1) continue;
        k = tdfs(j, k, head, next, post, stk);
    }
    free(head); free(next); free(stk);
    return LDLT_OK;
}

/* not used — symbolic.c builds full L pattern instead */
ldlt_status ldlt_col_counts(int32_t n, const int32_t *Ap, const int32_t *Ai,
                            const int32_t *parent, const int32_t *post,
                            int32_t *colcount) {
    (void)Ap; (void)Ai; (void)parent; (void)post;
    for (int32_t i = 0; i < n; ++i) colcount[i] = 0;
    return LDLT_OK;
}
