#ifndef LDLT_INTERNAL_H
#define LDLT_INTERNAL_H

#include "ldlt/ldlt.h"
#include <setjmp.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#ifndef LDLT_USE_SYSTEM_BLAS
#define LDLT_USE_SYSTEM_BLAS 0
#endif
#ifndef LDLT_HAVE_AMD
#define LDLT_HAVE_AMD 0
#endif
#ifndef LDLT_HAVE_METIS
#define LDLT_HAVE_METIS 0
#endif
#ifndef LDLT_HAVE_OPENMP
#define LDLT_HAVE_OPENMP 0
#endif
#ifndef LDLT_USE_POSIX_SIGNALS
#define LDLT_USE_POSIX_SIGNALS 0
#endif

#ifndef LDLT_THREAD_LOCAL
#if defined(_MSC_VER)
#define LDLT_THREAD_LOCAL __declspec(thread)
#elif defined(__STDC_VERSION__) && __STDC_VERSION__ >= 201112L && !defined(__STDC_NO_THREADS__)
#define LDLT_THREAD_LOCAL _Thread_local
#elif defined(__GNUC__) || defined(__clang__)
#define LDLT_THREAD_LOCAL __thread
#else
#define LDLT_THREAD_LOCAL
#endif
#endif

/* ---- ldlt_symbolic / ldlt_numeric internals ---- */

typedef struct {
    int32_t first_col;     /* global col index of first col in supernode */
    int32_t width;         /* k: number of cols */
    int32_t nrows_below;   /* m: rows strictly below the diagonal block */
    int32_t row_off;       /* offset into rows_pool of length nrows_below */
    int32_t *sn_perm;      /* reserved for local supernode permutations */
} ldlt_super;

struct ldlt_symbolic {
    int32_t n;
    int32_t nsuper;
    int32_t perm_identity;
    int32_t nlevels;
    int32_t max_level_width;
    int32_t *perm;         /* length n: perm[k] = original col placed at position k */
    int32_t *iperm;        /* inverse */
    int32_t *parent;       /* etree of permuted matrix, length n */
    int32_t *postorder;    /* length n */
    int32_t *col_to_super; /* length n */
    ldlt_super *super;     /* length nsuper */
    int32_t *rows_pool;    /* concatenated rows_below across supernodes */
    int64_t  rows_pool_len;
    int32_t *super_parent; /* supernodal etree, length nsuper; -1 for roots */
    int32_t *child_ptr;    /* children of each supernode, length nsuper+1 */
    int32_t *child_idx;    /* child_ptr[nsuper] entries */
    int32_t *edge_ptr;     /* child-to-parent row maps, length nsuper+1 */
    int32_t *edge_pos;     /* for child s, parent-front row positions of rows_below(s) */
    /* For each supernode s, list of supernode indices d that update s (left-looking). */
    int32_t *upd_ptr;      /* length nsuper+1 */
    int32_t *upd_idx;      /* upd_ptr[nsuper] entries */
    /* Sizes: per-supernode L panel size = (width+nrows_below)*width doubles. */
    int64_t *panel_off;    /* length nsuper+1; offsets into a packed L array (in doubles) */
    int64_t *d_off;        /* length nsuper+1; offsets into D array (in doubles) */
    int64_t  lnz;
    double   flops;
};

struct ldlt_numeric {
    const ldlt_symbolic *S; /* not owned */
    double *L;              /* packed panels; size = panel_off[nsuper] */
    double *D;              /* size = d_off[nsuper] = n */
    int32_t pivoted;        /* non-zero when using dynamic 1x1/2x2 pivots */
    int32_t npivots;        /* number of pivot blocks */
    int32_t n_piv_vars;     /* total pivoted variables, should equal S->n */
    int32_t *piv_size;      /* length <= n: block sizes, 1 or 2 */
    int32_t *piv_var_off;   /* length <= n+1: offsets into piv_vars */
    int32_t *piv_vars;      /* length n: permuted-space variable ids */
    int64_t *piv_d_off;     /* length <= n+1: offsets into piv_D */
    double  *piv_D;         /* packed D blocks: [d] or [d00,d10,d11] */
    int64_t *piv_l_row_off; /* length <= n+1: offsets into piv_l_rows */
    int64_t *piv_l_val_off; /* length <= n+1: offsets into piv_l_vals */
    int32_t *piv_l_rows;    /* permuted-space row ids for each block's L21 */
    double  *piv_l_vals;    /* row-major L21 values per block */
};

/* ---- helpers ---- */
typedef struct ldlt_error_trap {
    jmp_buf env;
    ldlt_status status;
    struct ldlt_error_trap *prev;
} ldlt_error_trap;

void ldlt_push_error_trap(ldlt_error_trap *trap);
void ldlt_pop_error_trap(ldlt_error_trap *trap);
void ldlt_raise_status(ldlt_status status);
void *ldlt_xmalloc(size_t bytes);
void *ldlt_xcalloc(size_t n, size_t sz);
void *ldlt_xrealloc(void *p, size_t bytes);
uint64_t ldlt_available_memory_bytes(void);
double ldlt_wall_time_seconds(void);

/* etree.c */
ldlt_status ldlt_etree(int32_t n, const int32_t *Ap, const int32_t *Ai, int32_t *parent);
ldlt_status ldlt_postorder(int32_t n, const int32_t *parent, int32_t *post);
ldlt_status ldlt_col_counts(int32_t n, const int32_t *Ap, const int32_t *Ai,
                            const int32_t *parent, const int32_t *post,
                            int32_t *colcount);

/* perm.c */
ldlt_status ldlt_permute_lower(int32_t n,
    const int32_t *Ap, const int32_t *Ai, const double *Ax,
    const int32_t *perm, const int32_t *iperm,
    int32_t **Bp_out, int32_t **Bi_out, double **Bx_out);
ldlt_status ldlt_permute_lower_pattern(int32_t n,
    const int32_t *Ap, const int32_t *Ai,
    const int32_t *perm, const int32_t *iperm,
    int32_t **Bp_out, int32_t **Bi_out);
void ldlt_apply_perm(int32_t n, const int32_t *perm, const double *x, double *y);
void ldlt_apply_iperm(int32_t n, const int32_t *perm, const double *x, double *y);

/* ordering.c */
ldlt_status ldlt_compute_ordering(int32_t n, const int32_t *Ap, const int32_t *Ai,
                                  const ldlt_options *opt, int32_t *perm);

/* amd_stub.c */
ldlt_status ldlt_amd_stub(int32_t n, const int32_t *Ap, const int32_t *Ai, int32_t *perm);

/* symbolic.c */
ldlt_status ldlt_build_symbolic(int32_t n, const int32_t *Ap, const int32_t *Ai,
                                const ldlt_options *opt, ldlt_symbolic *S);

/* factor.c / factor_kernels.c */
ldlt_status ldlt_factor_supernode(const struct ldlt_symbolic *S, struct ldlt_numeric *N,
                                  int32_t s, double *workspace, int32_t *rowmap,
                                  const ldlt_options *opt,
                                  const int32_t *Ap, const int32_t *Ai, const double *Ax);

ldlt_status ldlt_dense_ldl_nopiv(double *F, int32_t ldF, int32_t R, int32_t k,
                                 double *D, double tiny, int *indef);
ldlt_status ldlt_dense_front_nopiv(double *F, int32_t ldF, int32_t R, int32_t k,
                                   double *D, double *work, double tiny, int *indef);

/* io_mm.c (used by tests) */
ldlt_status ldlt_read_mm_symmetric(const char *path, int32_t *n_out,
                                   int32_t **Ap_out, int32_t **Ai_out,
                                   double **Ax_out);

#endif
