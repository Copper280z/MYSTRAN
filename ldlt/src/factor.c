/* Left-looking supernodal LDL^T numeric factorization.
   Parallelism: level-set scheduling. Each level is processed with #pragma omp parallel for;
   levels execute serially with an implicit barrier between them. */

#include "internal.h"
#include "blas_shim.h"
#include <math.h>
#include <stdio.h>

#if LDLT_HAVE_OPENMP
#include <omp.h>
#endif

typedef struct {
    int32_t level;
    int64_t nodes;
    int64_t lapack_try;
    int64_t lapack_accept;
    int64_t fallback_info;
    int64_t fallback_ipiv;
    int64_t fallback_l21;
    int64_t scalar_nodes;
    int64_t scalar_1x1;
    int64_t scalar_2x2;
    int64_t pivots;
    int64_t delayed_vars;
    int32_t max_R;
    int32_t max_cand;
    int32_t max_noff;
    double t_total;
    double t_assemble;
    double t_lapack;
    double t_scalar;
    double t_scalar_select;
    double t_scalar_update;
    double t_store;
} pivot_profile_level;

typedef struct {
    int enabled;
    double t_start;
    pivot_profile_level cur;
#if LDLT_HAVE_OPENMP
    omp_lock_t lock;
#endif
} pivot_profile;

typedef struct {
    int32_t R;
    int32_t cand;
    int32_t noff_initial;
    int32_t delayed_vars;
    int64_t lapack_try;
    int64_t lapack_accept;
    int64_t fallback_info;
    int64_t fallback_ipiv;
    int64_t fallback_l21;
    int64_t scalar_used;
    int64_t scalar_1x1;
    int64_t scalar_2x2;
    int64_t pivots;
    double t_total;
    double t_assemble;
    double t_lapack;
    double t_scalar;
    double t_scalar_select;
    double t_scalar_update;
    double t_store;
} pivot_node_profile;

static double prof_now_sec(void)
{
    return ldlt_wall_time_seconds();
}

static void pivot_profile_init(pivot_profile *P)
{
    memset(P, 0, sizeof(*P));
    P->enabled = getenv("LDLT_PROFILE_PIVOT") != NULL;
    if (!P->enabled) return;
    P->t_start = prof_now_sec();
#if LDLT_HAVE_OPENMP
    omp_init_lock(&P->lock);
#endif
    fprintf(stderr, "ldlt profile: pivot factorization profiling enabled\n");
}

static void pivot_profile_destroy(pivot_profile *P)
{
#if LDLT_HAVE_OPENMP
    if (P->enabled) omp_destroy_lock(&P->lock);
#else
    (void)P;
#endif
}

static void pivot_profile_begin_level(pivot_profile *P, int32_t level)
{
    if (!P || !P->enabled) return;
    memset(&P->cur, 0, sizeof(P->cur));
    P->cur.level = level;
}

static void pivot_profile_add_node(pivot_profile *P, const pivot_node_profile *N)
{
    if (!P || !P->enabled) return;
#if LDLT_HAVE_OPENMP
    omp_set_lock(&P->lock);
#endif
    pivot_profile_level *L = &P->cur;
    L->nodes++;
    L->lapack_try += N->lapack_try;
    L->lapack_accept += N->lapack_accept;
    L->fallback_info += N->fallback_info;
    L->fallback_ipiv += N->fallback_ipiv;
    L->fallback_l21 += N->fallback_l21;
    L->scalar_nodes += N->scalar_used;
    L->scalar_1x1 += N->scalar_1x1;
    L->scalar_2x2 += N->scalar_2x2;
    L->pivots += N->pivots;
    L->delayed_vars += N->delayed_vars;
    if (N->R > L->max_R) L->max_R = N->R;
    if (N->cand > L->max_cand) L->max_cand = N->cand;
    if (N->noff_initial > L->max_noff) L->max_noff = N->noff_initial;
    L->t_total += N->t_total;
    L->t_assemble += N->t_assemble;
    L->t_lapack += N->t_lapack;
    L->t_scalar += N->t_scalar;
    L->t_scalar_select += N->t_scalar_select;
    L->t_scalar_update += N->t_scalar_update;
    L->t_store += N->t_store;
#if LDLT_HAVE_OPENMP
    omp_unset_lock(&P->lock);
#endif
}

static void pivot_profile_end_level(pivot_profile *P)
{
    if (!P || !P->enabled) return;
    const pivot_profile_level *L = &P->cur;
    double elapsed = prof_now_sec() - P->t_start;
    fprintf(stderr,
            "ldlt profile: level=%d elapsed=%.3fs nodes=%lld pivots=%lld "
            "lapack=%lld accept=%lld scalar=%lld scalar_piv(1=%lld 2=%lld) "
            "fb(info=%lld ipiv=%lld l21=%lld) "
            "delayed=%lld max(R=%d cand=%d noff=%d) "
            "time(total=%.3f asm=%.3f lapack=%.3f scalar=%.3f "
            "select=%.3f update=%.3f store=%.3f)s\n",
            (int)L->level, elapsed, (long long)L->nodes, (long long)L->pivots,
            (long long)L->lapack_try, (long long)L->lapack_accept,
            (long long)L->scalar_nodes, (long long)L->scalar_1x1,
            (long long)L->scalar_2x2, (long long)L->fallback_info,
            (long long)L->fallback_ipiv, (long long)L->fallback_l21,
            (long long)L->delayed_vars, (int)L->max_R, (int)L->max_cand,
            (int)L->max_noff, L->t_total, L->t_assemble, L->t_lapack,
            L->t_scalar, L->t_scalar_select, L->t_scalar_update, L->t_store);
}

/* LAPACK dsytrf: blocked Bunch-Kaufman for symmetric indefinite matrices.
   Uses Fortran calling convention (trailing underscore, column-major). */
#if LDLT_USE_SYSTEM_BLAS
/* LAPACK dsytrf: blocked Bunch-Kaufman.
   Provide a forward declaration via the standard Fortran ABI. */
extern void dsytrf_(char *uplo, int *n, double *a, int *lda,
                    int *ipiv, double *work, int *lwork, int *info);
static void call_dsytrf(int n, double *A, int lda, int *ipiv,
                        double *work, int lwork, int *info)
{
    char uplo = 'L';
    dsytrf_(&uplo, &n, A, &lda, ipiv, work, &lwork, info);
}
#else
#define LDLT_NO_LAPACK 1
#endif

/* Scatter the permuted matrix's lower-triangle entries for cols [fc, fc+k) into F. */
static int scatter_B(int32_t fc, int32_t k,
                      const int32_t *Bp, const int32_t *Bi, const double *Bx,
                      const int32_t *rowmap,
                      double *F, int32_t ldF)
{
    int missed = 0;
    for (int32_t lc = 0; lc < k; ++lc) {
        int32_t col = fc + lc;
        for (int32_t p = Bp[col]; p < Bp[col+1]; ++p) {
            int32_t i = Bi[p];
            int32_t ipos = rowmap[i];
            if (ipos < 0) { missed++; continue; }
            F[ipos + (size_t)lc*ldF] += Bx[p];
        }
    }
    return missed;
}

/* Apply update from supernode d to s.
   Specializes the common kd=1 (width-1 supernode) case to a scalar rank-1 update,
   which is the dominant case for AMD-ordered SPD matrices and avoids 1x1 BLAS overhead. */
static void apply_update(const ldlt_symbolic *S, const ldlt_numeric *N,
                         int32_t d, int32_t s_width,
                         const int32_t *rowmap,
                         double *F, int32_t ldF,
                         double *scratch)
{
    const ldlt_super *Sd = &S->super[d];
    int32_t kd = Sd->width;
    int32_t md = Sd->nrows_below;
    const double *Ld = N->L + S->panel_off[d];     /* (kd+md) x kd col-major */
    const double *Dd = N->D + S->d_off[d];
    int32_t ldLd = kd + md;

    const int32_t *rb = S->rows_pool + Sd->row_off;

    /* Collect active rows of d's below-diag panel that fall in s's index set. */
    int32_t *top_src = (int32_t*)scratch;
    int32_t *off_src = top_src + md;
    int32_t *top_pos = off_src + md;
    int32_t *off_pos = top_pos + md;

    int32_t ntop = 0, noff = 0;
    for (int32_t t = 0; t < md; ++t) {
        int32_t r = rb[t];
        int32_t pos = rowmap[r];
        if (pos < 0) continue;
        if (pos < s_width) { top_src[ntop] = t; top_pos[ntop] = pos; ntop++; }
        else               { off_src[noff] = t; off_pos[noff] = pos; noff++; }
    }
    if (ntop == 0) return;

    /* ---- Fast path: kd == 1 (width-1 contributing supernode) — pure rank-1 update. ---- */
    if (kd == 1) {
        double dscale = Dd[0];
        /* Gather Ld below-diag entries at active rows */
        double *vt = (double*)(off_pos + md);            /* ntop doubles */
        double *vo = vt + ntop;                           /* noff doubles */
        for (int32_t i = 0; i < ntop; ++i) vt[i] = Ld[kd + top_src[i]];
        for (int32_t i = 0; i < noff; ++i) vo[i] = Ld[kd + off_src[i]];
        /* T11[i,j] = vt[i] * dscale * vt[j]; subtract from F[top_pos[i], top_pos[j]] */
        for (int32_t j = 0; j < ntop; ++j) {
            double aj = dscale * vt[j];
            int32_t fc = top_pos[j];
            for (int32_t i = 0; i < ntop; ++i) {
                F[top_pos[i] + (size_t)fc*ldF] -= vt[i] * aj;
            }
            for (int32_t i = 0; i < noff; ++i) {
                F[off_pos[i] + (size_t)fc*ldF] -= vo[i] * aj;
            }
        }
        return;
    }

    if (kd <= 8) {
        for (int32_t j = 0; j < ntop; ++j) {
            int32_t srcj = top_src[j];
            int32_t fc = top_pos[j];
            double wj[8];
            for (int32_t c = 0; c < kd; ++c)
                wj[c] = Dd[c] * Ld[(kd + srcj) + (size_t)c*ldLd];
            for (int32_t i = 0; i < ntop; ++i) {
                int32_t srci = top_src[i];
                double acc = 0.0;
                for (int32_t c = 0; c < kd; ++c)
                    acc += Ld[(kd + srci) + (size_t)c*ldLd] * wj[c];
                F[top_pos[i] + (size_t)fc*ldF] -= acc;
            }
            for (int32_t i = 0; i < noff; ++i) {
                int32_t srci = off_src[i];
                double acc = 0.0;
                for (int32_t c = 0; c < kd; ++c)
                    acc += Ld[(kd + srci) + (size_t)c*ldLd] * wj[c];
                F[off_pos[i] + (size_t)fc*ldF] -= acc;
            }
        }
        return;
    }

    /* ---- General path via dgemm ---- */
    double *Ubuf = (double*)(off_pos + md);
    double *Utop = Ubuf;                              /* ntop x kd */
    double *Uoff = Ubuf + (size_t)ntop * kd;          /* noff x kd */
    for (int32_t t = 0; t < ntop; ++t) {
        int32_t srct = top_src[t];
        for (int32_t c = 0; c < kd; ++c)
            Utop[t + (size_t)c*ntop] = Ld[(kd + srct) + (size_t)c*ldLd];
    }
    for (int32_t t = 0; t < noff; ++t) {
        int32_t srct = off_src[t];
        for (int32_t c = 0; c < kd; ++c)
            Uoff[t + (size_t)c*noff] = Ld[(kd + srct) + (size_t)c*ldLd];
    }
    double *Wtop = Uoff + (size_t)noff * kd;          /* ntop x kd */
    for (int32_t c = 0; c < kd; ++c) {
        double dc = Dd[c];
        for (int32_t t = 0; t < ntop; ++t)
            Wtop[t + (size_t)c*ntop] = Utop[t + (size_t)c*ntop] * dc;
    }
    double *T11 = Wtop + (size_t)ntop * kd;           /* ntop x ntop */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans,
                ntop, ntop, kd, 1.0, Wtop, ntop, Utop, ntop, 0.0, T11, ntop);
    double *T21 = T11 + (size_t)ntop * ntop;          /* noff x ntop */
    if (noff > 0) {
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans,
                    noff, ntop, kd, 1.0, Uoff, noff, Wtop, ntop, 0.0, T21, noff);
    }
    for (int32_t j = 0; j < ntop; ++j) {
        int32_t fc = top_pos[j];
        for (int32_t i = 0; i < ntop; ++i)
            F[top_pos[i] + (size_t)fc*ldF] -= T11[i + (size_t)j*ntop];
        for (int32_t i = 0; i < noff; ++i)
            F[off_pos[i] + (size_t)fc*ldF] -= T21[i + (size_t)j*noff];
    }
}

ldlt_status ldlt_factor_supernode(const ldlt_symbolic *S, ldlt_numeric *N,
                                  int32_t s, double *workspace, int32_t *rowmap,
                                  const ldlt_options *opt,
                                  const int32_t *Bp, const int32_t *Bi, const double *Bx)
{
    const ldlt_super *Ss = &S->super[s];
    int32_t fc = Ss->first_col;
    int32_t k = Ss->width;
    int32_t m = Ss->nrows_below;
    int32_t R = k + m;
    double *F = N->L + S->panel_off[s];
    int32_t ldF = R;

    for (int32_t i = 0; i < k; ++i) rowmap[fc + i] = i;
    const int32_t *rb = S->rows_pool + Ss->row_off;
    for (int32_t i = 0; i < m; ++i) rowmap[rb[i]] = k + i;

    for (int64_t i = 0; i < (int64_t)R * k; ++i) F[i] = 0.0;

    scatter_B(fc, k, Bp, Bi, Bx, rowmap, F, ldF);

    int32_t up_lo = S->upd_ptr[s], up_hi = S->upd_ptr[s+1];
    for (int32_t u = up_lo; u < up_hi; ++u) {
        int32_t d = S->upd_idx[u];
        apply_update(S, N, d, k, rowmap, F, ldF, workspace);
    }

    int indef_local = 0;
    double tiny = 1e-300;
    ldlt_status st = ldlt_dense_ldl_nopiv(F, ldF, R, k, N->D + S->d_off[s], tiny, &indef_local);
    if (st != LDLT_OK) return st;
    if (indef_local && opt && opt->require_pd) return LDLT_ERR_INDEFINITE;

    for (int32_t i = 0; i < k; ++i) rowmap[fc + i] = -1;
    for (int32_t i = 0; i < m; ++i) rowmap[rb[i]] = -1;

    return LDLT_OK;
}

static int32_t edge_run_end(const int32_t *emap, int32_t n, int32_t start)
{
    int32_t end = start + 1;
    while (end < n && emap[end] == emap[end - 1] + 1) end++;
    return end;
}

static size_t packed_lower_col(int32_t n, int32_t j)
{
    return (size_t)j * (size_t)n - ((size_t)j * (size_t)(j - 1)) / 2;
}

/* Inlined adds: per-column lengths in extend-add are typically small (often
   < 32), so BLAS call overhead dominated. Plain loops autovectorize and beat
   cblas_daxpy here. */
static void add_contiguous_lower_block(double *F, int32_t ldF, const double *C, int32_t nC,
                                       int32_t f0, int32_t c0, int32_t len)
{
    for (int32_t j = 0; j < len; ++j) {
        const double *src = C + packed_lower_col(nC, c0 + j);
        double *dst = F + (f0 + j) + (size_t)(f0 + j)*ldF;
        int32_t n = len - j;
        for (int32_t i = 0; i < n; ++i) dst[i] += src[i];
    }
}

static void add_contiguous_rect_block(double *F, int32_t ldF, const double *C, int32_t nC,
                                      int32_t fi0, int32_t fj0,
                                      int32_t ci0, int32_t cj0,
                                      int32_t ni, int32_t nj)
{
    for (int32_t j = 0; j < nj; ++j) {
        int32_t col = cj0 + j;
        const double *src = C + packed_lower_col(nC, col) + (ci0 - col);
        double *dst = F + fi0 + (size_t)(fj0 + j)*ldF;
        for (int32_t i = 0; i < ni; ++i) dst[i] += src[i];
    }
}

static void extend_add_contribution(double *F, int32_t ldF,
                                    const double *C, int32_t nC,
                                    const int32_t *emap, int32_t n)
{
    for (int32_t j0 = 0; j0 < n; ) {
        int32_t j1 = edge_run_end(emap, n, j0);
        int32_t fj0 = emap[j0];

        add_contiguous_lower_block(F, ldF, C, nC, fj0, j0, j1 - j0);

        for (int32_t i0 = j1; i0 < n; ) {
            int32_t i1 = edge_run_end(emap, n, i0);
            add_contiguous_rect_block(F, ldF, C, nC,
                                      emap[i0], fj0, i0, j0, i1 - i0, j1 - j0);
            i0 = i1;
        }
        j0 = j1;
    }
}

typedef struct {
    double *front;
    size_t front_cap;
} front_workspace;

static double *front_workspace_get(front_workspace *ws, size_t need)
{
    if (ws->front_cap < need) {
        ws->front = (double*)ldlt_xrealloc(ws->front, need * sizeof(double));
        ws->front_cap = need;
    }
    return ws->front;
}

static void zero_front_lower(double *F, int32_t R)
{
    if (R < 128) {
        memset(F, 0, (size_t)R * R * sizeof(double));
        return;
    }
    for (int32_t j = 0; j < R; ++j)
        memset(F + j + (size_t)j*R, 0, (size_t)(R - j) * sizeof(double));
}

static ldlt_status factor_front_node(const ldlt_symbolic *S, ldlt_numeric *N,
                                     int32_t s, double **contrib,
                                     int32_t *rowmap, double *work,
                                     front_workspace *front_ws,
                                     const ldlt_options *opt,
                                     const int32_t *Bp, const int32_t *Bi,
                                     const double *Bx)
{
    const ldlt_super *Ss = &S->super[s];
    int32_t fc = Ss->first_col;
    int32_t k = Ss->width;
    int32_t m = Ss->nrows_below;
    int32_t R = k + m;
    const int32_t *rb = S->rows_pool + Ss->row_off;

    double *F = front_workspace_get(front_ws, (size_t)R * R);
    zero_front_lower(F, R);
    for (int32_t i = 0; i < k; ++i) rowmap[fc + i] = i;
    for (int32_t i = 0; i < m; ++i) rowmap[rb[i]] = k + i;

    scatter_B(fc, k, Bp, Bi, Bx, rowmap, F, R);

    for (int32_t cp = S->child_ptr[s]; cp < S->child_ptr[s+1]; ++cp) {
        int32_t c = S->child_idx[cp];
        const ldlt_super *Sc = &S->super[c];
        int32_t cm = Sc->nrows_below;
        double *C = contrib[c];
        const int32_t *emap = S->edge_pos + S->edge_ptr[c];
        if (C) {
            extend_add_contribution(F, R, C, cm, emap, cm);
            free(C);
            contrib[c] = NULL;
        }
    }

    int indef_local = 0;
    ldlt_status st = ldlt_dense_front_nopiv(F, R, R, k, N->D + S->d_off[s],
                                            work, 1e-300, &indef_local);
    if (st == LDLT_OK && indef_local && opt && opt->require_pd)
        st = LDLT_ERR_INDEFINITE;

    if (st == LDLT_OK) {
        double *L = N->L + S->panel_off[s];
        for (int32_t j = 0; j < k; ++j)
            memcpy(L + (size_t)j*R, F + (size_t)j*R, (size_t)R * sizeof(double));
        if (m > 0) {
            double *C = (double*)ldlt_xmalloc(((size_t)m * (size_t)(m + 1) / 2) *
                                              sizeof(double));
            for (int32_t j = 0; j < m; ++j)
                memcpy(C + packed_lower_col(m, j), F + k + j + (size_t)(k+j)*R,
                       (size_t)(m - j) * sizeof(double));
            contrib[s] = C;
        }
    }

    for (int32_t i = 0; i < k; ++i) rowmap[fc + i] = -1;
    for (int32_t i = 0; i < m; ++i) rowmap[rb[i]] = -1;
    return st;
}

/* Compute supernode levels in the assembly tree. Children are level 0; a parent
   is one level above the deepest child. */
static void compute_levels(const ldlt_symbolic *S,
                           int32_t **level_ptr_out, int32_t **level_idx_out,
                           int32_t *nlevels_out)
{
    int32_t nsuper = S->nsuper;
    int32_t *lvl = (int32_t*)ldlt_xcalloc((size_t)(nsuper>0?nsuper:1), sizeof(int32_t));
    int32_t maxlvl = 0;
    for (int32_t s = 0; s < nsuper; ++s) {
        int32_t p = S->super_parent ? S->super_parent[s] : -1;
        if (p >= 0 && lvl[p] < lvl[s] + 1) {
            lvl[p] = lvl[s] + 1;
            if (lvl[p] > maxlvl) maxlvl = lvl[p];
        }
    }
    int32_t nlevels = maxlvl + 1;
    int32_t *lp = (int32_t*)ldlt_xcalloc((size_t)nlevels+1, sizeof(int32_t));
    for (int32_t s = 0; s < nsuper; ++s) lp[lvl[s]+1]++;
    for (int32_t i = 0; i < nlevels; ++i) lp[i+1] += lp[i];
    int32_t *li = (int32_t*)ldlt_xmalloc((size_t)(nsuper>0?nsuper:1) * sizeof(int32_t));
    int32_t *cur = (int32_t*)ldlt_xmalloc((size_t)nlevels * sizeof(int32_t));
    memcpy(cur, lp, (size_t)nlevels * sizeof(int32_t));
    for (int32_t s = 0; s < nsuper; ++s) li[cur[lvl[s]]++] = s;
    free(cur); free(lvl);
    *level_ptr_out = lp; *level_idx_out = li; *nlevels_out = nlevels;
}

typedef struct {
    int32_t nvars;
    int32_t *vars;
    double *A; /* packed lower, in vars order */
} pivot_contrib;

typedef struct {
    ldlt_numeric *N;
    int64_t d_cap;
    int64_t l_row_cap;
    int64_t l_val_cap;
    int indefinite;
} pivot_build;

static void pivot_contrib_free(pivot_contrib *C)
{
    if (!C) return;
    free(C->vars);
    free(C->A);
    free(C);
}

static void pivot_numeric_init(ldlt_numeric *N, int32_t n)
{
    int32_t nalloc = n > 0 ? n : 1;
    N->pivoted = 1;
    N->piv_size = (int32_t*)ldlt_xcalloc((size_t)nalloc, sizeof(int32_t));
    N->piv_var_off = (int32_t*)ldlt_xcalloc((size_t)nalloc + 1, sizeof(int32_t));
    N->piv_vars = (int32_t*)ldlt_xcalloc((size_t)nalloc, sizeof(int32_t));
    N->piv_d_off = (int64_t*)ldlt_xcalloc((size_t)nalloc + 1, sizeof(int64_t));
    N->piv_l_row_off = (int64_t*)ldlt_xcalloc((size_t)nalloc + 1, sizeof(int64_t));
    N->piv_l_val_off = (int64_t*)ldlt_xcalloc((size_t)nalloc + 1, sizeof(int64_t));
    N->piv_D = (double*)ldlt_xcalloc((size_t)(2 * nalloc), sizeof(double));
    N->piv_l_rows = (int32_t*)ldlt_xcalloc((size_t)nalloc, sizeof(int32_t));
    N->piv_l_vals = (double*)ldlt_xcalloc((size_t)nalloc, sizeof(double));
}

static void ensure_i32_cap(int32_t **p, int64_t *cap, int64_t need)
{
    if (*cap >= need) return;
    int64_t ncap = *cap > 0 ? *cap : 1;
    while (ncap < need) ncap *= 2;
    *p = (int32_t*)ldlt_xrealloc(*p, (size_t)ncap * sizeof(int32_t));
    *cap = ncap;
}

static void ensure_double_cap(double **p, int64_t *cap, int64_t need)
{
    if (*cap >= need) return;
    int64_t ncap = *cap > 0 ? *cap : 1;
    while (ncap < need) ncap *= 2;
    *p = (double*)ldlt_xrealloc(*p, (size_t)ncap * sizeof(double));
    *cap = ncap;
}

static void append_pivot_block(pivot_build *B, int32_t bs, const int32_t *pvars,
                               const double *Dblk, int32_t nrows,
                               const int32_t *rows, const double *lvals)
{
    ldlt_numeric *N = B->N;
    int32_t b = N->npivots;
    int32_t vo = N->piv_var_off[b];
    int64_t doff = N->piv_d_off[b];
    int64_t roff = N->piv_l_row_off[b];
    int64_t vloff = N->piv_l_val_off[b];
    int32_t dlen = (bs == 1) ? 1 : 3;
    int64_t vlen = (int64_t)nrows * bs;

    ensure_double_cap(&N->piv_D, &B->d_cap, doff + dlen);
    ensure_i32_cap(&N->piv_l_rows, &B->l_row_cap, roff + nrows);
    ensure_double_cap(&N->piv_l_vals, &B->l_val_cap, vloff + vlen);

    N->piv_size[b] = bs;
    for (int32_t i = 0; i < bs; ++i) N->piv_vars[vo + i] = pvars[i];
    for (int32_t i = 0; i < dlen; ++i) N->piv_D[doff + i] = Dblk[i];
    for (int32_t i = 0; i < nrows; ++i) N->piv_l_rows[roff + i] = rows[i];
    for (int64_t i = 0; i < vlen; ++i) N->piv_l_vals[vloff + i] = lvals[i];

    if (bs == 1) {
        if (Dblk[0] <= 0.0) B->indefinite = 1;
    } else {
        double det = Dblk[0] * Dblk[2] - Dblk[1] * Dblk[1];
        if (!(Dblk[0] > 0.0 && det > 0.0)) B->indefinite = 1;
    }

    N->npivots++;
    N->n_piv_vars += bs;
    N->piv_var_off[b + 1] = vo + bs;
    N->piv_d_off[b + 1] = doff + dlen;
    N->piv_l_row_off[b + 1] = roff + nrows;
    N->piv_l_val_off[b + 1] = vloff + vlen;
}

static void append_pivot_block_locked(pivot_build *B, int32_t bs,
                                      const int32_t *pvars,
                                      const double *Dblk, int32_t nrows,
                                      const int32_t *rows,
                                      const double *lvals,
                                      void *append_lock_ptr)
{
#if LDLT_HAVE_OPENMP
    omp_lock_t *lock = (omp_lock_t*)append_lock_ptr;
    if (lock) omp_set_lock(lock);
#else
    (void)append_lock_ptr;
#endif
    append_pivot_block(B, bs, pvars, Dblk, nrows, rows, lvals);
#if LDLT_HAVE_OPENMP
    if (lock) omp_unset_lock(lock);
#endif
}

static void dense_add_sym(double *F, int32_t ldF, int32_t i, int32_t j, double x)
{
    F[i + (size_t)j*ldF] += x;
    if (i != j) F[j + (size_t)i*ldF] += x;
}

static void swap_front_pos(double *F, int32_t R, int32_t *vars, int32_t a, int32_t b)
{
    if (a == b) return;
    int32_t tv = vars[a];
    vars[a] = vars[b];
    vars[b] = tv;
    for (int32_t j = 0; j < R; ++j) {
        double t = F[a + (size_t)j*R];
        F[a + (size_t)j*R] = F[b + (size_t)j*R];
        F[b + (size_t)j*R] = t;
    }
    for (int32_t i = 0; i < R; ++i) {
        double t = F[i + (size_t)a*R];
        F[i + (size_t)a*R] = F[i + (size_t)b*R];
        F[i + (size_t)b*R] = t;
    }
}

static void bring_pivot_pair_to_front(double *F, int32_t R, int32_t *vars,
                                      int32_t elim, int32_t *p, int32_t *q,
                                      int32_t bs)
{
    if (*p != elim) {
        int32_t oldp = *p;
        swap_front_pos(F, R, vars, oldp, elim);
        if (*q == elim) *q = oldp;
        else if (*q == oldp) *q = elim;
        *p = elim;
    }
    if (bs == 2 && *q != elim + 1) {
        swap_front_pos(F, R, vars, *q, elim + 1);
        *q = elim + 1;
    }
}

static ldlt_status select_pivot_block(const double *F, int32_t R,
                                      int32_t elim, int32_t cand_end,
                                      double alpha, double tiny,
                                      int32_t *piv0, int32_t *piv1,
                                      int32_t *bs)
{
    *bs = 0;
    *piv0 = *piv1 = -1;

    /* Fast path: classical Bunch-Kaufman first probe — try the leading column.
       If |F[elim,elim]| dominates its column off-diagonals by alpha, use a 1x1
       pivot at elim without scanning the full trailing matrix. This is the
       common case for well-scaled matrices and turns pivot selection from
       O(R^2) per step into O(R). */
    {
        double diag = fabs(F[elim + (size_t)elim*R]);
        double colmax = 0.0;
        const double *col = F + (size_t)elim * R;
        int finite = 1;
        for (int32_t i = elim + 1; i < R; ++i) {
            double a = col[i];
            if (!isfinite(a)) { finite = 0; break; }
            double aa = fabs(a);
            if (aa > colmax) colmax = aa;
        }
        if (finite && isfinite(diag) && diag > tiny &&
            (colmax == 0.0 || diag >= alpha * colmax)) {
            *piv0 = elim;
            *bs = 1;
            return LDLT_OK;
        }
    }

    /* Common case: find a stable 1x1 pivot by scanning candidate columns only.
       Avoid the full trailing R x R scale scan unless we need the 2x2 search. */
    for (int32_t j = elim; j < cand_end; ++j) {
        double diag = fabs(F[j + (size_t)j*R]);
        double colmax = 0.0;
        for (int32_t i = elim; i < R; ++i) {
            if (i == j) continue;
            double fij = F[i + (size_t)j*R];
            if (!isfinite(fij)) return LDLT_ERR_SINGULAR;
            double a = fabs(fij);
            if (a > colmax) colmax = a;
        }
        if (isfinite(diag) && diag > tiny &&
            (colmax <= tiny || diag >= alpha * colmax)) {
            *piv0 = j;
            *bs = 1;
            return LDLT_OK;
        }
    }

    double scale = 0.0;
    for (int32_t j = elim; j < R; ++j) {
        for (int32_t i = elim; i < R; ++i) {
            double a = F[i + (size_t)j*R];
            if (!isfinite(a)) return LDLT_ERR_SINGULAR;
            double aa = fabs(a);
            if (aa > scale) scale = aa;
        }
    }
    double eps = tiny * (scale > 1.0 ? scale : 1.0);

    double best = 0.0;
    int32_t bi = -1, bj = -1;
    double max_mult = 1.0 / alpha;
    for (int32_t j = elim; j < cand_end; ++j) {
        for (int32_t i = j + 1; i < cand_end; ++i) {
            double a = fabs(F[i + (size_t)j*R]);
            if (a <= best || a <= eps) continue;

            double d00 = F[j + (size_t)j*R];
            double d10 = F[i + (size_t)j*R];
            double d11 = F[i + (size_t)i*R];
            double det = d00 * d11 - d10 * d10;
            double dscale = fabs(d00);
            if (fabs(d10) > dscale) dscale = fabs(d10);
            if (fabs(d11) > dscale) dscale = fabs(d11);
            if (dscale < 1.0) dscale = 1.0;
            if (fabs(det) <= eps * dscale) continue;

            int stable = 1;
            for (int32_t r = elim; r < R && stable; ++r) {
                if (r == i || r == j) continue;
                double a0 = F[r + (size_t)j*R];
                double a1 = F[r + (size_t)i*R];
                double l0 = (a0 * d11 - a1 * d10) / det;
                double l1 = (-a0 * d10 + a1 * d00) / det;
                if (fabs(l0) > max_mult || fabs(l1) > max_mult)
                    stable = 0;
            }
            if (stable) {
                best = a;
                bi = i;
                bj = j;
            }
        }
    }
    if (best > 0.0) {
        *piv0 = bj;
        *piv1 = bi;
        *bs = 2;
    }
    return LDLT_OK;
}

static void pivot_debug_front(const char *why, int32_t s, const double *F, int32_t R,
                              int32_t elim, int32_t cand)
{
    if (!getenv("LDLT_DEBUG_PIVOT")) return;
    double maxdiag = 0.0, maxoff = 0.0;
    int32_t diag_i = -1, off_i = -1, off_j = -1;
    for (int32_t j = elim; j < cand; ++j) {
        double d = fabs(F[j + (size_t)j*R]);
        if (d > maxdiag) {
            maxdiag = d;
            diag_i = j;
        }
        for (int32_t i = j + 1; i < cand; ++i) {
            double a = fabs(F[i + (size_t)j*R]);
            if (a > maxoff) {
                maxoff = a;
                off_i = i;
                off_j = j;
            }
        }
    }
    fprintf(stderr,
            "ldlt pivot debug: %s s=%d R=%d elim=%d cand=%d rem=%d "
            "maxdiag=%g at %d maxoff=%g at (%d,%d)\n",
            why, s, R, elim, cand, R - elim,
            maxdiag, diag_i, maxoff, off_i, off_j);
}

static int add_front_var(int32_t v, int32_t *map, int32_t *vars,
                         int32_t *nvars, int32_t *touched, int32_t *ntouched)
{
    if (map[v] >= 0) return 0;
    map[v] = *nvars;
    vars[(*nvars)++] = v;
    touched[(*ntouched)++] = v;
    return 1;
}

static void scatter_pivot_front(int32_t fc, int32_t k,
                                const int32_t *Bp, const int32_t *Bi,
                                const double *Bx, const int32_t *map,
                                double *F, int32_t R)
{
    for (int32_t lc = 0; lc < k; ++lc) {
        int32_t col = fc + lc;
        int32_t jpos = map[col];
        if (jpos < 0) continue;
        for (int32_t p = Bp[col]; p < Bp[col+1]; ++p) {
            int32_t ipos = map[Bi[p]];
            if (ipos < 0) continue;
            dense_add_sym(F, R, ipos, jpos, Bx[p]);
        }
    }
}

static void add_child_contrib_to_front(const pivot_contrib *C, const int32_t *map,
                                       double *F, int32_t R)
{
    if (!C) return;
    for (int32_t j = 0; j < C->nvars; ++j) {
        int32_t fj = map[C->vars[j]];
        if (fj < 0) continue;
        const double *col = C->A + packed_lower_col(C->nvars, j);
        for (int32_t i = j; i < C->nvars; ++i) {
            int32_t fi = map[C->vars[i]];
            if (fi < 0) continue;
            dense_add_sym(F, R, fi, fj, col[i - j]);
        }
    }
}

typedef struct {
    double  *F;        int64_t F_cap;
    double  *Fsave;    int64_t Fsave_cap;   /* original front for LAPACK fallback */
    int     *ipiv;     int64_t ipiv_cap;     /* LAPACK IPIV */
    double  *dswork;   int64_t dswork_cap;   /* LAPACK dsytrf workspace */
    double  *A21;      int64_t A21_cap;      /* L21 buffer, noff x cand col-major */
    double  *V;        int64_t V_cap;        /* V = A21 * D^{-1}, noff x cand col-major */
    double  *lvals;    int64_t lvals_cap;    /* L column values, row-major (t*bs+r) */
    double  *l_cm;     int64_t l_cm_cap;     /* scalar path L buffer, column-major */
    double  *w_cm;     int64_t w_cm_cap;     /* scalar path D*L' buffer, column-major */
    int32_t *rows;     int64_t rows_cap;     /* row IDs buffer */
    int32_t *gather;   int64_t gather_cap;
    int32_t *touched;  int64_t touched_cap;
    int32_t *vars;     int64_t vars_cap;
} pivot_workspace;

static void pivot_workspace_free(pivot_workspace *W)
{
    free(W->F); free(W->Fsave); free(W->ipiv); free(W->dswork);
    free(W->A21); free(W->V); free(W->lvals); free(W->l_cm); free(W->w_cm);
    free(W->rows);
    free(W->gather); free(W->touched); free(W->vars);
    memset(W, 0, sizeof(*W));
}

static void *grow_buf(void *p, int64_t *cap, int64_t need, size_t elem)
{
    if (*cap >= need) return p;
    int64_t nc = *cap > 0 ? *cap : 1;
    while (nc < need) nc *= 2;
    *cap = nc;
    return ldlt_xrealloc(p, (size_t)nc * elem);
}

static ldlt_status factor_pivot_node(const ldlt_symbolic *S, pivot_build *B,
                                     int32_t s, pivot_contrib **contrib,
                                     int32_t *map, pivot_workspace *W,
                                     double alpha, double tiny,
                                     const int32_t *Bp, const int32_t *Bi,
                                     const double *Bx,
                                     void *append_lock_ptr,
                                     pivot_profile *profile)
{
    pivot_node_profile prof;
    memset(&prof, 0, sizeof(prof));
    double t_node0 = profile && profile->enabled ? prof_now_sec() : 0.0;
    double t0 = t_node0;
    const ldlt_super *Ss = &S->super[s];
    int32_t fc = Ss->first_col;
    int32_t k = Ss->width;
    int32_t m = Ss->nrows_below;
    const int32_t *rb = S->rows_pool + Ss->row_off;

    int32_t cap = k + m;
    for (int32_t cp = S->child_ptr[s]; cp < S->child_ptr[s+1]; ++cp) {
        int32_t c = S->child_idx[cp];
        if (contrib[c]) cap += contrib[c]->nvars;
    }
    if (cap < 1) cap = 1;
    W->gather  = (int32_t*)grow_buf(W->gather,  &W->gather_cap,  cap, sizeof(int32_t));
    W->touched = (int32_t*)grow_buf(W->touched, &W->touched_cap, cap, sizeof(int32_t));
    int32_t *gather = W->gather;
    int32_t *touched = W->touched;
    int32_t ngather = 0, ntouched = 0;

    for (int32_t i = 0; i < k; ++i)
        add_front_var(fc + i, map, gather, &ngather, touched, &ntouched);
    for (int32_t i = 0; i < m; ++i)
        add_front_var(rb[i], map, gather, &ngather, touched, &ntouched);
    for (int32_t cp = S->child_ptr[s]; cp < S->child_ptr[s+1]; ++cp) {
        int32_t c = S->child_idx[cp];
        pivot_contrib *C = contrib[c];
        if (!C) continue;
        for (int32_t i = 0; i < C->nvars; ++i)
            add_front_var(C->vars[i], map, gather, &ngather, touched, &ntouched);
    }

    int64_t need_vars = ngather > 0 ? ngather : 1;
    W->vars = (int32_t*)grow_buf(W->vars, &W->vars_cap, need_vars, sizeof(int32_t));
    int32_t *vars = W->vars;
    int32_t cand = 0, pos = 0;
    for (int32_t i = 0; i < ngather; ++i) {
        int32_t v = gather[i];
        if (S->col_to_super[v] <= s) vars[cand++] = v;
    }
    pos = cand;
    for (int32_t i = 0; i < ngather; ++i) {
        int32_t v = gather[i];
        if (S->col_to_super[v] > s) vars[pos++] = v;
    }
    for (int32_t i = 0; i < ngather; ++i) map[vars[i]] = i;

    int32_t R = ngather;
    int64_t F_need = (int64_t)(R > 0 ? R : 1) * (int64_t)(R > 0 ? R : 1);
    W->F = (double*)grow_buf(W->F, &W->F_cap, F_need, sizeof(double));
    double *F = W->F;
    memset(F, 0, (size_t)F_need * sizeof(double));

    scatter_pivot_front(fc, k, Bp, Bi, Bx, map, F, R);
    for (int32_t cp = S->child_ptr[s]; cp < S->child_ptr[s+1]; ++cp) {
        int32_t c = S->child_idx[cp];
        add_child_contrib_to_front(contrib[c], map, F, R);
        pivot_contrib_free(contrib[c]);
        contrib[c] = NULL;
    }
    if (profile && profile->enabled) {
        prof.t_assemble += prof_now_sec() - t0;
        prof.R = R;
        prof.cand = cand;
        prof.noff_initial = R - cand;
    }

    ldlt_status st = LDLT_OK;
    int32_t noff = R - cand; /* Schur complement size */
    int32_t rem_start = cand;

    if (cand == 0) goto store_contribution;

#ifndef LDLT_NO_LAPACK
    {
    /* --- LAPACK dsytrf path: blocked Bunch-Kaufman (BLAS3 performance) ---
       Accepted only for no-swap 1x1 pivot sequences. When off-front rows are
       present, explicitly check L21 multipliers before committing the update. */
        int n_cand = (int)cand;
        int lda_f  = (int)R;
        if (profile && profile->enabled) {
            prof.lapack_try = 1;
            t0 = prof_now_sec();
        }

        /* Grow workspace buffers. */
        W->ipiv = (int*)grow_buf(W->ipiv, &W->ipiv_cap, cand, sizeof(int));
        W->Fsave = (double*)grow_buf(W->Fsave, &W->Fsave_cap, F_need, sizeof(double));
        memcpy(W->Fsave, F, (size_t)F_need * sizeof(double));

        /* Workspace query. */
        int lwork = -1;
        double lwork_d = 0.0;
        int info = 0;
        call_dsytrf(n_cand, F, lda_f, W->ipiv, &lwork_d, lwork, &info);
        lwork = (lwork_d > 0) ? (int)lwork_d : 4 * n_cand;
        W->dswork = (double*)grow_buf(W->dswork, &W->dswork_cap, lwork, sizeof(double));

        /* Factorize: F[0:cand, 0:cand] = P * L * D * L^T * P^T (in-place). */
        call_dsytrf(n_cand, F, lda_f, W->ipiv, W->dswork, lwork, &info);
        if (profile && profile->enabled) prof.t_lapack += prof_now_sec() - t0;
        if (info != 0) {
            if (profile && profile->enabled) prof.fallback_info = 1;
            memcpy(F, W->Fsave, (size_t)F_need * sizeof(double));
            goto scalar_path;
        }

        int trivial_ipiv = 1;
        for (int32_t ki = 0; ki < cand; ++ki) {
            if (W->ipiv[ki] != (int)ki + 1) {
                trivial_ipiv = 0;
                break;
            }
        }
        if (!trivial_ipiv) {
            if (profile && profile->enabled) prof.fallback_ipiv = 1;
            memcpy(F, W->Fsave, (size_t)F_need * sizeof(double));
            goto scalar_path;
        }

        /* Apply IPIV permutation to vars[0:cand].
           LAPACK stores 1-indexed IPIV. For 'L':
             ipiv[k]>0: 1x1 pivot, rows k and ipiv[k]-1 (0-indexed) were swapped.
             ipiv[k]<0: 2x2 pivot at k,k+1; rows k+1 and -ipiv[k]-1 (0-indexed) were swapped. */
        for (int32_t ki = 0; ki < cand; ) {
            int iv = W->ipiv[ki];
            if (iv > 0) {
                int kp = iv - 1;
                if (kp != ki) { int32_t t = vars[ki]; vars[ki] = vars[kp]; vars[kp] = t; }
                ki++;
            } else {
                int kp = (-iv) - 1;
                if (kp != ki + 1) { int32_t t = vars[ki+1]; vars[ki+1] = vars[kp]; vars[kp] = t; }
                ki += 2;
            }
        }

        if (noff > 0) {
            /* Copy A21 = F[cand:R, 0:cand] to a compact buffer (noff x cand col-major).
               dsytrf only touches F[0:cand, 0:cand]; F[cand:R, 0:cand] is unchanged. */
            int64_t a21_need = (int64_t)noff * cand;
            W->A21 = (double*)grow_buf(W->A21, &W->A21_cap, a21_need > 0 ? a21_need : 1, sizeof(double));
            double *A21 = W->A21;
            for (int32_t j = 0; j < cand; ++j)
                for (int32_t i = 0; i < noff; ++i)
                    A21[i + (size_t)j*noff] = F[(cand+i) + (size_t)j*R];

            /* Apply column permutation from IPIV to A21. */
            for (int32_t ki = 0; ki < cand; ) {
                int iv = W->ipiv[ki];
                if (iv > 0) {
                    int kp = iv - 1;
                    if (kp != ki) {
                        for (int32_t i = 0; i < noff; ++i) {
                            double t = A21[i + (size_t)ki*noff];
                            A21[i + (size_t)ki*noff] = A21[i + (size_t)kp*noff];
                            A21[i + (size_t)kp*noff] = t;
                        }
                    }
                    ki++;
                } else {
                    int kp = (-iv) - 1;
                    if (kp != ki + 1) {
                        for (int32_t i = 0; i < noff; ++i) {
                            double t = A21[i + (size_t)(ki+1)*noff];
                            A21[i + (size_t)(ki+1)*noff] = A21[i + (size_t)kp*noff];
                            A21[i + (size_t)kp*noff] = t;
                        }
                    }
                    ki += 2;
                }
            }

            /* DTRSM: A21 = A21 * L^{-T}  (L from dsytrf, unit lower triangular, lda=R in F).
               For 2x2 pivot blocks dsytrf stores D[k+1,k] in F[k+1,k]; that is NOT an L
               element. Zero those slots before DTRSM (saving values), restore after. */
            for (int32_t ki = 0; ki < cand; ) {
                if (W->ipiv[ki] > 0) { ki++; }
                else { W->dswork[ki] = F[ki+1 + (size_t)ki*R]; F[ki+1 + (size_t)ki*R] = 0.0; ki += 2; }
            }
            cblas_dtrsm(CblasColMajor, CblasRight, CblasLower, CblasTrans, CblasUnit,
                        noff, cand, 1.0, F, R, A21, noff);
            for (int32_t ki = 0; ki < cand; ) {
                if (W->ipiv[ki] > 0) { ki++; }
                else { F[ki+1 + (size_t)ki*R] = W->dswork[ki]; ki += 2; }
            }

            /* Form V = A21 * D^{-1}  (scale/transform columns by D^{-1} blocks). */
            W->V = (double*)grow_buf(W->V, &W->V_cap, a21_need > 0 ? a21_need : 1, sizeof(double));
            double *V = W->V;
            memcpy(V, A21, (size_t)a21_need * sizeof(double));
            for (int32_t ki = 0; ki < cand; ) {
                int iv = W->ipiv[ki];
                if (iv > 0) {
                    double d = F[ki + (size_t)ki*R];
                    if (d == 0.0) { st = LDLT_ERR_SINGULAR; goto cleanup; }
                    double inv = 1.0 / d;
                    for (int32_t i = 0; i < noff; ++i) V[i + (size_t)ki*noff] *= inv;
                    ki++;
                } else {
                    double d00 = F[ki   + (size_t)ki*R];
                    double d10 = F[ki+1 + (size_t)ki*R];
                    double d11 = F[ki+1 + (size_t)(ki+1)*R];
                    double det = d00*d11 - d10*d10;
                    if (det == 0.0) { st = LDLT_ERR_SINGULAR; goto cleanup; }
                    for (int32_t i = 0; i < noff; ++i) {
                        double c0 = V[i + (size_t)ki*noff];
                        double c1 = V[i + (size_t)(ki+1)*noff];
                        V[i + (size_t)ki*noff]     = ( c0*d11 - c1*d10) / det;
                        V[i + (size_t)(ki+1)*noff] = (-c0*d10 + c1*d00) / det;
                    }
                    ki += 2;
                }
            }

            if (noff > 0) {
                double max_mult = 1.0 / alpha;
                int stable_l21 = 1;
                for (int32_t j = 0; j < cand && stable_l21; ++j) {
                    for (int32_t i = 0; i < noff; ++i) {
                        double lij = V[i + (size_t)j*noff];
                        if (!isfinite(lij) || fabs(lij) > max_mult) {
                            stable_l21 = 0;
                            break;
                        }
                    }
                }
                if (!stable_l21) {
                    if (profile && profile->enabled) prof.fallback_l21 = 1;
                    memcpy(F, W->Fsave, (size_t)F_need * sizeof(double));
                    goto scalar_path;
                }
            }

            /* Update Schur: F[cand:R, cand:R] -= A21 * V^T  (DGEMM, both A21 and V are noff x cand). */
            double *F22 = F + cand + (size_t)cand * R;
            cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans,
                        noff, noff, cand, -1.0, A21, noff, V, noff, 1.0, F22, R);
        }

        /* Extract pivot blocks from IPIV and append to the shared build B.
           Acquire the lock once for all blocks of this supernode so that the
           expensive dsytrf/DTRSM/DGEMM computation stays fully parallel.
           For each pivot k, the L column has:
             - L11 rows: vars[k+bs .. cand-1] with values from F (L stored by dsytrf)
             - L21 rows: vars[cand .. R-1]   with values from V (= A21*D^{-1} = L21)  */
        int64_t max_lvals = (int64_t)R * 2 + 1;
        W->lvals = (double*)grow_buf(W->lvals, &W->lvals_cap, max_lvals, sizeof(double));
        W->rows  = (int32_t*)grow_buf(W->rows,  &W->rows_cap,  R > 0 ? R : 1, sizeof(int32_t));
        double  *lvals_buf = W->lvals;
        int32_t *rows_buf  = W->rows;
        double  *Vf = (noff > 0) ? W->V : NULL;

#if LDLT_HAVE_OPENMP
        {
            omp_lock_t *lock = (omp_lock_t*)append_lock_ptr;
            if (lock) omp_set_lock(lock);
#endif
        for (int32_t ki = 0; ki < cand; ) {
            int32_t bs;
            double Dblk[3];
            if (W->ipiv[ki] > 0) {
                bs = 1;
                Dblk[0] = F[ki + (size_t)ki*R];
            } else {
                bs = 2;
                Dblk[0] = F[ki   + (size_t)ki*R];
                Dblk[1] = F[ki+1 + (size_t)ki*R];
                Dblk[2] = F[ki+1 + (size_t)(ki+1)*R];
            }
            int32_t n_l11 = cand - ki - bs;
            int32_t nrows_total = n_l11 + noff;

            /* L11 rows (within the cand block, below position ki+bs). */
            for (int32_t t = 0; t < n_l11; ++t) {
                rows_buf[t] = vars[ki + bs + t];
                for (int32_t r = 0; r < bs; ++r)
                    lvals_buf[t*bs + r] = F[(ki+bs+t) + (size_t)(ki+r)*R];
            }
            /* L21 rows (below the cand block): use V = A21*D^{-1} = L21, not A21. */
            for (int32_t t = 0; t < noff; ++t) {
                rows_buf[n_l11 + t] = vars[cand + t];
                for (int32_t r = 0; r < bs; ++r)
                    lvals_buf[(n_l11+t)*bs + r] = Vf[t + (size_t)(ki+r)*noff];
            }
            append_pivot_block(B, bs, vars + ki, Dblk, nrows_total, rows_buf, lvals_buf);
            if (profile && profile->enabled) prof.pivots++;
            ki += bs;
        }
#if LDLT_HAVE_OPENMP
            if (lock) omp_unset_lock(lock);
        }
#endif
        rem_start = cand;
        noff = R - rem_start;
        if (profile && profile->enabled) prof.lapack_accept = 1;
        goto store_contribution;
    }
#endif /* !LDLT_NO_LAPACK */

#ifndef LDLT_NO_LAPACK
scalar_path:
#endif
    /* --- Scalar Bunch-Kaufman ---
       Used either when LAPACK is unavailable or when off-front rows are present,
       because the pivot stability test has to include the whole assembled front. */
    {
        if (profile && profile->enabled) {
            prof.scalar_used = 1;
            t0 = prof_now_sec();
        }
        int32_t elim = 0;
        W->lvals = (double*)grow_buf(W->lvals, &W->lvals_cap,
                                     (int64_t)(R > 0 ? R : 1) * 2, sizeof(double));
        W->l_cm = (double*)grow_buf(W->l_cm, &W->l_cm_cap,
                                    (int64_t)(R > 0 ? R : 1) * 2, sizeof(double));
        W->w_cm = (double*)grow_buf(W->w_cm, &W->w_cm_cap,
                                    (int64_t)(R > 0 ? R : 1) * 2, sizeof(double));
        double *lvals = W->lvals;
        double *l_cm = W->l_cm;
        double *w_cm = W->w_cm;
        const int32_t blas_thresh = 64;
        while (elim < cand) {
            int32_t p = -1, q = -1, bs = 0;
            double t_select0 = profile && profile->enabled ? prof_now_sec() : 0.0;
            st = select_pivot_block(F, R, elim, cand, alpha, tiny, &p, &q, &bs);
            if (profile && profile->enabled)
                prof.t_scalar_select += prof_now_sec() - t_select0;
            if (st != LDLT_OK) break;
            if (bs == 0) break;
            bring_pivot_pair_to_front(F, R, vars, elim, &p, &q, bs);
            int32_t nrows = R - elim - bs;
            double Dblk[3] = {0.0, 0.0, 0.0};
            double t_update0 = profile && profile->enabled ? prof_now_sec() : 0.0;
            if (bs == 1) {
                if (profile && profile->enabled) prof.scalar_1x1++;
                double d = F[elim + (size_t)elim*R];
                if (!isfinite(d) || fabs(d) < tiny) {
                    pivot_debug_front("singular 1x1 pivot", s, F, R, elim, cand);
                    st = LDLT_ERR_SINGULAR;
                    break;
                }
                Dblk[0] = d;
                double inv = 1.0 / d;
                const double *acol = F + (size_t)elim * R + (elim + 1);
                for (int32_t t = 0; t < nrows; ++t) lvals[t] = acol[t] * inv;
                if (nrows >= blas_thresh) {
                    for (int32_t t = 0; t < nrows; ++t) {
                        l_cm[t] = lvals[t];
                        w_cm[t] = d * lvals[t];
                    }
                    double *F22 = F + (size_t)(elim + 1) * R + (elim + 1);
                    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans,
                                nrows, nrows, 1, -1.0,
                                l_cm, nrows, w_cm, nrows, 1.0, F22, R);
                } else {
                    for (int32_t jj = 0; jj < nrows; ++jj) {
                        int32_t j = elim + 1 + jj;
                        double dlj = d * lvals[jj];
                        double *Fcol = F + (size_t)j * R;
                        for (int32_t ii = jj; ii < nrows; ++ii) {
                            int32_t i = elim + 1 + ii;
                            double upd = lvals[ii] * dlj;
                            Fcol[i] -= upd;
                            if (i != j) F[j + (size_t)i*R] -= upd;
                        }
                    }
                }
                append_pivot_block_locked(B, 1, vars + elim, Dblk, nrows,
                                          vars + elim + 1, lvals,
                                          append_lock_ptr);
                if (profile && profile->enabled) prof.pivots++;
                elim += 1;
            } else {
                if (profile && profile->enabled) prof.scalar_2x2++;
                double d00 = F[elim + (size_t)elim*R];
                double d10 = F[(elim+1) + (size_t)elim*R];
                double d11 = F[(elim+1) + (size_t)(elim+1)*R];
                double det = d00*d11 - d10*d10;
                double dscale = fabs(d00) > fabs(d11) ? fabs(d00) : fabs(d11);
                if (fabs(d10) > dscale) dscale = fabs(d10);
                if (dscale < 1.0) dscale = 1.0;
                if (!isfinite(det) || fabs(det) <= tiny * dscale) {
                    pivot_debug_front("singular 2x2 pivot", s, F, R, elim, cand);
                    st = LDLT_ERR_SINGULAR;
                    break;
                }
                Dblk[0] = d00; Dblk[1] = d10; Dblk[2] = d11;
                const double *col0 = F + (size_t)elim * R + (elim + 2);
                const double *col1 = F + (size_t)(elim+1) * R + (elim + 2);
                for (int32_t t = 0; t < nrows; ++t) {
                    double a0 = col0[t], a1 = col1[t];
                    lvals[2*t]   = (a0*d11 - a1*d10) / det;
                    lvals[2*t+1] = (-a0*d10 + a1*d00) / det;
                }
                if (nrows >= blas_thresh) {
                    for (int32_t t = 0; t < nrows; ++t) {
                        double l0 = lvals[2*t];
                        double l1 = lvals[2*t + 1];
                        l_cm[t]         = l0;
                        l_cm[t + nrows] = l1;
                        w_cm[t]         = d00 * l0 + d10 * l1;
                        w_cm[t + nrows] = d10 * l0 + d11 * l1;
                    }
                    double *F22 = F + (size_t)(elim + 2) * R + (elim + 2);
                    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans,
                                nrows, nrows, 2, -1.0,
                                l_cm, nrows, w_cm, nrows, 1.0, F22, R);
                } else {
                    for (int32_t jj = 0; jj < nrows; ++jj) {
                        int32_t j = elim + 2 + jj;
                        double lj0 = lvals[2*jj], lj1 = lvals[2*jj+1];
                        double wj0 = d00*lj0 + d10*lj1, wj1 = d10*lj0 + d11*lj1;
                        double *Fcol = F + (size_t)j * R;
                        for (int32_t ii = jj; ii < nrows; ++ii) {
                            int32_t i = elim + 2 + ii;
                            double upd = lvals[2*ii]*wj0 + lvals[2*ii+1]*wj1;
                            Fcol[i] -= upd;
                            if (i != j) F[j + (size_t)i*R] -= upd;
                        }
                    }
                }
                append_pivot_block_locked(B, 2, vars + elim, Dblk, nrows,
                                          vars + elim + 2, lvals,
                                          append_lock_ptr);
                if (profile && profile->enabled) prof.pivots++;
                elim += 2;
            }
            if (profile && profile->enabled)
                prof.t_scalar_update += prof_now_sec() - t_update0;
        }
        if (profile && profile->enabled) prof.t_scalar += prof_now_sec() - t0;
        rem_start = elim;
        noff = R - rem_start;
        if (st == LDLT_OK && S->super_parent[s] < 0 && noff != 0) {
            pivot_debug_front("root has delayed variables", s, F, R, elim, cand);
            st = LDLT_ERR_SINGULAR;
        }
    }
#ifndef LDLT_NO_LAPACK
    goto store_contribution;

cleanup:
    ;
#endif

store_contribution:
    if (profile && profile->enabled) t0 = prof_now_sec();
    if (st == LDLT_OK) {
        int32_t parent = S->super_parent[s];
        if (noff > 0 && parent >= 0) {
            pivot_contrib *C = (pivot_contrib*)ldlt_xcalloc(1, sizeof(*C));
            C->nvars = noff;
            C->vars  = (int32_t*)ldlt_xmalloc((size_t)noff * sizeof(int32_t));
            C->A     = (double*)ldlt_xmalloc(((size_t)noff * (size_t)(noff+1) / 2) * sizeof(double));
            for (int32_t i = 0; i < noff; ++i) C->vars[i] = vars[rem_start + i];
            for (int32_t j = 0; j < noff; ++j) {
                double *dst = C->A + packed_lower_col(noff, j);
                for (int32_t i = j; i < noff; ++i)
                    dst[i - j] = F[(rem_start+i) + (size_t)(rem_start+j)*R];
            }
            contrib[s] = C;
        } else if (noff > 0 && parent < 0) {
            /* Root node should have no remainder. */
            st = LDLT_ERR_SINGULAR;
        }
    }
    if (profile && profile->enabled) {
        prof.delayed_vars = noff;
        prof.t_store += prof_now_sec() - t0;
    }

    for (int32_t i = 0; i < ntouched; ++i) map[touched[i]] = -1;
    if (profile && profile->enabled) {
        prof.t_total = prof_now_sec() - t_node0;
        pivot_profile_add_node(profile, &prof);
    }
    return st;
}

static ldlt_status ldlt_factorize_pivoted(const ldlt_symbolic *S,
                                          const int32_t *Ap, const int32_t *Ai,
                                          const double *Ax,
                                          const ldlt_options *opt,
                                          ldlt_numeric **out)
{
    double alpha = opt->pivot_threshold > 0.0
        ? opt->pivot_threshold
        : (1.0 + sqrt(17.0)) / 8.0;
    if (!(alpha > 0.0 && alpha <= 1.0)) return LDLT_ERR_INPUT;

    int32_t *Bp = NULL, *Bi = NULL;
    double *Bx = NULL;
    ldlt_status pst = ldlt_permute_lower(S->n, Ap, Ai, Ax, S->perm, S->iperm,
                                         &Bp, &Bi, &Bx);
    if (pst != LDLT_OK) return pst;

    ldlt_numeric *N = (ldlt_numeric*)ldlt_xcalloc(1, sizeof(*N));
    N->S = S;
    pivot_numeric_init(N, S->n);

    pivot_build B;
    memset(&B, 0, sizeof(B));
    B.N = N;
    B.d_cap = (int64_t)(2 * (S->n > 0 ? S->n : 1));
    B.l_row_cap = (int64_t)(S->n > 0 ? S->n : 1);
    B.l_val_cap = (int64_t)(S->n > 0 ? S->n : 1);

    int32_t *level_ptr = NULL, *level_idx = NULL, nlevels = 0;
    compute_levels(S, &level_ptr, &level_idx, &nlevels);
    pivot_contrib **contrib = (pivot_contrib**)ldlt_xcalloc(
        (size_t)(S->nsuper > 0 ? S->nsuper : 1), sizeof(pivot_contrib*));

    int32_t n = S->n;

#if LDLT_HAVE_OPENMP
    int nth = opt->nthreads > 0 ? opt->nthreads : omp_get_max_threads();
    if (nth < 1) nth = 1;
#else
    int nth = 1;
#endif

    /* Per-thread maps and workspaces; pivot_build B is shared and protected by a lock. */
    int32_t **maps = (int32_t**)ldlt_xmalloc((size_t)nth * sizeof(int32_t*));
    pivot_workspace *Ws = (pivot_workspace*)ldlt_xcalloc((size_t)nth, sizeof(pivot_workspace));
    for (int t = 0; t < nth; ++t) {
        maps[t] = (int32_t*)ldlt_xmalloc((size_t)(n > 0 ? n : 1) * sizeof(int32_t));
        for (int32_t i = 0; i < n; ++i) maps[t][i] = -1;
    }

#if LDLT_HAVE_OPENMP
    omp_lock_t append_lock;
    omp_init_lock(&append_lock);
#endif
    pivot_profile profile;
    pivot_profile_init(&profile);

    ldlt_status st = LDLT_OK;
    for (int32_t L = 0; L < nlevels && st == LDLT_OK; ++L) {
        int32_t lo = level_ptr[L], hi = level_ptr[L+1];
        pivot_profile_begin_level(&profile, L);
#if LDLT_HAVE_OPENMP
        int32_t width = hi - lo;
        if (nth > 1 && width > 1) {
            #pragma omp parallel for num_threads(nth) schedule(dynamic)
            for (int32_t idx = lo; idx < hi; ++idx) {
                if (st != LDLT_OK) continue;
                int tid = omp_get_thread_num();
                int32_t s = level_idx[idx];
                /* factor_pivot_node computes locally, then calls append_pivot_block
                   under the lock (inside factor_pivot_node for now: we pass the lock). */
                ldlt_status lst = factor_pivot_node(S, &B, s, contrib,
                                                    maps[tid], &Ws[tid],
                                                    alpha, 1e-300, Bp, Bi, Bx,
                                                    &append_lock, &profile);
                if (lst != LDLT_OK) {
                    #pragma omp atomic write
                    st = lst;
                }
            }
        } else
#endif
        {
            for (int32_t idx = lo; idx < hi && st == LDLT_OK; ++idx) {
                int32_t s = level_idx[idx];
                st = factor_pivot_node(S, &B, s, contrib,
                                       maps[0], &Ws[0],
                                       alpha, 1e-300, Bp, Bi, Bx,
                                       NULL, &profile);
            }
        }
        pivot_profile_end_level(&profile);
    }

    pivot_profile_destroy(&profile);

#if LDLT_HAVE_OPENMP
    omp_destroy_lock(&append_lock);
#endif

    for (int t = 0; t < nth; ++t) { free(maps[t]); pivot_workspace_free(&Ws[t]); }
    free(maps);
    free(Ws);

    for (int32_t s = 0; s < S->nsuper; ++s) pivot_contrib_free(contrib[s]);
    free(contrib);
    free(level_ptr);
    free(level_idx);
    free(Bp);
    free(Bi);
    free(Bx);

    if (st == LDLT_OK && N->n_piv_vars != S->n) st = LDLT_ERR_SINGULAR;
    if (st == LDLT_OK && opt->require_pd && B.indefinite) st = LDLT_ERR_INDEFINITE;
    if (st != LDLT_OK) {
        ldlt_free_numeric(N);
        return st;
    }
    *out = N;
    return LDLT_OK;
}

ldlt_status ldlt_factorize(const ldlt_symbolic *S,
                           const int32_t *Ap, const int32_t *Ai, const double *Ax,
                           const ldlt_options *opt, ldlt_numeric **out)
{
    if (!S || !out) return LDLT_ERR_INPUT;
    *out = NULL;
    ldlt_options local; ldlt_options_default(&local);
    if (opt) local = *opt;

    ldlt_numeric * volatile Ntrap = NULL;
    ldlt_error_trap trap;
    ldlt_push_error_trap(&trap);
    if (setjmp(trap.env) != 0) {
        ldlt_pop_error_trap(&trap);
        ldlt_free_numeric((ldlt_numeric*)Ntrap);
        return trap.status == LDLT_OK ? LDLT_ERR_NOMEM : trap.status;
    }

    if (local.pivot == LDLT_PIVOT_BUNCH_KAUFMAN)
    {
        ldlt_status st = ldlt_factorize_pivoted(S, Ap, Ai, Ax, &local, out);
        ldlt_pop_error_trap(&trap);
        return st;
    }
    if (local.pivot != LDLT_PIVOT_NONE) {
        ldlt_pop_error_trap(&trap);
        return LDLT_ERR_UNSUPPORTED;
    }

    /* Pre-permute A once into B = P A P^T (lower CSC) so the per-supernode scatter
       is a simple column scan with no permutation indirection. */
    int32_t *Bp = NULL, *Bi = NULL; double *Bx = NULL;
    ldlt_status pst = ldlt_permute_lower(S->n, Ap, Ai, Ax, S->perm, S->iperm,
                                         &Bp, &Bi, &Bx);
    if (pst != LDLT_OK) {
        ldlt_pop_error_trap(&trap);
        return pst;
    }

    ldlt_numeric *N = (ldlt_numeric*)ldlt_xcalloc(1, sizeof(*N));
    Ntrap = N;
    N->S = S;
    int64_t Lsz = S->panel_off[S->nsuper];
    if (Lsz < 1) Lsz = 1;
    if ((uint64_t)Lsz > (uint64_t)SIZE_MAX / sizeof(double)) {
        free(Bp); free(Bi); free(Bx);
        ldlt_free_numeric(N);
        ldlt_pop_error_trap(&trap);
        return LDLT_ERR_NOMEM;
    }
    uint64_t avail = ldlt_available_memory_bytes();
    uint64_t lbytes = (uint64_t)Lsz * (uint64_t)sizeof(double);
    if (avail > 0 && (long double)lbytes > (long double)avail * 0.90L) {
        fprintf(stderr,
                "ldlt: numeric factor storage requires %.2f GB, available %.2f GB"
                " -- attempting anyway.\n",
                (double)lbytes / 1073741824.0,
                (double)avail / 1073741824.0);
    }
    N->L = (double*)ldlt_xcalloc((size_t)Lsz, sizeof(double));
    N->D = (double*)ldlt_xcalloc((size_t)(S->n>0?S->n:1), sizeof(double));

    int32_t max_md = 0, max_kd = 0;
    for (int32_t d = 0; d < S->nsuper; ++d) {
        if (S->super[d].nrows_below > max_md) max_md = S->super[d].nrows_below;
        if (S->super[d].width > max_kd) max_kd = S->super[d].width;
    }
    if (max_md < 1) max_md = 1;
    if (max_kd < 1) max_kd = 1;
    if ((uint64_t)max_md > (uint64_t)SIZE_MAX / (uint64_t)max_kd) {
        free(Bp); free(Bi); free(Bx);
        ldlt_free_numeric(N);
        ldlt_pop_error_trap(&trap);
        return LDLT_ERR_NOMEM;
    }
    size_t work_doubles = (size_t)max_md * max_kd;
    size_t blocked_work = (size_t)(max_md + max_kd) * 32;
    if (blocked_work > work_doubles) work_doubles = blocked_work;
    if (work_doubles < 1) work_doubles = 1;

    int32_t *level_ptr = NULL, *level_idx = NULL, nlevels = 0;
    compute_levels(S, &level_ptr, &level_idx, &nlevels);

    int32_t n = S->n;
    ldlt_status global_status = LDLT_OK;

#if LDLT_HAVE_OPENMP
    int nth = local.nthreads > 0 ? local.nthreads : omp_get_max_threads();
    if (nth < 1) nth = 1;
#else
    int nth = 1;
#endif
    int32_t **rowmaps = (int32_t**)ldlt_xmalloc((size_t)nth * sizeof(int32_t*));
    double **works = (double**)ldlt_xmalloc((size_t)nth * sizeof(double*));
    front_workspace *fronts = (front_workspace*)ldlt_xcalloc((size_t)nth, sizeof(front_workspace));
    double **contrib = (double**)ldlt_xcalloc((size_t)(S->nsuper>0?S->nsuper:1),
                                             sizeof(double*));
    for (int t = 0; t < nth; ++t) {
        rowmaps[t]  = (int32_t*)ldlt_xmalloc((size_t)(n>0?n:1) * sizeof(int32_t));
        for (int32_t i = 0; i < n; ++i) rowmaps[t][i] = -1;
        works[t] = (double*)ldlt_xmalloc(work_doubles * sizeof(double));
    }

    for (int32_t L = 0; L < nlevels && global_status == LDLT_OK; ++L) {
        int32_t lo = level_ptr[L], hi = level_ptr[L+1];
#if LDLT_HAVE_OPENMP
        if (nth > 1 && hi - lo > 1) {
            #pragma omp parallel for num_threads(nth) schedule(dynamic)
            for (int32_t idx = lo; idx < hi; ++idx) {
                int32_t s = level_idx[idx];
#if LDLT_HAVE_OPENMP
                int tid = omp_get_thread_num();
#else
                int tid = 0;
#endif
                ldlt_status st = factor_front_node(S, N, s, contrib, rowmaps[tid],
                                                   works[tid], &fronts[tid],
                                                   &local, Bp, Bi, Bx);
                if (st != LDLT_OK) {
#if LDLT_HAVE_OPENMP
                    #pragma omp atomic write
#endif
                    global_status = st;
                }
            }
        } else
#endif
        {
            for (int32_t idx = lo; idx < hi; ++idx) {
                int32_t s = level_idx[idx];
                ldlt_status st = factor_front_node(S, N, s, contrib, rowmaps[0],
                                                   works[0], &fronts[0],
                                                   &local, Bp, Bi, Bx);
                if (st != LDLT_OK) {
                    global_status = st;
                    break;
                }
            }
        }
    }

    for (int t = 0; t < nth; ++t) { free(rowmaps[t]); free(works[t]); free(fronts[t].front); }
    for (int32_t s = 0; s < S->nsuper; ++s) free(contrib[s]);
    free(contrib);
    free(fronts);
    free(rowmaps); free(works);
    free(level_ptr); free(level_idx);
    free(Bp); free(Bi); free(Bx);

    if (global_status != LDLT_OK) {
        ldlt_free_numeric(N);
        ldlt_pop_error_trap(&trap);
        return global_status;
    }
    *out = N;
    Ntrap = NULL;
    ldlt_pop_error_trap(&trap);
    return LDLT_OK;
}
