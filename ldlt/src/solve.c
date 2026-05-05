/* Triangular solves with the supernodal factor: L y = Pb, D z = y, L^T x = z, b = P^T x.
   Serial in v1; uses BLAS only for the dense per-supernode operations. */

#include "internal.h"
#include "blas_shim.h"

#ifndef LDLT_SMALL_SOLVE_MAX
#define LDLT_SMALL_SOLVE_MAX 2
#endif

static void small_forward_solve(int32_t k, int32_t nrhs, const double *L, int32_t ldL,
                                double *X, int32_t ldX)
{
    for (int32_t r = 0; r < nrhs; ++r) {
        double *x = X + (size_t)r*ldX;
        for (int32_t j = 0; j < k; ++j) {
            double xj = x[j];
            for (int32_t i = j + 1; i < k; ++i)
                x[i] -= L[i + (size_t)j*ldL] * xj;
        }
    }
}

static void small_backward_solve(int32_t k, int32_t nrhs, const double *L, int32_t ldL,
                                 double *X, int32_t ldX)
{
    for (int32_t r = 0; r < nrhs; ++r) {
        double *x = X + (size_t)r*ldX;
        for (int32_t j = k - 1; j >= 0; --j) {
            double xj = x[j];
            for (int32_t i = 0; i < j; ++i)
                x[i] -= L[j + (size_t)i*ldL] * xj;
        }
    }
}

static void small_apply_l21(int32_t k, int32_t m, int32_t nrhs, const double *L21,
                            int32_t ldL, const int32_t *rb,
                            const double *Xtop, double *X, int32_t ldX)
{
    for (int32_t r = 0; r < nrhs; ++r) {
        const double *xt = Xtop + (size_t)r*ldX;
        for (int32_t t = 0; t < m; ++t) {
            double acc = 0.0;
            for (int32_t j = 0; j < k; ++j)
                acc += L21[t + (size_t)j*ldL] * xt[j];
            X[rb[t] + (size_t)r*ldX] -= acc;
        }
    }
}

static void small_apply_l21_trans(int32_t k, int32_t m, int32_t nrhs, const double *L21,
                                  int32_t ldL, const int32_t *rb,
                                  double *Xtop, const double *X, int32_t ldX)
{
    for (int32_t r = 0; r < nrhs; ++r) {
        double *xt = Xtop + (size_t)r*ldX;
        for (int32_t j = 0; j < k; ++j) {
            double acc = 0.0;
            for (int32_t t = 0; t < m; ++t)
                acc += L21[t + (size_t)j*ldL] * X[rb[t] + (size_t)r*ldX];
            xt[j] -= acc;
        }
    }
}

static ldlt_status solve_pivoted(const ldlt_numeric *N, int32_t nrhs,
                                 double *B, int32_t ldb)
{
    const ldlt_symbolic *S = N->S;
    int32_t n = S->n;
    if (ldb < n) return LDLT_ERR_INPUT;
    if (nrhs == 0 || n == 0) return LDLT_OK;
    if (N->n_piv_vars != n) return LDLT_ERR_SINGULAR;

    int identity = S->perm_identity;
    double *Xbuf = NULL;
    double *X = B;
    int32_t ldX = ldb;
    if (!identity || ldb != n) {
        Xbuf = (double*)ldlt_xmalloc((size_t)n * (size_t)nrhs * sizeof(double));
        for (int32_t r = 0; r < nrhs; ++r)
            for (int32_t i = 0; i < n; ++i)
                Xbuf[i + (size_t)r*n] = B[S->perm[i] + (size_t)r*ldb];
        X = Xbuf;
        ldX = n;
    }

    for (int32_t b = 0; b < N->npivots; ++b) {
        int32_t bs = N->piv_size[b];
        int32_t vo = N->piv_var_off[b];
        int64_t ro = N->piv_l_row_off[b];
        int64_t rhi = N->piv_l_row_off[b + 1];
        int64_t vlo = N->piv_l_val_off[b];
        const int32_t *pvars = N->piv_vars + vo;
        const int32_t *rows = N->piv_l_rows + ro;
        const double *vals = N->piv_l_vals + vlo;
        int32_t nrows = (int32_t)(rhi - ro);
        for (int32_t r = 0; r < nrhs; ++r) {
            double *x = X + (size_t)r*ldX;
            if (bs == 1) {
                double x0 = x[pvars[0]];
                for (int32_t t = 0; t < nrows; ++t)
                    x[rows[t]] -= vals[t] * x0;
            } else {
                double x0 = x[pvars[0]];
                double x1 = x[pvars[1]];
                for (int32_t t = 0; t < nrows; ++t)
                    x[rows[t]] -= vals[2*t] * x0 + vals[2*t + 1] * x1;
            }
        }
    }

    for (int32_t b = 0; b < N->npivots; ++b) {
        int32_t bs = N->piv_size[b];
        int32_t vo = N->piv_var_off[b];
        int64_t doff = N->piv_d_off[b];
        const int32_t *pvars = N->piv_vars + vo;
        const double *D = N->piv_D + doff;
        for (int32_t r = 0; r < nrhs; ++r) {
            double *x = X + (size_t)r*ldX;
            if (bs == 1) {
                if (D[0] == 0.0) {
                    free(Xbuf);
                    return LDLT_ERR_SINGULAR;
                }
                x[pvars[0]] /= D[0];
            } else {
                double det = D[0] * D[2] - D[1] * D[1];
                if (det == 0.0) {
                    free(Xbuf);
                    return LDLT_ERR_SINGULAR;
                }
                double y0 = x[pvars[0]];
                double y1 = x[pvars[1]];
                x[pvars[0]] = (D[2] * y0 - D[1] * y1) / det;
                x[pvars[1]] = (-D[1] * y0 + D[0] * y1) / det;
            }
        }
    }

    for (int32_t b = N->npivots - 1; b >= 0; --b) {
        int32_t bs = N->piv_size[b];
        int32_t vo = N->piv_var_off[b];
        int64_t ro = N->piv_l_row_off[b];
        int64_t rhi = N->piv_l_row_off[b + 1];
        int64_t vlo = N->piv_l_val_off[b];
        const int32_t *pvars = N->piv_vars + vo;
        const int32_t *rows = N->piv_l_rows + ro;
        const double *vals = N->piv_l_vals + vlo;
        int32_t nrows = (int32_t)(rhi - ro);
        for (int32_t r = 0; r < nrhs; ++r) {
            double *x = X + (size_t)r*ldX;
            if (bs == 1) {
                double acc = 0.0;
                for (int32_t t = 0; t < nrows; ++t)
                    acc += vals[t] * x[rows[t]];
                x[pvars[0]] -= acc;
            } else {
                double acc0 = 0.0, acc1 = 0.0;
                for (int32_t t = 0; t < nrows; ++t) {
                    double xr = x[rows[t]];
                    acc0 += vals[2*t] * xr;
                    acc1 += vals[2*t + 1] * xr;
                }
                x[pvars[0]] -= acc0;
                x[pvars[1]] -= acc1;
            }
        }
    }

    if (Xbuf) {
        for (int32_t r = 0; r < nrhs; ++r)
            for (int32_t i = 0; i < n; ++i)
                B[S->perm[i] + (size_t)r*ldb] = Xbuf[i + (size_t)r*n];
        free(Xbuf);
    }
    return LDLT_OK;
}

ldlt_status ldlt_solve(const ldlt_numeric *N, int32_t nrhs, double *B, int32_t ldb)
{
    if (!N || !B || nrhs < 0) return LDLT_ERR_INPUT;
    const ldlt_symbolic *S = N->S;
    int32_t n = S->n;
    if (ldb < n) return LDLT_ERR_INPUT;
    if (nrhs == 0 || n == 0) return LDLT_OK;

    ldlt_error_trap trap;
    ldlt_push_error_trap(&trap);
    if (setjmp(trap.env) != 0) {
        ldlt_pop_error_trap(&trap);
        return trap.status == LDLT_OK ? LDLT_ERR_NOMEM : trap.status;
    }

    if (N->pivoted) {
        ldlt_status st = solve_pivoted(N, nrhs, B, ldb);
        ldlt_pop_error_trap(&trap);
        return st;
    }

    int identity = S->perm_identity;

    double *Xbuf = NULL;
    double *X = B;
    int32_t ldX = ldb;
    if (!identity || ldb != n) {
        Xbuf = (double*)ldlt_xmalloc((size_t)n * (size_t)nrhs * sizeof(double));
        for (int32_t r = 0; r < nrhs; ++r)
            for (int32_t i = 0; i < n; ++i)
                Xbuf[i + (size_t)r*n] = B[S->perm[i] + (size_t)r*ldb];
        X = Xbuf;
        ldX = n;
    }

    int32_t max_m = 0;
    for (int32_t s = 0; s < S->nsuper; ++s)
        if (S->super[s].nrows_below > max_m) max_m = S->super[s].nrows_below;
    double *work = NULL;
    if (max_m > 0)
        work = (double*)ldlt_xmalloc((size_t)max_m * (size_t)nrhs * sizeof(double));

    /* Forward: solve L y = X for each rhs column */
    for (int32_t s = 0; s < S->nsuper; ++s) {
        const ldlt_super *Ss = &S->super[s];
        int32_t fc = Ss->first_col, k = Ss->width, m = Ss->nrows_below;
        int32_t R = k + m;
        const double *L = N->L + S->panel_off[s];
        if (k <= LDLT_SMALL_SOLVE_MAX) {
            small_forward_solve(k, nrhs, L, R, X + fc, ldX);
        } else {
            cblas_dtrsm(CblasColMajor, CblasLeft, CblasLower, CblasNoTrans, CblasUnit,
                        k, nrhs, 1.0, L, R, X + fc, ldX);
        }
        if (m > 0) {
            const int32_t *rb = S->rows_pool + Ss->row_off;
            if (k == 1) {
                const double *l21 = L + 1;
                for (int32_t r = 0; r < nrhs; ++r) {
                    double x0 = X[fc + (size_t)r*ldX];
                    for (int32_t t = 0; t < m; ++t)
                        X[rb[t] + (size_t)r*ldX] -= l21[t] * x0;
                }
            } else if (k <= LDLT_SMALL_SOLVE_MAX) {
                small_apply_l21(k, m, nrhs, L + k, R, rb, X + fc, X, ldX);
            } else {
                cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans,
                            m, nrhs, k, 1.0, L + k, R, X + fc, ldX, 0.0, work, m);
                for (int32_t r = 0; r < nrhs; ++r)
                    for (int32_t t = 0; t < m; ++t)
                        X[rb[t] + (size_t)r*ldX] -= work[t + (size_t)r*m];
            }
        }
    }

    /* Backward: first apply D^{-1}, then solve L^T x = z in reverse order. */
    for (int32_t s = S->nsuper - 1; s >= 0; --s) {
        const ldlt_super *Ss = &S->super[s];
        int32_t fc = Ss->first_col, k = Ss->width, m = Ss->nrows_below;
        int32_t R = k + m;
        const double *L = N->L + S->panel_off[s];
        const double *D = N->D + S->d_off[s];
        for (int32_t j = 0; j < k; ++j) {
            double inv = 1.0 / D[j];
            for (int32_t r = 0; r < nrhs; ++r)
                X[(fc+j) + (size_t)r*ldX] *= inv;
        }
        if (m > 0) {
            const int32_t *rb = S->rows_pool + Ss->row_off;
            if (k == 1) {
                const double *l21 = L + 1;
                for (int32_t r = 0; r < nrhs; ++r) {
                    double acc = 0.0;
                    for (int32_t t = 0; t < m; ++t)
                        acc += l21[t] * X[rb[t] + (size_t)r*ldX];
                    X[fc + (size_t)r*ldX] -= acc;
                }
            } else if (k <= LDLT_SMALL_SOLVE_MAX) {
                small_apply_l21_trans(k, m, nrhs, L + k, R, rb, X + fc, X, ldX);
            } else {
                for (int32_t r = 0; r < nrhs; ++r)
                    for (int32_t t = 0; t < m; ++t)
                        work[t + (size_t)r*m] = X[rb[t] + (size_t)r*ldX];
                cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans,
                            k, nrhs, m, -1.0, L + k, R, work, m, 1.0, X + fc, ldX);
            }
        }
        if (k <= LDLT_SMALL_SOLVE_MAX) {
            small_backward_solve(k, nrhs, L, R, X + fc, ldX);
        } else {
            cblas_dtrsm(CblasColMajor, CblasLeft, CblasLower, CblasTrans, CblasUnit,
                        k, nrhs, 1.0, L, R, X + fc, ldX);
        }
    }

    if (Xbuf) {
        for (int32_t r = 0; r < nrhs; ++r)
            for (int32_t i = 0; i < n; ++i)
                B[S->perm[i] + (size_t)r*ldb] = Xbuf[i + (size_t)r*n];
        free(Xbuf);
    }
    free(work);
    ldlt_pop_error_trap(&trap);
    return LDLT_OK;
}
