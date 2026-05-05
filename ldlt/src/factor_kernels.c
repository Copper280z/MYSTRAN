/* Dense diagonal-block LDL^T (no pivot) and helpers.
   Designed so the call site can be swapped with a Bunch-Kaufman variant later
   without changing the rest of factor.c. */

#include "internal.h"
#include "blas_shim.h"
#include <math.h>

/* In-place factor of column-major (R x k) panel F into:
   - strict-lower of F[:k,:k] = L11 (unit diagonal, implicit)
   - F[k:,:k] = L21
   - D[0..k) = diagonal entries
   - F[j,j] is set to 0 to mark the implicit unit diagonal.
   Returns LDLT_ERR_INDEFINITE if any |D[j]| < tiny. */
ldlt_status ldlt_dense_ldl_nopiv(double *F, int32_t ldF, int32_t R, int32_t k,
                                 double *D, double tiny, int *indef)
{
    if (indef) *indef = 0;
    for (int32_t j = 0; j < k; ++j) {
        double dj = F[j + (size_t)j*ldF];
        if (!isfinite(dj) || fabs(dj) < tiny) return LDLT_ERR_INDEFINITE;
        if (dj < 0.0 && indef) *indef = 1;
        D[j] = dj;
        for (int32_t i = j+1; i < R; ++i) F[i + (size_t)j*ldF] /= dj;
        for (int32_t c = j+1; c < k; ++c) {
            double ljc = F[c + (size_t)j*ldF];
            double scale = dj * ljc;
            for (int32_t i = c; i < R; ++i)
                F[i + (size_t)c*ldF] -= F[i + (size_t)j*ldF] * scale;
        }
        F[j + (size_t)j*ldF] = 0.0; /* implicit unit diagonal */
    }
    return LDLT_OK;
}

ldlt_status ldlt_dense_front_nopiv(double *F, int32_t ldF, int32_t R, int32_t k,
                                   double *D, double *work, double tiny, int *indef)
{
    if (indef) *indef = 0;
    if (k <= 0) return LDLT_OK;

    const int32_t nb = 32;
    if (k < nb) {
        ldlt_status st = ldlt_dense_ldl_nopiv(F, ldF, R, k, D, tiny, indef);
        if (st != LDLT_OK) return st;

        int32_t m = R - k;
        if (m <= 0) return LDLT_OK;

        const double *L21 = F + k;
        double *W = work;
        int positive = 1;
        for (int32_t c = 0; c < k; ++c) if (D[c] <= 0.0) positive = 0;
        for (int32_t c = 0; c < k; ++c) {
            double dc = positive ? sqrt(D[c]) : D[c];
            for (int32_t i = 0; i < m; ++i)
                W[i + (size_t)c*m] = L21[i + (size_t)c*ldF] * dc;
        }

        double *F22 = F + k + (size_t)k*ldF;
        if (positive) {
            cblas_dsyrk(CblasColMajor, CblasLower, CblasNoTrans,
                        m, k, -1.0, W, m, 1.0, F22, ldF);
        } else {
            cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans,
                        m, m, k, -1.0, W, m, L21, ldF, 1.0, F22, ldF);
        }
        return LDLT_OK;
    }

    for (int32_t b = 0; b < k; b += nb) {
        int32_t jb = k - b;
        if (jb > nb) jb = nb;
        int32_t Rb = R - b;
        int indef_block = 0;
        ldlt_status st = ldlt_dense_ldl_nopiv(F + b + (size_t)b*ldF, ldF, Rb, jb,
                                              D + b, tiny, &indef_block);
        if (st != LDLT_OK) return st;
        if (indef_block && indef) *indef = 1;

        int32_t m = Rb - jb;
        if (m <= 0) continue;
        const double *L21 = F + (b + jb) + (size_t)b*ldF;
        double *W = work;
        int positive = 1;
        for (int32_t c = 0; c < jb; ++c) if (D[b + c] <= 0.0) positive = 0;
        for (int32_t c = 0; c < jb; ++c) {
            double dc = positive ? sqrt(D[b + c]) : D[b + c];
            for (int32_t i = 0; i < m; ++i)
                W[i + (size_t)c*m] = L21[i + (size_t)c*ldF] * dc;
        }
        double *F22 = F + (b + jb) + (size_t)(b + jb)*ldF;
        if (positive) {
            cblas_dsyrk(CblasColMajor, CblasLower, CblasNoTrans,
                        m, jb, -1.0, W, m, 1.0, F22, ldF);
        } else {
            cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans,
                        m, m, jb, -1.0, W, m, L21, ldF, 1.0, F22, ldF);
        }
    }
    return LDLT_OK;
}
