#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "cholmod.h"

typedef long long int fptr;

typedef struct {
    cholmod_common common;
    cholmod_factor *factor;
} factors_t;

void c_fortran_dcholmod_factor_(int *n, int *nnz, double *values, int *rowind,
                                int *colptr, fptr *f_factors, int *info)
{
    factors_t *factors = NULL;
    cholmod_sparse *A = NULL;
    int *ai = NULL;
    int *ap = NULL;
    double *ax = NULL;
    int i;

    *info = 0;
    *f_factors = 0;

    factors = (factors_t *) malloc(sizeof(factors_t));
    if (factors == NULL) {
        *info = -1;
        return;
    }

    memset(factors, 0, sizeof(*factors));
    cholmod_start(&factors->common);
    factors->common.print = 0;

    A = cholmod_allocate_sparse((size_t) *n, (size_t) *n, (size_t) *nnz,
                                1, 1, 1, CHOLMOD_REAL, &factors->common);
    if (A == NULL) {
        *info = -2;
        cholmod_finish(&factors->common);
        free(factors);
        return;
    }

    ai = (int *) A->i;
    ap = (int *) A->p;
    ax = (double *) A->x;

    for (i = 0; i < *nnz; ++i) {
        ai[i] = rowind[i] - 1;
        ax[i] = values[i];
    }
    for (i = 0; i <= *n; ++i) {
        ap[i] = colptr[i] - 1;
    }

    factors->factor = cholmod_analyze(A, &factors->common);
    if (factors->factor == NULL) {
        *info = (factors->common.status != CHOLMOD_OK) ? -factors->common.status : -3;
        cholmod_free_sparse(&A, &factors->common);
        cholmod_finish(&factors->common);
        free(factors);
        return;
    }

    if (!cholmod_factorize(A, factors->factor, &factors->common)) {
        if (factors->common.status == CHOLMOD_NOT_POSDEF) {
            *info = (int) factors->factor->minor + 1;
        } else if (factors->common.status != CHOLMOD_OK) {
            *info = -factors->common.status;
        } else {
            *info = -4;
        }
        cholmod_free_factor(&factors->factor, &factors->common);
        cholmod_free_sparse(&A, &factors->common);
        cholmod_finish(&factors->common);
        free(factors);
        return;
    }

    cholmod_free_sparse(&A, &factors->common);
    *f_factors = (fptr) factors;
}

void c_fortran_dcholmod_solve_(int *n, double *rhs, fptr *f_factors, int *info)
{
    factors_t *factors;
    cholmod_dense *B = NULL;
    cholmod_dense *X = NULL;
    double *bx;
    double *xx;
    int i;

    *info = 0;
    if (*f_factors == 0) {
        *info = -1;
        return;
    }

    factors = (factors_t *) *f_factors;

    B = cholmod_allocate_dense((size_t) *n, 1, (size_t) *n, CHOLMOD_REAL, &factors->common);
    if (B == NULL) {
        *info = -2;
        return;
    }

    bx = (double *) B->x;
    for (i = 0; i < *n; ++i) {
        bx[i] = rhs[i];
    }

    X = cholmod_solve(CHOLMOD_A, factors->factor, B, &factors->common);
    if (X == NULL) {
        *info = (factors->common.status != CHOLMOD_OK) ? -factors->common.status : -3;
        cholmod_free_dense(&B, &factors->common);
        return;
    }

    xx = (double *) X->x;
    for (i = 0; i < *n; ++i) {
        rhs[i] = xx[i];
    }

    cholmod_free_dense(&X, &factors->common);
    cholmod_free_dense(&B, &factors->common);
}

void c_fortran_dcholmod_free_(fptr *f_factors, int *info)
{
    factors_t *factors;

    *info = 0;
    if (*f_factors == 0) {
        return;
    }

    factors = (factors_t *) *f_factors;
    cholmod_free_factor(&factors->factor, &factors->common);
    cholmod_finish(&factors->common);
    free(factors);
    *f_factors = 0;
}
