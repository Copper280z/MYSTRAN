/*
 * Fortran-callable interface to the experimental ldlt solver.
 *
 * Drop-in replacement for c_fortran_qdldl.c / c_fortran_cholmod.c: it
 * implements c_fortran_qdldl_ and c_qdldl_get_neg_count_ so the Fortran
 * solver path does not need to change.
 *
 * Input matrix convention:
 *   The Fortran caller passes a symmetric matrix in CRS format:
 *     colptr = I_MATIN (Fortran 1-based row pointers)
 *     rowind = J_MATIN (Fortran 1-based column indices)
 *   This routine converts to 0-based lower-triangular CSC, which is the
 *   storage format expected by ldlt_analyze/ldlt_factorize.
 *
 * iopt:
 *   1 = SPD LDL^T factorization
 *   2 = solve (reuses stored factors)
 *   3 = free all storage
 */

#include <limits.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "internal.h"

/* OpenBLAS exposes openblas_set_num_threads() when linked directly.
 * Declare it weakly so the code works with other BLAS implementations too. */
#if LDLT_USE_SYSTEM_BLAS
extern void openblas_set_num_threads(int) __attribute__((weak));
extern int openblas_get_num_threads(void) __attribute__((weak));
#endif

static void blas_set_single_threaded(int *saved) {
#if LDLT_USE_SYSTEM_BLAS
  if (openblas_get_num_threads) {
    *saved = openblas_get_num_threads();
    openblas_set_num_threads(1);
  } else {
    *saved = -1;
  }
#else
  *saved = -1;
#endif
}

static void blas_restore_threads(int saved) {
#if LDLT_USE_SYSTEM_BLAS
  if (saved > 0 && openblas_set_num_threads)
    openblas_set_num_threads(saved);
#else
  (void)saved;
#endif
}

/* Pointer-sized integer to pass handles back to Fortran as INTEGER(DBL_LONG).
 */
typedef intptr_t fptr;

static int s_last_neg_count = 0;

typedef struct {
  ldlt_symbolic *S;
  ldlt_numeric *N;
} ldlt_factors_t;

static double now_sec(void) { return ldlt_wall_time_seconds(); }

static double step_done(const char *label, double t0) {
  double t1 = now_sec();
  printf("LDLT %-39s %8.3f s\n", label, t1 - t0);
  fflush(stdout);
  return t1;
}

static void *malloc_array(size_t n, size_t elem_size) {
  if (n == 0 || elem_size == 0) {
    n = 1;
    elem_size = 1;
  }
  if (n > SIZE_MAX / elem_size) {
    return NULL;
  }
  return malloc(n * elem_size);
}

/* Build lower-triangular CSC (row >= col) from symmetric 0-based CRS. */
static int build_lower_csc(int n, const int *row_ptr, const int *col_idx,
                           const double *values, int32_t *Ap, int32_t *Ai,
                           double *Ax) {
  memset(Ap, 0, (size_t)(n + 1) * sizeof(int32_t));
  for (int i = 0; i < n; ++i) {
    for (int k = row_ptr[i]; k < row_ptr[i + 1]; ++k) {
      int j = col_idx[k];
      if (i >= j) {
        Ap[j + 1]++;
      }
    }
  }
  for (int j = 0; j < n; ++j) {
    Ap[j + 1] += Ap[j];
  }

  int *cursor = (int *)malloc_array((size_t)n, sizeof(int));
  if (!cursor) {
    return -1;
  }
  for (int j = 0; j < n; ++j) {
    cursor[j] = Ap[j];
  }

  for (int i = 0; i < n; ++i) {
    for (int k = row_ptr[i]; k < row_ptr[i + 1]; ++k) {
      int j = col_idx[k];
      if (i >= j) {
        int pos = cursor[j]++;
        Ai[pos] = (int32_t)i;
        Ax[pos] = values[k];
      }
    }
  }

  free(cursor);
  return Ap[n];
}

static int validate_crs(int n, int nnz, const int *row_ptr,
                        const int *col_idx) {
  if (n < 0 || nnz < 0 || row_ptr[0] != 0 || row_ptr[n] != nnz) {
    return 0;
  }
  for (int i = 0; i < n; ++i) {
    if (row_ptr[i] > row_ptr[i + 1]) {
      return 0;
    }
  }
  for (int k = 0; k < nnz; ++k) {
    if (col_idx[k] < 0 || col_idx[k] >= n) {
      return 0;
    }
  }
  return 1;
}

void c_fortran_qdldl_(int *iopt, int *n, int *nnz, int *nrhs, double *values,
                      int *rowind, int *colptr, double *b, int *ldb,
                      fptr *f_factors, int *info) {
  *info = 0;

  if (*iopt == 1) {
#if !LDLT_USE_SYSTEM_BLAS
    printf("LDLT WARNING: This LDLT was built with the slow reference BLAS\n");
#endif
    int N = *n;
    int NNZ = *nnz;
    double t0 = now_sec();
    double t_step;

    printf("LDLT SPD factorization begin  (n=%d, nnz=%d)\n", N, NNZ);
#if !LDLT_USE_SYSTEM_BLAS
    printf("LDLT WARNING: built without a system CBLAS; using reference (slow) BLAS.\n"
           "LDLT          Install libopenblas-dev (or equivalent) and reconfigure.\n");
#endif
    fflush(stdout);

    if (N < 0 || NNZ < 0) {
      fprintf(stderr, "LDLT: invalid matrix dimensions (n=%d, nnz=%d).\n", N,
              NNZ);
      *info = -1;
      return;
    }

    int *row_ptr = (int *)malloc_array((size_t)N + 1, sizeof(int));
    int *col_idx = (int *)malloc_array((size_t)NNZ, sizeof(int));
    if (!row_ptr || !col_idx) {
      *info = -1;
      free(row_ptr);
      free(col_idx);
      return;
    }

    for (int i = 0; i <= N; ++i) {
      row_ptr[i] = colptr[i] - 1;
    }
    for (int k = 0; k < NNZ; ++k) {
      col_idx[k] = rowind[k] - 1;
    }

    if (!validate_crs(N, NNZ, row_ptr, col_idx)) {
      fprintf(stderr, "LDLT: invalid CRS input (n=%d, nnz=%d).\n", N, NNZ);
      *info = -1;
      free(row_ptr);
      free(col_idx);
      return;
    }

    int64_t nnz_lt64 = 0;
    for (int i = 0; i < N; ++i) {
      for (int k = row_ptr[i]; k < row_ptr[i + 1]; ++k) {
        if (i >= col_idx[k]) {
          nnz_lt64++;
          if (nnz_lt64 > INT32_MAX) {
            fprintf(stderr,
                    "LDLT: lower-triangular matrix exceeds 32-bit indexing.\n");
            *info = -1;
            free(row_ptr);
            free(col_idx);
            return;
          }
        }
      }
    }
    int nnz_lt = (int)nnz_lt64;

    int32_t *Ap = (int32_t *)malloc_array((size_t)N + 1, sizeof(int32_t));
    int32_t *Ai = (int32_t *)malloc_array((size_t)nnz_lt, sizeof(int32_t));
    double *Ax = (double *)malloc_array((size_t)nnz_lt, sizeof(double));
    if (!Ap || !Ai || !Ax) {
      *info = -1;
      free(row_ptr);
      free(col_idx);
      free(Ap);
      free(Ai);
      free(Ax);
      return;
    }

    nnz_lt = build_lower_csc(N, row_ptr, col_idx, values, Ap, Ai, Ax);
    free(row_ptr);
    free(col_idx);
    if (nnz_lt < 0) {
      *info = -1;
      free(Ap);
      free(Ai);
      free(Ax);
      return;
    }

    t_step = step_done("CRS extract lower-tri CSC", t0);

    ldlt_options opt;
    ldlt_options_default(&opt);
#if LDLT_HAVE_METIS
    opt.ordering = LDLT_ORDER_METIS;
    printf("LDLT ordering: METIS nested dissection\n");
#elif LDLT_HAVE_AMD
    opt.ordering = LDLT_ORDER_AMD;
    printf("LDLT ordering: AMD\n");
#else
    opt.ordering = LDLT_ORDER_NATURAL;
    printf("LDLT ordering: natural (no reordering)\n");
#endif
    opt.require_pd = 1;
    opt.pivot = LDLT_PIVOT_NONE;

    ldlt_symbolic *S = NULL;
    ldlt_status st = ldlt_analyze((int32_t)N, Ap, Ai, &opt, &S);
    if (st != LDLT_OK) {
      fprintf(stderr, "LDLT: analyze failed: %s\n", ldlt_status_str(st));
      *info = -1;
      free(Ap);
      free(Ai);
      free(Ax);
      return;
    }
    t_step = step_done("analyze (SPD ordering/symbolic)", t_step);

    int saved_blas_threads;
    blas_set_single_threaded(&saved_blas_threads);
    ldlt_numeric *Num = NULL;
    st = ldlt_factorize(S, Ap, Ai, Ax, &opt, &Num);
    blas_restore_threads(saved_blas_threads);
    free(Ap);
    free(Ai);
    free(Ax);
    if (st != LDLT_OK) {
      fprintf(stderr, "LDLT: factorize failed: %s\n", ldlt_status_str(st));
      *info = (st == LDLT_ERR_INDEFINITE || st == LDLT_ERR_SINGULAR) ? N : -1;
      ldlt_free_symbolic(S);
      return;
    }
    t_step = step_done("factorize (SPD LDL^T)", t_step);

    s_last_neg_count = 0;
    printf("LDLT SPD D diagonal -- positive=%d  negative=0  (n=%d)\n", N, N);
    printf("LDLT stats -- supernodes=%d  lnz=%lld  levels=%d  "
           "max_level_width=%d\n",
           (int)ldlt_nsuper(S), (long long)ldlt_lnz(S), (int)ldlt_nlevels(S),
           (int)ldlt_max_level_width(S));
    printf("LDLT factorization total elapsed  %8.3f s\n", now_sec() - t0);
    fflush(stdout);

    ldlt_factors_t *facs = (ldlt_factors_t *)malloc(sizeof(ldlt_factors_t));
    if (!facs) {
      *info = -1;
      ldlt_free_numeric(Num);
      ldlt_free_symbolic(S);
      return;
    }
    facs->S = S;
    facs->N = Num;
    *f_factors = (fptr)facs;

  } else if (*iopt == 2) {
    ldlt_factors_t *facs = (ldlt_factors_t *)*f_factors;
    int n_rhs = (*nrhs > 0) ? *nrhs : 1;
    int leading_dim = (*ldb > 0) ? *ldb : ldlt_n(facs->S);
    ldlt_status st =
        ldlt_solve(facs->N, (int32_t)n_rhs, b, (int32_t)leading_dim);
    if (st != LDLT_OK) {
      fprintf(stderr, "LDLT: solve failed: %s\n", ldlt_status_str(st));
      *info = -1;
    }

  } else if (*iopt == 3) {
    ldlt_factors_t *facs = (ldlt_factors_t *)*f_factors;
    if (facs) {
      ldlt_free_numeric(facs->N);
      ldlt_free_symbolic(facs->S);
      free(facs);
    }
    *f_factors = (fptr)0;

  } else {
    fprintf(stderr, "Invalid iopt=%d passed to c_fortran_qdldl()\n", *iopt);
    *info = -99;
  }
}

int c_qdldl_get_neg_count_(void) { return s_last_neg_count; }
