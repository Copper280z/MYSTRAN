/*
 * Fortran-callable interface to QDLDL (LDL^T sparse factorization).
 *
 * Matches the iopt-based interface of c_fortran_dgssv so callers can switch
 * between SuperLU and QDLDL by changing only SPARSE_FLAVOR in the BDF deck.
 *
 * Input matrix convention (same as c_fortran_dgssv):
 *   The Fortran caller passes a symmetric matrix in CRS format:
 *     colptr  = I_MATIN (Fortran 1-based row pointers)
 *     rowind  = J_MATIN (Fortran 1-based column indices)
 *   This routine extracts the upper-triangular part (col >= row) and converts
 *   to CSC for QDLDL.  Diagonal entries are included.
 *
 * METIS_NodeND fill-reducing ordering is applied when HAVE_METIS is defined.
 *
 * QDLDL_int is int64_t so L column-pointer arrays never overflow.
 * Intermediate arrays (CRS row indices, CSC row indices, METIS adjacency)
 * use plain int (int32) to keep peak memory low on large structural models.
 *
 * iopt:
 *   1 = LDL^T factorization
 *   2 = triangular solve (reuses stored factors)
 *   3 = free all storage
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "qdldl/include/qdldl.h"

#if (HAVE_METIS)
#include "metis.h"
#endif

/* Pointer-sized integer to pass handles back to Fortran as INTEGER(DBL_LONG) */
typedef long long int fptr;

static QDLDL_int s_last_neg_count = 0;

typedef struct {
  int         n;
  QDLDL_int   nnz_L;
  QDLDL_int  *Lp;    /* int64 column pointers (sum = nnz_L, may exceed 2^31) */
  int        *Li;    /* int32 row indices (values < n < 2^31)                 */
  double     *Lx;
  double     *D;
  double     *Dinv;
  int        *perm;  /* fill-reducing permutation: new index = perm[old]      */
  int        *iperm; /* inverse permutation:       old index = iperm[new]     */
  QDLDL_int   pos_count;
} qdldl_factors_t;

static int validate_perm_pair(int n, const int *old_to_new,
                              const int *new_to_old) {
  int *seen = (int *)calloc((size_t)n, sizeof(int));
  if (!seen) return 0;

  for (int old = 0; old < n; old++) {
    int newidx = old_to_new[old];
    if (newidx < 0 || newidx >= n || seen[newidx] ||
        new_to_old[newidx] != old) {
      free(seen);
      return 0;
    }
    seen[newidx] = 1;
  }

  free(seen);
  return 1;
}

/* ============================================================================
 * Wall-clock timer helpers
 * ========================================================================= */
static double now_sec(void) {
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

static double qdldl_step_done(const char *label, double t0) {
  double t1 = now_sec();
  printf("QDLDL  %-36s  %8.3f s\n", label, t1 - t0);
  fflush(stdout);
  return t1;
}

/* ============================================================================
 * LDL^T triangular solve using int32 Li (avoids per-call int64 promotion).
 * Mirrors QDLDL_Lsolve / QDLDL_Ltsolve / QDLDL_solve but takes int* Li.
 * ========================================================================= */
static void my_Lsolve(int n, const QDLDL_int *Lp, const int *Li,
                      const double *Lx, double *x) {
  for (int i = 0; i < n; i++) {
    double val = x[i];
    for (QDLDL_int j = Lp[i]; j < Lp[i + 1]; j++)
      x[Li[j]] -= Lx[j] * val;
  }
}

static void my_Ltsolve(int n, const QDLDL_int *Lp, const int *Li,
                       const double *Lx, double *x) {
  for (int i = n - 1; i >= 0; i--) {
    double val = x[i];
    for (QDLDL_int j = Lp[i]; j < Lp[i + 1]; j++)
      val -= Lx[j] * x[Li[j]];
    x[i] = val;
  }
}

static void my_solve(int n, const QDLDL_int *Lp, const int *Li,
                     const double *Lx, const double *Dinv, double *x) {
  my_Lsolve(n, Lp, Li, Lx, x);
  for (int i = 0; i < n; i++) x[i] *= Dinv[i];
  my_Ltsolve(n, Lp, Li, Lx, x);
}

/* ============================================================================
 * Insertion sort for small arrays — used to sort row indices within each
 * column of the permuted CSC.  Fast when column lengths are small (O(log n)
 * with nested-dissection ordering).
 * ========================================================================= */
static void isort2(int *keys, double *vals, int len) {
  for (int i = 1; i < len; i++) {
    int    k = keys[i];
    double v = vals[i];
    int    j = i - 1;
    while (j >= 0 && keys[j] > k) {
      keys[j + 1] = keys[j];
      vals[j + 1] = vals[j];
      j--;
    }
    keys[j + 1] = k;
    vals[j + 1] = v;
  }
}

/* ============================================================================
 * Build upper-triangular (row <= col) CSC from a symmetric CRS matrix.
 *
 * Inputs (0-based):
 *   row_ptr[0..N]   CRS row pointers
 *   col_idx[0..NNZ] CRS column indices
 *   values [0..NNZ] CRS values
 *
 * Outputs (allocated by caller with sizes n+1, nnz_ut, nnz_ut):
 *   Ap  — int32 CSC column pointers
 *   Ai  — int32 CSC row indices (sorted ascending within each column)
 *   Ax  — double values
 *
 * Returns nnz_ut (number of upper-tri + diagonal entries).
 * ========================================================================= */
static int build_upper_csc(int N, const int *row_ptr, const int *col_idx,
                            const double *values, int *Ap, int *Ai,
                            double *Ax) {
  /* Pass 1: count entries per column. */
  memset(Ap, 0, (size_t)(N + 1) * sizeof(int));
  for (int i = 0; i < N; i++)
    for (int k = row_ptr[i]; k < row_ptr[i + 1]; k++)
      if (col_idx[k] >= i)
        Ap[col_idx[k] + 1]++;

  for (int j = 0; j < N; j++)
    Ap[j + 1] += Ap[j];

  int nnz_ut = Ap[N];

  /* Pass 2: fill Ai/Ax.  Because i iterates 0..N-1, row indices are inserted
   * in strictly increasing order within each column — no sort needed.         */
  int *cursor = (int *)malloc((size_t)N * sizeof(int));
  if (!cursor) return -1;
  memcpy(cursor, Ap, (size_t)N * sizeof(int));

  for (int i = 0; i < N; i++) {
    for (int k = row_ptr[i]; k < row_ptr[i + 1]; k++) {
      int j = col_idx[k];
      if (j >= i) {
        int pos = cursor[j]++;
        Ai[pos] = i;
        Ax[pos] = values[k];
      }
    }
  }
  free(cursor);
  return nnz_ut;
}

/* ============================================================================
 * Main Fortran-callable entry point.
 * ========================================================================= */
void c_fortran_qdldl_(int *iopt, int *n, int *nnz, int *nrhs, double *values,
                      int *rowind, int *colptr, double *b, int *ldb,
                      fptr *f_factors, int *info) {
  *info = 0;

  /* -----------------------------------------------------------------------
   * iopt == 1: LDL^T factorization
   * --------------------------------------------------------------------- */
  if (*iopt == 1) {

    int    N   = *n;
    int    NNZ = *nnz;
    double t0, t_step;

    printf("QDLDL  factorization begin  (n=%d, nnz=%d)\n", N, NNZ);
    fflush(stdout);
    t0 = now_sec();

    /* --- 1. Convert 1-based Fortran CRS to 0-based C ------------------- */
    int *row_ptr = (int *)malloc((size_t)(N + 1) * sizeof(int));
    int *col_idx = (int *)malloc((size_t)NNZ * sizeof(int));
    if (!row_ptr || !col_idx) { *info = -1; free(row_ptr); free(col_idx); return; }

    for (int i = 0; i <= N; i++) row_ptr[i] = colptr[i] - 1;
    for (int k = 0; k < NNZ; k++) col_idx[k] = rowind[k] - 1;

    /* --- 2. Build upper-tri CSC directly from CRS (no triplet array) --- */
    /* First count to allocate. */
    int nnz_ut = 0;
    for (int i = 0; i < N; i++)
      for (int k = row_ptr[i]; k < row_ptr[i + 1]; k++)
        if (col_idx[k] >= i) nnz_ut++;

    int    *orig_Ap = (int    *)malloc((size_t)(N + 1) * sizeof(int));
    int    *orig_Ai = (int    *)malloc((size_t)nnz_ut  * sizeof(int));
    double *orig_Ax = (double *)malloc((size_t)nnz_ut  * sizeof(double));
    double *orig_diag = (double *)calloc((size_t)N, sizeof(double));

    if (!orig_Ap || !orig_Ai || !orig_Ax || !orig_diag) {
      *info = -1;
      free(row_ptr); free(col_idx);
      free(orig_Ap); free(orig_Ai); free(orig_Ax); free(orig_diag);
      return;
    }

    nnz_ut = build_upper_csc(N, row_ptr, col_idx, values,
                              orig_Ap, orig_Ai, orig_Ax);
    if (nnz_ut < 0) {
      *info = -1;
      free(row_ptr); free(col_idx);
      free(orig_Ap); free(orig_Ai); free(orig_Ax); free(orig_diag);
      return;
    }

    /* Collect diagonal for MAXRATIO. */
    for (int j = 0; j < N; j++)
      for (int p = orig_Ap[j]; p < orig_Ap[j + 1]; p++)
        if (orig_Ai[p] == j) { orig_diag[j] = orig_Ax[p]; break; }

    free(row_ptr);
    free(col_idx);

    t_step = qdldl_step_done("CRS extract upper-tri + diagonal", t0);

    /* --- 3. METIS fill-reducing ordering ------------------------------- */
    int *perm  = (int *)malloc((size_t)N * sizeof(int));
    int *iperm = (int *)malloc((size_t)N * sizeof(int));
    if (!perm || !iperm) {
      *info = -1;
      free(orig_Ap); free(orig_Ai); free(orig_Ax); free(orig_diag);
      free(perm); free(iperm);
      return;
    }

#if (HAVE_METIS)
    {
      /* Build symmetric adjacency graph from the upper-tri CSC.
       * Off-diagonal upper-tri entry (i,j) with i<j contributes edges i↔j. */
      int nnz_offdiag = nnz_ut - N; /* subtract diagonal entries */
      int *deg   = (int *)calloc((size_t)N, sizeof(int));
      int *xadj  = (int *)malloc((size_t)(N + 1) * sizeof(int));
      int *adjncy = (int *)malloc((size_t)(2 * nnz_offdiag) * sizeof(int));
      if (!deg || !xadj || !adjncy) {
        *info = -1;
        free(deg); free(xadj); free(adjncy);
        free(orig_Ap); free(orig_Ai); free(orig_Ax); free(orig_diag);
        free(perm); free(iperm);
        return;
      }

      for (int j = 0; j < N; j++)
        for (int p = orig_Ap[j]; p < orig_Ap[j + 1]; p++) {
          int i = orig_Ai[p];
          if (i < j) { deg[i]++; deg[j]++; }
        }

      xadj[0] = 0;
      for (int i = 0; i < N; i++) xadj[i + 1] = xadj[i] + deg[i];
      memset(deg, 0, (size_t)N * sizeof(int));

      for (int j = 0; j < N; j++)
        for (int p = orig_Ap[j]; p < orig_Ap[j + 1]; p++) {
          int i = orig_Ai[p];
          if (i < j) {
            adjncy[xadj[i] + deg[i]++] = j;
            adjncy[xadj[j] + deg[j]++] = i;
          }
        }
      free(deg);

      t_step = qdldl_step_done("METIS adjacency graph build", t_step);

      idx_t nv = (idx_t)N;
      int metis_status = METIS_NodeND(&nv, (idx_t *)xadj, (idx_t *)adjncy,
                                      NULL, NULL, (idx_t *)perm,
                                      (idx_t *)iperm);
      free(xadj);
      free(adjncy);

      if (metis_status != METIS_OK) {
        fprintf(stderr, "QDLDL: METIS_NodeND failed with status %d.\n",
                metis_status);
        *info = 1;
        free(orig_Ap); free(orig_Ai); free(orig_Ax); free(orig_diag);
        free(perm); free(iperm);
        return;
      }

      /* METIS returns perm[new] = old and iperm[old] = new.  The rest of
       * this wrapper uses perm[old] = new and iperm[new] = old. */
      int *tmp_perm = perm;
      perm = iperm;
      iperm = tmp_perm;

      if (!validate_perm_pair(N, perm, iperm)) {
        fprintf(stderr, "QDLDL: METIS_NodeND returned an invalid permutation.\n");
        *info = 1;
        free(orig_Ap); free(orig_Ai); free(orig_Ax); free(orig_diag);
        free(perm); free(iperm);
        return;
      }

      t_step = qdldl_step_done("METIS_NodeND ordering", t_step);
    }
#else
    printf("QDLDL  WARNING: METIS not available — using identity ordering.\n");
    printf("QDLDL           Fill-in will be large for structural meshes.\n");
    fflush(stdout);
    for (int i = 0; i < N; i++) { perm[i] = i; iperm[i] = i; }
#endif

    /* --- 4. Build permuted upper-tri CSC --------------------------------
     * For each entry (row=i, col=j) in original CSC, the permuted entry is
     * (perm[i], perm[j]).  Keep upper-tri: if perm[i] > perm[j] swap them.
     *
     * Two-pass approach: count into permuted columns, then fill.
     * Rows within each permuted column are not sorted after filling, so we
     * sort each column individually — fast since column lengths are small
     * (O(log n) with nested-dissection ordering).                          */

    int *perm_Ap_int = (int *)calloc((size_t)(N + 1), sizeof(int));
    if (!perm_Ap_int) {
      *info = -1;
      free(orig_Ap); free(orig_Ai); free(orig_Ax); free(orig_diag);
      free(perm); free(iperm);
      return;
    }

    for (int j = 0; j < N; j++)
      for (int p = orig_Ap[j]; p < orig_Ap[j + 1]; p++) {
        int r = perm[orig_Ai[p]], c = perm[j];
        if (r > c) { int t = r; r = c; c = t; }
        perm_Ap_int[c + 1]++;
      }
    for (int j = 0; j < N; j++)
      perm_Ap_int[j + 1] += perm_Ap_int[j];

    /* perm_Ap_int[N] == nnz_ut (permutation is bijective) */
    int    *perm_Ai = (int    *)malloc((size_t)nnz_ut * sizeof(int));
    double *perm_Ax = (double *)malloc((size_t)nnz_ut * sizeof(double));
    int    *cursor  = (int    *)malloc((size_t)N       * sizeof(int));
    if (!perm_Ai || !perm_Ax || !cursor) {
      *info = -1;
      free(perm_Ap_int); free(perm_Ai); free(perm_Ax); free(cursor);
      free(orig_Ap); free(orig_Ai); free(orig_Ax); free(orig_diag);
      free(perm); free(iperm);
      return;
    }
    memcpy(cursor, perm_Ap_int, (size_t)N * sizeof(int));

    for (int j = 0; j < N; j++)
      for (int p = orig_Ap[j]; p < orig_Ap[j + 1]; p++) {
        int    r = perm[orig_Ai[p]], c = perm[j];
        double v = orig_Ax[p];
        if (r > c) { int t = r; r = c; c = t; }
        int pos = cursor[c]++;
        perm_Ai[pos] = r;
        perm_Ax[pos] = v;
      }
    free(cursor);
    free(orig_Ap); free(orig_Ai); free(orig_Ax);

    /* Sort row indices within each permuted column. */
    for (int j = 0; j < N; j++) {
      int start = perm_Ap_int[j];
      int len   = perm_Ap_int[j + 1] - start;
      if (len > 1) isort2(perm_Ai + start, perm_Ax + start, len);
    }

    t_step = qdldl_step_done("permute + sort to CSC", t_step);

    /* --- 5. Convert perm_Ap_int to QDLDL_int for QDLDL API ------------- */
    QDLDL_int *Ap = (QDLDL_int *)malloc((size_t)(N + 1) * sizeof(QDLDL_int));
    QDLDL_int *Ai = (QDLDL_int *)malloc((size_t)nnz_ut  * sizeof(QDLDL_int));
    double    *Ax = (double    *)malloc((size_t)nnz_ut  * sizeof(double));
    if (!Ap || !Ai || !Ax) {
      *info = -1;
      free(Ap); free(Ai); free(Ax);
      free(perm_Ap_int); free(perm_Ai); free(perm_Ax);
      free(orig_diag); free(perm); free(iperm);
      return;
    }
    for (int j = 0; j <= N; j++) Ap[j] = (QDLDL_int)perm_Ap_int[j];
    for (int k = 0; k < nnz_ut; k++) {
      Ai[k] = (QDLDL_int)perm_Ai[k];
      Ax[k] = perm_Ax[k];
    }
    free(perm_Ap_int); free(perm_Ai); free(perm_Ax);

    /* --- 6. Compute elimination tree ----------------------------------- */
    QDLDL_int *etree = (QDLDL_int *)malloc((size_t)N * sizeof(QDLDL_int));
    QDLDL_int *Lnz   = (QDLDL_int *)malloc((size_t)N * sizeof(QDLDL_int));
    QDLDL_int *work  = (QDLDL_int *)malloc((size_t)N * sizeof(QDLDL_int));
    if (!etree || !Lnz || !work) {
      *info = -1;
      free(etree); free(Lnz); free(work);
      free(Ap); free(Ai); free(Ax);
      free(orig_diag); free(perm); free(iperm);
      return;
    }

    QDLDL_int sumLnz = QDLDL_etree(N, Ap, Ai, work, Lnz, etree);
    free(work);

    if (sumLnz < 0) {
      if (sumLnz == -1)
        fprintf(stderr,
                "QDLDL_etree failed: empty column or lower-triangular entry "
                "(n=%d, nnz=%d).\n", N, NNZ);
      else
        fprintf(stderr,
                "QDLDL_etree failed: L nonzero count overflows int64 "
                "(n=%d, nnz=%d).  Matrix too large for QDLDL.\n", N, NNZ);
      *info = 1;
      free(etree); free(Lnz);
      free(Ap); free(Ai); free(Ax);
      free(perm); free(iperm); free(orig_diag);
      return;
    }

    t_step = qdldl_step_done("QDLDL_etree (elim. tree)", t_step);
    printf("QDLDL  L nonzeros (predicted) = %lld\n", (long long)sumLnz);
    fflush(stdout);

    /* --- 7. Allocate L and factor workspace ---------------------------- */
    QDLDL_int *Lp     = (QDLDL_int *)malloc((size_t)(N + 1) * sizeof(QDLDL_int));
    int       *Li_int = (int       *)malloc((size_t)sumLnz   * sizeof(int));
    double    *Lx     = (double    *)malloc((size_t)sumLnz   * sizeof(double));
    double    *D      = (double    *)malloc((size_t)N        * sizeof(double));
    double    *Dinv   = (double    *)malloc((size_t)N        * sizeof(double));

    QDLDL_bool *bwork  = (QDLDL_bool *)malloc((size_t)N        * sizeof(QDLDL_bool));
    QDLDL_int  *iwork3 = (QDLDL_int  *)malloc((size_t)(3 * N)  * sizeof(QDLDL_int));
    double     *fwork  = (double     *)malloc((size_t)N        * sizeof(double));

    if (!Lp || !Li_int || !Lx || !D || !Dinv || !bwork || !iwork3 || !fwork) {
      fprintf(stderr,
              "QDLDL: failed to allocate factor storage (sumLnz=%lld, n=%d)\n",
              (long long)sumLnz, N);
      *info = -1;
      free(Lp); free(Li_int); free(Lx); free(D); free(Dinv);
      free(bwork); free(iwork3); free(fwork);
      free(etree); free(Lnz);
      free(Ap); free(Ai); free(Ax);
      free(perm); free(iperm); free(orig_diag);
      return;
    }

    t_step = qdldl_step_done("allocate L/D/workspace", t_step);

    Lp[0] = 0;
    QDLDL_int pos_count = QDLDL_factor(N, Ap, Ai, Ax, Lp, Li_int, Lx,
                                       D, Dinv, Lnz, etree, bwork, iwork3, fwork);

    free(etree); free(Lnz); free(bwork); free(iwork3); free(fwork);
    free(Ap); free(Ai); free(Ax);

    if (pos_count < 0) {
      fprintf(stderr, "QDLDL_factor: zero pivot — matrix not quasidefinite.\n");
      *info = N;
      free(Lp); free(Li_int); free(Lx); free(D); free(Dinv);
      free(perm); free(iperm); free(orig_diag);
      return;
    }

    t_step = qdldl_step_done("QDLDL_factor (LDL^T)", t_step);

    QDLDL_int neg_count = (QDLDL_int)N - pos_count;
    s_last_neg_count = neg_count;

    printf("QDLDL  L nonzeros (actual)    = %lld\n", (long long)Lp[N]);
    printf("QDLDL  D diagonal — positive=%lld  negative=%lld  (n=%d)\n",
           (long long)pos_count, (long long)neg_count, N);
    fflush(stdout);

    /* --- 8. Compute and print MAXRATIO --------------------------------- */
    double maxratio    = -1.0;
    int    maxratio_dof = 0;
    for (int i = 0; i < N; i++) {
      double di = D[perm[i]];
      if (di != 0.0) {
        double ratio = fabs(orig_diag[i] / di);
        if (ratio > maxratio) { maxratio = ratio; maxratio_dof = i + 1; }
      }
    }
    printf("QDLDL  max |A_ii/D_ii| (MAXRATIO) = %.6E  at DOF %d\n",
           maxratio, maxratio_dof);
    free(orig_diag);

    qdldl_step_done("MAXRATIO computation", t_step);
    printf("QDLDL  factorization total elapsed  %8.3f s\n", now_sec() - t0);
    fflush(stdout);

    /* --- 9. Store factors handle --------------------------------------- */
    qdldl_factors_t *facs = (qdldl_factors_t *)malloc(sizeof(qdldl_factors_t));
    if (!facs) { *info = -1; return; }
    facs->n         = N;
    facs->nnz_L     = sumLnz;
    facs->Lp        = Lp;
    facs->Li        = Li_int;
    facs->Lx        = Lx;
    facs->D         = D;
    facs->Dinv      = Dinv;
    facs->perm      = perm;
    facs->iperm     = iperm;
    facs->pos_count = pos_count;

    *f_factors = (fptr)facs;

  /* -----------------------------------------------------------------------
   * iopt == 2: triangular solve  (reuse stored LDL^T factors)
   * --------------------------------------------------------------------- */
  } else if (*iopt == 2) {

    qdldl_factors_t *facs = (qdldl_factors_t *)*f_factors;
    int N = facs->n;

    double *y = (double *)malloc((size_t)N * sizeof(double));
    if (!y) { *info = -1; return; }

    for (int i = 0; i < N; i++) y[facs->perm[i]] = b[i];

    my_solve(N, facs->Lp, facs->Li, facs->Lx, facs->Dinv, y);

    for (int i = 0; i < N; i++) b[facs->iperm[i]] = y[i];

    free(y);

  /* -----------------------------------------------------------------------
   * iopt == 3: free all storage
   * --------------------------------------------------------------------- */
  } else if (*iopt == 3) {

    qdldl_factors_t *facs = (qdldl_factors_t *)*f_factors;
    free(facs->Lp);
    free(facs->Li);
    free(facs->Lx);
    free(facs->D);
    free(facs->Dinv);
    free(facs->perm);
    free(facs->iperm);
    free(facs);
    *f_factors = (fptr)0;

  } else {
    fprintf(stderr, "Invalid iopt=%d passed to c_fortran_qdldl()\n", *iopt);
    *info = -99;
  }
}

/* ============================================================================
 * Return the number of negative D entries from the most recent factorization.
 * Called from Fortran as C_QDLDL_GET_NEG_COUNT() after iopt=1.
 * ========================================================================= */
int c_qdldl_get_neg_count_(void) { return (int)s_last_neg_count; }
