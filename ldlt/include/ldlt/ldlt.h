#ifndef LDLT_H
#define LDLT_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct ldlt_symbolic ldlt_symbolic;
typedef struct ldlt_numeric  ldlt_numeric;

typedef enum {
    LDLT_OK = 0,
    LDLT_ERR_INPUT,
    LDLT_ERR_NOMEM,
    LDLT_ERR_INDEFINITE,
    LDLT_ERR_SINGULAR,
    LDLT_ERR_UNSUPPORTED,
} ldlt_status;

typedef enum {
    LDLT_ORDER_NATURAL = 0,
    LDLT_ORDER_AMD,
    LDLT_ORDER_METIS,
    LDLT_ORDER_USER,
} ldlt_ordering;

typedef enum {
    LDLT_PIVOT_NONE = 0,
    LDLT_PIVOT_BUNCH_KAUFMAN, /* dynamic 1x1/2x2 symmetric pivoting */
} ldlt_pivot_strategy;

typedef struct {
    ldlt_ordering ordering;
    const int32_t *user_perm;       /* used when ordering == LDLT_ORDER_USER */
    int32_t  relax_supernodes;      /* amalgamation slack; 0 disables relaxation */
    int32_t  nthreads;              /* 0 => OpenMP default */
    ldlt_pivot_strategy pivot;
    double   pivot_threshold;       /* <=0 => default Bunch-Kaufman threshold */
    int      require_pd;            /* if non-zero, fail on non-positive D */
} ldlt_options;

void ldlt_options_default(ldlt_options *opt);

/* Input: lower triangle of symmetric matrix in CSC, sorted row indices per column. */
ldlt_status ldlt_analyze(int32_t n,
                         const int32_t *Ap, const int32_t *Ai,
                         const ldlt_options *opt,
                         ldlt_symbolic **out);

ldlt_status ldlt_factorize(const ldlt_symbolic *S,
                           const int32_t *Ap, const int32_t *Ai,
                           const double  *Ax,
                           const ldlt_options *opt,
                           ldlt_numeric **out);

/* B is column-major, n rows, nrhs cols, leading dim ldb. Solved in place. */
ldlt_status ldlt_solve(const ldlt_numeric *N,
                       int32_t nrhs, double *B, int32_t ldb);

/* Diagnostics. */
int32_t ldlt_n(const ldlt_symbolic *S);
int32_t ldlt_nsuper(const ldlt_symbolic *S);
int64_t ldlt_lnz(const ldlt_symbolic *S);
double  ldlt_flops(const ldlt_symbolic *S);
int32_t ldlt_nlevels(const ldlt_symbolic *S);
int32_t ldlt_max_level_width(const ldlt_symbolic *S);

void ldlt_free_symbolic(ldlt_symbolic *S);
void ldlt_free_numeric(ldlt_numeric *N);

const char *ldlt_status_str(ldlt_status s);

#ifdef __cplusplus
}
#endif
#endif /* LDLT_H */
