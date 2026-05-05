/* CHOLMOD reference benchmarks: supernodal LL^T and simplicial LDL^T.
   Built only when -DLDLT_HAVE_CHOLMOD=1 is set (cholmod found by meson). */

#include <benchmark/benchmark.h>

#ifdef LDLT_HAVE_CHOLMOD

#include <cholmod.h>
#include <vector>
#include <cmath>
#include <cstdint>

static void make_tridiag(int32_t n, std::vector<int32_t>& Ap, std::vector<int32_t>& Ai,
                         std::vector<double>& Ax) {
    Ap.assign(n+1, 0);
    for (int32_t j = 0; j < n; ++j) Ap[j+1] = Ap[j] + (j < n-1 ? 2 : 1);
    Ai.resize(Ap[n]); Ax.resize(Ap[n]);
    for (int32_t j = 0; j < n; ++j) {
        int32_t p = Ap[j];
        Ai[p] = j; Ax[p] = 4.0;
        if (j < n-1) { Ai[p+1] = j+1; Ax[p+1] = -1.0; }
    }
}

/* Build a CHOLMOD lower-triangle CSC sparse matrix from our tridiag.
   stype = -1 means symmetric, lower triangle stored. */
static cholmod_sparse *make_chm_tridiag(cholmod_common *c, int32_t n) {
    std::vector<int32_t> Ap, Ai; std::vector<double> Ax;
    make_tridiag(n, Ap, Ai, Ax);
    cholmod_sparse *A = cholmod_allocate_sparse(n, n, Ap[n], 1, 1, -1, CHOLMOD_REAL, c);
    int *Cp = (int*)A->p; int *Ci = (int*)A->i; double *Cx = (double*)A->x;
    for (int32_t j = 0; j <= n; ++j) Cp[j] = Ap[j];
    for (int32_t k = 0; k < Ap[n]; ++k) { Ci[k] = Ai[k]; Cx[k] = Ax[k]; }
    return A;
}

static double residual_tridiag(int32_t n, const double *X, int64_t ldX, int32_t nrhs) {
    double worst = 0.0;
    for (int32_t rhs = 0; rhs < nrhs; ++rhs) {
        const double *x = X + (int64_t)rhs * ldX;
        double rn = 0.0;
        for (int32_t i = 0; i < n; ++i) {
            double ax = 4.0 * x[i];
            if (i > 0) ax -= x[i-1];
            if (i + 1 < n) ax -= x[i+1];
            double r = ax - 1.0;
            rn += r * r;
        }
        double rel = std::sqrt(rn / (n > 0 ? double(n) : 1.0));
        if (rel > worst || !std::isfinite(rel)) worst = rel;
    }
    return worst;
}

static void record_chm_tridiag_residual(benchmark::State& st, cholmod_common *c,
                                        cholmod_factor *L, int32_t n, int32_t nrhs) {
    cholmod_dense *B = cholmod_zeros(n, nrhs, CHOLMOD_REAL, c);
    for (int64_t i = 0; i < (int64_t)n * nrhs; ++i) ((double*)B->x)[i] = 1.0;
    cholmod_dense *X = cholmod_solve(CHOLMOD_A, L, B, c);
    bool ok = X != nullptr && c->status == CHOLMOD_OK;
    if (ok) {
        double rel = residual_tridiag(n, (const double*)X->x, (int64_t)X->d, nrhs);
        if (!std::isfinite(rel) || rel > 1e-8) st.counters["bad_resid"] = rel;
    } else {
        st.SkipWithError("CHOLMOD solve failed");
    }
    if (X) cholmod_free_dense(&X, c);
    cholmod_free_dense(&B, c);
}

static void BM_CHOLMOD_Supernodal_LLT(benchmark::State& st) {
    int32_t n = st.range(0);
    cholmod_common c; cholmod_start(&c);
    c.final_ll = 1;             /* keep as LL^T */
    c.supernodal = CHOLMOD_SUPERNODAL;
    cholmod_sparse *A = make_chm_tridiag(&c, n);
    cholmod_factor *L = cholmod_analyze(A, &c);
    bool checked_residual = false;
    for (auto _ : st) {
        int ok = cholmod_factorize(A, L, &c);
        if (!ok || c.status != CHOLMOD_OK || L->minor != (size_t)n) {
            st.SkipWithError("CHOLMOD factorization failed");
            break;
        }
        if (!checked_residual) {
            st.PauseTiming();
            record_chm_tridiag_residual(st, &c, L, n, 1);
            checked_residual = true;
            st.ResumeTiming();
        }
    }
    cholmod_free_factor(&L, &c);
    cholmod_free_sparse(&A, &c);
    cholmod_finish(&c);
    st.SetItemsProcessed(int64_t(st.iterations()) * n);
}
BENCHMARK(BM_CHOLMOD_Supernodal_LLT)->Arg(100)->Arg(1000)->Arg(10000);

static void BM_CHOLMOD_Simplicial_LDLT(benchmark::State& st) {
    int32_t n = st.range(0);
    cholmod_common c; cholmod_start(&c);
    c.final_ll = 0;             /* LDL^T */
    c.supernodal = CHOLMOD_SIMPLICIAL;
    cholmod_sparse *A = make_chm_tridiag(&c, n);
    cholmod_factor *L = cholmod_analyze(A, &c);
    bool checked_residual = false;
    for (auto _ : st) {
        int ok = cholmod_factorize(A, L, &c);
        if (!ok || c.status != CHOLMOD_OK || L->minor != (size_t)n) {
            st.SkipWithError("CHOLMOD factorization failed");
            break;
        }
        if (!checked_residual) {
            st.PauseTiming();
            record_chm_tridiag_residual(st, &c, L, n, 1);
            checked_residual = true;
            st.ResumeTiming();
        }
    }
    cholmod_free_factor(&L, &c);
    cholmod_free_sparse(&A, &c);
    cholmod_finish(&c);
    st.SetItemsProcessed(int64_t(st.iterations()) * n);
}
BENCHMARK(BM_CHOLMOD_Simplicial_LDLT)->Arg(100)->Arg(1000)->Arg(10000);

static void BM_CHOLMOD_Solve(benchmark::State& st) {
    int32_t n = st.range(0); int32_t nrhs = st.range(1);
    cholmod_common c; cholmod_start(&c);
    c.final_ll = 0;
    c.supernodal = CHOLMOD_SIMPLICIAL;
    cholmod_sparse *A = make_chm_tridiag(&c, n);
    cholmod_factor *L = cholmod_analyze(A, &c);
    int ok = cholmod_factorize(A, L, &c);
    if (!ok || c.status != CHOLMOD_OK || L->minor != (size_t)n) {
        st.SkipWithError("CHOLMOD factorization failed");
        cholmod_free_factor(&L, &c);
        cholmod_free_sparse(&A, &c);
        cholmod_finish(&c);
        return;
    }
    cholmod_dense *B = cholmod_zeros(n, nrhs, CHOLMOD_REAL, &c);
    for (int i = 0; i < n*nrhs; ++i) ((double*)B->x)[i] = 1.0;
    bool checked_residual = false;
    for (auto _ : st) {
        cholmod_dense *X = cholmod_solve(CHOLMOD_A, L, B, &c);
        if (!X || c.status != CHOLMOD_OK) {
            st.SkipWithError("CHOLMOD solve failed");
            if (X) cholmod_free_dense(&X, &c);
            break;
        }
        if (!checked_residual) {
            st.PauseTiming();
            double rel = residual_tridiag(n, (const double*)X->x, (int64_t)X->d, nrhs);
            if (!std::isfinite(rel) || rel > 1e-8) st.counters["bad_resid"] = rel;
            checked_residual = true;
            st.ResumeTiming();
        }
        cholmod_free_dense(&X, &c);
    }
    cholmod_free_dense(&B, &c);
    cholmod_free_factor(&L, &c);
    cholmod_free_sparse(&A, &c);
    cholmod_finish(&c);
}
BENCHMARK(BM_CHOLMOD_Solve)->Args({1000,1})->Args({1000,8})->Args({10000,1});

#endif /* LDLT_HAVE_CHOLMOD */
