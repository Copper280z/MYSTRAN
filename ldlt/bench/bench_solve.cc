#include <benchmark/benchmark.h>
#include "ldlt/ldlt.h"
#include <cmath>
#include <vector>

static double residual_tridiag(int32_t n, const double *X, int32_t nrhs) {
    double worst = 0.0;
    for (int32_t rhs = 0; rhs < nrhs; ++rhs) {
        const double *x = X + (size_t)rhs * n;
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

static void BM_Solve_Tridiag(benchmark::State& st) {
    int32_t n = st.range(0);
    int32_t nrhs = st.range(1);
    std::vector<int32_t> Ap(n+1, 0), Ai;
    for (int32_t j = 0; j < n; ++j) Ap[j+1] = Ap[j] + (j < n-1 ? 2 : 1);
    Ai.resize(Ap[n]); std::vector<double> Ax(Ap[n]);
    for (int32_t j = 0; j < n; ++j) {
        int32_t p = Ap[j]; Ai[p]=j; Ax[p]=4.0;
        if (j<n-1) { Ai[p+1]=j+1; Ax[p+1]=-1.0; }
    }
    ldlt_options opt; ldlt_options_default(&opt);
    opt.ordering = LDLT_ORDER_NATURAL;
    ldlt_symbolic *S=nullptr;
    ldlt_analyze(n, Ap.data(), Ai.data(), &opt, &S);
    ldlt_numeric *N=nullptr;
    ldlt_status rc = ldlt_factorize(S, Ap.data(), Ai.data(), Ax.data(), &opt, &N);
    if (rc != LDLT_OK) {
        st.SkipWithError(ldlt_status_str(rc));
        ldlt_free_symbolic(S);
        return;
    }
    std::vector<double> B(n*nrhs, 1.0);
    bool checked_residual = false;
    for (auto _ : st) {
        auto X = B;
        rc = ldlt_solve(N, nrhs, X.data(), n);
        if (rc != LDLT_OK) { st.SkipWithError(ldlt_status_str(rc)); break; }
        if (!checked_residual) {
            st.PauseTiming();
            double rel = residual_tridiag(n, X.data(), nrhs);
            if (!std::isfinite(rel) || rel > 1e-8) st.counters["bad_resid"] = rel;
            checked_residual = true;
            st.ResumeTiming();
        }
        benchmark::DoNotOptimize(X);
    }
    ldlt_free_numeric(N); ldlt_free_symbolic(S);
}
BENCHMARK(BM_Solve_Tridiag)->Args({1000,1})->Args({1000,8})->Args({10000,1});
