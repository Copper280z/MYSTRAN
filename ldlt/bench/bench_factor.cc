#include <benchmark/benchmark.h>
#include "ldlt/ldlt.h"
#include <cmath>
#include <vector>

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

static double residual_tridiag(int32_t n, const std::vector<double>& x) {
    double rn = 0.0;
    for (int32_t i = 0; i < n; ++i) {
        double ax = 4.0 * x[i];
        if (i > 0) ax -= x[i-1];
        if (i + 1 < n) ax -= x[i+1];
        double r = ax - 1.0;
        rn += r * r;
    }
    return std::sqrt(rn / (n > 0 ? double(n) : 1.0));
}

static void record_tridiag_residual(benchmark::State& st, const ldlt_numeric *N,
                                    int32_t n) {
    std::vector<double> x(n, 1.0);
    ldlt_status rc = ldlt_solve(N, 1, x.data(), n);
    if (rc != LDLT_OK) {
        st.SkipWithError(ldlt_status_str(rc));
        return;
    }
    double rel = residual_tridiag(n, x);
    if (!std::isfinite(rel) || rel > 1e-8) {
        st.counters["bad_resid"] = rel;
    }
}

static void BM_Factor_Tridiag(benchmark::State& st) {
    int32_t n = st.range(0);
    std::vector<int32_t> Ap, Ai; std::vector<double> Ax;
    make_tridiag(n, Ap, Ai, Ax);
    ldlt_options opt; ldlt_options_default(&opt);
    opt.ordering = LDLT_ORDER_NATURAL;
    ldlt_symbolic *S=nullptr;
    ldlt_analyze(n, Ap.data(), Ai.data(), &opt, &S);
    bool checked_residual = false;
    for (auto _ : st) {
        ldlt_numeric *N=nullptr;
        ldlt_status rc = ldlt_factorize(S, Ap.data(), Ai.data(), Ax.data(), &opt, &N);
        if (rc != LDLT_OK) { st.SkipWithError(ldlt_status_str(rc)); break; }
        if (!checked_residual) {
            st.PauseTiming();
            record_tridiag_residual(st, N, n);
            checked_residual = true;
            st.ResumeTiming();
        }
        ldlt_free_numeric(N);
    }
    ldlt_free_symbolic(S);
    st.SetItemsProcessed(int64_t(st.iterations()) * n);
}
BENCHMARK(BM_Factor_Tridiag)->Arg(100)->Arg(1000)->Arg(10000);
