/* Optional SuiteSparse fixtures: enabled when LDLT_FIXTURE_DIR is set at build time
   and a manifest of .mtx files has been fetched by matrices/fetch_suitesparse.py. */

#include <gtest/gtest.h>
#include "ldlt/ldlt.h"
#include <cstdlib>
#include <cstdio>
#include <filesystem>
#include <vector>
#include <string>
#include <random>
#include <cmath>

extern "C" {
struct LDLT_INT { int dummy; };
ldlt_status ldlt_read_mm_symmetric(const char *path, int32_t *n_out,
                                   int32_t **Ap_out, int32_t **Ai_out, double **Ax_out);
}

static std::string fixture_dir() {
    if (const char *e = std::getenv("LDLT_FIXTURE_DIR")) return e;
#ifdef LDLT_FIXTURE_DIR_DEFAULT
    return LDLT_FIXTURE_DIR_DEFAULT;
#else
    return "matrices";
#endif
}

class SuiteSparseFixture : public ::testing::TestWithParam<std::string> {};

TEST_P(SuiteSparseFixture, AnalyzeFactorSolve) {
    namespace fs = std::filesystem;
    fs::path p = fs::path(fixture_dir()) / GetParam();
    if (!fs::exists(p)) GTEST_SKIP() << "missing fixture: " << p;

    int32_t n=0, *Ap=nullptr, *Ai=nullptr; double *Ax=nullptr;
    ASSERT_EQ(ldlt_read_mm_symmetric(p.string().c_str(), &n, &Ap, &Ai, &Ax), LDLT_OK);

    ldlt_options opt; ldlt_options_default(&opt);
    opt.ordering = LDLT_ORDER_AMD; /* falls back to identity if AMD unavailable */
    ldlt_symbolic *S=nullptr;
    ASSERT_EQ(ldlt_analyze(n, Ap, Ai, &opt, &S), LDLT_OK);
    ldlt_numeric *N=nullptr;
    auto fst = ldlt_factorize(S, Ap, Ai, Ax, &opt, &N);
    if (fst == LDLT_ERR_INDEFINITE) {
        GTEST_SKIP() << "matrix needs pivoting; this fixture test uses the no-pivot path";
    }
    ASSERT_EQ(fst, LDLT_OK);

    std::mt19937 rng(42);
    std::uniform_real_distribution<double> U(-1.0, 1.0);
    std::vector<double> b(n), x(n);
    for (auto &v : b) v = U(rng);
    x = b;
    ASSERT_EQ(ldlt_solve(N, 1, x.data(), n), LDLT_OK);

    /* residual ||Ax - b|| / ||b|| */
    std::vector<double> Axv(n, 0.0);
    for (int32_t j = 0; j < n; ++j)
        for (int32_t kk = Ap[j]; kk < Ap[j+1]; ++kk) {
            int32_t i = Ai[kk]; double v = Ax[kk];
            Axv[i] += v * x[j];
            if (i != j) Axv[j] += v * x[i];
        }
    double rn=0, bn=0;
    for (int32_t i = 0; i < n; ++i) { double r = Axv[i]-b[i]; rn += r*r; bn += b[i]*b[i]; }
    EXPECT_LT(std::sqrt(rn/(bn>0?bn:1.0)), 1e-8);

    ldlt_free_numeric(N); ldlt_free_symbolic(S);
    free(Ap); free(Ai); free(Ax);
}

INSTANTIATE_TEST_SUITE_P(SmallSPD, SuiteSparseFixture,
    ::testing::Values("bcsstk01.mtx"));
