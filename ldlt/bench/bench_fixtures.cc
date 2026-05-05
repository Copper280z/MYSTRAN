/* SuiteSparse fixture benchmarks: each fixture run against ours, CHOLMOD-supernodal,
   and CHOLMOD-simplicial. Fixture path resolved from $LDLT_FIXTURE_DIR or "matrices/". */

#include <benchmark/benchmark.h>
#include "bench_fixtures.h"
#include "ldlt/ldlt.h"
#include <algorithm>
#include <cctype>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <chrono>
#include <filesystem>
#include <iostream>
#include <mutex>
#include <string>
#include <unordered_map>
#include <vector>

extern "C" ldlt_status ldlt_read_mm_symmetric(const char *path, int32_t *n_out,
    int32_t **Ap_out, int32_t **Ai_out, double **Ax_out);

#ifdef LDLT_HAVE_CHOLMOD
#include <cholmod.h>
#endif

struct Mtx { int32_t n=0; int32_t *Ap=nullptr, *Ai=nullptr; double *Ax=nullptr; };

static constexpr double kResidualTol = 1e-8;

static std::string fixture_path(const char *name) {
    std::filesystem::path raw(name);
    if (raw.is_absolute()) return raw.string();
    const char *dir = std::getenv("LDLT_FIXTURE_DIR");
    if (dir) return (std::filesystem::path(dir) / name).string();
    /* Search: current dir "matrices/", then source-relative "../matrices/". */
    for (const char *prefix : {"matrices", "../matrices"}) {
        auto p = std::filesystem::path(prefix) / name;
        if (std::filesystem::exists(p)) return p.string();
    }
    return (std::filesystem::path("matrices") / name).string();
}

static const Mtx *load(const char *name) {
    static std::mutex m;
    static std::unordered_map<std::string, Mtx> cache;
    std::lock_guard<std::mutex> g(m);
    auto it = cache.find(name);
    if (it != cache.end()) return it->second.n ? &it->second : nullptr;
    Mtx M;
    if (ldlt_read_mm_symmetric(fixture_path(name).c_str(), &M.n, &M.Ap, &M.Ai, &M.Ax) != LDLT_OK) {
        cache[name] = {};
        return nullptr;
    }
    auto [iter, _] = cache.emplace(name, M);
    return &iter->second;
}

static double residual_ones_rhs(const Mtx *M, const double *X, int64_t ldX, int32_t nrhs) {
    double worst = 0.0;
    std::vector<double> r(M->n);
    for (int32_t rhs = 0; rhs < nrhs; ++rhs) {
        std::fill(r.begin(), r.end(), -1.0);
        const double *x = X + (int64_t)rhs * ldX;
        for (int32_t j = 0; j < M->n; ++j) {
            for (int32_t p = M->Ap[j]; p < M->Ap[j+1]; ++p) {
                int32_t i = M->Ai[p];
                double v = M->Ax[p];
                r[i] += v * x[j];
                if (i != j) r[j] += v * x[i];
            }
        }
        double rn = 0.0;
        for (double v : r) rn += v * v;
        double rel = std::sqrt(rn / (M->n > 0 ? double(M->n) : 1.0));
        if (rel > worst || !std::isfinite(rel)) worst = rel;
    }
    return worst;
}

static void report_residual_if_bad(benchmark::State &st, double rel) {
    st.counters["resid"] = rel;
    if (std::isfinite(rel) && rel <= kResidualTol) return;
    st.counters["bad_resid"] = rel;
}

static void record_ours_residual(benchmark::State &st, const Mtx *M,
                                 const ldlt_numeric *N, int32_t nrhs) {
    std::vector<double> X((size_t)M->n * nrhs, 1.0);
    ldlt_status rc = ldlt_solve(N, nrhs, X.data(), M->n);
    if (rc != LDLT_OK) {
        st.SkipWithError(ldlt_status_str(rc));
        return;
    }
    report_residual_if_bad(st, residual_ones_rhs(M, X.data(), M->n, nrhs));
}

/* ---- ours ---- */
static void run_ours(benchmark::State &st, const char *name) {
    const Mtx *M = load(name);
    if (!M) { st.SkipWithError("fixture missing — run python3 matrices/fetch_suitesparse.py"); return; }
    ldlt_options opt; ldlt_options_default(&opt);
    opt.ordering = LDLT_ORDER_METIS;
    ldlt_symbolic *S = nullptr;
    if (ldlt_analyze(M->n, M->Ap, M->Ai, &opt, &S) != LDLT_OK) {
        st.SkipWithError("analyze failed"); return;
    }
    bool checked_residual = false;
    for (auto _ : st) {
        ldlt_numeric *N = nullptr;
        ldlt_status rc = ldlt_factorize(S, M->Ap, M->Ai, M->Ax, &opt, &N);
        if (rc != LDLT_OK) { st.SkipWithError(ldlt_status_str(rc)); break; }
        if (!checked_residual) {
            st.PauseTiming();
            record_ours_residual(st, M, N, 1);
            checked_residual = true;
            st.ResumeTiming();
        }
        ldlt_free_numeric(N);
    }
    st.SetItemsProcessed(int64_t(st.iterations()) * M->n);
    st.counters["nnz"] = double(M->Ap[M->n]);
    st.counters["lnz_est"] = double(ldlt_lnz(S));
    st.counters["nlevels"] = double(ldlt_nlevels(S));
    st.counters["max_level_width"] = double(ldlt_max_level_width(S));
    ldlt_free_symbolic(S);
}

static void run_ours_indef(benchmark::State &st, const char *name) {
    const Mtx *M = load(name);
    if (!M) { st.SkipWithError("fixture missing — run python3 matrices/fetch_suitesparse.py"); return; }
    ldlt_options opt; ldlt_options_default(&opt);
    opt.ordering = LDLT_ORDER_METIS;
    opt.pivot = LDLT_PIVOT_BUNCH_KAUFMAN;
    ldlt_symbolic *S = nullptr;
    if (ldlt_analyze(M->n, M->Ap, M->Ai, &opt, &S) != LDLT_OK) {
        st.SkipWithError("analyze failed"); return;
    }
    bool checked_residual = false;
    for (auto _ : st) {
        ldlt_numeric *N = nullptr;
        ldlt_status rc = ldlt_factorize(S, M->Ap, M->Ai, M->Ax, &opt, &N);
        if (rc != LDLT_OK) { st.SkipWithError(ldlt_status_str(rc)); break; }
        if (!checked_residual) {
            st.PauseTiming();
            record_ours_residual(st, M, N, 1);
            checked_residual = true;
            st.ResumeTiming();
        }
        ldlt_free_numeric(N);
    }
    st.SetItemsProcessed(int64_t(st.iterations()) * M->n);
    st.counters["nnz"] = double(M->Ap[M->n]);
    st.counters["lnz_est"] = double(ldlt_lnz(S));
    st.counters["nlevels"] = double(ldlt_nlevels(S));
    st.counters["max_level_width"] = double(ldlt_max_level_width(S));
    ldlt_free_symbolic(S);
}

static void run_ours_solve(benchmark::State &st, const char *name, int32_t nrhs) {
    const Mtx *M = load(name);
    if (!M) { st.SkipWithError("fixture missing — run python3 matrices/fetch_suitesparse.py"); return; }
    ldlt_options opt; ldlt_options_default(&opt);
    opt.ordering = LDLT_ORDER_METIS;
    ldlt_symbolic *S = nullptr;
    if (ldlt_analyze(M->n, M->Ap, M->Ai, &opt, &S) != LDLT_OK) {
        st.SkipWithError("analyze failed"); return;
    }
    ldlt_numeric *N = nullptr;
    ldlt_status rc = ldlt_factorize(S, M->Ap, M->Ai, M->Ax, &opt, &N);
    if (rc != LDLT_OK) {
        st.SkipWithError(ldlt_status_str(rc));
        ldlt_free_symbolic(S);
        return;
    }

    std::vector<double> B((size_t)M->n * nrhs, 1.0);
    bool checked_residual = false;
    for (auto _ : st) {
        std::vector<double> X = B;
        rc = ldlt_solve(N, nrhs, X.data(), M->n);
        if (rc != LDLT_OK) { st.SkipWithError(ldlt_status_str(rc)); break; }
        if (!checked_residual) {
            st.PauseTiming();
            report_residual_if_bad(st, residual_ones_rhs(M, X.data(), M->n, nrhs));
            checked_residual = true;
            st.ResumeTiming();
        }
        benchmark::DoNotOptimize(X.data());
    }
    st.SetItemsProcessed(int64_t(st.iterations()) * M->n * nrhs);
    ldlt_free_numeric(N);
    ldlt_free_symbolic(S);
}

#ifdef LDLT_HAVE_CHOLMOD
static cholmod_sparse *to_chm(cholmod_common *c, const Mtx *M) {
    cholmod_sparse *A = cholmod_allocate_sparse(M->n, M->n, M->Ap[M->n], 1, 1, -1, CHOLMOD_REAL, c);
    int *Cp = (int*)A->p; int *Ci = (int*)A->i; double *Cx = (double*)A->x;
    for (int32_t j = 0; j <= M->n; ++j) Cp[j] = M->Ap[j];
    for (int32_t k = 0; k < M->Ap[M->n]; ++k) { Ci[k] = M->Ai[k]; Cx[k] = M->Ax[k]; }
    return A;
}

/* Force METIS so fill-in matches our run_ours (which uses LDLT_ORDER_METIS).
   Without this, CHOLMOD's default analyze picks among AMD/METIS by trial fill,
   which can give a different L and a misleading head-to-head. */
static void chm_force_metis(cholmod_common *c) {
    c->nmethods = 1;
    c->method[0].ordering = CHOLMOD_METIS;
    c->postorder = 1;
}

static bool chm_factor_ok(cholmod_common *c, const cholmod_factor *L, int32_t n) {
    return c->status == CHOLMOD_OK && L && L->minor == (size_t)n;
}

static double elapsed_ms(std::chrono::steady_clock::time_point t0,
                         std::chrono::steady_clock::time_point t1) {
    return std::chrono::duration<double, std::milli>(t1 - t0).count();
}

static void record_chm_residual(benchmark::State &st, const Mtx *M,
                                cholmod_common *c, cholmod_factor *L,
                                int32_t nrhs) {
    cholmod_dense *B = cholmod_zeros(M->n, nrhs, CHOLMOD_REAL, c);
    for (int64_t i = 0; i < (int64_t)M->n * nrhs; ++i) ((double*)B->x)[i] = 1.0;
    cholmod_dense *X = cholmod_solve(CHOLMOD_A, L, B, c);
    bool ok = X != nullptr && c->status == CHOLMOD_OK;
    if (ok) {
        double rel = residual_ones_rhs(M, (const double*)X->x, (int64_t)X->d, nrhs);
        report_residual_if_bad(st, rel);
    } else {
        st.SkipWithError("CHOLMOD solve failed");
    }
    if (X) cholmod_free_dense(&X, c);
    cholmod_free_dense(&B, c);
}

static void run_chm(benchmark::State &st, const char *name, int supernodal) {
    const Mtx *M = load(name);
    if (!M) { st.SkipWithError("fixture missing"); return; }
    cholmod_common c; cholmod_start(&c);
    c.final_ll = supernodal ? 1 : 0;
    c.supernodal = supernodal ? CHOLMOD_SUPERNODAL : CHOLMOD_SIMPLICIAL;
    chm_force_metis(&c);
    cholmod_sparse *A = to_chm(&c, M);
    cholmod_factor *L = cholmod_analyze(A, &c);
    bool checked_residual = false;
    for (auto _ : st) {
        int ok = cholmod_factorize(A, L, &c);
        if (!ok || !chm_factor_ok(&c, L, M->n)) {
            st.SkipWithError("CHOLMOD factorization failed");
            break;
        }
        if (!checked_residual) {
            st.PauseTiming();
            record_chm_residual(st, M, &c, L, 1);
            checked_residual = true;
            st.ResumeTiming();
        }
    }
    st.SetItemsProcessed(int64_t(st.iterations()) * M->n);
    st.counters["chm_lnz"] = double(L->xsize ? L->xsize : L->nzmax);
    cholmod_free_factor(&L, &c);
    cholmod_free_sparse(&A, &c);
    cholmod_finish(&c);
}

static void run_chm_solve(benchmark::State &st, const char *name, int32_t nrhs) {
    const Mtx *M = load(name);
    if (!M) { st.SkipWithError("fixture missing"); return; }
    cholmod_common c; cholmod_start(&c);
    c.final_ll = 1;
    c.supernodal = CHOLMOD_SUPERNODAL;
    cholmod_sparse *A = to_chm(&c, M);
    cholmod_factor *L = cholmod_analyze(A, &c);
    int ok = cholmod_factorize(A, L, &c);
    if (!ok || !chm_factor_ok(&c, L, M->n)) {
        st.SkipWithError("CHOLMOD factorization failed");
        cholmod_free_factor(&L, &c);
        cholmod_free_sparse(&A, &c);
        cholmod_finish(&c);
        return;
    }
    cholmod_dense *B = cholmod_zeros(M->n, nrhs, CHOLMOD_REAL, &c);
    for (int64_t i = 0; i < (int64_t)M->n * nrhs; ++i) ((double*)B->x)[i] = 1.0;
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
            report_residual_if_bad(st, residual_ones_rhs(M, (const double*)X->x,
                                                        (int64_t)X->d, nrhs));
            checked_residual = true;
            st.ResumeTiming();
        }
        benchmark::DoNotOptimize(X);
        cholmod_free_dense(&X, &c);
    }
    st.SetItemsProcessed(int64_t(st.iterations()) * M->n * nrhs);
    cholmod_free_dense(&B, &c);
    cholmod_free_factor(&L, &c);
    cholmod_free_sparse(&A, &c);
    cholmod_finish(&c);
}
#endif

static void run_ours_e2e(benchmark::State &st, const char *name, int32_t nrhs);
#ifdef LDLT_HAVE_CHOLMOD
static void run_chm_e2e(benchmark::State &st, const char *name, int32_t nrhs);
#endif

namespace {
enum class fixture_class {
    spd,
    indef,
};

struct curated_fixture {
    const char *label;
    const char *file_stem;
    fixture_class klass;
    bool register_solve;
    bool register_e2e;
};

constexpr curated_fixture kCuratedFixtures[] = {
    {"bcsstk01", "bcsstk01", fixture_class::spd, false, false},
    {"bcsstk14", "bcsstk14", fixture_class::spd, false, false},
    {"msc00726", "msc00726", fixture_class::spd, false, false},
    {"pwtk", "pwtk", fixture_class::spd, false, true},
    {"ct20stif", "ct20stif", fixture_class::spd, false, false},
    {"Dubcova2", "Dubcova2", fixture_class::spd, true, true},
    {"bmwcra_1", "bmwcra_1", fixture_class::spd, false, false},
    {"inline_1", "inline_1", fixture_class::spd, false, false},
    {"bundle_adj", "bundle_adj", fixture_class::spd, false, false},
    {"oilpan", "oilpan", fixture_class::spd, false, true},
    {"pde2961", "pde2961", fixture_class::indef, false, false},
    {"dynamicSoaringProblem_4", "dynamicSoaringProblem_4", fixture_class::indef, false, false},
    {"c_27", "c-27", fixture_class::indef, false, false},
    {"analytics", "analytics", fixture_class::indef, false, false},
    {"c_big", "c-big", fixture_class::indef, false, false},
    {"dawson5", "dawson5", fixture_class::indef, false, false},
    {"bmw3_2", "bmw3_2", fixture_class::indef, false, false},
};

static void register_curated_fixture(const curated_fixture &spec, ldlt_bench_fixture_mode mode) {
    const std::string full = std::string(spec.file_stem) + ".mtx";
    const bool use_pivot = mode == ldlt_bench_fixture_mode::pivot;
    auto register_ours = [&](const std::string &name) {
        benchmark::RegisterBenchmark(name.c_str(),
            [full, use_pivot](benchmark::State &st) {
                if (use_pivot) run_ours_indef(st, full.c_str());
                else run_ours(st, full.c_str());
            })->Unit(benchmark::kMillisecond);
    };

    if (spec.klass == fixture_class::spd) {
        register_ours("Ours_" + std::string(spec.label));
#ifdef LDLT_HAVE_CHOLMOD
        benchmark::RegisterBenchmark(("CHM_Super_" + std::string(spec.label)).c_str(),
            [full](benchmark::State &st) { run_chm(st, full.c_str(), 1); })
            ->Unit(benchmark::kMillisecond);
#endif
        if (spec.register_e2e) {
            benchmark::RegisterBenchmark(("Ours_E2E_" + std::string(spec.label) + "_1").c_str(),
                [full](benchmark::State &st) { run_ours_e2e(st, full.c_str(), 1); })
                ->Unit(benchmark::kMillisecond);
#ifdef LDLT_HAVE_CHOLMOD
            benchmark::RegisterBenchmark(("CHM_E2E_" + std::string(spec.label) + "_1").c_str(),
                [full](benchmark::State &st) { run_chm_e2e(st, full.c_str(), 1); })
                ->Unit(benchmark::kMillisecond);
#endif
        }
        if (spec.register_solve) {
            benchmark::RegisterBenchmark(("Ours_Solve_" + std::string(spec.label) + "_1").c_str(),
                [full](benchmark::State &st) { run_ours_solve(st, full.c_str(), 1); })
                ->Unit(benchmark::kMillisecond);
            benchmark::RegisterBenchmark(("Ours_Solve_" + std::string(spec.label) + "_8").c_str(),
                [full](benchmark::State &st) { run_ours_solve(st, full.c_str(), 8); })
                ->Unit(benchmark::kMillisecond);
#ifdef LDLT_HAVE_CHOLMOD
            benchmark::RegisterBenchmark(("CHM_Super_Solve_" + std::string(spec.label) + "_1").c_str(),
                [full](benchmark::State &st) { run_chm_solve(st, full.c_str(), 1); })
                ->Unit(benchmark::kMillisecond);
            benchmark::RegisterBenchmark(("CHM_Super_Solve_" + std::string(spec.label) + "_8").c_str(),
                [full](benchmark::State &st) { run_chm_solve(st, full.c_str(), 8); })
                ->Unit(benchmark::kMillisecond);
#endif
        }
    } else {
        register_ours("Ours_Indef_" + std::string(spec.label));
#ifdef LDLT_HAVE_CHOLMOD
        benchmark::RegisterBenchmark(("CHM_Simpl_Indef_" + std::string(spec.label)).c_str(),
            [full](benchmark::State &st) { run_chm(st, full.c_str(), 0); })
            ->Unit(benchmark::kMillisecond);
#endif
    }
}
}  // namespace

void ldlt_register_curated_fixture_benchmarks(ldlt_bench_fixture_mode mode) {
    for (const curated_fixture &spec : kCuratedFixtures) {
        register_curated_fixture(spec, mode);
    }
}

/* ---------- End-to-end: analyze + factorize + solve, single iteration ---------- */
static void run_ours_e2e(benchmark::State &st, const char *name, int32_t nrhs) {
    const Mtx *M = load(name);
    if (!M) { st.SkipWithError("fixture missing"); return; }
    std::vector<double> B((size_t)M->n * nrhs, 1.0);
    bool checked_residual = false;
    double analyze_ms = 0.0, factor_ms = 0.0, solve_ms = 0.0;
    int64_t samples = 0;
    for (auto _ : st) {
        ldlt_options opt; ldlt_options_default(&opt);
        opt.ordering = LDLT_ORDER_METIS;
        ldlt_symbolic *S = nullptr;
        auto t0 = std::chrono::steady_clock::now();
        ldlt_status rc = ldlt_analyze(M->n, M->Ap, M->Ai, &opt, &S);
        auto t1 = std::chrono::steady_clock::now();
        if (rc != LDLT_OK) { st.SkipWithError("analyze failed"); break; }
        ldlt_numeric *N = nullptr;
        auto t2 = std::chrono::steady_clock::now();
        rc = ldlt_factorize(S, M->Ap, M->Ai, M->Ax, &opt, &N);
        auto t3 = std::chrono::steady_clock::now();
        if (rc != LDLT_OK) { st.SkipWithError(ldlt_status_str(rc)); ldlt_free_symbolic(S); break; }
        std::vector<double> X = B;
        auto t4 = std::chrono::steady_clock::now();
        rc = ldlt_solve(N, nrhs, X.data(), M->n);
        auto t5 = std::chrono::steady_clock::now();
        if (rc != LDLT_OK) { st.SkipWithError(ldlt_status_str(rc)); ldlt_free_numeric(N); ldlt_free_symbolic(S); break; }
        analyze_ms += elapsed_ms(t0, t1);
        factor_ms += elapsed_ms(t2, t3);
        solve_ms += elapsed_ms(t4, t5);
        samples++;
        if (!checked_residual) {
            st.PauseTiming();
            report_residual_if_bad(st, residual_ones_rhs(M, X.data(), M->n, nrhs));
            checked_residual = true;
            st.ResumeTiming();
        }
        benchmark::DoNotOptimize(X.data());
        ldlt_free_numeric(N);
        ldlt_free_symbolic(S);
    }
    st.SetItemsProcessed(int64_t(st.iterations()) * M->n);
    if (samples > 0) {
        st.counters["analyze_ms"] = benchmark::Counter(analyze_ms / samples);
        st.counters["factor_ms"] = benchmark::Counter(factor_ms / samples);
        st.counters["solve_ms"] = benchmark::Counter(solve_ms / samples);
    }
}

#ifdef LDLT_HAVE_CHOLMOD
static void run_chm_e2e(benchmark::State &st, const char *name, int32_t nrhs) {
    const Mtx *M = load(name);
    if (!M) { st.SkipWithError("fixture missing"); return; }
    bool checked_residual = false;
    double analyze_ms = 0.0, factor_ms = 0.0, solve_ms = 0.0;
    int64_t samples = 0;
    for (auto _ : st) {
        cholmod_common c; cholmod_start(&c);
        c.final_ll = 1;
        c.supernodal = CHOLMOD_SUPERNODAL;
        chm_force_metis(&c);
        cholmod_sparse *A = to_chm(&c, M);
        auto t0 = std::chrono::steady_clock::now();
        cholmod_factor *L = cholmod_analyze(A, &c);
        auto t1 = std::chrono::steady_clock::now();
        auto t2 = std::chrono::steady_clock::now();
        int ok = cholmod_factorize(A, L, &c);
        auto t3 = std::chrono::steady_clock::now();
        if (!ok || !chm_factor_ok(&c, L, M->n)) {
            st.SkipWithError("CHOLMOD factorization failed");
            cholmod_free_factor(&L, &c);
            cholmod_free_sparse(&A, &c);
            cholmod_finish(&c);
            break;
        }
        cholmod_dense *Bd = cholmod_zeros(M->n, nrhs, CHOLMOD_REAL, &c);
        for (int64_t i = 0; i < (int64_t)M->n * nrhs; ++i) ((double*)Bd->x)[i] = 1.0;
        auto t4 = std::chrono::steady_clock::now();
        cholmod_dense *X = cholmod_solve(CHOLMOD_A, L, Bd, &c);
        auto t5 = std::chrono::steady_clock::now();
        if (!X || c.status != CHOLMOD_OK) {
            st.SkipWithError("CHOLMOD solve failed");
            if (X) cholmod_free_dense(&X, &c);
            cholmod_free_dense(&Bd, &c);
            cholmod_free_factor(&L, &c);
            cholmod_free_sparse(&A, &c);
            cholmod_finish(&c);
            break;
        }
        analyze_ms += elapsed_ms(t0, t1);
        factor_ms += elapsed_ms(t2, t3);
        solve_ms += elapsed_ms(t4, t5);
        samples++;
        if (!checked_residual) {
            st.PauseTiming();
            report_residual_if_bad(st, residual_ones_rhs(M, (const double*)X->x,
                                                        (int64_t)X->d, nrhs));
            checked_residual = true;
            st.ResumeTiming();
        }
        benchmark::DoNotOptimize(X);
        cholmod_free_dense(&X, &c);
        cholmod_free_dense(&Bd, &c);
        cholmod_free_factor(&L, &c);
        cholmod_free_sparse(&A, &c);
        cholmod_finish(&c);
    }
    st.SetItemsProcessed(int64_t(st.iterations()) * M->n);
    if (samples > 0) {
        st.counters["analyze_ms"] = benchmark::Counter(analyze_ms / samples);
        st.counters["factor_ms"] = benchmark::Counter(factor_ms / samples);
        st.counters["solve_ms"] = benchmark::Counter(solve_ms / samples);
    }
}
#endif

#define REGISTER_E2E(NAME, NRHS) \
    static void Ours_E2E_##NAME##_##NRHS(benchmark::State &st) \
        { run_ours_e2e(st, #NAME ".mtx", NRHS); } \
    BENCHMARK(Ours_E2E_##NAME##_##NRHS)->Unit(benchmark::kMillisecond);
#ifdef LDLT_HAVE_CHOLMOD
#define REGISTER_E2E_CHM(NAME, NRHS) \
    static void CHM_E2E_##NAME##_##NRHS(benchmark::State &st) \
        { run_chm_e2e(st, #NAME ".mtx", NRHS); } \
    BENCHMARK(CHM_E2E_##NAME##_##NRHS)->Unit(benchmark::kMillisecond);
#else
#define REGISTER_E2E_CHM(NAME, NRHS)
#endif

REGISTER_E2E(Dubcova2, 1)
REGISTER_E2E_CHM(Dubcova2, 1)
REGISTER_E2E(oilpan, 1)
REGISTER_E2E_CHM(oilpan, 1)
REGISTER_E2E(pwtk, 1)
REGISTER_E2E_CHM(pwtk, 1)

static void Ours_Solve_Dubcova2_1(benchmark::State &st) { run_ours_solve(st, "Dubcova2.mtx", 1); }
static void Ours_Solve_Dubcova2_8(benchmark::State &st) { run_ours_solve(st, "Dubcova2.mtx", 8); }
BENCHMARK(Ours_Solve_Dubcova2_1)->Unit(benchmark::kMillisecond);
BENCHMARK(Ours_Solve_Dubcova2_8)->Unit(benchmark::kMillisecond);
#ifdef LDLT_HAVE_CHOLMOD
static void CHM_Super_Solve_Dubcova2_1(benchmark::State &st) { run_chm_solve(st, "Dubcova2.mtx", 1); }
static void CHM_Super_Solve_Dubcova2_8(benchmark::State &st) { run_chm_solve(st, "Dubcova2.mtx", 8); }
BENCHMARK(CHM_Super_Solve_Dubcova2_1)->Unit(benchmark::kMillisecond);
BENCHMARK(CHM_Super_Solve_Dubcova2_8)->Unit(benchmark::kMillisecond);
#endif

static std::string bench_safe_name(const std::filesystem::path &path) {
    std::string s = path.stem().string();
    if (s.empty()) s = path.filename().string();
    for (char &c : s) {
        if (!(std::isalnum((unsigned char)c) || c == '_' || c == '-')) c = '_';
    }
    return s;
}

void ldlt_register_matrix_dir_benchmarks(const std::string &dir) {
    std::filesystem::path root(dir);
    std::error_code ec;
    if (!std::filesystem::is_directory(root, ec)) {
        std::cerr << "ldlt_bench: --matrix_dir is not a directory: " << dir << "\n";
        return;
    }

    std::vector<std::filesystem::path> paths;
    for (const auto &entry : std::filesystem::recursive_directory_iterator(root, ec)) {
        if (ec) break;
        if (!entry.is_regular_file()) continue;
        std::filesystem::path p = entry.path();
        if (p.extension() == ".mtx") paths.push_back(std::filesystem::absolute(p));
    }
    std::sort(paths.begin(), paths.end());

    for (const auto &path : paths) {
        std::string full = path.string();
        std::string name = bench_safe_name(path);
        benchmark::RegisterBenchmark(("Ours_User_Factor_" + name).c_str(),
            [full](benchmark::State &st) { run_ours(st, full.c_str()); })
            ->Unit(benchmark::kMillisecond);
        benchmark::RegisterBenchmark(("Ours_User_E2E_" + name).c_str(),
            [full](benchmark::State &st) { run_ours_e2e(st, full.c_str(), 1); })
            ->Unit(benchmark::kMillisecond);
#ifdef LDLT_HAVE_CHOLMOD
        benchmark::RegisterBenchmark(("CHM_Super_User_Factor_" + name).c_str(),
            [full](benchmark::State &st) { run_chm(st, full.c_str(), 1); })
            ->Unit(benchmark::kMillisecond);
        benchmark::RegisterBenchmark(("CHM_Super_User_E2E_" + name).c_str(),
            [full](benchmark::State &st) { run_chm_e2e(st, full.c_str(), 1); })
            ->Unit(benchmark::kMillisecond);
#endif
    }
    if (paths.empty()) {
        std::cerr << "ldlt_bench: no .mtx files found under " << dir << "\n";
    }
}
