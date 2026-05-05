#include "bench_fixtures.h"

namespace {
struct AutoRegisterPivot {
    AutoRegisterPivot() { ldlt_register_curated_fixture_benchmarks(ldlt_bench_fixture_mode::pivot); }
};
AutoRegisterPivot g_auto_register_pivot;
}  // namespace
