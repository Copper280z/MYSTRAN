#include "bench_fixtures.h"

namespace {
struct AutoRegisterSpd {
    AutoRegisterSpd() { ldlt_register_curated_fixture_benchmarks(ldlt_bench_fixture_mode::spd); }
};
AutoRegisterSpd g_auto_register_spd;
}  // namespace
