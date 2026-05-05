#pragma once

#include <string>

enum class ldlt_bench_fixture_mode {
    spd,
    pivot,
};

void ldlt_register_curated_fixture_benchmarks(ldlt_bench_fixture_mode mode);
void ldlt_register_matrix_dir_benchmarks(const std::string &dir);
