#include <benchmark/benchmark.h>

#include <string>
#include <vector>

void ldlt_register_matrix_dir_benchmarks(const std::string &dir);

int main(int argc, char **argv) {
    std::vector<char*> bench_argv;
    bench_argv.reserve((size_t)argc);
    bench_argv.push_back(argv[0]);

    std::vector<std::string> matrix_dirs;
    const std::string flag = "--matrix_dir=";
    for (int i = 1; i < argc; ++i) {
        std::string arg(argv[i]);
        if (arg.rfind(flag, 0) == 0) {
            matrix_dirs.push_back(arg.substr(flag.size()));
        } else if (arg == "--matrix_dir" && i + 1 < argc) {
            matrix_dirs.push_back(argv[++i]);
        } else {
            bench_argv.push_back(argv[i]);
        }
    }

    int bench_argc = (int)bench_argv.size();
    benchmark::Initialize(&bench_argc, bench_argv.data());
    for (const std::string &dir : matrix_dirs)
        ldlt_register_matrix_dir_benchmarks(dir);
    benchmark::RunSpecifiedBenchmarks();
    benchmark::Shutdown();
    return 0;
}
