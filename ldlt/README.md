# ldlt - supernodal sparse LDL^T solver (C11)

Three-phase API: `ldlt_analyze` -> `ldlt_factorize` -> `ldlt_solve`. The default path is a
parallel left-looking supernodal factor via OpenMP tasks. Indefinite matrices can use
`LDLT_PIVOT_BUNCH_KAUFMAN`, which enables dynamic 1x1/2x2 symmetric pivoting with delayed
pivots.

## Meson Build

```sh
meson setup build
meson compile -C build
```

Tests and benchmarks are opt-in:

```sh
meson setup build-dev -Dwith_tests=true -Dwith_bench=true
meson compile -C build-dev
meson test -C build-dev
build-dev/bench/ldlt_bench_spd
build-dev/bench/ldlt_bench_pivot
```

Options (`meson configure build -Dkey=value`):
- `blas=auto|system|reference` - default auto (Accelerate on macOS, then cblas/openblas/blas; reference fallback).
- `ordering=auto|amd|metis|internal` - picks SuiteSparse AMD / METIS if available, else internal stub (identity).
- `with_openmp=true|false`.
- `with_fortran_interface=true|false` - builds `ldlt_fortran`, the QDLDL-compatible Fortran bridge.
- `with_tests=true|false`.
- `with_bench=true|false`.
- `large_matrices=true` - enable heavier SuiteSparse fixtures.

## CMake Build

```sh
cmake -S . -B build-cmake -DLDLT_BLAS=AUTO -DLDLT_ORDERING=AUTO
cmake --build build-cmake
```

Useful options:
- `LDLT_BLAS=AUTO|SYSTEM|REFERENCE`.
- `LDLT_ORDERING=AUTO|AMD|METIS|INTERNAL`.
- `LDLT_ENABLE_OPENMP=ON|OFF`.
- `LDLT_BUILD_FORTRAN_INTERFACE=ON|OFF`.
- `LDLT_BUILD_TESTS=ON|OFF`.
- `LDLT_BUILD_BENCH=ON|OFF`.

## Test fixtures

`python3 matrices/fetch_suitesparse.py` downloads .mtx files listed in `matrices/manifest.txt`.
Set `LDLT_FIXTURE_DIR=$PWD/matrices` so the parameterized SuiteSparse tests pick them up.
