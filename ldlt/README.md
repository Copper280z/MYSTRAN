# ldlt — supernodal sparse LDLᵀ solver (C23)

Three-phase API: `ldlt_analyze` → `ldlt_factorize` → `ldlt_solve`. The default path is a
parallel left-looking supernodal factor via OpenMP tasks. Indefinite matrices can use
`LDLT_PIVOT_BUNCH_KAUFMAN`, which enables dynamic 1x1/2x2 symmetric pivoting with delayed
pivots.

## Build

```sh
meson setup build
meson compile -C build
meson test -C build
build/bench/ldlt_bench_spd
build/bench/ldlt_bench_pivot
```

Options (`meson configure build -Dkey=value`):
- `blas=auto|system|reference` — default auto (Accelerate on macOS, then cblas/openblas/blas; reference fallback).
- `ordering=auto|amd|metis|internal` — picks SuiteSparse AMD / METIS if available, else internal stub (identity).
- `with_openmp=true|false`.
- `large_matrices=true` — enable heavier SuiteSparse fixtures.

## Test fixtures

`python3 matrices/fetch_suitesparse.py` downloads .mtx files listed in `matrices/manifest.txt`.
Set `LDLT_FIXTURE_DIR=$PWD/matrices` so the parameterized SuiteSparse tests pick them up.
