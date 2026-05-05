benchmark_dep = dependency('benchmark', required: false)
if not benchmark_dep.found()
  bench_proj = subproject('google-benchmark', required: false)
  if bench_proj.found()
    benchmark_dep = bench_proj.get_variable('google_benchmark_dep')
  endif
endif
if not benchmark_dep.found()
  warning('google-benchmark not available; bench disabled')
  subdir_done()
endif

cholmod_dep = dependency('cholmod', required: false)
bench_args = []
bench_deps = [ldlt_dep, benchmark_dep]
if cholmod_dep.found()
  bench_args += '-DLDLT_HAVE_CHOLMOD=1'
  bench_deps += cholmod_dep
endif

executable('ldlt_bench_pivot',
  files('bench_main.cc', 'bench_factor.cc', 'bench_solve.cc', 'bench_cholmod.cc', 'bench_fixtures.cc', 'bench_fixtures_pivot.cc'),
  dependencies: bench_deps,
  cpp_args: bench_args,
)

executable('ldlt_bench_spd',
  files('bench_main.cc', 'bench_factor.cc', 'bench_solve.cc', 'bench_cholmod.cc', 'bench_fixtures.cc', 'bench_fixtures_spd.cc'),
  dependencies: bench_deps,
  cpp_args: bench_args,
)
