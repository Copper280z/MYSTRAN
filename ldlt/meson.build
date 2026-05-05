project('ldlt', 'c', 'cpp',
  version: '0.1.0',
  default_options: [
    'c_std=c2x',
    'cpp_std=c++20',
    'warning_level=2',
    'b_ndebug=if-release',
  ],
  meson_version: '>=1.1.0',
)

cc = meson.get_compiler('c')
add_project_arguments('-D_POSIX_C_SOURCE=200809L', language: 'c')

# Tuning flags applied to all targets (release builds benefit most).
tuning_args = []
foreach a : ['-march=native', '-funroll-loops']
  if cc.has_argument(a)
    tuning_args += a
  endif
endforeach
if tuning_args.length() > 0
  add_project_arguments(tuning_args, language: 'c')
  add_project_arguments(tuning_args, language: 'cpp')
endif

# ---- BLAS ----
blas_opt = get_option('blas')
blas_dep = dependency('', required: false)
use_system_blas = false
if blas_opt != 'reference'
  # Prefer Accelerate on macOS, then generic blas/cblas.
  if host_machine.system() == 'darwin'
    blas_dep = dependency('Accelerate', required: false)
  endif
  if not blas_dep.found()
    blas_dep = dependency('cblas', required: false)
  endif
  if not blas_dep.found()
    blas_dep = dependency('openblas', required: false)
  endif
  if not blas_dep.found()
    blas_dep = dependency('blas', required: false)
  endif
  if blas_dep.found()
    use_system_blas = true
  elif blas_opt == 'system'
    error('blas=system but no system BLAS found')
  endif
endif

blas_defs = []
if use_system_blas
  blas_defs += '-DLDLT_USE_SYSTEM_BLAS=1'
  if host_machine.system() == 'darwin'
    blas_defs += '-DLDLT_BLAS_ACCELERATE=1'
  endif
else
  blas_defs += '-DLDLT_USE_SYSTEM_BLAS=0'
  message('Using internal reference BLAS (slow). Pass -Dblas=system to require system BLAS.')
endif

# ---- Ordering ----
ord_opt = get_option('ordering')
amd_dep = dependency('', required: false)
metis_dep = dependency('', required: false)
if ord_opt != 'internal'
  amd_dep = dependency('amd', required: false)
  metis_dep = dependency('metis', required: false)
  if not metis_dep.found()
    metis_lib = cc.find_library('metis', dirs: ['/opt/homebrew/lib', '/usr/local/lib'], required: false)
    metis_inc = include_directories('/opt/homebrew/include', is_system: true)
    if metis_lib.found() and cc.has_header('metis.h', args: ['-I/opt/homebrew/include'])
      metis_dep = declare_dependency(dependencies: metis_lib, include_directories: metis_inc)
    endif
  endif
  if ord_opt == 'amd' and not amd_dep.found()
    error('ordering=amd but SuiteSparse AMD not found')
  endif
  if ord_opt == 'metis' and not metis_dep.found()
    error('ordering=metis but METIS not found')
  endif
endif

ord_defs = []
if amd_dep.found()
  ord_defs += '-DLDLT_HAVE_AMD=1'
endif
if metis_dep.found()
  ord_defs += '-DLDLT_HAVE_METIS=1'
endif

# ---- OpenMP ----
omp_dep = dependency('', required: false)
if get_option('with_openmp')
  omp_dep = dependency('openmp', required: false)
endif
omp_defs = omp_dep.found() ? ['-DLDLT_HAVE_OPENMP=1'] : []

ldlt_inc = include_directories('include', 'src')

ldlt_sources = files(
  'src/blas_ref.c',
  'src/util.c',
  'src/perm.c',
  'src/etree.c',
  'src/symbolic.c',
  'src/ordering.c',
  'src/amd_stub.c',
  'src/factor.c',
  'src/factor_kernels.c',
  'src/solve.c',
  'src/io_mm.c',
)

ldlt_lib = static_library('ldlt',
  ldlt_sources,
  include_directories: ldlt_inc,
  c_args: blas_defs + ord_defs + omp_defs,
  dependencies: [blas_dep, amd_dep, metis_dep, omp_dep],
  install: true,
)

ldlt_dep = declare_dependency(
  link_with: ldlt_lib,
  include_directories: ldlt_inc,
  compile_args: blas_defs + ord_defs + omp_defs,
  dependencies: [blas_dep, amd_dep, metis_dep, omp_dep],
)

install_headers('include/ldlt/ldlt.h', subdir: 'ldlt')

if get_option('with_tests')
  subdir('tests')
endif
if get_option('with_bench')
  subdir('bench')
endif

summary({
  'BLAS': use_system_blas ? 'system' : 'reference',
  'AMD':  amd_dep.found(),
  'METIS': metis_dep.found(),
  'OpenMP': omp_dep.found(),
}, section: 'ldlt')
