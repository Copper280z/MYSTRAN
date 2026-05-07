project('ldlt', 'c', 'cpp',
  version: '0.1.0',
  default_options: [
    'c_std=c11',
    'cpp_std=c++20',
    'warning_level=2',
    'b_ndebug=if-release',
  ],
  meson_version: '>=1.1.0',
)

cc = meson.get_compiler('c')
posix_defs = []
if host_machine.system() != 'windows'
  posix_defs += '-DLDLT_USE_POSIX_SIGNALS=1'
  if host_machine.system() != 'darwin'
    add_project_arguments('-D_POSIX_C_SOURCE=200809L', language: 'c')
  endif
endif

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
blas_deps = []
use_system_blas = false
if blas_opt != 'reference'
  # Prefer Accelerate on macOS, then generic blas/cblas.
  if host_machine.system() == 'darwin'
    blas_dep = dependency('Accelerate', required: false)
    use_system_blas = blas_dep.found()
  endif
  foreach blas_name : ['cblas', 'openblas', 'blas']
    if not use_system_blas
      candidate_dep = dependency(blas_name, required: false)
      if candidate_dep.found() and cc.has_header('cblas.h', dependencies: candidate_dep) and cc.has_function('cblas_dgemm', prefix: '#include <cblas.h>', dependencies: candidate_dep)
        blas_dep = candidate_dep
        use_system_blas = true
      endif
    endif
  endforeach
  if use_system_blas
    blas_deps += blas_dep
    if host_machine.system() != 'darwin'
      lapack_dep = dependency('lapack', required: false)
      if lapack_dep.found()
        blas_deps += lapack_dep
      endif
    endif
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
ldlt_defs = blas_defs + ord_defs + omp_defs + posix_defs
ldlt_deps = blas_deps + [amd_dep, metis_dep, omp_dep]
ldlt_c_args = ldlt_defs
if use_system_blas and host_machine.system() == 'darwin'
  ldlt_c_args += '-Wno-deprecated-declarations'
endif

ldlt_inc = include_directories('include', 'src')

ldlt_sources = files(
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
if not use_system_blas
  ldlt_sources += files('src/blas_ref.c')
endif

ldlt_lib = static_library('ldlt',
  ldlt_sources,
  include_directories: ldlt_inc,
  c_args: ldlt_c_args,
  dependencies: ldlt_deps,
  install: true,
)

ldlt_dep = declare_dependency(
  link_with: ldlt_lib,
  include_directories: ldlt_inc,
  compile_args: ldlt_defs,
  dependencies: ldlt_deps,
)

install_headers('include/ldlt/ldlt.h', subdir: 'ldlt')

if get_option('with_fortran_interface')
  ldlt_fortran_lib = static_library('ldlt_fortran',
    files('fortran/c_fortran_ldlt.c'),
    include_directories: ldlt_inc,
    c_args: ldlt_c_args,
    install: true,
  )

  ldlt_fortran_dep = declare_dependency(
    link_with: [ldlt_fortran_lib, ldlt_lib],
    include_directories: ldlt_inc,
    compile_args: ldlt_defs,
    dependencies: ldlt_deps,
  )
endif

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
