gtest_dep = dependency('gtest', main: true, required: false)
if not gtest_dep.found()
  gtest_proj = subproject('gtest', required: false)
  if gtest_proj.found()
    gtest_dep = gtest_proj.get_variable('gtest_main_dep')
  endif
endif

if not gtest_dep.found()
  warning('GoogleTest not available; tests disabled')
  subdir_done()
endif

test_exe = executable('ldlt_tests',
  files('test_main.cc', 'test_basic.cc', 'test_solve.cc', 'test_matrices.cc'),
  dependencies: [ldlt_dep, gtest_dep],
)
test('ldlt-tests', test_exe, timeout: 120)
