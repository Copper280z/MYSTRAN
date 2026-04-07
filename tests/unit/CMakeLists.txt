# Unit tests for MYSTRAN utility subroutines.
#
# Each test is a standalone Fortran program that compiles only the module(s) it
# exercises plus minimal stubs for MYSTRAN's global module dependencies.  The
# stubs live in tests/unit/stubs/ and provide just enough symbols to satisfy the
# compiler/linker without pulling in the full MYSTRAN module tree.
#
# All tests are registered with CTest and labelled "unit" so they can be run in
# isolation with:
#   ctest --test-dir <build_dir> -L unit -V

cmake_minimum_required(VERSION 3.18)

set(MYSTRAN_SRC "${PROJECT_SOURCE_DIR}/Source")
set(UNIT_STUBS  "${CMAKE_CURRENT_SOURCE_DIR}/stubs")

# Use a separate Fortran module output directory so stub modules (SCONTR, IOUNT1,
# SUBR_BEGEND_LEVELS) don't collide with the real modules compiled for the main
# MYSTRAN binary.
set(UNIT_TEST_MOD_DIR "${CMAKE_BINARY_DIR}/tests/unit/mod")

# ============================================================
# Shared support library: stub modules + real low-level modules
# that the subroutines under test USE.
# ============================================================
add_library(unit_test_support STATIC
  # Real modules (no MYSTRAN-internal dependencies)
  ${MYSTRAN_SRC}/Modules/PENTIUM_II_KIND.f90
  ${MYSTRAN_SRC}/Modules/TIMDAT.f90
  ${MYSTRAN_SRC}/Modules/CONSTANTS_1.f90

  # Stub replacements for large global-state modules
  ${UNIT_STUBS}/stub_scontr.f90
  ${UNIT_STUBS}/stub_iount1.f90
  ${UNIT_STUBS}/stub_subr_begend_levels.f90

  # Real interface/USE_IF modules (only pull in OURTIM_Interface)
  ${MYSTRAN_SRC}/Interfaces/OURTIM_Interface.f90
  ${MYSTRAN_SRC}/USE_IFs/CROSS_USE_IFs.f90
  ${MYSTRAN_SRC}/USE_IFs/MATMULT_FFF_USE_IFs.f90
  ${MYSTRAN_SRC}/USE_IFs/MATMULT_FFF_T_USE_IFs.f90

  # No-op OURTIM subroutine — satisfies the linker; never called in unit tests
  # because WRT_LOG = 0 in the stub
  ${UNIT_STUBS}/stub_ourtim.f90
)

set_target_properties(unit_test_support PROPERTIES
  Fortran_MODULE_DIRECTORY ${UNIT_TEST_MOD_DIR}
)

target_include_directories(unit_test_support PUBLIC
  ${MYSTRAN_SRC}/INCLUDE
)

# ============================================================
# Subroutines under test — compiled against the stub support lib
# ============================================================
add_library(subroutines_under_test STATIC
  ${MYSTRAN_SRC}/UTIL/OUTER_PRODUCT.f90
  ${MYSTRAN_SRC}/UTIL/CROSS.f90
  ${MYSTRAN_SRC}/UTIL/MATMULT_FFF.f90
  ${MYSTRAN_SRC}/UTIL/MATMULT_FFF_T.f90
)

set_target_properties(subroutines_under_test PROPERTIES
  Fortran_MODULE_DIRECTORY ${UNIT_TEST_MOD_DIR}
)

target_link_libraries(subroutines_under_test PUBLIC unit_test_support)

# ============================================================
# Helper macro to register a unit test
# ============================================================
macro(add_unit_test test_name source_file)
  add_executable(${test_name} ${source_file})
  set_target_properties(${test_name} PROPERTIES
    Fortran_MODULE_DIRECTORY ${UNIT_TEST_MOD_DIR}
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests/unit"
  )
  target_link_libraries(${test_name} PRIVATE subroutines_under_test unit_test_support)
  add_test(NAME unit/${test_name} COMMAND ${test_name})
  set_tests_properties(unit/${test_name} PROPERTIES LABELS "unit")
endmacro()

# ============================================================
# Register tests
# ============================================================
add_unit_test(test_outer_product test_outer_product.f90)
add_unit_test(test_cross_product test_cross_product.f90)
add_unit_test(test_matmult       test_matmult.f90)
