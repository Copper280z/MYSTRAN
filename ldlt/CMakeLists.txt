cmake_minimum_required(VERSION 3.18)

project(ldlt
  VERSION 0.1.0
  LANGUAGES C CXX
)

include(GNUInstallDirs)
include(CheckFunctionExists)

option(LDLT_BUILD_TESTS "Build LDLT tests" OFF)
option(LDLT_BUILD_BENCH "Build LDLT benchmarks" OFF)
option(LDLT_BUILD_FORTRAN_INTERFACE "Build the Fortran-callable LDLT bridge" ON)
option(LDLT_ENABLE_OPENMP "Enable OpenMP when available" ON)

set(LDLT_BLAS "AUTO" CACHE STRING "BLAS backend: AUTO, SYSTEM, or REFERENCE")
set_property(CACHE LDLT_BLAS PROPERTY STRINGS AUTO SYSTEM REFERENCE)
set(LDLT_ORDERING "AUTO" CACHE STRING "Ordering backend: AUTO, AMD, METIS, or INTERNAL")
set_property(CACHE LDLT_ORDERING PROPERTY STRINGS AUTO AMD METIS INTERNAL)

set(ldlt_core_sources
  src/util.c
  src/perm.c
  src/etree.c
  src/symbolic.c
  src/ordering.c
  src/amd_stub.c
  src/factor.c
  src/factor_kernels.c
  src/solve.c
  src/io_mm.c
)

set(ldlt_use_system_blas OFF)
set(ldlt_blas_libraries)
set(ldlt_blas_include_dirs)
set(ldlt_blas_definitions LDLT_USE_SYSTEM_BLAS=0)

string(TOUPPER "${LDLT_BLAS}" LDLT_BLAS_UPPER)
if(NOT LDLT_BLAS_UPPER MATCHES "^(AUTO|SYSTEM|REFERENCE)$")
  message(FATAL_ERROR "LDLT_BLAS must be AUTO, SYSTEM, or REFERENCE")
endif()

if(NOT LDLT_BLAS_UPPER STREQUAL "REFERENCE")
  if(APPLE)
    find_library(LDLT_ACCELERATE_FRAMEWORK Accelerate)
    if(LDLT_ACCELERATE_FRAMEWORK)
      set(ldlt_use_system_blas ON)
      set(ldlt_blas_libraries "${LDLT_ACCELERATE_FRAMEWORK}")
      set(ldlt_blas_definitions LDLT_USE_SYSTEM_BLAS=1 LDLT_BLAS_ACCELERATE=1)
    endif()
  endif()

  if(NOT ldlt_use_system_blas)
    find_path(LDLT_CBLAS_INCLUDE_DIR NAMES cblas.h
      HINTS /usr/include /usr/local/include
      PATH_SUFFIXES x86_64-linux-gnu openblas openblas-pthread openblas-openmp)
    find_library(LDLT_CBLAS_LIBRARY NAMES openblas cblas
      HINTS /usr/lib /usr/local/lib
      PATH_SUFFIXES x86_64-linux-gnu openblas-pthread openblas-openmp)
    find_package(LAPACK QUIET)
    if(LDLT_CBLAS_INCLUDE_DIR AND LDLT_CBLAS_LIBRARY)
      set(_ldlt_old_required_includes "${CMAKE_REQUIRED_INCLUDES}")
      set(_ldlt_old_required_libraries "${CMAKE_REQUIRED_LIBRARIES}")
      set(CMAKE_REQUIRED_INCLUDES "${LDLT_CBLAS_INCLUDE_DIR}")
      set(CMAKE_REQUIRED_LIBRARIES "${LDLT_CBLAS_LIBRARY}")
      check_function_exists(cblas_dgemm LDLT_HAVE_CBLAS_DGEMM)
      set(CMAKE_REQUIRED_INCLUDES "${_ldlt_old_required_includes}")
      set(CMAKE_REQUIRED_LIBRARIES "${_ldlt_old_required_libraries}")
      unset(_ldlt_old_required_includes)
      unset(_ldlt_old_required_libraries)
      if(LDLT_HAVE_CBLAS_DGEMM)
        set(ldlt_use_system_blas ON)
        set(ldlt_blas_include_dirs "${LDLT_CBLAS_INCLUDE_DIR}")
        set(ldlt_blas_libraries "${LDLT_CBLAS_LIBRARY}")
        if(LAPACK_FOUND)
          list(APPEND ldlt_blas_libraries ${LAPACK_LIBRARIES})
        endif()
        set(ldlt_blas_definitions LDLT_USE_SYSTEM_BLAS=1)
        if(LDLT_CBLAS_LIBRARY MATCHES "[Oo]pen[Bb][Ll][Aa][Ss]")
          list(APPEND ldlt_blas_definitions LDLT_HAVE_OPENBLAS=1)
        endif()
      endif()
    endif()
  endif()

  if(NOT ldlt_use_system_blas AND LDLT_BLAS_UPPER STREQUAL "SYSTEM")
    message(FATAL_ERROR "LDLT_BLAS=SYSTEM requested, but no usable CBLAS/LAPACK backend was found")
  endif()
endif()

if(NOT ldlt_use_system_blas)
  list(APPEND ldlt_core_sources src/blas_ref.c)
endif()

string(TOUPPER "${LDLT_ORDERING}" LDLT_ORDERING_UPPER)
if(NOT LDLT_ORDERING_UPPER MATCHES "^(AUTO|AMD|METIS|INTERNAL)$")
  message(FATAL_ERROR "LDLT_ORDERING must be AUTO, AMD, METIS, or INTERNAL")
endif()

set(ldlt_ordering_libraries)
set(ldlt_ordering_include_dirs)
set(ldlt_ordering_definitions)

if(NOT LDLT_ORDERING_UPPER STREQUAL "INTERNAL")
  find_path(LDLT_AMD_INCLUDE_DIR NAMES amd.h)
  find_library(LDLT_AMD_LIBRARY NAMES amd)
  find_path(LDLT_METIS_INCLUDE_DIR NAMES metis.h)
  find_library(LDLT_METIS_LIBRARY NAMES metis)

  if((LDLT_ORDERING_UPPER STREQUAL "AUTO" OR LDLT_ORDERING_UPPER STREQUAL "AMD") AND
     LDLT_AMD_INCLUDE_DIR AND LDLT_AMD_LIBRARY)
    list(APPEND ldlt_ordering_include_dirs "${LDLT_AMD_INCLUDE_DIR}")
    list(APPEND ldlt_ordering_libraries "${LDLT_AMD_LIBRARY}")
    list(APPEND ldlt_ordering_definitions LDLT_HAVE_AMD=1)
  elseif(LDLT_ORDERING_UPPER STREQUAL "AMD")
    message(FATAL_ERROR "LDLT_ORDERING=AMD requested, but SuiteSparse AMD was not found")
  endif()

  if((LDLT_ORDERING_UPPER STREQUAL "AUTO" OR LDLT_ORDERING_UPPER STREQUAL "METIS") AND
     LDLT_METIS_INCLUDE_DIR AND LDLT_METIS_LIBRARY)
    list(APPEND ldlt_ordering_include_dirs "${LDLT_METIS_INCLUDE_DIR}")
    list(APPEND ldlt_ordering_libraries "${LDLT_METIS_LIBRARY}")
    list(APPEND ldlt_ordering_definitions LDLT_HAVE_METIS=1)
  elseif(LDLT_ORDERING_UPPER STREQUAL "METIS")
    message(FATAL_ERROR "LDLT_ORDERING=METIS requested, but METIS was not found")
  endif()
endif()

set(ldlt_openmp_libraries)
set(ldlt_openmp_definitions)
if(LDLT_ENABLE_OPENMP)
  find_package(OpenMP COMPONENTS C QUIET)
  if(OpenMP_C_FOUND)
    list(APPEND ldlt_openmp_libraries OpenMP::OpenMP_C)
    list(APPEND ldlt_openmp_definitions LDLT_HAVE_OPENMP=1)
  endif()
endif()

add_library(ldlt ${ldlt_core_sources})
target_compile_features(ldlt PUBLIC c_std_11)
target_include_directories(ldlt
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
  PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/src"
    ${ldlt_blas_include_dirs}
    ${ldlt_ordering_include_dirs}
)
target_compile_definitions(ldlt
  PUBLIC
    ${ldlt_blas_definitions}
    ${ldlt_ordering_definitions}
    ${ldlt_openmp_definitions}
)
if(UNIX AND NOT APPLE)
  target_compile_definitions(ldlt PRIVATE _POSIX_C_SOURCE=200809L)
endif()
if(UNIX)
  target_compile_definitions(ldlt PRIVATE LDLT_USE_POSIX_SIGNALS=1)
endif()
target_link_libraries(ldlt
  PUBLIC
    ${ldlt_blas_libraries}
    ${ldlt_ordering_libraries}
    ${ldlt_openmp_libraries}
)

if(APPLE AND ldlt_use_system_blas)
  set_source_files_properties(
    src/factor.c
    src/factor_kernels.c
    src/solve.c
    PROPERTIES COMPILE_OPTIONS "-Wno-deprecated-declarations"
  )
endif()

add_library(ldlt::ldlt ALIAS ldlt)

if(LDLT_BUILD_FORTRAN_INTERFACE)
  add_library(ldlt_fortran fortran/c_fortran_ldlt.c)
  target_compile_features(ldlt_fortran PUBLIC c_std_11)
  target_include_directories(ldlt_fortran
    PUBLIC
      $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
      $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
    PRIVATE
      "${CMAKE_CURRENT_SOURCE_DIR}/src"
  )
  target_link_libraries(ldlt_fortran PUBLIC ldlt)
  add_library(ldlt::fortran ALIAS ldlt_fortran)
endif()

install(TARGETS ldlt
  EXPORT ldltTargets
  ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
)
if(TARGET ldlt_fortran)
  install(TARGETS ldlt_fortran
    EXPORT ldltTargets
    ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
    LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
    RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
  )
endif()
install(FILES include/ldlt/ldlt.h DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/ldlt")
install(EXPORT ldltTargets
  NAMESPACE ldlt::
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/ldlt"
)

if(LDLT_BUILD_TESTS)
  find_package(GTest QUIET)
  if(GTest_FOUND)
    add_executable(ldlt_tests
      tests/test_main.cc
      tests/test_basic.cc
      tests/test_solve.cc
      tests/test_matrices.cc
    )
    target_link_libraries(ldlt_tests PRIVATE ldlt GTest::gtest_main)
    enable_testing()
    add_test(NAME ldlt-tests COMMAND ldlt_tests)
  else()
    message(WARNING "GTest not found; LDLT tests disabled")
  endif()
endif()

if(LDLT_BUILD_BENCH)
  find_package(benchmark QUIET)
  if(benchmark_FOUND)
    set(ldlt_bench_sources
      bench/bench_main.cc
      bench/bench_factor.cc
      bench/bench_solve.cc
      bench/bench_cholmod.cc
      bench/bench_fixtures.cc
    )
    add_executable(ldlt_bench_pivot ${ldlt_bench_sources} bench/bench_fixtures_pivot.cc)
    add_executable(ldlt_bench_spd ${ldlt_bench_sources} bench/bench_fixtures_spd.cc)
    target_link_libraries(ldlt_bench_pivot PRIVATE ldlt benchmark::benchmark)
    target_link_libraries(ldlt_bench_spd PRIVATE ldlt benchmark::benchmark)
  else()
    message(WARNING "google-benchmark not found; LDLT benchmarks disabled")
  endif()
endif()

message(STATUS "LDLT BLAS backend: ${LDLT_BLAS_UPPER} (system=${ldlt_use_system_blas})")
message(STATUS "LDLT ordering backend: ${LDLT_ORDERING_UPPER}")
