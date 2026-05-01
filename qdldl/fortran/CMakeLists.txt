# Fortran bindings for QDLDL.
#
# When included via add_subdirectory from the top-level qdldl CMakeLists,
# the `qdldlstatic` target and its include directories are already defined.
#
# To build standalone (outside the qdldl tree), point CMake at the qdldl
# install prefix so that find_package(qdldl) resolves `qdldlstatic`.

cmake_minimum_required(VERSION 3.18)

# Only declare a project when built as a standalone tree (not as subdir).
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  project(qdldl_fortran C)
  find_package(qdldl REQUIRED)
endif()

# ---- Fortran-callable C wrapper library ------------------------------------
add_library(qdldl_fortran STATIC c_fortran_qdldl.c)

target_include_directories(qdldl_fortran
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../include>
    $<INSTALL_INTERFACE:include>
)

target_link_libraries(qdldl_fortran PUBLIC qdldlstatic)

# ---- METIS fill-reducing ordering (optional) --------------------------------
# When built as a subdirectory of the MYSTRAN tree, HAVE_METIS and the `metis`
# target are already set by the parent.  When built standalone, search for a
# system-installed METIS via find_path / find_library.
option(QDLDL_FORTRAN_USE_METIS "Use METIS for fill-reducing ordering" ON)

if(QDLDL_FORTRAN_USE_METIS)
  if(DEFINED HAVE_METIS AND HAVE_METIS)
    # Parent build already found METIS; reuse its target.
    target_compile_definitions(qdldl_fortran PRIVATE HAVE_METIS=1)
    target_link_libraries(qdldl_fortran PUBLIC metis)
  else()
    find_path(METIS_INCLUDE_DIR  metis.h  PATH_SUFFIXES metis)
    find_library(METIS_LIBRARY   NAMES metis)

    if(METIS_INCLUDE_DIR AND METIS_LIBRARY)
      message(STATUS "Found system METIS: ${METIS_LIBRARY}")
      target_compile_definitions(qdldl_fortran PRIVATE HAVE_METIS=1)
      target_include_directories(qdldl_fortran PRIVATE ${METIS_INCLUDE_DIR})
      target_link_libraries(qdldl_fortran PUBLIC ${METIS_LIBRARY})
    else()
      message(STATUS "METIS not found — identity ordering will be used (set METIS_INCLUDE_DIR / METIS_LIBRARY to override)")
    endif()
  endif()
endif()

# ---- Fortran example -------------------------------------------------------
cmake_dependent_option(
  QDLDL_FORTRAN_EXAMPLE
  "Build the Fortran example program"
  ON
  "QDLDL_FORTRAN_BINDINGS" OFF
)

if(QDLDL_FORTRAN_EXAMPLE)
  enable_language(Fortran)
  add_executable(qdldl_fortran_example ../examples/fortran_example.f90)
  target_link_libraries(qdldl_fortran_example PRIVATE qdldl_fortran)
  set_target_properties(qdldl_fortran_example
    PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/out"
  )
endif()
