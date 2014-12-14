cmake_minimum_required (VERSION 2.6)
project (cint)
set(cint_VERSION "2.3.0")

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -O2 -DNDEBUG")
set(CMAKE_VERBOSE_MAKEFILE OFF)

#set(CMAKE_INCLUDE_CURRENT_DIR ON)

find_package(BLAS REQUIRED)

include_directories(
  ${PROJECT_SOURCE_DIR}/include
  ${PROJECT_SOURCE_DIR}/src
  ${PROJECT_BINARY_DIR}/src)

configure_file(
  "${PROJECT_SOURCE_DIR}/src/config.h.in"
  "${PROJECT_BINARY_DIR}/src/config.h")

set(cintSrc 
  src/c2f.c src/cart2sph.c src/cint1e.c src/cint2e.c src/cint_bas.c
  src/fblas.c src/g1e.c src/g2e.c src/misc.c src/optimizer.c
  src/rys_roots.c
  src/autocode/auto_intor1.c src/autocode/auto_intor2.c
  src/autocode/auto_intor3.c src/autocode/gaunt1.c src/autocode/grad1.c
  src/autocode/grad2.c src/autocode/grad3.c)


if (${ENABLE_STATIC})
  add_library(cintStatic STATIC ${cintSrc})
  set_target_properties(cintStatic PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}
    OUTPUT_NAME cint)

  install(TARGETS cintStatic DESTINATION lib)
  install(FILES ${PROJECT_SOURCE_DIR}/include/cint.h DESTINATION include)

else()

  add_library(cint SHARED ${cintSrc})

  set_target_properties(cint PROPERTIES
    VERSION ${cint_VERSION}
    LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})

  target_link_libraries(cint ${BLAS_LIBRARIES})

  install(TARGETS cint DESTINATION lib)
  install(FILES ${PROJECT_SOURCE_DIR}/include/cint.h DESTINATION include)

  if (${ENABLE_EXAMPLE})
    enable_language(Fortran)
    add_subdirectory(examples)

    find_package(PythonInterp)
    message(STATUS "Found python  ${PYTHON_EXECUTABLE}")
    if (PYTHONINTERP_FOUND)
      enable_testing()
      add_test(NAME cinttest
        COMMAND ${PYTHON_EXECUTABLE} ${PROJECT_SOURCE_DIR}/testsuite/test_cint.py)
    endif()
  endif()

endif()
