cmake_minimum_required (VERSION 3.15)
project (cint C)
set(cint_VERSION_MAJOR "6")
set(cint_VERSION_MINOR "1")
set(cint_VERSION_PATCH "3")
set(cint_VERSION_TWEAK "0")
set(cint_VERSION "${cint_VERSION_MAJOR}.${cint_VERSION_MINOR}.${cint_VERSION_PATCH}")
set(cint_SOVERSION "${cint_VERSION_MAJOR}")

if ("${CMAKE_BUILD_TYPE}" STREQUAL "")
  set(CMAKE_BUILD_TYPE RELWITHDEBINFO)
  message("Build type RELWITHDEBINFO")
else()
  message("Build type ${CMAKE_BUILD_TYPE}")
endif()
set(CMAKE_VERBOSE_MAKEFILE OFF)

option(BUILD_MARCH_NATIVE "gcc flag -march=native" off)
if (BUILD_MARCH_NATIVE)
  include(CheckCCompilerFlag)
  CHECK_C_COMPILER_FLAG("-march=native" COMPILER_SUPPORTS_MARCH_NATIVE)
  if(COMPILER_SUPPORTS_MARCH_NATIVE)
    if ("${CMAKE_C_COMPILER_ID}" STREQUAL "Intel")
      message("Add CFLAGS -march=native -unroll-aggressive")
      set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -unroll-aggressive -ipo")
      set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fno-math-errno")
    else()
      message("Add CFLAGS -march=native -ftree-vectorize")
      set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -march=native -funroll-loops -ftree-vectorize")
      set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fno-math-errno")
    endif()
  endif()
endif()

#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -O2 -DNDEBUG")

set(CMAKE_C_STANDARD 99)

include(GNUInstallDirs)
#set(CMAKE_INCLUDE_CURRENT_DIR ON)

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
find_package(QUADMATH)
if(QUADMATH_FOUND)
  set(HAVE_QUADMATH_H 1 CACHE INTERNAL "Have QUADMATH")
endif()

include(CheckFunctionExists)
set(CMAKE_REQUIRED_INCLUDES math.h)
set(CMAKE_REQUIRED_LIBRARIES m)
check_function_exists(expl HAVE_EXPL)
check_function_exists(sqrtl HAVE_SQRTL)
#set(CMAKE_REQUIRED_INCLUDES quadmath.h)
#set(CMAKE_REQUIRED_LIBRARIES quadmath)
#check_function_exists(fabsq HAVE_QUADMATH_H)


configure_file(
  "${PROJECT_SOURCE_DIR}/src/cint_config.h.in"
  "${PROJECT_BINARY_DIR}/src/cint_config.h")

configure_file(
  "${PROJECT_SOURCE_DIR}/include/cint.h.in"
  "${PROJECT_BINARY_DIR}/include/cint.h")

set(cintSrc
  src/c2f.c src/cart2sph.c src/cint1e.c src/cint2e.c src/cint_bas.c
  src/fblas.c src/g1e.c src/g2e.c src/misc.c src/optimizer.c
  src/fmt.c src/rys_wheeler.c src/eigh.c src/rys_roots.c src/find_roots.c
  src/cint2c2e.c src/g2c2e.c src/cint3c2e.c src/g3c2e.c
  src/cint3c1e.c src/g3c1e.c src/breit.c
  src/cint1e_a.c src/cint3c1e_a.c
  src/cint1e_grids.c src/g1e_grids.c
  src/autocode/breit1.c src/autocode/dkb.c src/autocode/gaunt1.c
  src/autocode/grad1.c src/autocode/grad2.c src/autocode/hess.c
  src/autocode/int3c1e.c src/autocode/int3c2e.c src/autocode/intor1.c
  src/autocode/intor2.c src/autocode/intor3.c src/autocode/intor4.c
  src/autocode/deriv3.c src/autocode/int1e_grids1.c
  src/autocode/deriv4.c src/autocode/lresc.c)

#option(WITH_RANGE_COULOMB "Range separated Coulomb" on)
if(WITH_RANGE_COULOMB)
# defined in cint_config.h
#  add_definitions(-DWITH_RANGE_COULOMB)
#  message("Enabled WITH_RANGE_COULOMB")
endif(WITH_RANGE_COULOMB)

if(WITH_POLYNOMIAL_FIT)
  set(cintSrc ${cintSrc} src/polyfits.c src/sr_rys_polyfits.c)
  add_definitions(-DWITH_POLYNOMIAL_FIT)
  message("Enabled WITH_POLYNOMIAL_FIT")
endif(WITH_POLYNOMIAL_FIT)

if(WITH_COULOMB_ERF)
  message("WITH_COULOMB_ERF is deprecated since v6.0")
endif(WITH_COULOMB_ERF)

if(WITH_F12)
  set(cintSrc ${cintSrc} src/g2e_f12.c src/stg_roots.c src/cint2e_f12.c)
  add_definitions(-DWITH_F12)
  message("Enabled WITH_F12")
endif(WITH_F12)

# Gaussian type geminal code are incorrect
if(WITH_GTG)
#  set(cintSrc ${cintSrc} src/g2e_gtg.c src/cint2e_gtg.c src/cint3c2e_gtg.c
#    src/cint2c2e_gtg.c)
#  add_definitions(-DWITH_GTG)
#  message("Enabled WITH_GTG. Note there are bugs in gtg type integrals")
  message("WITH_GTG is deprecated since v6.0")
endif(WITH_GTG)

if(WITH_4C1E)
  set(cintSrc ${cintSrc} src/cint4c1e.c src/g4c1e.c)
  message("Enabled WITH_4C1E. Note there are bugs in 4c1e integral functions")
endif(WITH_4C1E)

if(PYPZPX)
  add_definitions(-DPYPZPX)
  message("P orbitals convention (py, pz, px)")
else()
  message("P orbitals convention (px, py, pz)")
endif(PYPZPX)

if(MIN_EXPCUTOFF)
  message("Minimal cutoff in exponent for integral screening ${MIN_EXPCUTOFF}")
  add_definitions(-DMIN_EXPCUTOFF=${MIN_EXPCUTOFF})
endif()

if(KEEP_GOING)
  message("Do not trigger hard exit for numerical issues in Rys quadrature")
  add_definitions(-DKEEP_GOING)
endif()

option(WITH_FORTRAN "Fortran interface" on)
if(WITH_FORTRAN)
  add_definitions(-DWITH_FORTRAN)
  message("Enable Fortran interface")
else()
  message("Exclude Fortran interface")
endif(WITH_FORTRAN)

option(WITH_CINT2_INTERFACE "Old libcint (version 2) interface" on)
if(WITH_CINT2_INTERFACE)
  add_definitions(-DWITH_CINT2_INTERFACE)
  message("Enable old cint (version 2) interface")
else()
  message("Exclude old cint (version 2) interface")
endif(WITH_CINT2_INTERFACE)

option(BUILD_SHARED_LIBS "build shared libraries" 1)
option(ENABLE_EXAMPLE "build examples" 0)
option(ENABLE_TEST "build tests" 0)
option(ENABLE_STATIC "Enforce static library build" 0)
if(QUICK_TEST)
  set(RUN_QUICK_TEST --quick)
endif()

if(ENABLE_STATIC)
  set(BUILD_SHARED_LIBS 0)
endif()


set(CintHeaders
  ${PROJECT_SOURCE_DIR}/include/cint_funcs.h
  ${PROJECT_BINARY_DIR}/include/cint.h)


add_library(cint ${cintSrc})

target_include_directories(cint
  PUBLIC 
  $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
  INTERFACE
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
  PRIVATE
  ${PROJECT_BINARY_DIR}/src
  ${PROJECT_SOURCE_DIR}/src
)


set_target_properties(cint PROPERTIES
  VERSION ${cint_VERSION}
  SOVERSION ${cint_SOVERSION}
  LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})
if(QUADMATH_FOUND)
  target_link_libraries(cint quadmath)
endif()
target_link_libraries(cint "-lm")

# for consistency between FetchContent / add_subdirectory
# and find_package, so one can always link to cint::cint.
add_library(cint::cint ALIAS cint)

# install headers
install(FILES ${CintHeaders} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} COMPONENT "dev")

# install the library and associate the installed targets
# with export set "cint-targets"
install(TARGETS cint
        EXPORT cint-targets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

# generate cintTargets.cmake and prepend the namespace "cint::"
# to exported targets.
export(EXPORT cint-targets
    FILE ${CMAKE_CURRENT_BINARY_DIR}/cintTargets.cmake
        NAMESPACE cint::
)

# install cintTargets.cmake and associate this file with
# the export set "cint-targets"
install(EXPORT cint-targets
        FILE cintTargets.cmake
        NAMESPACE cint::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/cint
)

include(CMakePackageConfigHelpers)

# configure cintConfig.cmake
configure_package_config_file(${PROJECT_SOURCE_DIR}/cmake/cintConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/cintConfig.cmake
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/cint)

# generate cintConfigVersion.cmake
write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/cintConfigVersion.cmake"
    VERSION ${cint_VERSION}
    COMPATIBILITY SameMajorVersion
)

# and install both of these two files...
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/cintConfig.cmake"
              "${CMAKE_CURRENT_BINARY_DIR}/cintConfigVersion.cmake"
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/cint
)

# configure the pkg-config file...
configure_file(${PROJECT_SOURCE_DIR}/cmake/cint.pc.in
               ${CMAKE_CURRENT_BINARY_DIR}/cint.pc @ONLY)

# and install it.
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/cint.pc"
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig
)


if(ENABLE_EXAMPLE)
  enable_language(Fortran)
  find_package(OpenMP)
  if(OPENMP_FOUND)
    set(HAVE_OPENMP 1)
  else ()
    set(OpenMP_C_FLAGS " ")
  endif()
  add_subdirectory(examples)
endif()

if(BUILD_SHARED_LIBS AND ENABLE_TEST)
  find_package(PythonInterp)
  message(STATUS "Found python  ${PYTHON_EXECUTABLE}")
  if(PYTHONINTERP_FOUND)
    enable_testing()
    add_test(NAME cinttest
      COMMAND ${PYTHON_EXECUTABLE} ${PROJECT_SOURCE_DIR}/testsuite/test_cint.py ${RUN_QUICK_TEST})
    add_test(NAME cint3c2etest
      COMMAND ${PYTHON_EXECUTABLE} ${PROJECT_SOURCE_DIR}/testsuite/test_3c2e.py ${RUN_QUICK_TEST})
  endif()
endif()

include(CPackComponent)

set(CPACK_PACKAGE_CHECKSUM "SHA512")
set(CPACK_GENERATOR "DEB")
set(CPACK_DEBIAN_COMPRESSION_TYPE "xz")
set(CPACK_PACKAGE_VENDOR "Qiming @sunqm Sun")
set(CPACK_PACKAGE_CONTACT "https://github.com/sunqm/libcint")
set(CPACK_PACKAGE_HOMEPAGE_URL "https://github.com/sunqm/libcint")
set(CPACK_RESOURCE_FILE_LICENSE "${PROJECT_SOURCE_DIR}/LICENSE")
set(CPACK_RESOURCE_FILE_README "${PROJECT_SOURCE_DIR}/README.rst")
set(CPACK_PACKAGE_DESCRIPTION_FILE "${PROJECT_SOURCE_DIR}/doc/libcint.CPack.txt")

set(CPACK_PACKAGE_VERSION_MAJOR "${cint_VERSION_MAJOR}")
set(CPACK_PACKAGE_VERSION_MINOR "${cint_VERSION_MINOR}")
set(CPACK_PACKAGE_VERSION_PATCH "${cint_VERSION_PATCH}")
set(CPACK_PACKAGE_VERSION_TWEAK "${cint_VERSION_TWEAK}")

set(CPACK_COMPONENTS_GROUPING "IGNORE")

set(CPACK_MONOLITHIC_INSTALL OFF)
set(CPACK_DEB_COMPONENT_INSTALL ON)
set(CPACK_RPM_COMPONENT_INSTALL ON)
set(CPACK_NSIS_COMPONENT_INSTALL ON)

set("CPACK_DEBIAN_LIB_PACKAGE_NAME" "libcint")
set("CPACK_DEBIAN_DEV_PACKAGE_NAME" "${CPACK_DEBIAN_LIB_PACKAGE_NAME}-dev")

set(CPACK_DEBIAN_LIB_PACKAGE_DEPENDS "libopenblas0-pthread")
set(CPACK_DEBIAN_DEV_PACKAGE_DEPENDS "${CPACK_DEBIAN_LIB_PACKAGE_NAME},libopenblas-pthread-dev")
if(HAVE_QUADMATH)
  set(CPACK_DEBIAN_LIB_PACKAGE_DEPENDS "${CPACK_DEBIAN_LIB_PACKAGE_DEPENDS},libquadmath0")
endif()

set(CPACK_COMPONENTS_GROUPING "IGNORE")

include(CPack)
