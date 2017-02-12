find_path(QUADMATH_INCLUDE_DIR
  NAMES quadmath.h
  PATHS
  ${QUADMATH_ROOT_DIR}/include
  ${QUADMATH_ROOT_DIR}/usr/include
  ${QUADMATH_ROOT_DIR}/usr/local/include
)

find_library(QUADMATH_LIBRARY
  NAMES quadmath
  PATHS 
  ${QUADMATH_ROOT_DIR}
  ${QUADMATH_ROOT_DIR}/usr/lib
  ${QUADMATH_ROOT_DIR}/usr/lib64
  ${QUADMATH_ROOT_DIR}/usr/local/lib
  ${QUADMATH_ROOT_DIR}/lib
  ${QUADMATH_ROOT_DIR}/lib64
)

include(FindPackageHandleStandardArgs) 
find_package_handle_standard_args(QUADMATH
  FOUND_VAR
  QUADMATH_FOUND
  REQUIRED_VARS
  QUADMATH_INCLUDE_DIR
  QUADMATH_LIBRARY
)

mark_as_advanced(QUADMATH_INCLUDE_DIR, QUADMATH_LIBRARY)
