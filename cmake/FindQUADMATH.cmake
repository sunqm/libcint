find_library(QUADMATH_LIBRARY
  NAMES quadmath
  PATHS 
  /usr/lib /usr/lib64
  /usr/local/lib /usr/local/lib64
  /usr/x86_64-linux-gnu/*
  /usr/lib/gcc/x86_64-linux-gnu/*
  ${QUADMATH_ROOT_DIR}
)

include(FindPackageHandleStandardArgs) 
find_package_handle_standard_args(QUADMATH
  FOUND_VAR
  QUADMATH_FOUND
  REQUIRED_VARS
  QUADMATH_LIBRARY
)

mark_as_advanced(QUADMATH_INCLUDE_DIR, QUADMATH_LIBRARY)
