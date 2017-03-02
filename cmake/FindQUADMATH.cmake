set(QUADMATH_TEST_SOURCE
"
#include <quadmath.h>
int main() {
fabsq(1);
}
")

include(CheckCSourceCompiles)
set(SAFE_CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS}")
set(CMAKE_REQUIRED_LIBRARIES "quadmath")
set(CHECK_QUADMATH_ADD_LIBRARIES LINK_LIBRARIES quadmath)
check_c_source_compiles("${QUADMATH_TEST_SOURCE}" QUADMATH_FOUND)
set(CMAKE_REQUIRED_FLAGS "${SAFE_CMAKE_REQUIRED_FLAGS}")

