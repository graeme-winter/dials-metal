cmake_minimum_required(VERSION 3.20)
project(dials_metal_index LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

# No third-party dependencies, deliberately. The spot finder hand-rolled its
# msgpack writer for the same reason: this has to build on a laptop and on a
# beamline machine with nothing installed on it.
add_library(mxi
  src/linalg.cc
  src/geometry.cc
  src/predict.cc
)
target_include_directories(mxi PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_compile_options(mxi PRIVATE -Wall -Wextra -Wpedantic)

enable_testing()
add_executable(mxi_tests
  tests/main.cc
  tests/test_geometry.cc
  tests/test_real_geometry.cc
)
target_link_libraries(mxi_tests PRIVATE mxi)
target_compile_options(mxi_tests PRIVATE -Wall -Wextra)
add_test(NAME mxi_tests COMMAND mxi_tests)
