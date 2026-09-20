cmake_minimum_required(VERSION 3.20)
project(dials_metal LANGUAGES CXX)

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
  src/json.cc
  src/expt.cc
  src/refl.cc
  src/fft.cc
  src/index.cc
  src/refine.cc
  src/args.cc
  src/derivatives.cc
  src/shoebox.cc
  src/profile_model.cc
  src/profile_grid.cc
  src/mask.cc
)
target_include_directories(mxi PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_compile_options(mxi PRIVATE -Wall -Wextra -Wpedantic)

add_executable(mxi_index apps/mxi_index.cc)
target_link_libraries(mxi_index PRIVATE mxi)
target_compile_options(mxi_index PRIVATE -Wall -Wextra)

add_executable(mxi_refine apps/mxi_refine.cc)
target_link_libraries(mxi_refine PRIVATE mxi)
target_compile_options(mxi_refine PRIVATE -Wall -Wextra)

add_executable(mxi_residuals apps/mxi_residuals.cc)
target_link_libraries(mxi_residuals PRIVATE mxi)
target_compile_options(mxi_residuals PRIVATE -Wall -Wextra)

add_executable(mxi_profile apps/mxi_profile.cc)
target_link_libraries(mxi_profile PRIVATE mxi)
target_compile_options(mxi_profile PRIVATE -Wall -Wextra)

add_executable(mxi_mask apps/mxi_mask.cc)
target_link_libraries(mxi_mask PRIVATE mxi)
target_compile_options(mxi_mask PRIVATE -Wall -Wextra)

add_executable(mxi_grid apps/mxi_grid.cc)
target_link_libraries(mxi_grid PRIVATE mxi)
target_compile_options(mxi_grid PRIVATE -Wall -Wextra)

# The spot finder, which needs HDF5 and the bitshuffle submodule. Both are
# checked rather than required, so that a machine with neither still configures
# and still builds and tests everything else -- the clean-clone check has caught
# several real problems in this repository and it has to keep working on a
# machine that cannot build the whole of it.
#
# It is added unmodified. Its own CMakeLists calls find_package(HDF5 REQUIRED),
# which would abort the whole configure rather than skip one subdirectory, so
# the decision is made here before the subdirectory is entered.
find_package(HDF5 QUIET COMPONENTS C)
if(NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/spotfinder/third_party/bitshuffle/LICENSE")
    message(STATUS
            "spotfinder: bitshuffle submodule not checked out, skipping. "
            "git submodule update --init --recursive")
elseif(NOT HDF5_FOUND)
    message(STATUS "spotfinder: HDF5 not found, skipping")
else()
    message(STATUS "spotfinder: building")
    add_subdirectory(spotfinder)
endif()

enable_testing()
add_executable(mxi_tests
  tests/main.cc
  tests/test_geometry.cc
  tests/test_real_geometry.cc
  tests/test_multi_axis.cc
  tests/test_index.cc
  tests/test_refine.cc
  tests/test_args.cc
  tests/test_derivatives.cc
  tests/test_precision.cc
  tests/test_expt_format.cc
  tests/test_shoebox.cc
  tests/test_mask.cc
  tests/test_profile_model.cc
)
target_link_libraries(mxi_tests PRIVATE mxi)
target_compile_options(mxi_tests PRIVATE -Wall -Wextra)
add_test(NAME mxi_tests COMMAND mxi_tests)
