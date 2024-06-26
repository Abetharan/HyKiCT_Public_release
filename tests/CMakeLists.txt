set (CMAKE_CXX_STANDARD 17)
find_package(GTEST QUIET)
# NOTE: the upper case GTEST! CK
if(NOT GTEST_FOUND)
  # Download and unpack googletest at configure time
  # but only if needed! CK
  configure_file(${CMAKE_SOURCE_DIR}/cmake-modules/GoogleTest-CMakeLists.txt.in
    ${CMAKE_BINARY_DIR}/googletest-download/CMakeLists.txt)
  execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
    RESULT_VARIABLE result
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/googletest-download)
  if(result)
    message(FATAL_ERROR "CMake step for googletest failed: ${result}")
  endif()

  execute_process(COMMAND ${CMAKE_COMMAND} --build .
    RESULT_VARIABLE result
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/googletest-download)
  if(result)
    message(FATAL_ERROR "Build step for googletest failed: ${result}")
  endif()

  # Add googletest directly to our build. This defines
  # the gtest and gtest_main targets.
  add_subdirectory(${CMAKE_BINARY_DIR}/googletest-src
    ${CMAKE_BINARY_DIR}/googletest-build
    EXCLUDE_FROM_ALL)
endif()

set_target_properties(gmock PROPERTIES
    IMPORTED_LOCATION ${PROJECT_BINARY_DIR}/prefix/lib/${CMAKE_STATIC_LIBRARY_PREFIX}yaml-cpp${CMAKE_STATIC_LIBRARY_SUFFIX}
)
# Now simply link against gtest as needed. Eg
add_library(test_suite 
            ${HyKiCT_SOURCE_DIR}/tests/RadiationTests.cpp
            ${HyKiCT_SOURCE_DIR}/tests/RadiationTests.h 
            ${HyKiCT_SOURCE_DIR}/tests/ImplicitRadiationTests.cpp
            ${HyKiCT_SOURCE_DIR}/tests/ImplicitRadiationTests.h 
            ${HyKiCT_SOURCE_DIR}/tests/IntSourceTests.cpp
            ${HyKiCT_SOURCE_DIR}/tests/IntSourceTests.h 
            ${HyKiCT_SOURCE_DIR}/tests/GridDataTests.cpp
            ${HyKiCT_SOURCE_DIR}/tests/ExtSourceTests.cpp
            ${HyKiCT_SOURCE_DIR}/tests/LoadInTablesTests.cpp
            ${HyKiCT_SOURCE_DIR}/tests/FluidDynamicsTests.cpp
            )

target_link_libraries(test_suite gtest HYKICT::HyKiCT_lib)
add_executable(test_driver tests.cpp)
target_compile_features(test_driver PRIVATE cxx_std_17)
target_link_libraries(test_driver test_suite)

