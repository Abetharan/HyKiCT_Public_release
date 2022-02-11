# Works with 3.13 and tested through 3.15
cmake_minimum_required(VERSION 3.13...3.16)

# Project name and a few useful settings. Other commands can pick up the results
project(HyKiCT
    VERSION 2.0
    DESCRIPTION "A 1-D Lagrangian Gray Diffusive Rad-Hydro Transport Code Designed For Coupling with VFP Code"
    LANGUAGES CXX)

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake-modules" ${CMAKE_MODULE_PATH})
string(REGEX REPLACE "( -DNDEBUG$|-DNDEBUG )" "" CMAKE_CXX_FLAGS_BETATEST "${CMAKE_CXX_FLAGS_RELEASE}" ) 
string(REGEX REPLACE "( -DNDEBUG$|-DNDEBUG )" "" CMAKE_C_FLAGS_BETATEST "${CMAKE_C_FLAGS_RELEASE}" )    
string(REGEX REPLACE "-DNDEBUG " "" CMAKE_CXX_FLAGS_RELWITHDEBUG "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -DDEBUG" )
string(REGEX REPLACE "-DNDEBUG " "" CMAKE_C_FLAGS_RELWITHDEBUG "${CMAKE_C_FLAGS_RELWITHDEBINFO} -DDEBUG" )
#Ignore empty file errors 
set(CMAKE_CXX_ARCHIVE_FINISH "<CMAKE_RANLIB> -no_warning_for_no_symbols -c <TARGET>")

include_directories(include)
### Require out-of-source builds
file(TO_CMAKE_PATH "${PROJECT_BINARY_DIR}/CMakeLists.txt" LOC_PATH)
if(EXISTS "${LOC_PATH}")
    message(FATAL_ERROR "You cannot build in a source directory (or any directory with a CMakeLists.txt file). Please make a build subdirectory. Feel free to remove CMakeCache.txt and CMakeFiles.")
endif()

# Only do these if this is the main project, and not if it is included through add_subdirectory
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)

    # Optionally set things like CMAKE_CXX_STANDARD, CMAKE_POSITION_INDEPENDENT_CODE here
    
    # Let's ensure -std=c++xx instead of -std=g++xx
    set(CMAKE_CXX_EXTENSIONS OFF)

    # Let's nicely support folders in IDE's
    set_property(GLOBAL PROPERTY USE_FOLDERS ON)

    # find_package(PETSc REQUIRED)
    # message (STATUS "PETSC COMPILER ${PETSC_COMPILER}")
    # add_definitions (${PETSC_DEFINITIONS})
    # add_definitions (-g)
    # message (STATUS "PETSC_DEFINITIONS ${PETSC_DEFINITIONS}")
    find_package(yaml-cpp REQUIRED)
    if(${yaml-cpp_FOUND})
    	message(STATUS "FOUND YAML_CPP")
    endif()
    # include_directories(${PETSC_INCLUDES})
    # find_package(Boost COMPONENTS thread REQUIRED)
    # INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR} )
    # Testing only available if this is the main app
    # Note this needs to be done in the main CMakeLists
    # since it calls enable_testing, which must be in the
    # main CMakeLists.

    add_subdirectory (src)

    # The executable code is here
    add_subdirectory(apps)
        #[[ Docs only available if this is the main app
        #find_package(Doxygen)
        #if(Doxygen_FOUND)
        #add_subdirectory(docs)
        #else()
        #message(STATUS "Doxygen not found, not building docs")
        #endif()]]
endif()

# Testing only available if this is the main app
# Emergency override MODERN_CMAKE_BUILD_TESTING provided as well
enable_testing()
include(GoogleTest)
add_subdirectory(tests)

