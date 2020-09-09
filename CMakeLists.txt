# Start configuration
cmake_minimum_required(VERSION 2.8)
project(nca CXX)

# Default to Release build type
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Type of build" FORCE)
endif()
message(STATUS "-------- BUILD-TYPE: ${CMAKE_BUILD_TYPE} -------------")

# Use shared libraries
set(BUILD_SHARED_LIBS ON)

# Load TRIQS and Cpp2Py
find_package(TRIQS 2.0 REQUIRED)
find_package(Cpp2Py 1.5 REQUIRED)

# Default Install directory to TRIQS_ROOT if not given. Checks an absolute name is given.
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT OR (NOT IS_ABSOLUTE ${CMAKE_INSTALL_PREFIX}))
  message(STATUS "No install prefix given (or invalid). Defaulting to TRIQS_ROOT")
  set(CMAKE_INSTALL_PREFIX ${TRIQS_ROOT} CACHE PATH "default install path" FORCE)
endif()
message(STATUS "-------- CMAKE_INSTALL_PREFIX: ${CMAKE_INSTALL_PREFIX} -------------")

# Use git to detect NCA version and git hash
execute_process(COMMAND bash "-c" "${TRIQS_GIT_EXECUTABLE} describe | sed \"s/-\\([0-9]*\\)-.*/.\\1/g\""
 WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
 OUTPUT_VARIABLE   NCA_VERSION
 OUTPUT_STRIP_TRAILING_WHITESPACE
)
if(NOT NCA_VERSION)
  set(NCA_VERSION 1.0)
endif()
triqs_get_git_hash_of_source_dir(NCA_GIT_HASH)
MESSAGE(STATUS "NCA version : ${NCA_VERSION}")
message(STATUS "Git hash: ${NCA_GIT_HASH}")

# The std for all targets
add_compile_options(-std=c++17 $<$<CONFIG:Release>:-march=native>)

# Common to all targets (ctint, python, tests).
include_directories(c++)

# To be fixed!!
triqs_set_rpath_for_all_targets()

# Main c++ code
add_subdirectory(c++/nca)

# Python interface
if (${TRIQS_WITH_PYTHON_SUPPORT})
 add_subdirectory(python/nca)
endif()

# Tests
option(Build_Tests "Enable Tests" ON)
if (${Build_Tests})
 enable_testing()
 add_subdirectory(test)
endif()
