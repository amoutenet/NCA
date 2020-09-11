# Copy h5 files to binary dir
FILE(GLOB all_h5_files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.h5)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/${all_h5_files} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

macro(add_test_defs name)
 add_executable(${name} ${name}.cpp)
 target_link_libraries(${name} PRIVATE nca_c gtest)
 add_test(NAME ${name} COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${name})
 #triqs_add_cpp_test(${name})
endmacro()

add_test_defs(anderson)

# Perf
set(perf_sources "perf.cpp")
add_executable(perf ${perf_sources})
target_link_libraries(perf nca_c ${TRIQS_LIBRARY_ALL})
include_directories(${TRIQS_INCLUDE_ALL})
triqs_set_rpath_for_target(perf)
