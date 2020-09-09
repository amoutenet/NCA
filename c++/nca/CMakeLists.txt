FILE(GLOB_RECURSE cpp_sources RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)

add_library(nca_c ${cpp_sources})

target_link_libraries(nca_c PUBLIC triqs)

install(TARGETS nca_c DESTINATION lib)
