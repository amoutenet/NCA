# Build the python module
add_cpp2py_module(nca)

target_link_libraries(nca nca_c)

# We need to include the convertes.hxx files
target_include_directories(nca PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

# Configure the version
configure_file(version.py.in version.py)

# All Python file. Copy them in the build dir to have a complete package for the tests.
set(PYTHON_SOURCES __init__.py)

foreach(f ${PYTHON_SOURCES})
 configure_file(${f} ${f} COPYONLY)
endforeach()

# Install python module to proper location
set(PYTHON_LIB_DEST ${CPP2PY_PYTHON_LIB_DEST_ROOT}/nca)
install(TARGETS nca DESTINATION ${PYTHON_LIB_DEST})
install(FILES ${PYTHON_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/version.py DESTINATION ${PYTHON_LIB_DEST})
