find_package(PythonLibs)

pybind11_add_module(pymrcpp SHARED
    pymrcpp.cpp
    pyTmpEverything.cpp
    PyBoundingBox.cpp
    pyProject.cpp
)

target_include_directories(pymrcpp PRIVATE
    ${PYTHON_INCLUDE_DIRS}
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/external/pybind11/include/
    )

target_link_libraries(pymrcpp PRIVATE
    mrcpp
    )
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/test_python.py ${CMAKE_CURRENT_BINARY_DIR}/test_python.py)
