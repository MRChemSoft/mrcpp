find_package(PythonLibs)

pybind11_add_module(pymrcpp SHARED
    pymrcpp.cpp
    PyBoundingBox.cpp
    pyProject.cpp
    pyOperators.cpp
    pyMethods.cpp
    pyBases.cpp
    pyFundamentalTypes.cpp
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
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/test_pymrcpp3d.py ${CMAKE_CURRENT_BINARY_DIR}/test_pymrcpp3d.py)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/test_pymrcpp2d.py ${CMAKE_CURRENT_BINARY_DIR}/test_pymrcpp2d.py)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/test_pymrcpp1d.py ${CMAKE_CURRENT_BINARY_DIR}/test_pymrcpp1d.py)
