find_package(PythonLibs)

pybind11_add_module(pymrcpp SHARED
    pymrcpp.cpp
    pyBases.cpp
    #pyFundamentalTypes.cpp
    #pyGauss.cpp
    pyMethods.cpp
    pyOperators.cpp
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
