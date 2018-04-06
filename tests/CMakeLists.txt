add_executable(mrcpp-tests unit_tests.cpp)
add_subdirectory(mrcpp)

target_include_directories(mrcpp-tests PRIVATE
    ${CMAKE_SOURCE_DIR}/src/
    ${CMAKE_SOURCE_DIR}/src/mwcore
    ${CMAKE_SOURCE_DIR}/src/mwtrees
    ${CMAKE_SOURCE_DIR}/src/mwbuilders
    ${CMAKE_SOURCE_DIR}/src/mwfunctions
    ${CMAKE_SOURCE_DIR}/src/mwoperators
    ${CMAKE_SOURCE_DIR}/src/mwutils
    ${CMAKE_SOURCE_DIR}/src/gaussians
    )

target_link_libraries(mrcpp-tests
    mrcpp
    )
