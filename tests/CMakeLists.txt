add_executable(mrcpp-tests unit_tests.cpp)
add_subdirectory(mrcpp)

target_include_directories(mrcpp-tests PRIVATE
    ${CMAKE_SOURCE_DIR}/src/
    ${CMAKE_SOURCE_DIR}/include/
    )

target_link_libraries(mrcpp-tests
    mrcpp
    )
