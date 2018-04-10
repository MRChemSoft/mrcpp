add_executable(mrcpp-tests unit_tests.cpp)
add_subdirectory(mrcpp)

target_include_directories(mrcpp-tests PRIVATE
    ${CMAKE_SOURCE_DIR}/src/
    )

target_link_libraries(mrcpp-tests
    mrcpp
    )
