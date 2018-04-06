include_directories(${CMAKE_BINARY_DIR})
include_directories(${CMAKE_SOURCE_DIR}/include)

add_subdirectory(mrcpp)

add_executable(unit-tests unit_tests.cpp $<TARGET_OBJECTS:mrcpp-tests>)

target_link_libraries(unit-tests mrcpp)
