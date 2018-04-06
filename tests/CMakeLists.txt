include_directories(${CMAKE_BINARY_DIR})
include_directories(${CMAKE_SOURCE_DIR}/src)
link_directories(${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})

add_subdirectory(mrcpp)

add_executable(unit-tests unit_tests.cpp $<TARGET_OBJECTS:mrcpp-tests>)

target_link_libraries(unit-tests mrcpp)
