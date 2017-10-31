include(DefaultPathsMRCPP)
include_directories(${PROJECT_BINARY_DIR})
include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/external/catch)

add_subdirectory(mrcpp)

add_executable(unit-tests unit_tests.cpp $<TARGET_OBJECTS:mrcpp-tests>)

target_link_libraries(unit-tests mrcpp)
