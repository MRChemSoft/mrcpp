# Prepare "Catch" library for other executables
set(CATCH_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/external/catch)
add_library(Catch INTERFACE)
target_include_directories(Catch INTERFACE ${CATCH_INCLUDE_DIR})

add_executable(mrcpp-tests unit_tests.cpp)
add_subdirectory(mwbuilders)
add_subdirectory(mwcore)
add_subdirectory(mwfunctions)
add_subdirectory(mwoperators)
add_subdirectory(mwtrees)

target_include_directories(mrcpp-tests PRIVATE
    ${CMAKE_SOURCE_DIR}/src/
    ${CMAKE_SOURCE_DIR}/include/
    )

target_link_libraries(mrcpp-tests
    mrcpp
    Catch
    )
