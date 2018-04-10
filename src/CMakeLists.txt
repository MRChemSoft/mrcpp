add_library(mrcpp STATIC)

add_subdirectory(mwcore)
add_subdirectory(mwutils)
add_subdirectory(mwtrees)
add_subdirectory(mwbuilders)
add_subdirectory(mwfunctions)
add_subdirectory(mwoperators)
add_subdirectory(mwparallel)

target_include_directories(mrcpp PRIVATE
    ${CMAKE_BINARY_DIR}
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/include
    )

install(TARGETS mrcpp ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})

