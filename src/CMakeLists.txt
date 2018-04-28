add_library(mrcpp STATIC "")

add_subdirectory(core)
add_subdirectory(utils)
add_subdirectory(trees)
add_subdirectory(treebuilders)
add_subdirectory(functions)
add_subdirectory(operators)
add_subdirectory(mwparallel)

target_include_directories(mrcpp PRIVATE
    ${CMAKE_BINARY_DIR}
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/include
    )

install(TARGETS mrcpp ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})

