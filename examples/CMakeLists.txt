list(APPEND _examples
    tree_cleaner
    projection
    derivative
    addition
    poisson
    scf
    mpi_shared_tree
    mpi_send_tree
    mpi_matrix
    )

foreach(_example ${_examples})
    add_executable(${_example}.x ${_example}.cpp)
    target_include_directories(${_example}.x PRIVATE
        ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}
        )
    target_link_libraries(${_example}.x
        mrcpp
        ${BLAS_LIBRARIES}
        )
endforeach()

