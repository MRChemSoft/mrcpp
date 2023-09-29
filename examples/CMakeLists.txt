list(APPEND _examples
    tree_cleaner
    projection
    derivative
    addition
    multiplication
    poisson
    scf
    mpi_shared_tree
    mpi_send_tree
    mpi_matrix
	schrodinger_semigroup1d
    )

foreach(_example IN LISTS _examples)
    add_executable(${_example} ${_example}.cpp)
    target_link_libraries(${_example}
        mrcpp
        ${BLAS_LIBRARIES}
        )
endforeach()
