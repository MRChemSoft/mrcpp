include(DefaultPathsMRCPP)
include_directories (${CMAKE_SOURCE_DIR}/src)

link_directories (${CMAKE_INSTALL_PREFIX}/lib)

add_executable(projection projection.cpp)
add_executable(addition addition.cpp)
add_executable(multiplication multiplication.cpp)
add_executable(tree_cleaner tree_cleaner.cpp)
add_executable(derivative derivative.cpp)
add_executable(poisson poisson.cpp)
add_executable(scf scf.cpp)

target_link_libraries(projection mrcpp ${BLAS_LIBRARIES})
target_link_libraries(addition mrcpp ${BLAS_LIBRARIES})
target_link_libraries(multiplication mrcpp ${BLAS_LIBRARIES})
target_link_libraries(tree_cleaner mrcpp ${BLAS_LIBRARIES})
target_link_libraries(derivative mrcpp ${BLAS_LIBRARIES})
target_link_libraries(poisson mrcpp ${BLAS_LIBRARIES})
target_link_libraries(scf mrcpp ${BLAS_LIBRARIES})

