link_directories(${CMAKE_INSTALL_PREFIX}/lib)

if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/mrcpp.cpp)
    add_executable(mrcpp-pilot mrcpp.cpp)
    target_include_directories(mrcpp-pilot PRIVATE
        ${CMAKE_BINARY_DIR}
        ${CMAKE_SOURCE_DIR}/src
        )
    target_link_libraries(mrcpp-pilot
        mrcpp
        ${BLAS_LIBRARIES}
        )

endif()

