link_directories(${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR})

if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/mrcpp.cpp)
    add_executable(mrcpp-pilot mrcpp.cpp)
    target_include_directories(mrcpp-pilot PRIVATE
        ${CMAKE_BINARY_DIR}
        ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}
        )
    target_link_libraries(mrcpp-pilot
        mrcpp
        ${BLAS_LIBRARIES}
        )

endif()

