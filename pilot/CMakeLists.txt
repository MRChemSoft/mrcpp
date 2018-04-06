if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/mrcpp.cpp)
    add_executable(mrcpp-pilot.x mrcpp.cpp)
    target_include_directories(mrcpp-pilot.x PRIVATE
        ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}
        )
    target_link_libraries(mrcpp-pilot.x
        mrcpp
        ${BLAS_LIBRARIES}
        )
endif()

