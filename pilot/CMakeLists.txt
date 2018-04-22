if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/mrcpp.cpp)
    add_executable(mrcpp-pilot mrcpp.cpp)
    target_link_libraries(mrcpp-pilot
        mrcpp
        ${BLAS_LIBRARIES}
        )
endif()
