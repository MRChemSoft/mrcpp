include(DefaultPathsMRCPP)
include_directories (${CMAKE_SOURCE_DIR}/src)

link_directories (${CMAKE_INSTALL_PREFIX}/lib)

set (MRCPP_EXECUTABLE ${CMAKE_BINARY_DIR}/pilot/mrcpp-pilot.x)
set (MRCPP_DRIVER_DEBUG ON)

if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/mrcpp.cpp)
    add_executable(mrcpp-pilot.x mrcpp.cpp)
    target_link_libraries(mrcpp-pilot.x mrcpp ${BLAS_LIBRARIES})
endif()

