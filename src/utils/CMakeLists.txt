target_sources(mrcpp-objlib PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/MathUtils.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Plotter.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Printer.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Timer.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/mrcpp_mpi.cpp
)

list(APPEND utils_h
    ${CMAKE_CURRENT_SOURCE_DIR}/MathUtils.h
    ${CMAKE_CURRENT_SOURCE_DIR}/Plotter.h
    ${CMAKE_CURRENT_SOURCE_DIR}/Printer.h
    ${CMAKE_CURRENT_SOURCE_DIR}/Timer.h
    ${CMAKE_CURRENT_SOURCE_DIR}/mrcpp_mpi.h
    ${CMAKE_CURRENT_SOURCE_DIR}/mrcpp_omp.h
    )

# Sets install directory for all the headers in the list
foreach(_header IN LISTS utils_h)
 file(
    COPY ${_header}
    DESTINATION ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}/utils
    )
 install(
    FILES ${_header}
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}/utils
    )
endforeach()
