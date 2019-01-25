target_sources(mrcpp
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/math_utils.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Plotter.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Printer.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Timer.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/mpi_utils.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/details.cpp
  )

list(APPEND utils_h
    ${CMAKE_CURRENT_SOURCE_DIR}/Plotter.h
    ${CMAKE_CURRENT_SOURCE_DIR}/Printer.h
    ${CMAKE_CURRENT_SOURCE_DIR}/Timer.h
    ${CMAKE_CURRENT_SOURCE_DIR}/mpi_utils.h
    ${CMAKE_CURRENT_SOURCE_DIR}/omp_utils.h
    ${CMAKE_CURRENT_SOURCE_DIR}/math_utils.h
    ${CMAKE_CURRENT_SOURCE_DIR}/details.h
    )

# Sets install directory for all the headers in the list
foreach(_header IN LISTS utils_h)
 install(
   FILES ${_header}
   DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}/utils
   )
endforeach()
