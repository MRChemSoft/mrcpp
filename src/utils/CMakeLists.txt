target_sources(mrcpp
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/details.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/math_utils.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/mpi_utils.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/omp_utils.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/periodic_utils.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Plotter.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Printer.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Timer.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/tree_utils.cpp
  )

get_filename_component(_dirname ${CMAKE_CURRENT_LIST_DIR} NAME)

list(APPEND ${_dirname}_h
  ${CMAKE_CURRENT_SOURCE_DIR}/details.h
  ${CMAKE_CURRENT_SOURCE_DIR}/math_utils.h
  ${CMAKE_CURRENT_SOURCE_DIR}/mpi_utils.h
  ${CMAKE_CURRENT_SOURCE_DIR}/omp_utils.h
  ${CMAKE_CURRENT_SOURCE_DIR}/periodic_utils.h
  ${CMAKE_CURRENT_SOURCE_DIR}/Plotter.h
  ${CMAKE_CURRENT_SOURCE_DIR}/Printer.h
  ${CMAKE_CURRENT_SOURCE_DIR}/Timer.h
  ${CMAKE_CURRENT_SOURCE_DIR}/tree_utils.h
  )

# Sets install directory for all the headers in the list
file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}/${_dirname})
foreach(_header IN LISTS ${_dirname}_h)
  configure_file(
    ${_header}
    ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}/${_dirname}
    COPYONLY
    )
  install(
    FILES ${_header}
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}/utils
    )
endforeach()
