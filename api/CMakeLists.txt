list(APPEND api_list
    MWFunctions
    MWOperators
    Gaussians
    Parallel
    Plotter
    Printer
    Timer
    )

# Sets install directory for all the headers in the list
foreach(_header ${api_list})
  file(
    COPY ${_header}
    DESTINATION ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}
    )
  install(
    FILES ${_header}
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}
    )
endforeach()

