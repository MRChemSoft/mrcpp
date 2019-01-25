list(APPEND include_list
    constants.h
    eigen_disable_warnings.h
    macros.h
    mrcpp_declarations.h
    )

# Sets install directory for all the headers in the list
foreach(_header IN LISTS include_list)
  file(
    COPY ${_header}
    DESTINATION ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}
    )
  install(
    FILES ${_header}
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}
    )
endforeach()
