list(APPEND api_h
  ${CMAKE_CURRENT_LIST_DIR}/Gaussians
  ${CMAKE_CURRENT_LIST_DIR}/MWFunctions
  ${CMAKE_CURRENT_LIST_DIR}/MWOperators
  ${CMAKE_CURRENT_LIST_DIR}/Parallel
  ${CMAKE_CURRENT_LIST_DIR}/Plotter
  ${CMAKE_CURRENT_LIST_DIR}/Printer
  ${CMAKE_CURRENT_LIST_DIR}/Timer
  ${CMAKE_CURRENT_LIST_DIR}/constants.h
  ${CMAKE_CURRENT_LIST_DIR}/macros.h
  ${CMAKE_CURRENT_LIST_DIR}/mrcpp_declarations.h
  )

# Sets install directory for all the headers in the list
foreach(_header IN LISTS api_h)
  configure_file(
    ${_header}
    ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}
    COPYONLY
    )
endforeach()

# Append to mrcpp_headers list.
# These will be installed as the public headers of the target
list(APPEND mrcpp_headers "${api_h}")
set(mrcpp_headers "${mrcpp_headers}" PARENT_SCOPE)
