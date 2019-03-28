target_sources(mrcpp
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/BoysFunction.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/function_utils.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/LegendrePoly.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Polynomial.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/RepresentableFunction.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/GaussExp.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/GaussFunc.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/GaussPoly.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Gaussian.cpp
  )

get_filename_component(_dirname ${CMAKE_CURRENT_LIST_DIR} NAME)

list(APPEND ${_dirname}_h
  ${CMAKE_CURRENT_SOURCE_DIR}/AnalyticFunction.h
  ${CMAKE_CURRENT_SOURCE_DIR}/BoysFunction.h
  ${CMAKE_CURRENT_SOURCE_DIR}/function_utils.h
  ${CMAKE_CURRENT_SOURCE_DIR}/LegendrePoly.h
  ${CMAKE_CURRENT_SOURCE_DIR}/Polynomial.h
  ${CMAKE_CURRENT_SOURCE_DIR}/RepresentableFunction.h
  ${CMAKE_CURRENT_SOURCE_DIR}/GaussExp.h
  ${CMAKE_CURRENT_SOURCE_DIR}/GaussFunc.h
  ${CMAKE_CURRENT_SOURCE_DIR}/GaussPoly.h
  ${CMAKE_CURRENT_SOURCE_DIR}/Gaussian.h
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
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}/functions
    )
endforeach()
