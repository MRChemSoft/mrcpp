target_sources(mrcpp
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/ABGVOperator.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/BSOperator.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/ConvolutionOperator.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/CartesianConvolution.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/DerivativeConvolution.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/HeatOperator.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/HelmholtzKernel.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/HelmholtzOperator.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/IdentityConvolution.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/MWOperator.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/OperatorStatistics.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/PHOperator.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/PoissonKernel.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/PoissonOperator.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/TimeEvolutionOperator.cpp
  )

get_filename_component(_dirname ${CMAKE_CURRENT_LIST_DIR} NAME)

list(APPEND ${_dirname}_h
  ${CMAKE_CURRENT_SOURCE_DIR}/ABGVOperator.h
  ${CMAKE_CURRENT_SOURCE_DIR}/CartesianConvolution.h
  ${CMAKE_CURRENT_SOURCE_DIR}/ConvolutionOperator.h
  ${CMAKE_CURRENT_SOURCE_DIR}/BSOperator.h
  ${CMAKE_CURRENT_SOURCE_DIR}/DerivativeConvolution.h
  ${CMAKE_CURRENT_SOURCE_DIR}/DerivativeKernel.h
  ${CMAKE_CURRENT_SOURCE_DIR}/DerivativeOperator.h
  ${CMAKE_CURRENT_SOURCE_DIR}/HeatKernel.h
  ${CMAKE_CURRENT_SOURCE_DIR}/HeatOperator.h
  ${CMAKE_CURRENT_SOURCE_DIR}/HelmholtzKernel.h
  ${CMAKE_CURRENT_SOURCE_DIR}/HelmholtzOperator.h
  ${CMAKE_CURRENT_SOURCE_DIR}/IdentityConvolution.h
  ${CMAKE_CURRENT_SOURCE_DIR}/IdentityKernel.h
  ${CMAKE_CURRENT_SOURCE_DIR}/MWOperator.h
  ${CMAKE_CURRENT_SOURCE_DIR}/OperatorState.h
  ${CMAKE_CURRENT_SOURCE_DIR}/OperatorStatistics.h
  ${CMAKE_CURRENT_SOURCE_DIR}/PHOperator.h
  ${CMAKE_CURRENT_SOURCE_DIR}/PoissonKernel.h
  ${CMAKE_CURRENT_SOURCE_DIR}/PoissonOperator.h
  ${CMAKE_CURRENT_SOURCE_DIR}/TimeEvolutionOperator.h
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
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}/operators
    )
endforeach()
