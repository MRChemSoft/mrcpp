target_sources(mrcpp
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/ABGVOperator.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/ConvolutionOperator.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/HelmholtzKernel.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/HelmholtzOperator.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/MWOperator.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/OperatorStatistics.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/PHOperator.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/PoissonKernel.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/PoissonOperator.cpp
  )

list(APPEND operators_h
  ${CMAKE_CURRENT_SOURCE_DIR}/ABGVOperator.h
  ${CMAKE_CURRENT_SOURCE_DIR}/ConvolutionOperator.h
  ${CMAKE_CURRENT_SOURCE_DIR}/DerivativeConvolution.h
  ${CMAKE_CURRENT_SOURCE_DIR}/DerivativeKernel.h
  ${CMAKE_CURRENT_SOURCE_DIR}/DerivativeOperator.h
  ${CMAKE_CURRENT_SOURCE_DIR}/GreensKernel.h
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
  )

# Sets install directory for all the headers in the list
foreach(_header IN LISTS operators_h)
 file(
    COPY ${_header}
    DESTINATION ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}/operators
    )
 install(
   FILES ${_header}
   DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}/operators
   )
endforeach()

