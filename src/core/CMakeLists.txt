target_sources(mrcpp
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/CrossCorrelation.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/CrossCorrelationCache.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/FilterCache.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/GaussQuadrature.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/InterpolatingBasis.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/LegendreBasis.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/MWFilter.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/ObjectCache.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/QuadratureCache.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/ScalingBasis.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/SchrodingerEvolution_CrossCorrelation.cpp
  )

get_filename_component(_dirname ${CMAKE_CURRENT_LIST_DIR} NAME)

list(APPEND ${_dirname}_h
  ${CMAKE_CURRENT_SOURCE_DIR}/CrossCorrelation.h
  ${CMAKE_CURRENT_SOURCE_DIR}/CrossCorrelationCache.h
  ${CMAKE_CURRENT_SOURCE_DIR}/FilterCache.h
  ${CMAKE_CURRENT_SOURCE_DIR}/GaussQuadrature.h
  ${CMAKE_CURRENT_SOURCE_DIR}/InterpolatingBasis.h
  ${CMAKE_CURRENT_SOURCE_DIR}/LegendreBasis.h
  ${CMAKE_CURRENT_SOURCE_DIR}/MWFilter.h
  ${CMAKE_CURRENT_SOURCE_DIR}/ObjectCache.h
  ${CMAKE_CURRENT_SOURCE_DIR}/QuadratureCache.h
  ${CMAKE_CURRENT_SOURCE_DIR}/ScalingBasis.h
  ${CMAKE_CURRENT_SOURCE_DIR}/ScalingCache.h
  ${CMAKE_CURRENT_SOURCE_DIR}/SchrodingerEvolution_CrossCorrelation.h
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
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}/core
    )
endforeach()
