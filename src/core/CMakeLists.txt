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
  )

list(APPEND core_h
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
  )

# Sets install directory for all the headers in the list
foreach(_header IN LISTS core_h)
 install(
   FILES ${_header}
   DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}/core
   )
endforeach()
