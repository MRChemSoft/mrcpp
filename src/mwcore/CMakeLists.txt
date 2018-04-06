target_sources(mrcpp PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/CrossCorrelation.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/CrossCorrelationCache.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/FilterCache.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/GaussQuadrature.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/MWFilter.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/ObjectCache.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/QuadratureCache.cpp
)

list(APPEND mwcore_h
    ${CMAKE_CURRENT_SOURCE_DIR}/CrossCorrelation.h
    ${CMAKE_CURRENT_SOURCE_DIR}/CrossCorrelationCache.h
    ${CMAKE_CURRENT_SOURCE_DIR}/FilterCache.h
    ${CMAKE_CURRENT_SOURCE_DIR}/GaussQuadrature.h
    ${CMAKE_CURRENT_SOURCE_DIR}/MWFilter.h
    ${CMAKE_CURRENT_SOURCE_DIR}/ObjectCache.h
    ${CMAKE_CURRENT_SOURCE_DIR}/QuadratureCache.h
    ${CMAKE_CURRENT_SOURCE_DIR}/ScalingCache.h
    )                                                       
                                                                                 
# Sets install directory for all the headers in the list                         
foreach(_header ${mwcore_h})                                                 
 file(
    COPY ${_header}
    DESTINATION ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}/src
    )
 install(
    FILES ${_header}
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}/src
    )
endforeach() 

