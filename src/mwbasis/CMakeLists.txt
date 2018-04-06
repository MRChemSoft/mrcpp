target_sources(mrcpp PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/InterpolatingBasis.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/LegendreBasis.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/LegendrePoly.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Polynomial.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/ScalingBasis.cpp
)

list(APPEND mwbasis_h
    ${CMAKE_CURRENT_SOURCE_DIR}/InterpolatingBasis.h
    ${CMAKE_CURRENT_SOURCE_DIR}/LegendreBasis.h
    ${CMAKE_CURRENT_SOURCE_DIR}/LegendrePoly.h
    ${CMAKE_CURRENT_SOURCE_DIR}/Polynomial.h
    ${CMAKE_CURRENT_SOURCE_DIR}/ScalingBasis.h
    )                                                       
                                                                                 
# Sets install directory for all the headers in the list                         
foreach(_header ${mwbasis_h})                                                 
 file(
    COPY ${_header}
    DESTINATION ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}/src
    )
 install(
    FILES ${_header}
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}/src
    )
endforeach() 

