target_sources(mrcpp PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/LegendrePoly.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Polynomial.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/RepresentableFunction.cpp
)

list(APPEND mwfunctions_h
    ${CMAKE_CURRENT_SOURCE_DIR}/AnalyticFunction.h
    ${CMAKE_CURRENT_SOURCE_DIR}/LegendrePoly.h
    ${CMAKE_CURRENT_SOURCE_DIR}/Polynomial.h
    ${CMAKE_CURRENT_SOURCE_DIR}/RepresentableFunction.h
    )                                                       
                                                                                 
# Sets install directory for all the headers in the list                         
foreach(_header ${mwfunctions_h})                                                 
 file(
    COPY ${_header}
    DESTINATION ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}/src
    )
 install(
    FILES ${_header}
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}/src
    )
endforeach() 

