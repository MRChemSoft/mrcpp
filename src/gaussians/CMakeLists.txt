target_sources(mrcpp PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/GaussExp.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/GaussFunc.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/GaussPoly.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Gaussian.cpp
)

list(APPEND gaussians_h
    ${CMAKE_CURRENT_SOURCE_DIR}/GaussExp.h
    ${CMAKE_CURRENT_SOURCE_DIR}/GaussFunc.h
    ${CMAKE_CURRENT_SOURCE_DIR}/GaussPoly.h
    ${CMAKE_CURRENT_SOURCE_DIR}/Gaussian.h
    )                                                       
                                                                                 
# Sets install directory for all the headers in the list                         
foreach(_header ${gaussians_h})                                                 
 file(
    COPY ${_header}
    DESTINATION ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}/src
    )
 install(
    FILES ${_header}
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}/src
    )
endforeach() 

