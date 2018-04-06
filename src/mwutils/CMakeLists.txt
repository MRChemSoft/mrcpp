target_sources(mrcpp PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/BoysFunction.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/MathUtils.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Plotter.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Printer.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Timer.cpp
)

list(APPEND mwutils_h
    ${CMAKE_CURRENT_SOURCE_DIR}/BoysFunction.h
    ${CMAKE_CURRENT_SOURCE_DIR}/MathUtils.h
    ${CMAKE_CURRENT_SOURCE_DIR}/Plotter.h
    ${CMAKE_CURRENT_SOURCE_DIR}/Printer.h
    ${CMAKE_CURRENT_SOURCE_DIR}/Timer.h
    )                                                       
                                                                                 
# Sets install directory for all the headers in the list                         
foreach(_header ${mwutils_h})                                                 
 file(
    COPY ${_header}
    DESTINATION ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}/src
    )
 install(
    FILES ${_header}
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}/src
    )
endforeach() 

