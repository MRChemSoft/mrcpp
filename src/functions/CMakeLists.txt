target_sources(mrcpp-objlib PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/BoysFunction.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/LegendrePoly.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Polynomial.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/RepresentableFunction.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/GaussExp.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/GaussFunc.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/GaussPoly.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Gaussian.cpp
)

list(APPEND functions_h
    ${CMAKE_CURRENT_SOURCE_DIR}/BoysFunction.h
    ${CMAKE_CURRENT_SOURCE_DIR}/AnalyticFunction.h
    ${CMAKE_CURRENT_SOURCE_DIR}/LegendrePoly.h
    ${CMAKE_CURRENT_SOURCE_DIR}/Polynomial.h
    ${CMAKE_CURRENT_SOURCE_DIR}/RepresentableFunction.h
    ${CMAKE_CURRENT_SOURCE_DIR}/GaussExp.h
    ${CMAKE_CURRENT_SOURCE_DIR}/GaussFunc.h
    ${CMAKE_CURRENT_SOURCE_DIR}/GaussPoly.h
    ${CMAKE_CURRENT_SOURCE_DIR}/Gaussian.h
    )

# Sets install directory for all the headers in the list
foreach(_header IN LISTS functions_h)
 file(
    COPY ${_header}
    DESTINATION ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}/functions
    )
 install(
    FILES ${_header}
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}/functions
    )
endforeach()
