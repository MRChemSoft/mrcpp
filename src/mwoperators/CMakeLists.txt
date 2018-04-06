target_sources(mrcpp PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/ABGVOperator.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/BandWidth.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/ConvolutionOperator.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/HelmholtzKernel.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/HelmholtzOperator.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/MWOperator.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/OperatorStatistics.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/PHOperator.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/PoissonKernel.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/PoissonOperator.cpp
)

list(APPEND mwoperators_h
    ABGVOperator.h
    BandWidth.h
    ConvolutionOperator.h
    DerivativeConvolution.h
    DerivativeKernel.h
    DerivativeOperator.h
    GreensKernel.h
    HelmholtzKernel.h
    HelmholtzOperator.h
    IdentityConvolution.h
    IdentityKernel.h
    MWOperator.h
    OperatorState.h
    OperatorStatistics.h
    PHOperator.h
    PoissonKernel.h
    PoissonOperator.h
    )                                                       
                                                                                 
# Sets install directory for all the headers in the list                         
foreach(_header ${mwoperators_h})                                                 
 file(
    COPY ${_header}
    DESTINATION ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}/src
    )
 install(
    FILES ${_header}
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}/src
    )
endforeach() 

