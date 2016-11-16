
add_library(mwoperators STATIC
    ABGVCalculator.cpp
    ConvolutionOperator.cpp
    CrossCorrelationCalculator.cpp
    HelmholtzKernel.cpp
    MWOperator.cpp
    PHCalculator.cpp
    PoissonKernel.cpp
)

target_link_libraries(mwoperators
    mrcpp
    mwanalytic
    mwcore 
    mwtrees
    mwbuilders
)

install(TARGETS mwoperators
        ARCHIVE DESTINATION lib
        CONFIGURATIONS Debug-devel Release-devel
        )

