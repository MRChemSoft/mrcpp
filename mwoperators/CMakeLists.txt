
add_library(mwoperators STATIC
    PHCalculator.cpp
    ABGVCalculator.cpp
    ConvolutionOperator.cpp
    HelmholtzKernel.cpp
    PoissonKernel.cpp
    MWOperator.cpp
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

