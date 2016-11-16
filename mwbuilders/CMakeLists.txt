
add_library(mwbuilders STATIC 
    ConvolutionCalculator.cpp
    DerivativeCalculator.cpp
    ProjectionCalculator.cpp
    TreeBuilder.cpp
)

target_link_libraries(mwbuilders
    mrcpp
    mwanalytic
    mwcore
    mwtrees
    mwoperators
)

install(TARGETS mwbuilders
        ARCHIVE DESTINATION lib
        CONFIGURATIONS Debug-devel Release-devel
        )

