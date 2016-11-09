
add_library(mwbuilders STATIC 
    MWFDCalculator.cpp
    FiniteDifferenceOperator.cpp
    HarrisonDerivativeCalculator.cpp
    CrossCorrelationCalculator.cpp
    DerivativeCalculator.cpp
    OperApplicationCalculator.cpp
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

