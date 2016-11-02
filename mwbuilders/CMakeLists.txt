
add_library(mwbuilders STATIC 
    CrossCorrelationCalculator.cpp
    DerivativeCalculator.cpp
    GridCleaner.cpp
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

