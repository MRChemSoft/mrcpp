
add_library(mwcore STATIC 
    CrossCorrelation.cpp
    CrossCorrelationCache.cpp
    FilterCache.cpp
    GaussQuadrature.cpp
    InterpolatingBasis.cpp
    LegendreBasis.cpp
    MWFilter.cpp
    QuadratureCache.cpp
)

target_link_libraries(mwcore
    mrcpp
    mwanalytic
    mwtrees
    mwoperators
    mwbuilders
)

install(TARGETS mwcore
        ARCHIVE DESTINATION lib
        CONFIGURATIONS Debug-devel Release-devel
        )

