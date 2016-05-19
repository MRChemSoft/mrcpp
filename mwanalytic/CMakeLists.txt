
add_library(mwanalytic STATIC 
    GaussExp.cpp
    GaussFunc.cpp
    GaussPoly.cpp
    Gaussian.cpp
    LegendrePoly.cpp
    Polynomial.cpp
)

target_link_libraries(mwanalytic
    mrcpp
    mwcore
    mwtrees
    mwoperators
    mwbuilders
)

install(TARGETS mwanalytic
        ARCHIVE DESTINATION lib
        CONFIGURATIONS Debug-devel Release-devel
        )

