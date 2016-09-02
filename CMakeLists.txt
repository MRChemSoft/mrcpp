include_directories(mwbuilders)
include_directories(mwoperators)
include_directories(mwtrees)
include_directories(mwcore)
include_directories(mwanalytic)

add_subdirectory(mwbuilders)
add_subdirectory(mwoperators)
add_subdirectory(mwtrees)
add_subdirectory(mwcore)
add_subdirectory(mwanalytic)

add_library(mrcpp STATIC 
    BoundingBox.cpp
    HilbertPath.cpp
    MathUtils.cpp
    NodeBox.cpp
    OperatorStatistics.cpp
    RepresentableFunction.cpp
    TelePrompter.cpp
    SerialTree.cpp
    parallel.cpp
    Plot.cpp
    MPI.h
    MPI.cpp
)

target_link_libraries(mrcpp
    mwanalytic
    mwcore 
    mwtrees
    mwoperators
    mwbuilders
    rt
)

install(TARGETS mrcpp
        ARCHIVE DESTINATION lib
        CONFIGURATIONS Debug-devel Release-devel
        )

