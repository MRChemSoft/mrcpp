
add_library(mwtrees STATIC 
    FunctionNode.cpp
    FunctionTree.cpp
    MWNode.cpp
    MWTree.cpp
    OperatorNode.cpp
    OperatorTree.cpp
    OperatorTreeVector.cpp
    ProjectedNode.cpp
)

target_link_libraries(mwtrees
    mrcpp
    mwanalytic
    mwcore
    mwoperators
    mwbuilders
)

install(TARGETS mwtrees
        ARCHIVE DESTINATION lib
        CONFIGURATIONS Debug-devel Release-devel
        )

