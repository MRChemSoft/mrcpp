add_library(mrcpp STATIC 
    BoundingBox.cpp
    FunctionNode.cpp
    FunctionTree.cpp
    FunctionTreeVector.cpp
    GenNode.cpp
    HilbertPath.cpp
    MWNode.cpp
    MWTree.cpp
    MultiResolutionAnalysis.cpp
    NodeBox.cpp
    ProjectedNode.cpp
    RepresentableFunction.cpp
    SerialFunctionTree.cpp
    SerialOperatorTree.cpp
    SerialTree.cpp
    TreeIterator.cpp
    parallel_mpi.cpp
)

add_subdirectory(gaussians)
add_subdirectory(mwcore)
add_subdirectory(mwutils)
add_subdirectory(mwoperators)
add_subdirectory(mwbuilders)

target_include_directories(mrcpp PRIVATE
    ${CMAKE_BINARY_DIR}
    ${CMAKE_SOURCE_DIR}/src/gaussians
    ${CMAKE_SOURCE_DIR}/src/mwcore
    ${CMAKE_SOURCE_DIR}/src/mwbasis
    ${CMAKE_SOURCE_DIR}/src/mwutils
    ${CMAKE_SOURCE_DIR}/src/mwoperators
    ${CMAKE_SOURCE_DIR}/src/mwbuilders
    )

install(TARGETS mrcpp ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})


list(APPEND headers_list 
    AnalyticFunction.h
    BoundingBox.h
    FunctionNode.h
    FunctionTree.h
    FunctionTreeVector.h
    GenNode.h
    HilbertIterator.h
    HilbertPath.h
    LebesgueIterator.h
    MWNode.h
    MWTree.h
    MultiResolutionAnalysis.h
    NodeBox.h
    NodeIndex.h
    ProjectedNode.h
    RepresentableFunction.h
    SerialFunctionTree.h
    SerialOperatorTree.h
    SerialTree.h
    TreeIterator.h
    constants.h
    eigen_disable_warnings.h
    macros.h
    mrcpp_declarations.h
    parallel_mpi.h
    parallel_omp.h
    )                                                       
                                                                                 
# Sets install directory for all the headers in the list                         
foreach(_header ${headers_list})                                                 
 file(
    COPY ${_header}
    DESTINATION ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}/src
    )
 install(
    FILES ${_header}
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}/src
    )
endforeach() 

