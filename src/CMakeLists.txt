add_library(mrcpp STATIC 
    BoundingBox.cpp
    BoysFunction.cpp
    CopyAdaptor.cpp
    CrossCorrelation.cpp
    CrossCorrelationCache.cpp
    CrossCorrelationCalculator.cpp
    FilterCache.cpp
    FunctionNode.cpp
    FunctionTree.cpp
    FunctionTreeVector.cpp
    GaussQuadrature.cpp
    GenNode.cpp
    grid.cpp
    HilbertPath.cpp
    InterpolatingBasis.cpp
    LegendreBasis.cpp
    LegendrePoly.cpp
    MWFilter.cpp
    MWNode.cpp
    MWTree.cpp
    project.cpp
    add.cpp
    multiply.cpp
    MathUtils.cpp
    MultiResolutionAnalysis.cpp
    NodeBox.cpp
    ObjectCache.cpp
    Polynomial.cpp
    ProjectedNode.cpp
    ProjectionCalculator.cpp
    QuadratureCache.cpp
    RepresentableFunction.cpp
    ScalingBasis.cpp
    SerialFunctionTree.cpp
    SerialOperatorTree.cpp
    SerialTree.cpp
    TreeBuilder.cpp
    TreeIterator.cpp
    parallel_mpi.cpp
)

add_subdirectory(gaussians)
add_subdirectory(mwutils)
add_subdirectory(mwoperators)

target_include_directories(mrcpp PRIVATE
    ${CMAKE_BINARY_DIR}
    ${CMAKE_SOURCE_DIR}/src/gaussians
    ${CMAKE_SOURCE_DIR}/src/mwutils
    ${CMAKE_SOURCE_DIR}/src/mwoperators
    )

install(TARGETS mrcpp ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})


list(APPEND headers_list 
    AdditionCalculator.h
    AnalyticAdaptor.h
    AnalyticFunction.h
    BoundingBox.h
    BoysFunction.h
    CopyAdaptor.h
    CrossCorrelation.h
    CrossCorrelationCache.h
    CrossCorrelationCalculator.h
    DefaultCalculator.h
    FilterCache.h
    FunctionNode.h
    FunctionTree.h
    FunctionTreeVector.h
    GaussQuadrature.h
    GenNode.h
    HilbertIterator.h
    HilbertPath.h
    InterpolatingBasis.h
    LebesgueIterator.h
    LegendreBasis.h
    LegendrePoly.h
    MWFilter.h
    MWNode.h
    MWTree.h
    MathUtils.h
    MultiResolutionAnalysis.h
    MultiplicationCalculator.h
    NodeBox.h
    NodeIndex.h
    ObjectCache.h
    Polynomial.h
    ProjectedNode.h
    ProjectionCalculator.h
    QuadratureCache.h
    RepresentableFunction.h
    ScalingBasis.h
    ScalingCache.h
    SerialFunctionTree.h
    SerialOperatorTree.h
    SerialTree.h
    TreeAdaptor.h
    TreeBuilder.h
    TreeCalculator.h
    TreeIterator.h
    WaveletAdaptor.h
    add.h
    constants.h
    eigen_disable_warnings.h
    grid.h
    macros.h
    mrcpp_declarations.h
    multiply.h
    parallel_mpi.h
    parallel_omp.h
    project.h
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

