include_directories(mwarithmetics)

set(mrcpp_SOURCES
    BoundingBox.cpp
    CrossCorrelation.cpp
    CrossCorrelationCache.cpp
    CrossCorrelationCalculator.cpp
    ConvolutionOperator.cpp
    FilterCache.cpp
    FunctionNode.cpp
    FunctionTree.cpp
    GaussExp.cpp
    GaussFunc.cpp
    GaussPoly.cpp
    GaussQuadrature.cpp
    Gaussian.cpp
    GenNode.cpp
    GridCleaner.cpp
    HelmholtzKernel.cpp
    HilbertPath.cpp
    InterpolatingBasis.cpp
    LegendreBasis.cpp
    LegendrePoly.cpp
    MWFilter.cpp
    MWNode.cpp
    MWTree.cpp
    MWOperator.cpp
    MathUtils.cpp
    NodeBox.cpp
    OperApplicationCalculator.cpp
    OperatorNode.cpp
    OperatorStatistics.cpp
    OperatorTree.cpp
    OperatorTreeVector.cpp
    PoissonKernel.cpp
    Polynomial.cpp
    ProjectedNode.cpp
    ProjectionCalculator.cpp
    QuadratureCache.cpp
    RepresentableFunction.cpp
    TelePrompter.cpp
    TreeAllocator.cpp
    TreeBuilder.cpp
    WaveletAdaptor.cpp
    parallel.cpp
)

add_library(mrcpp STATIC ${mrcpp_SOURCES})

install(TARGETS mrcpp
        ARCHIVE DESTINATION lib
        CONFIGURATIONS Debug-devel Release-devel
        )

