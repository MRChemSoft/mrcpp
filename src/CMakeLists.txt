include(DefaultPathsMRCPP)

add_library(mrcpp STATIC 
    ABGVCalculator.cpp
    BoundingBox.cpp
    ConvolutionCalculator.cpp
    ConvolutionOperator.cpp
    CrossCorrelation.cpp
    CrossCorrelationCache.cpp
    CrossCorrelationCalculator.cpp
    DerivativeCalculator.cpp
    FilterCache.cpp
    FunctionNode.cpp
    FunctionTree.cpp
    GaussExp.cpp
    GaussFunc.cpp
    GaussPoly.cpp
    GaussQuadrature.cpp
    Gaussian.cpp
    HelmholtzKernel.cpp
    HilbertPath.cpp
    InterpolatingBasis.cpp
    LegendreBasis.cpp
    LegendrePoly.cpp
    MWFilter.cpp
    MWNode.cpp
    MWOperator.cpp
    MWTree.cpp
    MathUtils.cpp
    NodeBox.cpp
    OperatorNode.cpp
    OperatorStatistics.cpp
    OperatorTree.cpp
    PHCalculator.cpp
    Plotter.cpp
    PoissonKernel.cpp
    Polynomial.cpp
    ProjectedNode.cpp
    ProjectionCalculator.cpp
    QuadratureCache.cpp
    RepresentableFunction.cpp
    SerialFunctionTree.cpp
    SerialOperatorTree.cpp
    SerialTree.cpp
    Printer.cpp
    TreeBuilder.cpp
    parallel.cpp
)

install(TARGETS mrcpp
        ARCHIVE DESTINATION lib
        CONFIGURATIONS Debug-devel Release-devel
        )
