
set(mrcpp_SOURCES
    AnalyticAdaptor.h
    BoundingBox.h
    CopyAdaptor.h
    DefaultCalculator.h
    FilterCache.h
    FunctionNode.h
    FunctionTree.h
    GaussExp.h
    GaussFunc.h
    GaussPoly.h
    GaussQuadrature.h
    Gaussian.h
    GenNode.h
    GridCleaner.h
    GridGenerator.h
    HilbertIterator.h
    HilbertPath.h
    InterpolatingBasis.h
    LebesgueIterator.h
    LegendreBasis.h
    LegendrePoly.h
    MWAdder.h
    MWFilter.h
    MWMultiplier.h
    MWNode.h
    MWProjector.h
    MWTree.h
    MathUtils.h
    MultiResolutionAnalysis.h
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
    TelePrompter.h
    TreeAdaptor.h
    TreeAllocator.h
    TreeBuilder.h
    TreeCalculator.h
    TreeIterator.h
    WaveletAdaptor.h
    parallel.h
    BoundingBox.cpp
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
    HilbertPath.cpp
    InterpolatingBasis.cpp
    LegendreBasis.cpp
    LegendrePoly.cpp
    MWFilter.cpp
    MWNode.cpp
    MWTree.cpp
    MathUtils.cpp
    NodeBox.cpp
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

