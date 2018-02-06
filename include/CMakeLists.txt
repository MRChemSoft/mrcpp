list(APPEND headers_list 
    ABGVCalculator.h
    ABGVOperator.h
    AdditionCalculator.h
    AnalyticAdaptor.h
    AnalyticFunction.h
    BandWidth.h
    BandWidthAdaptor.h
    BoundingBox.h
    BoysFunction.h
    ConvolutionCalculator.h
    ConvolutionOperator.h
    CopyAdaptor.h
    CrossCorrelation.h
    CrossCorrelationCache.h
    CrossCorrelationCalculator.h
    DefaultCalculator.h
    DerivativeCalculator.h
    DerivativeConvolution.h
    DerivativeKernel.h
    DerivativeOperator.h
    FilterCache.h
    FunctionNode.h
    FunctionTree.h
    FunctionTreeVector.h
    GaussExp.h
    GaussFunc.h
    GaussPoly.h
    GaussQuadrature.h
    Gaussian.h
    GenNode.h
    GreensKernel.h
    HelmholtzKernel.h
    HelmholtzOperator.h
    HilbertIterator.h
    HilbertPath.h
    HydrogenicFunction.h
    IdentityConvolution.h
    IdentityKernel.h
    InterpolatingBasis.h
    LebesgueIterator.h
    LegendreBasis.h
    LegendrePoly.h
    MWFilter.h
    MWNode.h
    MWOperator.h
    MWTree.h
    MathUtils.h
    MultiResolutionAnalysis.h
    MultiplicationCalculator.h
    NodeBox.h
    NodeIndex.h
    ObjectCache.h
    OperatorAdaptor.h
    OperatorNode.h
    OperatorState.h
    OperatorStatistics.h
    OperatorTree.h
    PHCalculator.h
    PHOperator.h
    Plotter.h
    PoissonKernel.h
    PoissonOperator.h
    Polynomial.h
    Printer.h
    ProjectedNode.h
    ProjectionCalculator.h
    QuadratureCache.h
    RepresentableFunction.h
    ScalingBasis.h
    ScalingCache.h
    SerialFunctionTree.h
    SerialOperatorTree.h
    SerialTree.h
    Timer.h
    TreeAdaptor.h
    TreeBuilder.h
    TreeCalculator.h
    TreeIterator.h
    WaveletAdaptor.h
    add.h
    apply.h
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
   file(COPY ${_header} DESTINATION ${CMAKE_BINARY_DIR}/include)
   install(FILES ${_header} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
endforeach() 

list(APPEND api_list
    MRCPP/MWFunctions
    MRCPP/MWOperators
    MRCPP/Gaussians
    MRCPP/Parallel
    MRCPP/Plotter
    MRCPP/Printer
    MRCPP/Timer
    )

# Sets install directory for all the headers in the list                         
foreach(_header ${api_list})                                                 
   file(COPY ${_header} DESTINATION ${CMAKE_BINARY_DIR}/include/MRCPP)
   install(FILES ${_header} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/MRCPP)
endforeach() 
