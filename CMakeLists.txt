
set(mrcpp_SOURCES
            parallel.cpp
            parallel.h
            constants.h
            macros.h
            MathUtils.cpp
            MathUtils.h
            BoundingBox.cpp
            BoundingBox.h
            NodeBox.cpp
            NodeBox.h
            MultiResolutionAnalysis.h
            MWFilter.cpp
            MWFilter.h
            HilbertPath.cpp
            HilbertPath.h
            NodeIndex.h
            TelePrompter.cpp
            TelePrompter.h
            GaussQuadrature.cpp
            GaussQuadrature.h
            ObjectCache.h
            FilterCache.cpp
            FilterCache.h
            ScalingCache.h
            QuadratureCache.cpp
            QuadratureCache.h
            RepresentableFunction.cpp
            RepresentableFunction.h
            LegendrePoly.cpp
            LegendrePoly.h
            Polynomial.cpp
            Polynomial.h
            ScalingBasis.h
            LegendreBasis.h
            InterpolatingBasis.cpp
            InterpolatingBasis.h
            MWTree.cpp
            MWTree.h
            MWNode.cpp
            MWNode.h
            FunctionTree.cpp
            FunctionTree.h
            TreeAllocator.cpp
            TreeAllocator.h
            FunctionNode.cpp
            FunctionNode.h
            ProjectedNode.cpp
            ProjectedNode.h
            GenNode.cpp
            GenNode.h
            Gaussian.cpp
            Gaussian.h
            GaussPoly.cpp
            GaussPoly.h
            GaussFunc.cpp
            GaussFunc.h
            GaussExp.cpp
            GaussExp.h
            TreeBuilder.cpp
            TreeBuilder.h
            MWProjector.h
            MWAdder.h
            MWMultiplier.h
            GridCleaner.h
            GridGenerator.h
            TreeAdaptor.h
            AnalyticAdaptor.h
            CopyAdaptor.h
            WaveletAdaptor.cpp
            WaveletAdaptor.h
            TreeCalculator.h
            ProjectionCalculator.cpp
            ProjectionCalculator.h
        )

add_library(mrcpp STATIC ${mrcpp_SOURCES})

install(TARGETS mrcpp
        ARCHIVE DESTINATION lib
        CONFIGURATIONS Debug-devel Release-devel
        )

