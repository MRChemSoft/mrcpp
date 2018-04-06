#pragma once

#include <vector>
#include <set>

namespace mrcpp {

class Timer;
class Printer;
template<int D> class Plotter;

template <int D> class Gaussian;
template <int D> class GaussFunc;
template <int D> class GaussPoly;
template <int D> class GaussExp;

template <int D> class BoundingBox;
template <int D> class NodeBox;
template <int D> class NodeIndex;
template <int D> class NodeIndexComp;

class SharedMemory;
class ScalingBasis;
class LegendreBasis;
class InterpolatingBasis;

template <int D> class RepresentableFunction;
template <int D> class MultiResolutionAnalysis;

template <int D> class MWTree;
template <int D> class FunctionTree;
template <int D> class FunctionTreeVector;
class OperatorTree;
template <int D> class SerialTree;
template <int D> class SerialFunctionTree;
class SerialOperatorTree;

template <int D> class MWNode;
template <int D> class FunctionNode;
template <int D> class ProjectedNode;
template <int D> class GenNode;
class OperatorNode;

template <int D> class IdentityConvolution;
template <int D> class DerivativeConvolution;
template <int D> class ConvolutionOperator;
class PoissonOperator;
class HelmholtzOperator;
template<int D> class DerivativeOperator;
template<int D> class ABGVOperator;
template<int D> class PHOperator;

class GreensKernel;
class IdentityKernel;
class DerivativeKernel;
class PoissonKernel;
class HelmholtzKernel;

template <int D> class TreeBuilder;
template <int D> class TreeCalculator;
template <int D> class DefaultCalculator;
template <int D> class ProjectionCalculator;
template <int D> class AdditionCalculator;
template <int D> class MultiplicationCalculator;
template <int D> class ConvolutionCalculator;
template <int D> class DerivativeCalculator;
class CrossCorrelationCalculator;

template <int D> class TreeAdaptor;
template <int D> class AnalyticAdaptor;
template <int D> class WaveletAdaptor;
template <int D> class CopyAdaptor;

template <int D> class TreeIterator;
template <int D> class LebesgueIterator;
template <int D> class HilbertIterator;
template <int D> class HilbertPath;
template <int D> class IteratorNode;

class BandWidth;
template <int D> class OperatorState;

#define OperatorTreeVector std::vector<OperatorTree *>
#define MWNodeVector std::vector<MWNode<D> *>
#define NodeIndexSet std::set<const NodeIndex<D> *, NodeIndexComp<D> >

}
