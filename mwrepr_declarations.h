#ifndef MWREPR_DECLARATIONS_H_
#define MWREPR_DECLARATIONS_H_

#include <vector>
#include <set>

template <int D> class BoundingBox;
template <int D> class NodeBox;
template <int D> class NodeIndex;
template <int D> class NodeIndexComp;

template <int D> class RepresentableFunction;
template <int D> class MultiResolutionAnalysis;

template <int D> class MWTree;
template <int D> class MWNode;

template <int D> class FunctionTree;
template <int D> class FunctionNode;
template <int D> class ProjectedNode;
template <int D> class GenNode;

template <int D> class TreeBuilder;
template <int D> class MWProjector;
template <int D> class MWAdder;
template <int D> class MWMultiplier;
template <int D> class MWOperator;
template <int D> class GridGenerator;
template <int D> class GridCleaner;

template <int D> class TreeCalculator;
template <int D> class DefaultCalculator;
template <int D> class ProjectionCalculator;
template <int D> class AdditionCalculator;
template <int D> class MultiplicationCalculator;
template <int D> class IdentityCalculator;

template <int D> class TreeAdaptor;
template <int D> class AnalyticAdaptor;
template <int D> class WaveletAdaptor;
template <int D> class CopyAdaptor;

template <int D> class TreeIterator;
template <int D> class LebesgueIterator;
template <int D> class HilbertIterator;
template <int D> class HilbertPath;

#define MWNodeVector std::vector<MWNode<D> *>
#define NodeIndexSet std::set<const NodeIndex<D> *, NodeIndexComp<D> >

#endif /* MWREPR_DECLARATIONS_H_*/
