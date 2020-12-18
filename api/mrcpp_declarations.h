/*
 * MRCPP, a numerical library based on multiresolution analysis and
 * the multiwavelet basis which provide low-scaling algorithms as well as
 * rigorous error control in numerical computations.
 * Copyright (C) 2020 Stig Rune Jensen, Jonas Juselius, Luca Frediani and contributors.
 *
 * This file is part of MRCPP.
 *
 * MRCPP is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * MRCPP is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with MRCPP.  If not, see <https://www.gnu.org/licenses/>.
 *
 * For information on the complete list of contributors to MRCPP, see:
 * <https://mrcpp.readthedocs.io/>
 */

#pragma once

#include <array>
#include <functional>
#include <set>
#include <vector>

namespace mrcpp {

class Timer;
class Printer;
template <int D> class Plotter;

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
class OperatorTree;

template <int D> class NodeAllocator;
template <int D> class ProjectedNodeAllocator;
template <int D> class GenNodeAllocator;
class OperatorNodeAllocator;

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
template <int D> class DerivativeOperator;
template <int D> class ABGVOperator;
template <int D> class PHOperator;

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

using OperatorTreeVector = std::vector<OperatorTree *>;
template <int D> using Coord = std::array<double, D>;
template <int D> using MWNodeVector = std::vector<MWNode<D> *>;

template <typename T, typename U> using FMap_ = std::function<T(U)>;
typedef FMap_<double, double> FMap;

} // namespace mrcpp
