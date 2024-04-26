/*
 * MRCPP, a numerical library based on multiresolution analysis and
 * the multiwavelet basis which provide low-scaling algorithms as well as
 * rigorous error control in numerical computations.
 * Copyright (C) 2021 Stig Rune Jensen, Jonas Juselius, Luca Frediani and contributors.
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
template <int D, typename T = double> class Plotter;

template <int D, typename T = double> class Gaussian;
template <int D, typename T = double> class GaussFunc;
template <int D, typename T = double> class GaussPoly;
template <int D, typename T = double> class GaussExp;

template <int D> class BoundingBox;
template <int D, typename T = double> class NodeBox;
template <int D> class NodeIndex;
template <int D> class NodeIndexComp;

template <typename T = double> class SharedMemory;
class ScalingBasis;
class LegendreBasis;
class InterpolatingBasis;

template <int D, typename T = double> class RepresentableFunction;
template <int D> class MultiResolutionAnalysis;

template <int D, typename T = double> class MWTree;
template <int D, typename T = double> class FunctionTree;
class OperatorTree;

template <int D, typename T = double> class NodeAllocator;

template <int D, typename T = double> class MWNode;
template <int D, typename T = double> class FunctionNode;
class OperatorNode;

template <int D> class IdentityConvolution;
template <int D> class DerivativeConvolution;
template <int D> class ConvolutionOperator;
class CartesianConvolution;
class PoissonOperator;
class HelmholtzOperator;
template <int D> class DerivativeOperator;
template <int D> class ABGVOperator;
template <int D> class PHOperator;

template <int D> class IdentityKernel;
template <int D> class DerivativeKernel;
class PoissonKernel;
class HelmholtzKernel;

template <int D, typename T = double> class TreeBuilder;
template <int D, typename T = double> class TreeCalculator;
template <int D, typename T = double> class DefaultCalculator;
template <int D, typename T = double> class ProjectionCalculator;
template <int D, typename T = double> class AdditionCalculator;
template <int D, typename T = double> class MultiplicationCalculator;
template <int D, typename T = double> class ConvolutionCalculator;
template <int D, typename T = double> class DerivativeCalculator;
class CrossCorrelationCalculator;

template <int D, typename T = double> class TreeAdaptor;
template <int D, typename T = double> class AnalyticAdaptor;
template <int D, typename T = double> class WaveletAdaptor;
template <int D, typename T = double> class CopyAdaptor;

template <int D, typename T = double> class TreeIterator;
template <int D, typename T = double> class IteratorNode;

class BandWidth;
template <int D, typename T = double> class OperatorState;

template <int D> using Coord = std::array<double, D>;
template <int D, typename T = double> using MWNodeVector = std::vector<MWNode<D, T> *>;

template <typename T = double, typename U = double> using FMap = std::function<T(U)>;

} // namespace mrcpp
