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

/**
 * @file PHOperator.cpp
 * @brief Implementation of a derivative operator assembled via PHCalculator.
 *
 * @details
 * This module builds a single-component multiwavelet operator that approximates
 * a spatial derivative of order 1 or 2. Construction proceeds by:
 *  1) creating a @ref PHCalculator tailored to the current scaling basis and
 *     requested derivative order,
 *  2) expanding to an @ref OperatorTree with a @ref TreeBuilder and a simple
 *     @ref BandWidthAdaptor (bandwidth = 1),
 *  3) transforming/caching the operator for efficient application.
 *
 * The class derives from @ref DerivativeOperator and uses the MRA’s root scale
 * by default. The operator is stored as a single separable component and exposed
 * through the common @ref MWOperator interface.
 */

#include "PHOperator.h"
#include "treebuilders/BandWidthAdaptor.h"
#include "treebuilders/PHCalculator.h"
#include "treebuilders/TreeBuilder.h"
#include "utils/Printer.h"
#include "utils/Timer.h"

namespace mrcpp {

/**
 * @brief Construct a PH-based derivative operator.
 *
 * @tparam D Spatial dimension (1, 2, or 3).
 * @param mra   MultiResolutionAnalysis defining the domain and basis.
 * @param order Derivative order (supported: 1 or 2).
 *
 * @details
 * Initializes the base @ref DerivativeOperator at the MRA root scale and
 * triggers internal assembly via @ref initialize(). The resulting expansion
 * contains a single operator block (rank-1 in the separable sense).
 */
template <int D>
PHOperator<D>::PHOperator(const MultiResolutionAnalysis<D> &mra, int order)
        : DerivativeOperator<D>(mra, mra.getRootScale(), -10) {
    this->order = order;
    initialize();
}

/**
 * @brief Assemble the operator tree for the requested derivative order.
 *
 * @details
 * - Creates a @ref PHCalculator using the MRA’s scaling basis and the stored
 *   derivative order.
 * - Uses a @ref BandWidthAdaptor with bandwidth 1 and the MRA’s maximum scale.
 * - Builds an @ref OperatorTree with @ref TreeBuilder, computes its squared
 *   norm, and prepares the node cache for application.
 * - Stores the built tree as a single raw term and initializes the operator
 *   expansion with @ref initOperExp(1).
 */
template <int D> void PHOperator<D>::initialize() {
    auto o_mra = this->getOperatorMRA();

    TreeBuilder<2> builder;

    auto &basis = this->MRA.getScalingBasis();
    PHCalculator calculator(basis, this->order);

    int bw = 1; // Operator bandwidth
    int max_scale = this->MRA.getMaxScale();
    BandWidthAdaptor adaptor(bw, max_scale);

    auto o_tree = std::make_unique<OperatorTree>(o_mra, MachineZero);
    builder.build(*o_tree, calculator, adaptor, -1);

    Timer trans_t;
    o_tree->calcSquareNorm();
    o_tree->setupOperNodeCache();
    print::time(10, "Time transform", trans_t);
    print::separator(10, ' ');

    this->raw_exp.push_back(std::move(o_tree));
    this->initOperExp(1);
}

/* Explicit template instantiations */
template class PHOperator<1>;
template class PHOperator<2>;
template class PHOperator<3>;

} // namespace mrcpp