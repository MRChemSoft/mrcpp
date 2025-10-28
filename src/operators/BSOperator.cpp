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
 * @file BSOperator.cpp
 * @brief Assembly of a smooth multiresolution derivative operator (“BS” operator).
 *
 * ## Purpose
 * Build and cache a derivative operator in the multiresolution basis of a given
 * #mrcpp::MultiResolutionAnalysis. This operator is intended for **smooth**
 * functions; for non-smooth or discontinuous data, prefer #mrcpp::ABGVOperator.
 *
 * ## What happens here
 * 1. The constructor stores the requested derivative **order** (1, 2, or 3) and
 *    delegates to `initialize()`.
 * 2. `initialize()`:
 *    - chooses a (small) **bandwidth** (nearest-neighbor coupling, `bw = 1`),
 *    - creates a #mrcpp::BSCalculator that provides the local operator blocks
 *      in the MRA scaling basis for the selected derivative order,
 *    - wraps assembly with a #mrcpp::BandWidthAdaptor to enforce sparsity
 *      across all scales,
 *    - uses #mrcpp::TreeBuilder to assemble an #mrcpp::OperatorTree,
 *    - finalizes the operator (computes norms, builds per-node caches),
 *    - registers the resulting tree in the base #mrcpp::DerivativeOperator,
 *      and initializes the internal operator expansion for fast application.
 *
 * ## Notes
 * - The chosen bandwidth (`bw = 1`) yields a compact stencil (nearest neighbors).
 * - `calcSquareNorm()` is invoked once to precompute norms; this can aid later
 *   conditioning/thresholding steps that use these norms.
 * - `setupOperNodeCache()` prepares per-node data needed for efficient
 *   application of the operator during transforms/apply calls.
 *
 * ## Performance/usage
 * After construction, applying the operator to MR coefficient vectors is cheap
 * and can be repeated many times. The build cost is paid once per (MRA, order).
 */

#include "BSOperator.h"
#include "treebuilders/BSCalculator.h"
#include "treebuilders/BandWidthAdaptor.h"
#include "treebuilders/TreeBuilder.h"
#include "utils/Printer.h"
#include "utils/Timer.h"

namespace mrcpp {

/**
 * @brief Construct a smooth (“BS”) multiresolution derivative operator.
 *
 * @param mra   Multiresolution analysis that defines the domain, basis, and scales.
 * @param order Derivative order (supported: 1, 2, or 3).
 *
 * The operator is anchored at the MRA’s **root scale** (via the base
 * #mrcpp::DerivativeOperator constructor) and immediately assembled by
 * calling `initialize()`. The internal representation is stored as an
 * #mrcpp::OperatorTree and cached for fast application.
 */
template <int D>
BSOperator<D>::BSOperator(const MultiResolutionAnalysis<D> &mra, int order)
        : DerivativeOperator<D>(mra, mra.getRootScale()) {
    this->order = order;
    initialize();
}

/**
 * @brief Build and cache the “BS” derivative operator.
 *
 * **Assembly pipeline**
 * 1. Select operator bandwidth `bw = 1` (nearest-neighbor coupling).
 * 2. Query the operator MRA (`getOperatorMRA()`), which carries the scaling
 *    basis and max scale.
 * 3. Instantiate:
 *    - #mrcpp::BSCalculator with the scaling basis and the requested derivative
 *      order (generates local operator blocks),
 *    - #mrcpp::BandWidthAdaptor with `(bw, maxScale)` to enforce sparsity,
 *    - #mrcpp::TreeBuilder to assemble the global #mrcpp::OperatorTree.
 * 4. Build into a fresh `OperatorTree(oper_mra, MachineZero)`:
 *    - `MachineZero` is used as a numerical floor for tree entries.
 * 5. Finalize:
 *    - `calcSquareNorm()` precomputes norms (useful for later compression/metrics),
 *    - `setupOperNodeCache()` creates per-node caches for fast application.
 * 6. Store the assembled tree in `raw_exp` (owned by the base class) and call
 *    `initOperExp(1)` to finalize the expansion with a single raw operator.
 */
template <int D> void BSOperator<D>::initialize() {
    int bw = 1; // Operator bandwidth: nearest-neighbor coupling
    auto oper_mra = this->getOperatorMRA();

    TreeBuilder<2> builder;                               // 2: binary tree arity in 1D blocks
    BSCalculator    calculator(oper_mra.getScalingBasis(), this->order);
    BandWidthAdaptor adaptor(bw, oper_mra.getMaxScale()); // enforce sparsity across scales

    // Assemble the operator tree with numerical floor MachineZero
    auto o_tree = std::make_unique<OperatorTree>(oper_mra, MachineZero);
    builder.build(*o_tree, calculator, adaptor, -1 /* all levels */);

    // Finalize and cache per-node data for fast application
    Timer trans_t;
    o_tree->calcSquareNorm();          // precompute norms (once)
    o_tree->setupOperNodeCache();      // build caches for fast apply
    print::time(10, "Time transform", trans_t);
    print::separator(10, ' ');

    // Register this raw operator with the base class and initialize expansion
    this->raw_exp.push_back(std::move(o_tree));
    this->initOperExp(1);
}

// Explicit instantiations
template class BSOperator<1>;
template class BSOperator<2>;
template class BSOperator<3>;

} // namespace mrcpp
