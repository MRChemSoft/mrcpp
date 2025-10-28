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

#include "ABGVOperator.h"
#include "treebuilders/ABGVCalculator.h"
#include "treebuilders/BandWidthAdaptor.h"
#include "treebuilders/TreeBuilder.h"
#include "trees/OperatorTree.h"
#include "utils/Printer.h"
#include "utils/Timer.h"

namespace mrcpp {

/**
 * # ABGV finite-difference(-like) operator on an MRA
 *
 * This operator implements a family of first-derivative stencils controlled by two
 * boundary parameters \p a and \p b (see below). The operator is **assembled once**
 * as an `OperatorTree` in the multiresolution basis of the provided MRA and can then
 * be applied repeatedly to vectors/functions defined on the same MRA.
 *
 * ## Boundary parameters
 * The pair `(a,b)` selects a particular linear combination of forward/backward bias:
 *
 * - `a = 0.0`, `b = 0.0`  → strictly local “center” difference (bandwidth 0)
 * - `a = 0.5`, `b = 0.5`  → semi-local **central** difference (bandwidth 1)
 * - `a = 1.0`, `b = 0.0`  → semi-local **forward** difference (bandwidth 1)
 * - `a = 0.0`, `b = 1.0`  → semi-local **backward** difference (bandwidth 1)
 *
 * Any non-zero `a` or `b` increases the operator’s bandwidth to 1 (one-ring coupling
 * between neighboring nodes at each scale), which the `BandWidthAdaptor` enforces.
 *
 * @tparam D Spatial dimension (1, 2, or 3).
 * @param mra The multiresolution analysis that defines basis, scales, and domain.
 * @param a   Left boundary parameter controlling asymmetry at the “minus” side.
 * @param b   Right boundary parameter controlling asymmetry at the “plus” side.
 *
 * @note The operator is built at the **root scale** of the provided MRA, and its
 *       internal representation (raw expansion) is cached for later applications.
 */
template <int D>
ABGVOperator<D>::ABGVOperator(const MultiResolutionAnalysis<D> &mra, double a, double b)
        : DerivativeOperator<D>(mra, mra.getRootScale()) {
    initialize(a, b);
}

/**
 * @brief Internal construction routine: builds a bandwidth-adapted OperatorTree.
 *
 * Steps (high level):
 *  1. **Bandwidth decision** — if either \p a or \p b is non-zero, set bandwidth = 1,
 *     otherwise 0. This determines how many neighbor interactions the operator will keep.
 *  2. **Calculator** — instantiate `ABGVCalculator` with the MRA’s scaling basis and (a,b).
 *     The calculator knows how to evaluate local operator blocks (stencil entries) in
 *     the chosen basis.
 *  3. **Adaptor** — create a `BandWidthAdaptor(bw, maxScale)` to prune any far-off
 *     couplings beyond the requested bandwidth across all scales.
 *  4. **Tree build** — use `TreeBuilder<2>` (matrix builder) to assemble an `OperatorTree`
 *     from root to finest scale with tolerance `MachineZero` and adaptor-controlled sparsity.
 *  5. **Finalize** — trigger norm computation and set up an operator-node cache for fast
 *     application; then store the finished tree in `raw_exp` and initialize the expansion.
 *
 * @param a Left boundary parameter.
 * @param b Right boundary parameter.
 *
 * @details
 * - `calcSquareNorm()` performs a pass that also ensures the internal transform state is
 *   consistent (it may trigger lazy transforms). We time this step for diagnostics.
 * - `setupOperNodeCache()` precomputes/cache-friendly structures for repeated operator
 *   application (e.g., fast traversal, block reuse).
 * - `initOperExp(1)` finalizes the operator’s internal expansion (single component here).
 */
template <int D> void ABGVOperator<D>::initialize(double a, double b) {
    // --- (1) Decide operator bandwidth from boundary parameters -------------------------
    int bw = 0; // 0 = strictly local, 1 = nearest-neighbor coupling
    if (std::abs(a) > MachineZero) bw = 1;
    if (std::abs(b) > MachineZero) bw = 1;

    // --- (2) Access the operator MRA ----------------------------------------------------
    auto oper_mra = this->getOperatorMRA();

    // --- (3) Prepare builder, calculator, and bandwidth adaptor -------------------------
    TreeBuilder<2> builder;  // <2> means: building a 2-index object (matrix/operator)
    ABGVCalculator calculator(oper_mra.getScalingBasis(), a, b);
    BandWidthAdaptor adaptor(bw, oper_mra.getMaxScale());

    // --- (4) Assemble the operator tree -------------------------------------------------
    // MachineZero: force exact assembly within floating point epsilon (no thresholding).
    auto o_tree = std::make_unique<OperatorTree>(oper_mra, MachineZero);
    builder.build(*o_tree, calculator, adaptor, -1 /* build all scales */);

    // --- (5) Finalize and cache ---------------------------------------------------------
    Timer trans_t;
    o_tree->calcSquareNorm();        // also ensures internal transforms are ready
    o_tree->setupOperNodeCache();    // allocate and fill fast-access caches
    print::time(10, "Time transform", trans_t);
    print::separator(10, ' ');

    // Keep the assembled operator as our raw expansion and finalize its use
    this->raw_exp.push_back(std::move(o_tree));
    this->initOperExp(1);            // single-operator expansion component
}

// Explicit template instantiations for 1D, 2D, and 3D operators.
template class ABGVOperator<1>;
template class ABGVOperator<2>;
template class ABGVOperator<3>;

} // namespace mrcpp