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

#include "CartesianConvolution.h"

#include "core/InterpolatingBasis.h"
#include "core/LegendreBasis.h"

#include "functions/GaussExp.h"
#include "functions/Gaussian.h"

#include "treebuilders/CrossCorrelationCalculator.h"
#include "treebuilders/OperatorAdaptor.h"
#include "treebuilders/TreeBuilder.h"
#include "treebuilders/grid.h"
#include "treebuilders/project.h"

#include "trees/BandWidth.h"
#include "trees/FunctionTreeVector.h"
#include "trees/OperatorTree.h"

#include "utils/Printer.h"
#include "utils/Timer.h"
#include "utils/math_utils.h"

namespace mrcpp {

/**
 * @class CartesianConvolution
 * @brief 3D separable convolution operator assembled from a 1D Gaussian expansion.
 *
 * ### What this class builds
 * We construct a rank-`R` separable operator (with `R = kernel.size()`) that can
 * later be combined into Cartesian components (x, y, z). Internally we build
 * **three batches** of 1D operator trees from the same 1D Gaussian expansion,
 * corresponding to polynomial prefactors of degree 0, 1, and 2 multiplying the
 * Gaussian (i.e., monomials \f$x^0, x^1, x^2\f$ on the line). These three
 * batches are stored back-to-back in `raw_exp` and can be mapped onto the
 * (x, y, z) axes in any order using #setCartesianComponents.
 *
 * This pattern is useful for assembling vector kernels whose Cartesian
 * components differ only by the polynomial factor in each axis (e.g., kernels
 * proportional to \f$(1,\,x,\,x^2)\f$ times a radial Gaussian, or for forming
 * gradients / moments where different axes pick different polynomial orders).
 *
 * ### Precision knobs
 * - `prec` is the user‐requested overall build precision for the operator.
 * - We derive two internal precisions:
 *   - `o_prec = prec` for operator assembly;
 *   - `k_prec = prec / 10` for fitting/projection of the 1D kernel, slightly
 *     tighter so that the overall composition meets the requested tolerance.
 *
 * ### Memory layout of the built batches
 * After construction, `raw_exp` contains `3 * R` operator trees in this order:
 * ```
 * block 0: monomial power {0} for all R terms      (indices 0 ... R-1)
 * block 1: monomial power {1} for all R terms      (indices R ... 2R-1)
 * block 2: monomial power {2} for all R terms      (indices 2R ... 3R-1)
 * ```
 * The method #setCartesianComponents selects one of these three blocks per axis.
 */

/**
 * @brief Construct a Cartesian convolution operator on an MRA with a 1D Gaussian expansion.
 *
 * @param[in] mra     3D multiresolution analysis defining basis/domain/scales.
 * @param[in,out] kernel  1D Gaussian expansion \f$ \sum_{r=1}^R g_r(x) \f$ used to
 *                        generate the separable operator factors. Its length
 *                        determines the separation rank \f$R\f$.
 *                        **Note:** This function temporarily modifies the
 *                        monomial power of each Gaussian term and restores it
 *                        across the three assembly passes.
 * @param[in] prec    Target build precision for the operator.
 *
 * @details
 * **Assembly recipe (done three times):**
 *  1. For every term in the input 1D Gaussian expansion, set its monomial
 *     power to `{0}`, then call `initialize(...)` to build and append one
 *     operator tree per term (rank-`R` block).
 *  2. Repeat with monomial power `{1}` to build the second block (indices
 *     `R ... 2R-1`).
 *  3. Repeat with monomial power `{2}` to build the third block (indices
 *     `2R ... 3R-1`).
 *
 * After these three passes, we call `initOperExp(R)` to declare that downstream
 * separable composition will have rank \f$R\f$ (each axis picks one block).
 *
 * **Why powers {0,1,2}?**
 * Many Cartesian tensor kernels (e.g., derivatives, moments, or vector fields)
 * differ by low-order polynomial prefactors along each coordinate. Prebuilding
 * the families \f$\{0,1,2\}\f$ provides flexible combinations via
 * #setCartesianComponents without having to rebuild for each axis.
 */
CartesianConvolution::CartesianConvolution(const MultiResolutionAnalysis<3> &mra,
                                           GaussExp<1> &kernel,
                                           double prec)
        : ConvolutionOperator<3>(mra)
        , sep_rank(kernel.size()) {
    int oldlevel = Printer::setPrintLevel(0);

    // Configure precision: operator vs. kernel fit
    this->setBuildPrec(prec);
    auto o_prec = prec;          // Operator assembly precision
    auto k_prec = prec / 10.0;   // Kernel fitting precision (tighter on purpose)

    // --- Batch 0: monomial power {0} (constant prefactor) ---
    for (auto &k : kernel) k->setPow({0});
    this->initialize(kernel, k_prec, o_prec);

    // --- Batch 1: monomial power {1} (linear prefactor) ---
    for (auto &k : kernel) k->setPow({1});
    this->initialize(kernel, k_prec, o_prec);

    // --- Batch 2: monomial power {2} (quadratic prefactor) ---
    for (auto &k : kernel) k->setPow({2});
    this->initialize(kernel, k_prec, o_prec);

    // Tell the separable framework we will later combine per-axis using rank = sep_rank
    this->initOperExp(this->sep_rank);

    Printer::setPrintLevel(oldlevel);
}

/**
 * @brief Choose which prebuilt monomial block (0,1,2) to use for each Cartesian axis.
 *
 * @param[in] x  Block index used for the x-axis   (0 → power{0}, 1 → power{1}, 2 → power{2})
 * @param[in] y  Block index used for the y-axis   (same convention as above)
 * @param[in] z  Block index used for the z-axis   (same convention as above)
 *
 * @details
 * This function **does not** rebuild; it only wires the already constructed
 * 1D operator trees into the separable 3D operator slots. For separation rank
 * \f$R\f$, each block occupies a contiguous range of \f$R\f$ entries:
 *
 * - Block `x`: indices `[x*R, x*R + R - 1]` become the x-factors;
 * - Block `y`: indices `[y*R, y*R + R - 1]` become the y-factors;
 * - Block `z`: indices `[z*R, z*R + R - 1]` become the z-factors.
 *
 * You may reuse the same block on multiple axes if the physics warrants it
 * (e.g., isotropic components), or select different ones to form vector/tensor
 * kernels with distinct Cartesian prefactors.
 *
 * @warning Valid block indices are 0, 1, or 2. No bounds checking is performed here.
 */
void CartesianConvolution::setCartesianComponents(int x, int y, int z) {
    int x_shift = x * this->sep_rank;
    int y_shift = y * this->sep_rank;
    int z_shift = z * this->sep_rank;

    // Fill the separable operator slots (rank index i, axis 0/1/2) with the chosen blocks.
    for (int i = 0; i < this->sep_rank; i++) this->assign(i, 0, this->raw_exp[x_shift + i].get());
    for (int i = 0; i < this->sep_rank; i++) this->assign(i, 1, this->raw_exp[y_shift + i].get());
    for (int i = 0; i < this->sep_rank; i++) this->assign(i, 2, this->raw_exp[z_shift + i].get());
}

} // namespace mrcpp