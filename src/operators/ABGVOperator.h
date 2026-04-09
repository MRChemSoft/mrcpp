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

#include "DerivativeOperator.h"

namespace mrcpp {

/**
 * @class ABGVOperator
 * @brief Multiresolution first-derivative operator of Alpert–Beylkin–Gines–Vozovoi.
 *
 * This class builds a **first-order differential operator** in the
 * multiresolution (MR) basis defined by a given
 * #mrcpp::MultiResolutionAnalysis. The discrete representation follows
 * the construction in:
 *
 * - B. Alpert, G. Beylkin, D. Gines, and L. Vozovoi,
 *   *Adaptive Solution of Partial Differential Equations in Multiwavelet Bases*,
 *   J. Comput. Phys. **182** (2002) 149–190.
 *
 * ### When to use this operator
 * - **Recommended** for functions with **cusps, kinks, or discontinuities**,
 *   where strictly smooth (BS) operators tend to produce Gibbs-type artifacts.
 * - For **smooth** functions, prefer #mrcpp::BSOperator for slightly better
 *   accuracy/efficiency with smooth stencils.
 *
 * ### Boundary/stencil parameters \p a and \p b
 * The parameters `(a, b)` control the local stencil asymmetry at element
 * interfaces. Common choices:
 *
 * - `a = 0.0`, `b = 0.0`  → strictly local “center” rule (bandwidth 0)
 * - `a = 0.5`, `b = 0.5`  → semi-local **central** difference (bandwidth 1)
 * - `a = 1.0`, `b = 0.0`  → semi-local **forward** difference (bandwidth 1)
 * - `a = 0.0`, `b = 1.0`  → semi-local **backward** difference (bandwidth 1)
 *
 * Any non-zero `a` or `b` widens the coupling to nearest neighbors (bandwidth = 1)
 * across scales; this is enforced during assembly.
 *
 * ### Assembly and application
 * Internally, the operator is assembled once into an #mrcpp::OperatorTree
 * (stored in the base #mrcpp::DerivativeOperator). After construction,
 * applying the operator to MR coefficient vectors is cheap and can be done
 * repeatedly.
 *
 * @tparam D Spatial dimension (1, 2, or 3).
 */
template <int D> class ABGVOperator final : public DerivativeOperator<D> {
public:
    /**
     * @brief Construct the ABGV derivative operator on a given MRA.
     *
     * The constructor triggers an internal `initialize(a, b)` routine that:
     * 1. Decides the operator bandwidth from `(a, b)`.
     * 2. Builds the operator matrix blocks using the MRA’s scaling basis.
     * 3. Assembles an #mrcpp::OperatorTree with a bandwidth adaptor.
     * 4. Finalizes and caches the representation for fast application.
     *
     * @param mra Multiresolution analysis defining the domain, basis and scales.
     * @param a   Left-side boundary/stencil parameter (see class docs).
     * @param b   Right-side boundary/stencil parameter (see class docs).
     *
     * @note The operator is built at the MRA’s **root scale** and is valid for
     *       coefficient vectors defined on the same MRA.
     */
    ABGVOperator(const MultiResolutionAnalysis<D> &mra, double a, double b);

    ABGVOperator(const ABGVOperator &oper) = delete;
    ABGVOperator &operator=(const ABGVOperator &oper) = delete;

protected:
    /**
     * @brief Internal assembly routine (called by the constructor).
     *
     * Decides sparsity (bandwidth) from `(a, b)`, constructs the calculator
     * implementing the ABGV derivative in the given scaling basis, and uses a
     * `TreeBuilder` + `BandWidthAdaptor` to assemble and cache an operator tree.
     *
     * @param a Left boundary/stencil parameter.
     * @param b Right boundary/stencil parameter.
     */
    void initialize(double a, double b);
};

} // namespace mrcpp