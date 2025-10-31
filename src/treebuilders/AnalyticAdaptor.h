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

#include "MRCPP/constants.h"
#include "TreeAdaptor.h"

namespace mrcpp {

/**
 * @class AnalyticAdaptor
 * @brief Refinement policy that consults an analytic (representable) function.
 *
 * @tparam D Spatial dimension (1, 2, or 3).
 * @tparam T Coefficient/value type (e.g., `double`, `ComplexDouble`).
 *
 * @details
 * This adaptor requests refinement of a node when the provided analytic
 * function is **not yet visible** at the node's current scale and is **not
 * identically zero** on the node's cell. Concretely:
 *
 * - If `func->isVisibleAtScale(scale, kp1)` returns **true**, the node is
 *   considered sufficiently resolved at this scale → **do not split**.
 * - Else, if `func->isZeroOnInterval(lb, ub)` returns **true**, the function
 *   vanishes on the cell → **do not split**.
 * - Otherwise, the feature likely requires more resolution → **split**.
 *
 * The visibility test uses the node’s polynomial order `k+1` (via `getKp1()`)
 * as the quadrature/collocation count hint for the analytic oracle.
 *
 * ### Requirements on the analytic function
 * The `RepresentableFunction<D,T>` passed in must implement:
 * - `bool isVisibleAtScale(int scale, int nQuadPts) const;`
 * - `bool isZeroOnInterval(const double* lower, const double* upper) const;`
 *
 * ### Typical usage
 * @code{.cpp}
 * AnalyticFunction<3,double>  f(...);   // implements the required interface
 * AnalyticAdaptor<3,double>   adapt(f, mra.getMaxScale());
 * TreeBuilder<3,double>       builder;
 * DefaultCalculator<3,double> calc;
 * builder.build(tree, calc, adapt, -1); // maxIter: unbounded
 * @endcode
 */
template <int D, typename T>
class AnalyticAdaptor final : public TreeAdaptor<D, T> {
public:
    /**
     * @brief Construct an analytic-driven adaptor.
     * @param f  Analytic (representable) function used as refinement oracle.
     * @param ms Maximum allowed scale for splitting (forwarded to TreeAdaptor).
     */
    AnalyticAdaptor(const RepresentableFunction<D, T> &f, int ms)
            : TreeAdaptor<D, T>(ms)
            , func(&f) {}

private:
    /// Pointer to the refinement oracle (not owned).
    const RepresentableFunction<D, T> *func;

    /**
     * @brief Decide whether a node should be split.
     *
     * @param node Candidate node to test.
     * @return `true` if refinement is requested; `false` otherwise.
     *
     * @details
     * Uses the two-step logic described in the class documentation:
     *  1) skip split if visible at current scale,
     *  2) skip split if identically zero on the node's interval,
     *  3) otherwise split.
     */
    bool splitNode(const MWNode<D, T> &node) const override {
        int scale = node.getScale();
        int nQuadPts = node.getKp1();
        if (this->func->isVisibleAtScale(scale, nQuadPts)) return false;
        auto ub = node.getUpperBounds();
        auto lb = node.getLowerBounds();
        if (this->func->isZeroOnInterval(lb.data(), ub.data())) return false;
        return true;
    }
};

} // namespace mrcpp