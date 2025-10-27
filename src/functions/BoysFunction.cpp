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
 * @file BoysFunction.cpp
 *
 * @brief Numerically evaluates the Boys function
 * \f[
 *   F_n(x) \;=\; \int_{0}^{1} t^{2n}\,e^{-x\,t^2}\,dt
 * \f]
 * by projecting the integrand onto an adaptive multiresolution basis and then
 * integrating the resulting @ref FunctionTree.
 *
 * Design overview
 * ---------------
 * 1) The class derives from @ref RepresentableFunction in 1D so that it can be
 *    used wherever MRCPP expects a function object with `evalf`.
 * 2) Given an input coordinate `r`, we interpret `x = r[0]` and define the
 *    integrand
 *       g_x(t) = e^{-x t^2} · t^{2n},  t ∈ [0,1].
 * 3) We build a 1D @ref FunctionTree using an @ref MRA configured with an
 *    interpolating scaling basis (order 13 by default here), and call
 *    `project(prec, tree, f)` which adaptively refines the tree so that the
 *    projection error is below `prec`.
 * 4) Finally we call `tree.integrate()` to integrate the projected function on
 *    [0,1], which is the desired value F_n(x).
 *
 * Notes
 * -----
 * - The basis choice (`InterpolatingBasis(13)`) is a trade-off: sufficiently
 *   smooth to capture Gaussians well, while keeping stencil sizes reasonable.
 * - The adaptive projection concentrates resolution where the integrand has
 *   structure (e.g., for large x near t=0 the function is sharply peaked).
 * - Printing is muted during evaluation to keep the call side quiet.
 */

#include "BoysFunction.h"
#include "core/InterpolatingBasis.h"  // basis used in the MRA
#include "treebuilders/project.h"     // adaptive projection into a FunctionTree
#include "trees/FunctionTree.h"       // hierarchical representation + integrate()
#include "utils/Printer.h"

namespace mrcpp {

/**
 * @brief Construct a BoysFunction evaluator.
 *
 * @param n   Non-negative integer order in @f$F_n(x)@f$ (power @f$t^{2n}@f$).
 * @param p   Target projection precision (controls adaptive refinement).
 *
 * Internals:
 *  - `order` stores @p n.
 *  - `prec` stores the target accuracy threshold used by `project`.
 *  - `MRA` is initialised over a default 1D bounding box with an
 *    interpolating basis of order 13; this MRA is reused per evaluation.
 */
BoysFunction::BoysFunction(int n, double p)
        : RepresentableFunction<1, double>()
        , order(n)
        , prec(p)
        , MRA(BoundingBox<1>(), InterpolatingBasis(13)) {}

/**
 * @brief Evaluate @f$F_n(x)@f$ at the requested abscissa.
 *
 * @param r Coordinate container; the one and only component is @f$x=r[0]@f$.
 * @return  The value of @f$F_n(x)=\int_0^1 t^{2n} e^{-x t^2}\,dt@f$.
 *
 * Algorithm:
 *  1) Silence the printer and remember the old verbosity.
 *  2) Capture `x` and `n` and form a lambda `f(t)` representing the integrand
 *     on @f$t\in[0,1]@f$. We compute `t_2 = t^2`, `t_2n = (t^2)^n`, and return
 *     `exp(-x * t_2) * t_2n`. For `n=0`, `t_2n` is set to 1 for speed/stability.
 *  3) Build a fresh `FunctionTree<1,double>` bound to the stored `MRA`.
 *  4) Call `project(prec, tree, f)` to approximate `f` within the tolerance
 *     `prec` by adaptively refining nodes where needed.
 *  5) Call `tree.integrate()` to obtain the integral over [0,1].
 *  6) Restore the printer level and return the result.
 *
 * Accuracy remarks:
 *  - The achieved error depends on `prec`, the basis order, and the behaviour
 *    of the integrand (large x leads to rapid decay, which is well captured
 *    by the multiresolution approach).
 */
double BoysFunction::evalf(const Coord<1> &r) const {
    // Temporarily mute verbose output from the projection/integration.
    int oldlevel = Printer::setPrintLevel(0);

    int n = this->order;
    double x = r[0];

    // Integrand g_x(t) = exp(-x * t^2) * (t^2)^n over t in [0,1].
    // Written in terms of t^2 to reduce pow() evaluation count.
    auto f = [x, n](const Coord<1> &t) -> double {
        double t_2 = t[0] * t[0];
        double xt_2 = x * t_2;
        double t_2n = 1.0;
        if (n > 0) { t_2n = std::pow(t_2, n); }
        return std::exp(-xt_2) * t_2n;
    };

    // Build an adaptive representation of f on [0,1] and integrate it.
    FunctionTree<1, double> tree(this->MRA);
    mrcpp::project<1, double>(this->prec, tree, f);
    double result = tree.integrate();

    // Restore previous verbosity and return.
    Printer::setPrintLevel(oldlevel);
    return result;
}

} // namespace mrcpp