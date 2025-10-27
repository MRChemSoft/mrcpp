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

#include "RepresentableFunction.h"
#include "trees/MultiResolutionAnalysis.h"

namespace mrcpp {

/**
 * @class BoysFunction
 * @brief Adaptive multiresolution evaluator for the 1D Boys-type integral
 *        \f$F_n(x) = \int_{0}^{1} t^{2n}\,e^{-x\,t^2}\,dt\f$.
 *
 * What this class does
 * --------------------
 * This class provides an implementation of MRCPP's @ref RepresentableFunction
 * interface for the scalar function \f$F_n(x)\f$ of a single variable \f$x\f$.
 * Given an input abscissa \f$x\f$, `evalf()`:
 *  1. builds the integrand \f$g_x(t)=e^{-x t^2}\,t^{2n}\f$ on \f$t\in[0,1]\f$,
 *  2. projects it adaptively into a 1D multiresolution basis (using the
 *     `MultiResolutionAnalysis<1>` member),
 *  3. integrates the resulting @ref FunctionTree over the unit interval.
 *
 * Notes on conventions
 * --------------------
 * - In quantum-chemistry literature, the “Boys function” is often defined with
 *   an integral to \f$\infty\f$. Here it is the *unit-interval* variant
 *   \f$[0,1]\f$, consistent with the implementation in the corresponding .cpp.
 * - The basis family and order used by the `MRA` are chosen in the .cpp
 *   definition (currently an interpolating basis of order 13).
 *
 * Accuracy and performance
 * ------------------------
 * - The tolerance passed at construction (`prec`) controls the adaptive
 *   projection target. Smaller values yield higher accuracy at greater cost.
 * - The multiresolution approach concentrates degrees of freedom where the
 *   integrand has structure (e.g., near \f$t=0\f$ for large \f$x\f$).
 */
class BoysFunction final : public RepresentableFunction<1, double> {
public:
    /**
     * @brief Construct an evaluator for \f$F_n(x)\f$.
     *
     * @param n     Non-negative integer order in \f$F_n(x)\f$ (power \f$t^{2n}\f$).
     * @param prec  Target projection precision for the adaptive MRA
     *              (default \f$10^{-10}\f$).
     *
     * Implementation detail:
     * The `MRA` member is initialised in the .cpp with a default 1D bounding
     * box and a fixed scaling basis; this header does not constrain that choice.
     */
    BoysFunction(int n, double prec = 1.0e-10);

    /**
     * @brief Evaluate \f$F_n(x)\f$ at the given abscissa.
     *
     * @param r Coordinate container with a single component: \f$x = r[0]\f$.
     * @return  The numerical value of \f$F_n(x)\f$ obtained by adaptively
     *          projecting and integrating on \f$[0,1]\f$.
     *
     * Semantics:
     *  - Satisfies the @ref RepresentableFunction contract.
     *  - Internally constructs the integrand lambda and invokes the MRCPP
     *    `project` + `integrate` pipeline on the stored `MRA`.
     */
    double evalf(const Coord<1> &r) const override;

private:
    /** @brief Integer order \f$n\f$ in \f$F_n(x)\f$ (kept constant for the lifetime). */
    const int order;

    /** @brief Target projection tolerance for adaptive representation. */
    const double prec;

    /**
     * @brief Multiresolution context used to project/integrate the integrand.
     *
     * The concrete basis family and order are configured in the .cpp file.
     * The same `MRA` instance is reused across evaluations for efficiency.
     */
    MultiResolutionAnalysis<1> MRA;
};

} // namespace mrcpp