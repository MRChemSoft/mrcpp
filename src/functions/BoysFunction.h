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
 * @brief Adaptive multiresolution evaluator for the Boys function
 *        \f$ F_n(x) = \int_{0}^{1} t^{2n}\mathrm{e}^{-xt^2}\mathrm{d}t \f$,
 *        where \f$ x\ge0 \f$ and \f$ n\ge0 \f$
 */
class BoysFunction final : public RepresentableFunction<1, double> {
public:
    /**
     * @brief Construct an evaluator for \f$ F_n(x) \f$
     * @param n The order (\f$ \ge0 \f$) of the Boys function
     * @param prec Projection precision for the adaptive MRA (default \f$ 10^{-10} \f$)
     *
     * @details The `MRA` member is initialised in the .cpp with a default 1D bounding
     * box and a fixed scaling basis (currently an interpolating basis of order 13);
     * this header does not constrain that choice.
     */
    BoysFunction(int n, double prec = 1.0e-10);

    /**
     * @brief Evaluate \f$ F_n(x) \f$ at the given abscissa
     * @param r Coordinate container with a single component: \f$ x = r[0] \f$
     * @return  The numerical value of \f$ F_n(x) \f$
     */
    double evalf(const Coord<1> &r) const override;

private:
    const int order;

    const double prec;

    MultiResolutionAnalysis<1> MRA;
};

} // namespace mrcpp