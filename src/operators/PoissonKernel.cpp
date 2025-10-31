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

#include "PoissonKernel.h"

#include <cmath>

#include "functions/GaussFunc.h"
#include "utils/Printer.h"

namespace mrcpp {

PoissonKernel::PoissonKernel(double epsilon, double r_min, double r_max)
        : GaussExp<1>() {
    // Constructed on [rMin/rMax, 1.0], then rescaled to [rMin, rMax]
    double r0 = r_min / r_max;
    double r1 = r_max;

    // Choose t1, t2 so that tail contributions are below epsilon
    double t1 = 1.0L;
    while ((2.0 * t1 * std::exp(-t1)) > epsilon) t1 *= 1.1L;

    double t2 = 1.0L;
    while ((std::sqrt(t2) * std::exp(-t2) / r0) > epsilon) t2 *= 1.1L;

    // Truncation window [s1, s2] ensuring relative error ~ epsilon
    double s1 = -std::log(2.0 * t1);
    double s2 = std::log(t2 / (r0 * r0)) / 2.0;

    // Trapezoidal step size h determined from epsilon (empirical fit)
    double h = 1.0 / (0.2L - 0.47L * std::log10(epsilon));
    int n_exp = static_cast<int>(std::ceil((s2 - s1) / h) + 1);
    if (n_exp > MaxSepRank) MSG_ABORT("Maximum separation rank exceeded.");

    for (int i = 0; i < n_exp; i++) {
        double arg = s1 + h * i;
        double sinharg = std::sinh(arg);
        double cosharg = std::cosh(arg);
        double onepexp = 1.0 + std::exp(-sinharg);

        // Parameters before rescaling back to [r_min, r_max]
        double expo = 4.0L * (sinharg + std::log(onepexp)) * (sinharg + std::log(onepexp));
        double coef = h * (4.0L / root_pi) * cosharg / onepexp;

        // Rescale to physical interval
        expo *= 1.0 / (r1 * r1);
        coef *= 1.0 / r1;

        // Trapezoidal rule endpoint correction
        if (i == 0 || i == (n_exp - 1)) coef *= 1.0 / 2.0;

        GaussFunc<1> gFunc(expo, coef);
        this->append(gFunc);
    }
}

} // namespace mrcpp
