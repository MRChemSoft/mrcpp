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

#include "HelmholtzKernel.h"

#include <cmath>

#include "functions/GaussFunc.h"
#include "utils/Printer.h"

namespace mrcpp {

HelmholtzKernel::HelmholtzKernel(double mu, double epsilon, double r_min, double r_max)
        : GaussExp<1>() {
    const double r0 = r_min / r_max;
    const double r1 = r_max;
    const double mu_tilde = mu * r1;

    const long double t = std::max((-2.5L * std::log(epsilon)), 5.0L);
    const double s1 = -std::log(4.0L * t / (mu_tilde * mu_tilde)) / 2.0L;
    const double s2 =  std::log(t / (r0 * r0)) / 2.0L;

    const double h = 1.0 / (0.20L - 0.47L * std::log10(epsilon));
    const int n_exp = static_cast<int>(std::ceil((s2 - s1) / h) + 1.0);

    if (n_exp > MaxSepRank) MSG_ABORT("Maximum separation rank exceeded.");

    for (int i = 0; i < n_exp; ++i) {
        const double s = s1 + h * i;

        const double temp  = -2.0 * s;
        const double temp2 = - (mu_tilde * mu_tilde) * std::exp(temp) / 4.0 + s;

        double beta  = h * (2.0 / root_pi) * std::exp(temp2);
        double alpha = std::exp(2.0L * s);

        alpha *= 1.0 / (r1 * r1);
        beta  *= 1.0 / r1;

        if (i == 0 || i == (n_exp - 1)) beta *= 0.5;

        GaussFunc<1> gFunc(alpha, beta);
        this->append(gFunc);
    }
}

} // namespace mrcpp