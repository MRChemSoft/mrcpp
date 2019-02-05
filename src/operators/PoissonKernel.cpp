/*
 * MRCPP, a numerical library based on multiresolution analysis and
 * the multiwavelet basis which provide low-scaling algorithms as well as
 * rigorous error control in numerical computations.
 * Copyright (C) 2019 Stig Rune Jensen, Jonas Juselius, Luca Frediani and contributors.
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

/*
 *
 *
 *  \date Jul 7, 2009
 *  \author Jonas Juselius <jonas.juselius@uit.no> \n
 *          CTCC, University of Tromsø
 *
 * \breif
 */

#include "PoissonKernel.h"

#include <cmath>

#include "MRCPP/constants.h"

#include "functions/GaussFunc.h"
#include "utils/Printer.h"

namespace mrcpp {

/** generate an approximation of the 3d poisson kernel expanded in
 * gaussian functions this routine assumes that the expansion be centered
 */
void PoissonKernel::initializeKernel() {
    // Constructed on [rMin/rMax, 1.0], and then rescaled to [rMin,rMax]
    double r0 = this->rMin / this->rMax;
    double r1 = this->rMax;

    double t1 = 1.0L;
    while ((2 * t1 * std::exp(-t1)) > this->epsilon) { t1 *= 1.1L; }
    double t2 = 1.0L;
    while ((std::sqrt(t2) * std::exp(-t2) / r0) > this->epsilon) { t2 *= 1.1L; }

    // Set the truncation limits s1,s2 of the integral (integrate over [s1,s2])
    // for achieving relative error epsilon
    double s1 = -std::log(2 * t1);
    double s2 = std::log(t2 / (r0 * r0)) / 2;

    // Now, set the step size h for use in the trapezoidal rule for given MU
    double h = 1 / (0.2L - 0.47L * std::log10(this->epsilon));
    int n_exp = (int)std::ceil((s2 - s1) / h) + 1;

    if (n_exp > MaxSepRank) MSG_FATAL("Maximum separation rank exceeded.");

    for (int i = 0; i < n_exp; i++) {
        double arg = s1 + h * i;
        double sinharg = std::sinh(arg);
        double cosharg = std::cosh(arg);
        double onepexp = 1.0 + std::exp(-sinharg);

        double alpha = 4.0L * (sinharg + std::log(onepexp)) * (sinharg + std::log(onepexp));
        double beta = h * (4.0L / root_pi) * cosharg / onepexp;

        alpha *= 1.0 / (r1 * r1);
        beta *= 1.0 / r1;
        if (i == 0 or i == (n_exp - 1)) { beta *= 1.0 / 2.0; }

        GaussFunc<1> gFunc(alpha, beta);
        this->append(gFunc);
    }
    this->calcSquareNorm();
}

std::ostream &PoissonKernel::print(std::ostream &o) const {
    o << " PoissonKernel: " << std::endl;
    o << " epsilon:  " << this->epsilon << std::endl;
    o << " rMin:     " << this->rMin << std::endl;
    o << " rMax:     " << this->rMax << std::endl;
    return o;
}

} // namespace mrcpp
