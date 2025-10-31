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

#include "function_utils.h"

namespace mrcpp {

namespace function_utils {
double ObaraSaika_ab(int power_a, int power_b,
                     double pos_a, double pos_b,
                     double expo_a, double expo_b);
} // namespace function_utils

template <int D>
double function_utils::calc_overlap(const GaussFunc<D> &a, const GaussFunc<D> &b) {
    double S = 1.0;
    for (int d = 0; d < D; d++) {
        S *= ObaraSaika_ab(
            a.getPower()[d], b.getPower()[d],
            a.getPos()[d],   b.getPos()[d],
            a.getExp()[d],   b.getExp()[d]
        );
    }
    S *= a.getCoef() * b.getCoef();
    return S;
}

double function_utils::ObaraSaika_ab(int power_a, int power_b,
                                     double pos_a, double pos_b,
                                     double expo_a, double expo_b) {
    int i, j;
    double expo_p, mu, pos_p, x_ab, x_pa, x_pb, s_00;

    double s_coeff[64];

    expo_p = expo_a + expo_b;
    mu     = expo_a * expo_b / (expo_a + expo_b);
    pos_p  = (expo_a * pos_a + expo_b * pos_b) / expo_p;
    x_ab   = pos_a - pos_b;
    x_pa   = pos_p - pos_a;
    x_pb   = pos_p - pos_b;

    s_00 = pi / expo_p;
    s_00 = std::sqrt(s_00) * std::exp(-mu * x_ab * x_ab);

    s_coeff[0] = s_00;
    s_coeff[1] = x_pb * s_00;

    j = 1;
    while (j < power_b) {
        s_coeff[j + 1] = x_pb * s_coeff[j] + j * s_coeff[j - 1] / (2.0 * expo_p);
        j++;
    }

    s_coeff[j + 1] = s_coeff[j] - x_ab * s_coeff[j - 1];
    s_coeff[j + 2] = x_pa * s_coeff[j] + j * s_coeff[j - 1] / (2.0 * expo_p);

    i = 1;
    while (i < power_a) {
        int i_l = j + 2 * i + 1;
        int i_r = j + 2 * i + 2;

        s_coeff[i_l] = s_coeff[i_l - 1] - x_ab * s_coeff[i_l - 2];
        s_coeff[i_r] = x_pa * s_coeff[i_r - 2] + (j * s_coeff[i_r - 3] + i * s_coeff[i_r - 4]) / (2.0 * expo_p);

        i++;
    }

    return s_coeff[power_b + 2 * power_a];
}

template double function_utils::calc_overlap<1>(const GaussFunc<1> &a, const GaussFunc<1> &b);
template double function_utils::calc_overlap<2>(const GaussFunc<2> &a, const GaussFunc<2> &b);
template double function_utils::calc_overlap<3>(const GaussFunc<3> &a, const GaussFunc<3> &b);

} // namespace mrcpp