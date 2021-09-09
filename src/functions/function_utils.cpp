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
#include "utils/details.h"

#include "GaussExp.h"
#include "Gaussian.h"
#include "utils/Printer.h"
#include "utils/math_utils.h"

#include <memory>
#include <numeric>

namespace mrcpp {

namespace function_utils {
double ObaraSaika_ab(int power_a, int power_b, double pos_a, double pos_b, double expo_a, double expo_b);
} // namespace function_utils

/** @brief Generates a GaussExp that is semi-periodic around a unit-cell
 *
 * @returns Semi-periodic version of a Gaussian around a unit-cell
 * @param[in] inp: A Gaussian function that is to be made semi-periodic
 * @param[in] period: The period of the unit cell
 * @param[in] nStdDev: Number of standard diviations covered in each direction. Default 4.0
 *
 * @details nStdDev = 1, 2, 3 and 4 ensures atleast 68.27%, 95.45%, 99.73% and 99.99% of the
 * integral is conserved with respect to the integration limits.
 *
 */
template <int D> GaussExp<D> function_utils::periodify(const Gaussian<D> &inp, const std::array<double, D> &period, double nStdDev) {
    GaussExp<D> gauss_exp;
    auto pos_vec = std::vector<Coord<D>>();

    auto x_std = nStdDev * inp.getMaximumStandardDiviation();

    // This lambda function  calculates the number of neighbooring cells
    // requred to keep atleast x_stds of the integral conserved in the
    // unit-cell.
    auto neighbooring_cells = [period, x_std](auto pos) {
        auto needed_cells_vec = std::vector<int>();
        for (auto i = 0; i < D; i++) {
            auto upper_bound = pos[i] + x_std;
            auto lower_bound = pos[i] - x_std;
            // number of cells upp and down relative to the center of the Gaussian
            needed_cells_vec.push_back(std::ceil(upper_bound / period[i]));
        }

        return *std::max_element(needed_cells_vec.begin(), needed_cells_vec.end());
    };

    // Finding starting position
    auto startpos = inp.getPos();

    for (auto d = 0; d < D; d++) {
        startpos[d] = std::fmod(startpos[d], period[d]);
        if (startpos[d] < 0) startpos[d] += period[d];
    }

    auto nr_cells_upp_and_down = neighbooring_cells(startpos);
    for (auto d = 0; d < D; d++) { startpos[d] -= nr_cells_upp_and_down * period[d]; }

    auto tmp_pos = startpos;
    std::vector<double> v(2 * nr_cells_upp_and_down + 1);
    std::iota(v.begin(), v.end(), 0.0);
    auto cart = math_utils::cartesian_product(v, D);
    for (auto &c : cart) {
        for (auto i = 0; i < D; i++) c[i] *= period[i];
    }
    // Shift coordinates
    for (auto &c : cart) std::transform(c.begin(), c.end(), tmp_pos.begin(), c.begin(), std::plus<double>());
    // Go from vector to mrcpp::Coord
    for (auto &c : cart) {
        mrcpp::Coord<D> pos;
        std::copy_n(c.begin(), D, pos.begin());
        pos_vec.push_back(pos);
    }

    for (auto &pos : pos_vec) {
        auto *gauss = inp.copy();
        gauss->setPos(pos);
        gauss_exp.append(*gauss);
        delete gauss;
    }

    return gauss_exp;
}

template <int D> GaussExp<D> function_utils::periodify(const GaussExp<D> &gexp, const std::array<double, D> &period, double nStdDev) {
    GaussExp<D> out_exp;
    for (const auto &gauss : gexp) {
        auto periodic_gauss = periodify<D>(*gauss, period, nStdDev);
        out_exp.append(periodic_gauss);
    }
    return out_exp;
}

template <int D> double function_utils::calc_overlap(const GaussFunc<D> &a, const GaussFunc<D> &b) {
    double S = 1.0;
    for (int d = 0; d < D; d++) {
        S *= ObaraSaika_ab(a.getPower()[d], b.getPower()[d], a.getPos()[d], b.getPos()[d], a.getExp()[d], b.getExp()[d]);
    }
    S *= a.getCoef() * b.getCoef();
    return S;
}

/**  Compute the monodimensional overlap integral between two
 gaussian distributions by means of the Obara-Saika recursiive
 scheme

 \f[ S_{ij} = \int_{-\infty}^{+\infty} \,\mathrm{d} x
 (x-x_a)^{p_a}
 (x-x_b)^{p_b}
 e^{-c_a (x-x_a)^2}
 e^{-c_b (x-x_b)^2}\f]

 @param power_a \f$ p_a     \f$
 @param power_b \f$ p_b     \f$
 @param pos_a   \f$ x_a     \f$
 @param pos_b   \f$ x_b     \f$
 @param expo_a  \f$ c_a \f$
 @param expo_b  \f$ c_b \f$

 */
double function_utils::ObaraSaika_ab(int power_a, int power_b, double pos_a, double pos_b, double expo_a, double expo_b) {
    int i, j;
    double expo_p, mu, pos_p, x_ab, x_pa, x_pb, s_00;
    /* The highest angular momentum combination is l=20 for a and b
     * simulatnelusly */
    double s_coeff[64];

    //	if (out_of_bounds(power_a, 0, MAX_GAUSS_POWER) ||
    //		out_of_bounds(power_b, 0, MAX_GAUSS_POWER)
    //		) {
    //		PRINT_FUNC_NAME;
    //		INVALID_ARG_EXIT;
    //	}

    /* initialization of a hell of a lot of coefficients.... */
    expo_p = expo_a + expo_b;                           /* total exponent */
    mu = expo_a * expo_b / (expo_a + expo_b);           /* reduced exponent */
    pos_p = (expo_a * pos_a + expo_b * pos_b) / expo_p; /* center of charge */
    x_ab = pos_a - pos_b;                               /* X_{AB} */
    x_pa = pos_p - pos_a;                               /* X_{PA} */
    x_pb = pos_p - pos_b;                               /* X_{PB} */
    s_00 = pi / expo_p;
    s_00 = std::sqrt(s_00) * std::exp(-mu * x_ab * x_ab); /* overlap of two spherical gaussians */
    // int n_0j_coeff = 1 + power_b; /* n. of 0j coefficients needed */
    // int n_ij_coeff = 2 * power_a; /* n. of ij coefficients needed (i > 0) */

    /* we add 3 coeffs. to avoid a hell of a lot of if statements */
    /*    n_tot_coeff = n_0j_coeff + n_ij_coeff + 3;	*/
    /*    s_coeff = (double *) calloc(n_tot_coeff, sizeof(double));*/

    /* generate first two coefficients */
    s_coeff[0] = s_00;
    s_coeff[1] = x_pb * s_00;
    j = 1;
    /* generate the rest of the first row */
    while (j < power_b) {
        s_coeff[j + 1] = x_pb * s_coeff[j] + j * s_coeff[j - 1] / (2.0 * expo_p);
        j++;
    }
    /* generate the first two coefficients with i > 0 */
    s_coeff[j + 1] = s_coeff[j] - x_ab * s_coeff[j - 1];
    s_coeff[j + 2] = x_pa * s_coeff[j] + j * s_coeff[j - 1] / (2.0 * expo_p);
    i = 1;
    /* generate the remaining coefficients with i > 0 */
    while (i < power_a) {
        int i_l = j + 2 * i + 1;
        int i_r = j + 2 * i + 2;
        s_coeff[i_l] = s_coeff[i_l - 1] - x_ab * s_coeff[i_l - 2];
        s_coeff[i_r] = x_pa * s_coeff[i_r - 2] + (j * s_coeff[i_r - 3] + i * s_coeff[i_r - 4]) / (2.0 * expo_p);
        i++;
    }

    /*    free(s_coeff);*/
    return s_coeff[power_b + 2 * power_a];
}

template GaussExp<1> function_utils::periodify<1>(const Gaussian<1> &inp, const std::array<double, 1> &period, double nStdDev);
template GaussExp<2> function_utils::periodify<2>(const Gaussian<2> &inp, const std::array<double, 2> &period, double nStdDev);
template GaussExp<3> function_utils::periodify<3>(const Gaussian<3> &inp, const std::array<double, 3> &period, double nStdDev);
template GaussExp<1> function_utils::periodify<1>(const GaussExp<1> &inp, const std::array<double, 1> &period, double nStdDev);
template GaussExp<2> function_utils::periodify<2>(const GaussExp<2> &inp, const std::array<double, 2> &period, double nStdDev);
template GaussExp<3> function_utils::periodify<3>(const GaussExp<3> &inp, const std::array<double, 3> &period, double nStdDev);
template double function_utils::calc_overlap<1>(const GaussFunc<1> &a, const GaussFunc<1> &b);
template double function_utils::calc_overlap<2>(const GaussFunc<2> &a, const GaussFunc<2> &b);
template double function_utils::calc_overlap<3>(const GaussFunc<3> &a, const GaussFunc<3> &b);
} // namespace mrcpp
