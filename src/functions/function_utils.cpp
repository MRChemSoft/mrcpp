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
template <int D> std::shared_ptr<GaussExp<D>> function_utils::make_gaussian_periodic(const Gaussian<D> &inp, const std::array<double, D> &period, double nStdDev) {
    auto gauss_exp = std::make_shared<GaussExp<D>>();
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
        gauss_exp->append(*gauss);
        delete gauss;
    }

    return gauss_exp;
}

template std::shared_ptr<GaussExp<1>> function_utils::make_gaussian_periodic<1>(const Gaussian<1> &inp, const std::array<double, 1> &period, double nStdDev);
template std::shared_ptr<GaussExp<2>> function_utils::make_gaussian_periodic<2>(const Gaussian<2> &inp, const std::array<double, 2> &period, double nStdDev);
template std::shared_ptr<GaussExp<3>> function_utils::make_gaussian_periodic<3>(const Gaussian<3> &inp, const std::array<double, 3> &period, double nStdDev);
} // namespace mrcpp
