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

#include "function_utils.h"

#include "GaussExp.h"
#include "GaussFunc.h"
#include "utils/Printer.h"

namespace mrcpp {

/** @brief Projection of a GaussFunc centered at the Origin of a
 *  periodic unit cell into a MW representation
 *
 * @param[in, out] The GaussExp to be built using the input GaussFunc
 * @param[i] Input GaussFunc with center at the origin that is to be made periodic
 * @param[in] The scaling factor/period of the unit cell
 * @param[in] True/false periodicity of the unit cell, must be true
 *
 * This function takes a GaussFunc center at the Origin then makes it
 * into a GaussExp ensuring the entire function is captured within the
 * unit cell.
 *
 */

template <int D>
void function_utils::make_gaussian_periodic(GaussExp<D> &out,
                                            const GaussFunc<D> &inp,
                                            const std::array<double, D> &scaling_factor,
                                            const bool &periodic) {

    auto beta = inp.getCoef();
    auto alpha = inp.getExp();
    auto pos = inp.getPos();

    if (pos != Coord<D>{}) MSG_FATAL("Gaussian must be placed at the origin");
    if (not periodic) MSG_FATAL("The world has to be periodic");
    out.append(inp);
    if (D == 3 or D == 2) {
        for (auto i = 0; i < D; i++) {
            auto pos_0 = pos;
            pos_0[i] = scaling_factor[i];
            auto gauss_0 = GaussFunc<D>(alpha, beta, pos_0);
            out.append(gauss_0);
            if (D == 3) {
                auto pos_sf = scaling_factor;
                pos_sf[i] = 0.0;
                auto gauss_sf = GaussFunc<D>(alpha, beta, pos_sf);
                out.append(gauss_sf);
            }
        }
    }
    auto gauss_aaa = GaussFunc<D>(alpha, beta, scaling_factor); // Centered at corner diagonal to origin
    out.append(gauss_aaa);
}

template void function_utils::make_gaussian_periodic<1>(GaussExp<1> &out,
                                                        const GaussFunc<1> &inp,
                                                        const std::array<double, 1> &scaling_factor,
                                                        const bool &periodic);
template void function_utils::make_gaussian_periodic<2>(GaussExp<2> &out,
                                                        const GaussFunc<2> &inp,
                                                        const std::array<double, 2> &scaling_factor,
                                                        const bool &periodic);
template void function_utils::make_gaussian_periodic<3>(GaussExp<3> &out,
                                                        const GaussFunc<3> &inp,
                                                        const std::array<double, 3> &scaling_factor,
                                                        const bool &periodic);

} // namespace mrcpp
