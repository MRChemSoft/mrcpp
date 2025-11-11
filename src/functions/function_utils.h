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



#include "GaussExp.h"
#include "Gaussian.h"

namespace mrcpp {

// Forward declaration only: definition is provided in function_utils.cpp.
// Keeping this here avoids heavy includes and potential include cycles.
namespace function_utils {
template <int D>

/**  
 * @brief Compute the monodimensional overlap integral between two
 * gaussian distributions by means of the Obara-Saika recursive
 * scheme
 *
 * \f$ [ S_{ij} = \int_{-\infty}^{+\infty} \,\mathrm{d} x
 * (x-x_a)^{p_a}
 * (x-x_b)^{p_b}
 * e^{-c_a (x-x_a)^2}
 * e^{-c_b (x-x_b)^2} \f$
 *
 * @param power_a \f$ p_a     \f$
 * @param power_b \f$ p_b     \f$
 * @param pos_a   \f$ x_a     \f$
 * @param pos_b   \f$ x_b     \f$
 * @param expo_a  \f$ c_a \f$
 * @param expo_b  \f$ c_b \f$
 *
 * @return The value of the overlap integral as a d
 */
double function_utils::ObaraSaika_ab(int power_a, int power_b, double pos_a, double pos_b, double expo_a, double expo_b)

/**
 * @brief Compute the overlap integral between two Gaussian functions.
 * 
 * @param[in] a The first Gaussian function
 * @param[in] b The second Gaussian function
 *
 * @return The value of the overlap integral
 */
double calc_overlap(const GaussFunc<D> &a, const GaussFunc<D> &b);






} // namespace function_utils

} // namespace mrcpp