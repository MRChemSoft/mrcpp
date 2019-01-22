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

#include <algorithm>
#include <array>

namespace mrcpp {
namespace details {


/** @brief checks if all elements of an array of doubles are equal */
template<int D>
bool are_all_equal(const std::array<double, D> &exponent) {
        return std::all_of(exponent.begin(), exponent.end(),
               [ex = std::begin(exponent)](double i) {return i == *ex; });
}

/** @brief converts c_type arrays to std::arrays */
template<typename T, int D>
std::array<T, D>  convert_to_std_array(T *arr) {
    auto ret_arr = std::array<T, D>{};
    for (auto d = 0; d < D; d++) {
        ret_arr[d] = arr[d];
    }
    return ret_arr;
}

template bool are_all_equal<1>(const std::array<double, 1> &exponent);
template bool are_all_equal<2>(const std::array<double, 2> &exponent);
template bool are_all_equal<3>(const std::array<double, 3> &exponent);

template std::array<double, 1> convert_to_std_array<double, 1>(double *arr);
template std::array<double, 2> convert_to_std_array<double, 2>(double *arr);
template std::array<double, 3> convert_to_std_array<double, 3>(double *arr);

template std::array<int, 1> convert_to_std_array<int, 1>(int *arr);
template std::array<int, 2> convert_to_std_array<int, 2>(int *arr);
template std::array<int, 3> convert_to_std_array<int, 3>(int *arr);
} // namespace details
} // namespace mrcpp
