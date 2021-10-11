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

#include <algorithm>
#include <array>
#include <iostream>
#include <sstream>

namespace mrcpp {
namespace details {
bool directory_exists(std::string path);
std::string find_filters();
int get_memory_usage();
template <int D> bool are_all_equal(const std::array<double, D> &exponent);
template <typename T, size_t D> bool are_any(const std::array<T, D> &col, const T eq) {
    return std::any_of(col.cbegin(), col.cend(), [eq](const T &el) { return el == eq; });
};
template <typename T, int D> std::array<T, D> convert_to_std_array(T *arr);
template <typename T> auto stream_collection(const T &coll) -> std::string {
    std::ostringstream os;
    bool first = true;
    os << "[";
    for (auto elem : coll) {
        if (!first) os << ", ";
        os << elem;
        first = false;
    }
    os << "]";
    return os.str();
}
} // namespace details

template <typename T, size_t D> auto operator<<(std::ostream &os, const std::array<T, D> &coll) -> std::ostream & {
    return (os << details::stream_collection(coll));
}
} // namespace mrcpp
