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

#pragma once

#include <array>
#include <iostream>
#include <sstream>

namespace mrcpp {
namespace details {
int get_memory_usage();
template <int D> bool are_all_equal(const std::array<double, D> &exponent);
template <typename T, int D> std::array<T, D> convert_to_std_array(T *arr);
template <typename T> auto stream_collection(const T &coll) -> std::string {
    std::ostringstream os;
    os << "[";
    for (auto elem : coll) {
        os << elem;
        if (elem != coll.back()) os << ", ";
    }
    os << "]";
    return os.str();
}
} // namespace details

template <typename T> auto operator<<(std::ostream &os, const std::array<T, 1> &coll) -> std::ostream & {
    return (os << details::stream_collection(coll));
}

template <typename T> auto operator<<(std::ostream &os, const std::array<T, 2> &coll) -> std::ostream & {
    return (os << details::stream_collection(coll));
}

template <typename T> auto operator<<(std::ostream &os, const std::array<T, 3> &coll) -> std::ostream & {
    return (os << details::stream_collection(coll));
}
} // namespace mrcpp
