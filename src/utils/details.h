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
/**
 * @file
 * @brief Small cross-cutting utilities and helpers for MRCPP internals.
 *
 * This header declares:
 * - Filesystem and environment helpers (e.g., locating MW filter folders).
 * - Lightweight process information (Linux memory usage).
 * - Tiny array algorithms (equality/any checks, C-array → std::array conversion).
 * - Generic collection pretty-printer and a `std::array` stream operator.
 *
 * Most functions live in the `mrcpp::details` namespace to signal internal use.
 */

#include <algorithm>
#include <array>
#include <iostream>
#include <sstream>

namespace mrcpp {
/**
 * @namespace mrcpp::details
 * @brief Internal utilities; APIs may change without notice.
 */
namespace details {

/**
 * @brief Check whether a path refers to an existing directory.
 * @param path Path to check.
 * @return `true` if the path exists and is a directory, otherwise `false`.
 * @note Implementation typically uses `stat(2)` and is therefore
 *       POSIX-oriented.
 */
bool directory_exists(std::string path);

/**
 * @brief Locate the directory containing multiresolution filter files.
 *
 * The search strategy prefers an explicit environment override and then
 * compiled-in locations.
 *
 * @return Absolute/relative path to a directory with filter files.
 * @throws (implementation-defined) if no suitable directory is found.
 *
 * @details
 * The implementation checks (in order):
 *  1. Environment variable `MWFILTERS_DIR`, if set and points to a directory.
 *  2. Compiled-in source/install search paths (e.g., `mwfilters_source_dir()`,
 *     `mwfilters_install_dir()`).
 */
std::string find_filters();

/**
 * @brief Return the current process memory usage in kilobytes.
 * @return Resident (or data+stack) usage in kB, or a negative value on error.
 * @note Implemented via `/proc/self/statm`; available on Linux only.
 */
int get_memory_usage();

/**
 * @brief Check if all elements of a fixed-size array of doubles are equal.
 * @tparam D Array length.
 * @param exponent Input array.
 * @return `true` if all elements compare equal to the first element; otherwise `false`.
 * @warning Equality is tested with `==` (no tolerance).
 */
template <int D>
bool are_all_equal(const std::array<double, D> &exponent);

/**
 * @brief Test whether any element of an array equals a given value.
 * @tparam T Element type (must be equality comparable).
 * @tparam D Array length (compile-time).
 * @param col Array to scan.
 * @param eq  Value to compare against.
 * @return `true` if at least one element satisfies `element == eq`.
 * @complexity O(D).
 */
template <typename T, size_t D>
bool are_any(const std::array<T, D> &col, const T eq) {
    return std::any_of(col.cbegin(), col.cend(), [eq](const T &el) { return el == eq; });
};

/**
 * @brief Convert a C-style pointer to a fixed-size `std::array`.
 * @tparam T Element type.
 * @tparam D Number of elements to copy.
 * @param arr Pointer to at least `D` contiguous elements of type `T`.
 * @return `std::array<T, D>` with a shallow copy of the `D` elements.
 * @warning Caller is responsible for ensuring `arr` has at least `D` valid elements.
 */
template <typename T, int D>
std::array<T, D> convert_to_std_array(T *arr);

/**
 * @brief Render a collection to a compact bracketed string.
 * @tparam T A range/collection with range-for iteration and streamable elements.
 * @param coll Collection to print.
 * @return String like `"[e0, e1, ...]"`.
 * @note This utility underpins the `operator<<` overload for `std::array`.
 */
template <typename T>
auto stream_collection(const T &coll) -> std::string {
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

/**
 * @brief Stream insertion for `std::array`, producing a compact bracketed list.
 * @tparam T Element type (must be stream-insertable).
 * @tparam D Array length.
 * @param os Output stream.
 * @param coll Array to print.
 * @return Reference to @p os.
 * @sa mrcpp::details::stream_collection
 */
template <typename T, size_t D>
auto operator<<(std::ostream &os, const std::array<T, D> &coll) -> std::ostream & {
    return (os << details::stream_collection(coll));
}

} // namespace mrcpp