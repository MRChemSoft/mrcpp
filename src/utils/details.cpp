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

#include "details.h"

#include <algorithm>
#include <array>
#include <cstdlib>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>

#include "MRCPP/config.h"
#include "utils/Printer.h"

namespace mrcpp {
namespace details {
bool directory_exists(std::string path) {
    struct stat info;

    int statRC = stat(path.c_str(), &info);
    if (statRC != 0) {
        if (errno == ENOENT) { return 0; }  // something along the path does not exist
        if (errno == ENOTDIR) { return 0; } // something in path prefix is not a dir
        MSG_ABORT("Non-existent directory " + path);
    }

    return (info.st_mode & S_IFDIR) ? true : false;
}

std::string find_filters() {
    std::string filters;
    auto envvar = "";
    if (const char *env_p = std::getenv("MWFILTERS_DIR")) envvar = env_p;
    for (auto n : {envvar, mwfilters_source_dir(), mwfilters_install_dir()}) {
        if (details::directory_exists(n)) {
            filters = n;
            break;
        }
    }
    
    if (filters.empty()) {
        MSG_ABORT("Could not find a folder containing filters!");
    } else {
        return filters;
    }
}

// helper function: parse a string and returns the nth integer number
int get_val(char *line, int n) {
    char *p = line;
    int len = 0;
    for (int i = 0; i < n - 1; i++) {
        // jump over n-1 first numbers
        while (*p < '0' || *p > '9') p++;
        while (*p >= '0' && *p <= '9') p++;
    }
    while (*p < '0' || *p > '9') p++;
    char *s = p;
    while (*s >= '0' && *s <= '9') {
        s++;
        line[len] = p[len];
        len++;
    }
    p[len] = 0;
    return atoi(p);
}

/** @brief returns the current memory usage of this process, in kB */
int get_memory_usage() {
    int mem_val = -1;
    FILE *file = fopen("/proc/self/statm", "r");
    if (file != nullptr) {
        char line[80];
        while (fgets(line, 80, file) != nullptr) {
            mem_val = 4 * get_val(line, 6); // sixth number is data+stack in pages (4kB)
        }
        fclose(file);
    }
    return mem_val;
}

/** @brief checks if all elements of an array of doubles are equal */
template <int D> bool are_all_equal(const std::array<double, D> &exponent) {
    return std::all_of(exponent.begin(), exponent.end(), [ex = std::begin(exponent)](double i) { return i == *ex; });
}

/** @brief converts c_type arrays to std::arrays */
template <typename T, int D> std::array<T, D> convert_to_std_array(T *arr) {
    auto ret_arr = std::array<T, D>{};
    for (auto d = 0; d < D; d++) { ret_arr[d] = arr[d]; }
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
