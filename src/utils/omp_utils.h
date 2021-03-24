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

#define EIGEN_DONT_PARALLELIZE

#ifdef MRCPP_HAS_OMP
#include <omp.h>
#define mrcpp_get_max_threads() omp_get_max_threads()
#define mrcpp_get_num_threads() mrcpp::max_threads
#define mrcpp_get_thread_num() omp_get_thread_num()
#define MRCPP_INIT_OMP_LOCK() omp_init_lock(&this->omp_lock)
#define MRCPP_DESTROY_OMP_LOCK() omp_destroy_lock(&this->omp_lock)
#define MRCPP_SET_OMP_LOCK() omp_set_lock(&this->omp_lock)
#define MRCPP_UNSET_OMP_LOCK() omp_unset_lock(&this->omp_lock)
#define MRCPP_TEST_OMP_LOCK() omp_test_lock(&this->omp_lock)
#else
#define mrcpp_get_max_threads() 1
#define mrcpp_get_num_threads() 1
#define mrcpp_get_thread_num() 0
#define MRCPP_INIT_OMP_LOCK()
#define MRCPP_DESTROY_OMP_LOCK()
#define MRCPP_SET_OMP_LOCK()
#define MRCPP_UNSET_OMP_LOCK()
#define MRCPP_TEST_OMP_LOCK()
#endif

namespace mrcpp {
extern int max_threads;
void set_max_threads(int threads);
} // namespace mrcpp
