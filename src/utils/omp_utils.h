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
 * @brief OpenMP utilities and portability shims for MRCPP.
 *
 * This header centralizes MRCPP's interaction with OpenMP so client code can
 * compile and run both with and without OpenMP support. It provides:
 * - A consistent way to query the number of threads/rank of a thread.
 * - Lightweight lock helpers that compile away in non-OpenMP builds.
 * - A global cap on MRCPP-managed threads to avoid oversubscription with Eigen.
 *
 * @note Eigen is explicitly forced to single-threaded mode via
 *       `EIGEN_DONT_PARALLELIZE` to prevent nested parallelism when MRCPP
 *       runs OpenMP regions.
 */

/// Disable Eigen's internal multi-threading to avoid oversubscription.
#define EIGEN_DONT_PARALLELIZE

#ifdef MRCPP_HAS_OMP
  #include <omp.h>

  /**
   * @name Thread query helpers (OpenMP build)
   * @{
   */

  /// Maximum number of threads OpenMP may use for parallel regions.
  #define mrcpp_get_max_threads() omp_get_max_threads()

  /**
   * Number of threads MRCPP intends to use in parallel regions.
   *
   * @details This is capped by user/runtime policy via @ref mrcpp::set_max_threads
   * and may be lower than @c omp_get_max_threads() to respect node-level limits.
   */
  #define mrcpp_get_num_threads() mrcpp::max_threads

  /// Zero-based thread id within a running OpenMP parallel region.
  #define mrcpp_get_thread_num()  omp_get_thread_num()
  /** @} */

  /**
   * @name Lightweight lock helpers (OpenMP build)
   * @brief Macros that operate on a member `omp_lock_t omp_lock;`.
   * @{
   */
  #define MRCPP_INIT_OMP_LOCK()    omp_init_lock(&this->omp_lock)
  #define MRCPP_DESTROY_OMP_LOCK() omp_destroy_lock(&this->omp_lock)
  #define MRCPP_SET_OMP_LOCK()     omp_set_lock(&this->omp_lock)
  #define MRCPP_UNSET_OMP_LOCK()   omp_unset_lock(&this->omp_lock)
  #define MRCPP_TEST_OMP_LOCK()    omp_test_lock(&this->omp_lock)
  /** @} */

#else
  /**
   * @name Thread/query helpers (non-OpenMP build)
   * @brief Serial fallbacks so code compiles and runs without OpenMP.
   * @{
   */
  #define mrcpp_get_max_threads() 1   ///< Always 1 in serial builds.
  #define mrcpp_get_num_threads() 1   ///< Always 1 in serial builds.
  #define mrcpp_get_thread_num()  0   ///< Single thread has id 0.
  /** @} */

  /**
   * @name Lock helpers (non-OpenMP build)
   * @brief No-ops in serial builds.
   * @{
   */
  #define MRCPP_INIT_OMP_LOCK()
  #define MRCPP_DESTROY_OMP_LOCK()
  #define MRCPP_SET_OMP_LOCK()
  #define MRCPP_UNSET_OMP_LOCK()
  #define MRCPP_TEST_OMP_LOCK()
  /** @} */
#endif

namespace mrcpp {

/**
 * @brief Upper bound on threads MRCPP will request for OpenMP regions.
 *
 * @details This value is used by @c mrcpp_get_num_threads() and allows MRCPP
 * to honor node-level thread budgeting (e.g., when co-scheduled with MPI or
 * other threaded libraries). In non-OpenMP builds this remains 1.
 */
extern int max_threads;

/**
 * @brief Set the global thread cap used by MRCPP parallel regions.
 * @param threads Desired number of threads (clamped to at least 1).
 *
 * @note This does not change system-wide OpenMP settings; it only influences
 * MRCPP's internal use (e.g., via @c mrcpp_get_num_threads()).
 */
void set_max_threads(int threads);

} // namespace mrcpp