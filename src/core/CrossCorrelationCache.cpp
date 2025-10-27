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

/*
 *  File purpose (high level)
 *  -------------------------
 *  This file implements a small, thread-safe cache for CrossCorrelation
 *  objects, parameterized on the filter family (template parameter T).
 *
 *  Motivation:
 *    CrossCorrelation(order, type) loads two dense (K*K)×(2K) matrices
 *    from binary files. Loading them repeatedly is expensive. This cache
 *    stores one instance per (order, type) and returns references to it.
 *
 *  Template parameter T:
 *    - Must be one of the family tags (e.g. Interpol, Legendre).
 *    - The explicit instantiations at the end fix T to these two values.
 *
 *  Concurrency:
 *    - The cache uses MRCPP_SET_OMP_LOCK / MRCPP_UNSET_OMP_LOCK to guard the
 *      critical section that performs the initial load and insertion.
 *    - Once loaded, subsequent get() calls read from the cache without
 *      reloading (fast path).
 */

#include "CrossCorrelationCache.h"
#include "utils/Printer.h"

#include "MRCPP/constants.h"

using namespace Eigen;

namespace mrcpp {

/*
 * Constructor
 * -----------
 * Initialize the runtime 'type' field from the compile-time template
 * parameter T, and validate that it matches a known family.
 * If T is invalid, emit an error.
 */
template <int T> CrossCorrelationCache<T>::CrossCorrelationCache() {
    switch (T) {
        case (Interpol):
            type = Interpol;
            break;
        case (Legendre):
            type = Legendre;
            break;
        default:
            MSG_ERROR("Invalid CrossCorrelation type: " << T)
    }
}

/*
 * load(order)
 * -----------
 * Ensure that a CrossCorrelation for the given 'order' exists in the cache.
 * If not present, construct it and insert it with a memory budget hint.
 *
 * Steps:
 *  1) Acquire OpenMP lock (thread-safe insertion).
 *  2) If key 'order' is absent, allocate a new CrossCorrelation(order, type).
 *  3) Compute a crude memory footprint 'memo' for cache accounting:
 *       - getLMatrix().size() returns (#rows * #cols)
 *       - Multiply by 2 because we also store a Right matrix of same size
 *       - Multiply by sizeof(double) to get bytes
 *  4) Insert into the underlying ObjectCache keyed by 'order'.
 *  5) Release lock.
 */
template <int T> void CrossCorrelationCache<T>::load(int order) {
    MRCPP_SET_OMP_LOCK();
    if (not hasId(order)) {
        auto *ccc = new CrossCorrelation(order, type);
        int memo = ccc->getLMatrix().size() * 2 * sizeof(double);
        ObjectCache<CrossCorrelation>::load(order, ccc, memo);
    }
    MRCPP_UNSET_OMP_LOCK();
}

/*
 * get(order)
 * ----------
 * Return a reference to the cached CrossCorrelation for 'order'.
 * If missing, it will be loaded on-demand (calling load()).
 */
template <int T> CrossCorrelation &CrossCorrelationCache<T>::get(int order) {
    if (not hasId(order)) { load(order); }
    return ObjectCache<CrossCorrelation>::get(order);
}

/*
 * getLMatrix(order)
 * -----------------
 * Convenience accessor: returns a const reference to the Left matrix
 * for the requested 'order', auto-loading it if necessary.
 */
template <int T> const Eigen::MatrixXd &CrossCorrelationCache<T>::getLMatrix(int order) {
    if (not hasId(order)) { load(order); }
    return ObjectCache<CrossCorrelation>::get(order).getLMatrix();
}

/*
 * getRMatrix(order)
 * -----------------
 * Convenience accessor: returns a const reference to the Right matrix
 * for the requested 'order', auto-loading it if necessary.
 */
template <int T> const Eigen::MatrixXd &CrossCorrelationCache<T>::getRMatrix(int order) {
    if (not hasId(order)) { load(order); }
    return ObjectCache<CrossCorrelation>::get(order).getRMatrix();
}

/*
 * Explicit template instantiations
 * --------------------------------
 * Build concrete cache types for the known families:
 *   - CrossCorrelationCache<Interpol>
 *   - CrossCorrelationCache<Legendre>
 *
 * This ensures the compiler generates code for these two variants in
 * this translation unit, so users can link against them.
 */
template class CrossCorrelationCache<Interpol>;
template class CrossCorrelationCache<Legendre>;

} // namespace mrcpp