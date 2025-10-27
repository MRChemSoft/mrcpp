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
 *  This file implements a small, thread-safe cache for MWFilter objects
 *  (multiwavelet filter banks) keyed by polynomial order. The cache is
 *  parameterized by filter family via the template parameter T (e.g.,
 *  Interpol or Legendre).
 *
 *  Motivation:
 *    Loading/constructing MWFilter(order, type) may involve I/O and setup.
 *    Reusing the same filter for repeated calls is faster and reduces memory
 *    churn. This cache ensures a single instance per (order, type).
 *
 *  Concurrency:
 *    - Uses MRCPP_SET_OMP_LOCK / MRCPP_UNSET_OMP_LOCK to serialize the
 *      first-time construction and insertion into the cache.
 *    - After an entry exists, get() returns it without reloading.
 *
 *  Memory accounting:
 *    - A rough memory footprint (in bytes) is computed as
 *        f->getFilter().size() * sizeof(double)
 *      and passed to the base ObjectCache for bookkeeping/eviction policy.
 *
 *  Template parameter T:
 *    - Must be a valid family tag (Interpol or Legendre).
 *    - The explicit instantiations at the end of the file make sure code for
 *      these two variants is emitted by the compiler.
 */

#include "FilterCache.h"
#include "utils/Printer.h"

#include "MRCPP/constants.h"

using namespace Eigen;

namespace mrcpp {

/*
 * Constructor
 * -----------
 * Determine the runtime 'type' field from the compile-time template parameter T.
 * If T is not a recognized family, emit an error. Valid values are:
 *   - Interpol : interpolatory multiwavelet family
 *   - Legendre : Legendre multiwavelet family
 */
template <int T> FilterCache<T>::FilterCache() {
    switch (T) {
        case (Interpol):
            type = Interpol;
            break;
        case (Legendre):
            type = Legendre;
            break;
        default:
            MSG_ERROR("Invalid filter type: " << T);
    }
}

/*
 * load(order)
 * -----------
 * Ensure that an MWFilter for the given 'order' exists in the cache. If not,
 * construct it and insert it along with a memory estimate.
 *
 * Steps:
 *   1) Acquire OpenMP lock to prevent concurrent insertions.
 *   2) Check presence via hasId(order). If absent:
 *        - Allocate MWFilter(order, type).
 *        - Compute 'memo' as (#elements) * sizeof(double).
 *        - Insert into base ObjectCache keyed by 'order'.
 *   3) Release the lock.
 */
template <int T> void FilterCache<T>::load(int order) {
    MRCPP_SET_OMP_LOCK();
    if (not hasId(order)) {
        auto *f = new MWFilter(order, type);
        int memo = f->getFilter().size() * sizeof(double);
        ObjectCache<MWFilter>::load(order, f, memo);
    }
    MRCPP_UNSET_OMP_LOCK();
}

/*
 * get(order)
 * ----------
 * Retrieve a reference to the cached MWFilter for 'order'; if it doesn't
 * exist yet, load() is called lazily. The reference is owned by the cache.
 */
template <int T> MWFilter &FilterCache<T>::get(int order) {
    if (not hasId(order)) { load(order); }
    return ObjectCache<MWFilter>::get(order);
}

/*
 * getFilterMatrix(order)
 * ----------------------
 * Convenience accessor: returns a const reference to the underlying filter
 * matrix for the requested 'order'. Triggers lazy load if necessary.
 *
 * Notes:
 *  - MWFilter::getFilter() is expected to return an Eigen::MatrixXd (or
 *    compatible type) containing the filter taps laid out as used elsewhere
 *    in MRCPP.
 */
template <int T> const MatrixXd &FilterCache<T>::getFilterMatrix(int order) {
    if (not hasId(order)) { load(order); }
    return ObjectCache<MWFilter>::get(order).getFilter();
}

/*
 * Explicit template instantiations
 * --------------------------------
 * Instantiate the cache for the two standard families so clients can link
 * against these symbols without needing to compile this TU with their T.
 */
template class FilterCache<Interpol>;
template class FilterCache<Legendre>;

} // namespace mrcpp
