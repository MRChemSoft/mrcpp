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
 * Overview
 * --------
 * Implementation of the QuadratureCache singleton. This cache wraps
 * ObjectCache<GaussQuadrature> to manage Gauss–Legendre quadrature rules by
 * integer key `k` (the rule's order). It additionally tracks a *global*
 * integration domain [A,B] and a number of equal sub-intervals `intervals`
 * that should be applied to *all* cached rules.
 *
 * Responsibilities:
 *  - Provide default domain settings ([0,1], intervals=1).
 *  - Lazy-load GaussQuadrature objects for requested orders (load/get).
 *  - Propagate changes to [A,B] or `intervals` to any already-cached rules.
 *
 * Thread-safety:
 *  - First-time loads are guarded by MRCPP_SET_OMP_LOCK / MRCPP_UNSET_OMP_LOCK.
 *    Once an object is present in the cache, read access is lock-free.
 *
 * Memory accounting:
 *  - The `memo` passed to the base cache is a rough estimate: 2 * k * sizeof(double).
 *    (This is intentionally approximate and used only for coarse bookkeeping.)
 *
 * Notes:
 *  - Iteration over cached ids uses `for (int i = 0; i < getNObjs(); ++i)`.
 *    Slots may be empty (not all ids 0..high-water-mark are loaded).
 *  - Potential typo/bug (left as-is by request): in setIntervals(), the input
 *    validity check tests `this->intervals < 1` instead of `ivals < 1`.
 */

/*
 *
 *
 *  \date Jul 26, 2009
 *  \author Jonas Juselius <jonas.juselius@uit.no> \n
 *          CTCC, University of Tromsø
 *
 * \breif
 */

#include "QuadratureCache.h"
#include "utils/Printer.h"

namespace mrcpp {

/**
 * @brief Construct the cache with default domain and replication settings.
 *
 * Defaults:
 *  - A = 0.0, B = 1.0  → unit interval [0,1]
 *  - intervals = 1     → no subdivision (composite quadrature disabled)
 *
 * Actual GaussQuadrature objects are created lazily on demand in load().
 */
QuadratureCache::QuadratureCache() {
    this->A = 0.0;
    this->B = 1.0;
    this->intervals = 1;
}

/** @brief Trivial destructor; owned objects are freed by the base cache. */
QuadratureCache::~QuadratureCache() = default;

/**
 * @brief Ensure a GaussQuadrature of order k is present in the cache.
 *
 * Under the OMP lock:
 *  - If absent, allocate a new GaussQuadrature(k, A, B, intervals),
 *    compute a rough memory estimate, and insert it into ObjectCache.
 */
void QuadratureCache::load(int k) {
    MRCPP_SET_OMP_LOCK();
    if (not hasId(k)) {
        auto *gp = new GaussQuadrature(k, this->A, this->B, this->intervals);
        int memo = 2 * k * sizeof(double); // rough accounting only
        ObjectCache<GaussQuadrature>::load(k, gp, memo);
    }
    MRCPP_UNSET_OMP_LOCK();
}

/**
 * @brief Retrieve a reference to the cached quadrature of order k.
 *        Lazily loads it if not present yet.
 */
GaussQuadrature &QuadratureCache::get(int k) {
    if (not hasId(k)) { load(k); }
    return ObjectCache<GaussQuadrature>::get(k);
}

/**
 * @brief Update the global integration bounds to [a,b] and propagate the
 *        change to all already-cached GaussQuadrature objects.
 *
 * Behavior:
 *  - If the new bounds are effectively identical (within MachineZero), do nothing.
 *  - Otherwise, set A/B and iterate over existing ids; for each loaded entry,
 *    call .setBounds(a,b) so its scaled nodes/weights are rebuilt.
 */
void QuadratureCache::setBounds(double a, double b) {
    if (std::abs(this->A - a) < MachineZero and std::abs(this->B - b) < MachineZero) { return; }
    if (a >= b) { MSG_ERROR("Invalid Gauss interval, a > b."); }
    this->A = a;
    this->B = b;
    for (int i = 0; i < getNObjs(); i++) {
        if (hasId(i)) { ObjectCache<GaussQuadrature>::get(i).setBounds(a, b); }
    }
}

/**
 * @brief Update the global number of equal sub-intervals and propagate to all
 *        already-cached rules.
 *
 * Behavior:
 *  - If unchanged, return early.
 *  - Sanity check: intervals must be ≥ 1.
 *  - Iterate over existing ids; for each loaded entry, call .setIntervals(ivals).
 *
 * Note:
 *  - The input validity test uses `this->intervals < 1` (likely intended to be
 *    `ivals < 1`). Left unchanged intentionally.
 */
void QuadratureCache::setIntervals(int ivals) {
    if (ivals == this->intervals) { return; }
    if (this->intervals < 1) { MSG_ERROR("Invalid number of intervals, intervals < 1"); }
    for (int i = 0; i < getNObjs(); i++) {
        if (hasId(i)) { ObjectCache<GaussQuadrature>::get(i).setIntervals(ivals); }
    }
}

} // namespace mrcpp
