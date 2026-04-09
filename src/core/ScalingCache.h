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

#include "ObjectCache.h"

namespace mrcpp {

/**
 * @def getLegendreScalingCache(X)
 * @brief Convenience macro to bind a local reference @p X to the singleton
 *        ScalingCache specialized for LegendreBasis.
 *
 * Usage:
 *   getLegendreScalingCache(cache);
 *   auto& B = cache.get(order);
 */
#define getLegendreScalingCache(X) ScalingCache<LegendreBasis> &X = ScalingCache<LegendreBasis>::getInstance()

/**
 * @def getInterpolatingScalingCache(X)
 * @brief Convenience macro to bind a local reference @p X to the singleton
 *        ScalingCache specialized for InterpolatingBasis.
 *
 * Usage:
 *   getInterpolatingScalingCache(cache);
 *   auto& B = cache.get(order);
 */
#define getInterpolatingScalingCache(X)                                                                                \
    ScalingCache<InterpolatingBasis> &X = ScalingCache<InterpolatingBasis>::getInstance()

/**
 * @class ScalingCache
 * @tparam P  A concrete scaling-basis type (e.g., LegendreBasis, InterpolatingBasis).
 * @brief Thread-safe singleton cache for scaling bases keyed by polynomial order.
 *
 * Motivation
 * ----------
 * Constructing a scaling basis of order `k` (which internally prepares
 * polynomials, quadrature-derived maps, etc.) can be relatively expensive.
 * This cache guarantees that for a given template parameter P (basis family)
 * and a given order, exactly one instance is created and then reused.
 *
 * Design
 * ------
 * - Inherits from @ref ObjectCache<P>, which provides sparse indexed storage,
 *   memory accounting, and basic get/load/unload primitives.
 * - Singleton per `P` (Meyers singleton via getInstance()) so that all parts
 *   of the program share the same cache for the same basis family.
 * - Thread-safety: the first-time construction/insertion is protected by
 *   MRCPP_SET_OMP_LOCK / MRCPP_UNSET_OMP_LOCK. Reads after presence are fast.
 *
 * Memory accounting
 * -----------------
 * The `memo` value passed to ObjectCache is a *rough* byte estimate:
 *   memo ≈ 2 * (order+1)^2 * sizeof(double)
 * The constant factor “2” approximates two q×q matrices stored by a basis
 * (e.g., cvMap and vcMap), where q = order + 1. It is not an exact footprint,
 * but suffices for simple bookkeeping.
 */
template <class P> class ScalingCache final : public ObjectCache<P> {
public:
    /**
     * @brief Access the singleton instance for the template parameter P.
     *
     * One instance per concrete basis family exists process-wide.
     */
    static ScalingCache &getInstance() {
        static ScalingCache theScalingCache;
        return theScalingCache;
    }

    /**
     * @brief Ensure the basis of a given @p order is present in the cache.
     *
     * If absent, constructs a new P(order) under an OpenMP lock and inserts it
     * into the underlying ObjectCache with a rough memory estimate.
     */
    void load(int order) {
        MRCPP_SET_OMP_LOCK();
        if (not this->hasId(order)) {
            P *f = new P(order);
            int memo = 2 * SQUARE(order + 1) * sizeof(double); // approx
            ObjectCache<P>::load(order, f, memo);
        }
        MRCPP_UNSET_OMP_LOCK();
    }

    /**
     * @brief Retrieve the cached basis of a given @p order (lazy-loads if needed).
     * @return Reference to the basis object owned by the cache.
     */
    P &get(int order) {
        if (not this->hasId(order)) { load(order); }
        return ObjectCache<P>::get(order);
    }

private:
    /// Private constructor enforces the singleton pattern.
    ScalingCache() {}
    // Non-copyable / non-assignable.
    ScalingCache(const ScalingCache<P> &sc) = delete;
    ScalingCache<P> &operator=(const ScalingCache<P> &sc) = delete;
};

} // namespace mrcpp