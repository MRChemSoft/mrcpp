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
 * @brief Convenience shorthand to bind a local reference @p X to the singleton ScalingCache for LegendreBasis
 */
#define getLegendreScalingCache(X) ScalingCache<LegendreBasis> &X = ScalingCache<LegendreBasis>::getInstance()

/**
 * @def getInterpolatingScalingCache(X)
 * @brief Convenience shorthand to bind a local reference @p X to the singleton ScalingCache for InterpolatingBasis
 */
#define getInterpolatingScalingCache(X)                                                                                \
    ScalingCache<InterpolatingBasis> &X = ScalingCache<InterpolatingBasis>::getInstance()

/**
 * @class ScalingCache
 * @tparam P Concrete scaling-basis type to cache (e.g., LegendreBasis, InterpolatingBasis)
 * @brief Process-wide singleton cache for scaling-basis objects keyed by polynomial order
 *
 * @details
 * Constructing a scaling basis of order \f$ k \f$ (which internally builds orthonormal polynomials,
 * evaluates them at quadrature nodes, and assembles coefficient-to-value maps) can be expensive. This
 * cache guarantees that for a given basis family @p P and polynomial order, exactly one instance is
 * constructed and then reused for the lifetime of the program.
 *
 * One singleton exists per basis family (Meyers singleton via getInstance()). First-time construction
 * is protected by MRCPP_SET_OMP_LOCK / MRCPP_UNSET_OMP_LOCK; subsequent reads are lock-free. The memory
 * estimate passed to the underlying ObjectCache is \f$ 2(k+1)^2 \cdot \texttt{sizeof(double)} \f$, which
 * approximates the two \f$ q \times q \f$ matrices (cvMap and vcMap) stored by each basis. Copy and
 * assignment are deleted to enforce singleton semantics.
 *
 * @see ObjectCache for the base class providing sparse indexed storage and memory accounting
 * @see FilterCache for the analogous cache for MWFilter objects
 */
template <class P> class ScalingCache final : public ObjectCache<P> {
public:
    /**
     * @brief Access the singleton cache for basis family @p P
     *
     * @details The instance is created on first use and lives until program exit
     */
    static ScalingCache &getInstance() {
        static ScalingCache theScalingCache;
        return theScalingCache;
    }

    /**
     * @brief Ensure the basis of a given @p order is present in the cache
     * @param order Polynomial order \f$ k \geq 0 \f$ of the basis to load
     *
     * @details If absent, constructs a new P(@p order) under an OpenMP lock and inserts it into the
     * underlying ObjectCache with a memory estimate of \f$ 2(k+1)^2 \cdot \texttt{sizeof(double)} \f$
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
     * @brief Retrieve the cached basis of a given @p order, loading it on demand if absent
     * @param order Polynomial order \f$ k \geq 0 \f$
     * @return Reference to the basis object owned by the cache
     */
    P &get(int order) {
        if (not this->hasId(order)) { load(order); }
        return ObjectCache<P>::get(order);
    }

private:
    /**
     * @brief Private constructor enforces the singleton pattern
     */
    ScalingCache() {}
    // Non-copyable and non-assignable to maintain single instance semantics.
    ScalingCache(const ScalingCache<P> &sc) = delete;
    ScalingCache<P> &operator=(const ScalingCache<P> &sc) = delete;
};

} // namespace mrcpp