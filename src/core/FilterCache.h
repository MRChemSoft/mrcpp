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
 * FilterCache provides a process-wide cache for multiwavelet filter banks
 * (MWFilter) so that the same filter for a given polynomial order is created
 * and loaded exactly once and then reused. This avoids repeated I/O and setup.
 *
 * Design highlights:
 *  - There are different *families* of filters (e.g. Legendre vs Interpolating).
 *    We want caches for both, alive simultaneously. To achieve this, the
 *    concrete cache is a class template FilterCache<T>, where T encodes the
 *    family. Each T gets its own singleton instance.
 *
 *  - The cache is keyed by the *order* (polynomial order k). Loading an entry
 *    constructs MWFilter(order, type) and stores it internally for reuse.
 *
 *  - Thread-safety and the actual load/get logic are implemented in the .cpp
 *    using OpenMP locks (MRCPP_SET_OMP_LOCK / MRCPP_UNSET_OMP_LOCK).
 *
 * About this header:
 *  - Declares a tiny abstract façade (BaseFilterCache) to allow use via a
 *    non-templated base pointer/reference when the family is not known at
 *    compile time.
 *  - Declares the templated FilterCache<T> singleton with the minimal API:
 *      load(order), get(order), and getFilterMatrix(order).
 */

#pragma once

#include "MWFilter.h"
#include "ObjectCache.h"

#include <iostream>
#include <string>

namespace mrcpp {

/**
 * @def getFilterCache(T, X)
 * @brief Create a named reference @p X bound to the singleton FilterCache<T>.
 *
 * Usage:
 *   getFilterCache(Interpol, cache);
 *   const auto& H = cache.getFilterMatrix(order);
 *
 * @def getLegendreFilterCache(X)
 * @brief Convenience macro for FilterCache<Legendre>.
 *
 * @def getInterpolatingFilterCache(X)
 * @brief Convenience macro for FilterCache<Interpol>.
 */
#define getFilterCache(T, X) FilterCache<T> &X = FilterCache<T>::getInstance()
#define getLegendreFilterCache(X) FilterCache<Legendre> &X = FilterCache<Legendre>::getInstance()
#define getInterpolatingFilterCache(X) FilterCache<Interpol> &X = FilterCache<Interpol>::getInstance()

/**
 * @class BaseFilterCache
 * @brief Abstract façade over the templated filter cache.
 *
 * Rationale:
 *   Callers that do not know the filter family T at compile time can still
 *   interact with a cache through this non-templated interface. Concrete
 *   implementations are provided by FilterCache<T>.
 *
 * Notes:
 *   - Inherits from ObjectCache<MWFilter> to reuse generic cache plumbing.
 *   - Pure virtual methods delegate to the concrete implementation in
 *     FilterCache<T>.
 */
class BaseFilterCache : public ObjectCache<MWFilter> {
public:
    /// Ensure the filter for @p order exists in the cache (lazy load if needed).
    void load(int order) override = 0;

    /// Retrieve the cached MWFilter for @p order (loads it on demand).
    MWFilter &get(int order) override = 0;

    /// Convenience accessor: return the filter matrix (const) for @p order.
    virtual const Eigen::MatrixXd &getFilterMatrix(int order) = 0;
};

/**
 * @class FilterCache
 * @tparam T  Integer tag selecting the filter family (e.g., Interpol, Legendre).
 * @brief Singleton cache of MWFilter objects for a specific filter family.
 *
 * Key properties:
 *  - One singleton instance per family T (Meyers singleton via getInstance()).
 *  - API mirrors BaseFilterCache and ObjectCache<MWFilter>.
 *  - The constructor is private to enforce the singleton pattern.
 *  - Copy/assignment are deleted to prevent accidental duplication.
 *
 * Thread-safety:
 *  - The .cpp guards first-time loads with OpenMP locks.
 *  - Reads after an entry exists are lock-free through the base cache API.
 */
template <int T> class FilterCache final : public BaseFilterCache {
public:
    /**
     * @brief Access the singleton cache for the template family T.
     *
     * The instance is created on first use and lives until program exit.
     */
    static FilterCache &getInstance() {
        static FilterCache theFilterCache;
        return theFilterCache;
    }

    /// Ensure entry for @p order exists; loads it if missing (see .cpp).
    void load(int order) override;

    /// Retrieve the MWFilter for @p order; loads it if missing.
    MWFilter &get(int order) override;

    /// Convenience accessor returning a const reference to the filter matrix.
    const Eigen::MatrixXd &getFilterMatrix(int order) override;

protected:
    /**
     * @brief Runtime family/type code corresponding to template parameter T.
     *
     * Initialized in the private constructor; used to construct MWFilter(order, type).
     */
    int type;

private:
    /**
     * @brief Private constructor enforces the singleton pattern.
     *
     * Sets #type based on T and performs any minimal family-specific setup.
     * (Validation happens in the .cpp.)
     */
    FilterCache();

    // Non-copyable and non-assignable to maintain single instance semantics.
    FilterCache(FilterCache<T> const &fc) = delete;
    FilterCache &operator=(FilterCache<T> const &fc) = delete;
};

} // namespace mrcpp