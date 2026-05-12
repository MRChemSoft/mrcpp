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

#include "MWFilter.h"
#include "ObjectCache.h"

#include <iostream>
#include <string>

namespace mrcpp {

/**
 * @def getFilterCache(T, X)
 * @brief Bind a local reference @p X to the singleton FilterCache<T>
 *
 * @details Expands to <tt>FilterCache<T> &X = FilterCache<T>::getInstance()</tt>
 */

/**
 * @def getLegendreFilterCache(X)
 * @brief Convenience shorthand for @ref getFilterCache with the Legendre family
 */

/**
 * @def getInterpolatingFilterCache(X)
 * @brief Convenience shorthand for @ref getFilterCache with the Interpol family
 */
#define getFilterCache(T, X) FilterCache<T> &X = FilterCache<T>::getInstance()
#define getLegendreFilterCache(X) FilterCache<Legendre> &X = FilterCache<Legendre>::getInstance()
#define getInterpolatingFilterCache(X) FilterCache<Interpol> &X = FilterCache<Interpol>::getInstance()

/**
 * @class BaseFilterCache
 * @brief Non-templated abstract interface over the filter cache hierarchy
 *
 * @details Callers that do not know the filter family at compile time can interact with any FilterCache<T>
 * through this base class. Inherits from ObjectCache<MWFilter> for generic index-addressed storage,
 * and declares pure virtual load(), get(), and getFilterMatrix() for family-specific implementations.
 */
class BaseFilterCache : public ObjectCache<MWFilter> {
public:
    /// Ensure the filter for @p order exists in the cache (lazy load if needed)
    void load(int order) override = 0;

    /// Retrieve the cached MWFilter for @p order (loads it on demand)
    MWFilter &get(int order) override = 0;

    /// Convenience accessor: return the filter matrix (const) for @p order
    virtual const Eigen::MatrixXd &getFilterMatrix(int order) = 0;
};

/**
 * @class FilterCache
 * @tparam T Integer tag selecting the filter family (Interpol or Legendre)
 * @brief Process-wide singleton cache for MWFilter objects keyed by polynomial order
 *
 * @details One instance exists per family @p T (Meyers singleton via getInstance()). First-time loads
 * construct MWFilter(order, type) under an OpenMP lock; subsequent reads are lock-free. Copy and
 * assignment are deleted to enforce singleton semantics.
 *
 * @see MWFilter for the cached object type
 * @see BaseFilterCache for the non-templated interface
 */
template <int T> class FilterCache final : public BaseFilterCache {
public:
    /**
     * @brief Access the singleton cache for the template family T
     *
     * The instance is created on first use and lives until program exit
     */
    static FilterCache &getInstance() {
        static FilterCache theFilterCache;
        return theFilterCache;
    }

    /// Ensure entry for @p order exists; loads it if missing (see .cpp)
    void load(int order) override;

    /// Retrieve the MWFilter for @p order; loads it if missing
    MWFilter &get(int order) override;

    /// Convenience accessor returning a const reference to the filter matrix
    const Eigen::MatrixXd &getFilterMatrix(int order) override;

protected:
    int type; ///< Runtime family/type code derived from template parameter @p T; passed to MWFilter(order, type)

private:
    /**
     * @brief Private constructor enforces the singleton pattern
     *
     * @details Sets #type from @p T and performs minimal family-specific validation (see .cpp)
     */
    FilterCache();

    // Non-copyable and non-assignable to maintain single instance semantics.
    FilterCache(FilterCache<T> const &fc) = delete;
    FilterCache &operator=(FilterCache<T> const &fc) = delete;
};

} // namespace mrcpp