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

#include "CrossCorrelation.h"
#include "ObjectCache.h"

#include <iostream>
#include <string>

namespace mrcpp {

/**
 * @def getCrossCorrelationCache(T, X)
 * @brief Convenience macro to obtain a named reference to the singleton cache.
 *
 * Expands to:
 *   CrossCorrelationCache<T> &X = CrossCorrelationCache<T>::getInstance()
 *
 * Example:
 *   getCrossCorrelationCache(Interpol, ccc);
 *   const auto& L = ccc.getLMatrix(order);
 */
#define getCrossCorrelationCache(T, X) CrossCorrelationCache<T> &X = CrossCorrelationCache<T>::getInstance()

/**
 * @class CrossCorrelationCache
 * @brief Thread-safe cache for @ref CrossCorrelation objects, keyed by order.
 *
 * This cache avoids repeatedly loading the (potentially large) left/right
 * cross-correlation matrices from disk. One cache instance exists per filter
 * family, realized as a template parameter @p T (e.g., Interpol or Legendre).
 *
 * Design notes:
 *  - Singleton pattern (Meyers singleton) per @p T via getInstance().
 *  - Inherits from @ref ObjectCache<CrossCorrelation>, which provides the
 *    generic cache interface (load/get/hasId etc.).
 *  - Actual loading and synchronization details are implemented in the
 *    corresponding .cpp; OpenMP locks guard first-time insertions.
 *
 * @tparam T Filter family tag (int constant), e.g. Interpol or Legendre.
 */
template <int T> class CrossCorrelationCache final : public ObjectCache<CrossCorrelation> {
public:
    /**
     * @brief Access the unique cache instance for the template family @p T.
     *
     * Uses a function-local static (Meyers singleton). Thread-safe in C++11+.
     */
    static CrossCorrelationCache<T> &getInstance() {
        static CrossCorrelationCache<T> theCrossCorrelationCache;
        return theCrossCorrelationCache;
    }

    /**
     * @brief Ensure that the entry for @p order is present in the cache.
     *
     * If absent, constructs a new @ref CrossCorrelation(order, type) and
     * inserts it. See .cpp for locking and memory accounting.
     */
    void load(int order) override;

    /**
     * @brief Retrieve the cached @ref CrossCorrelation for @p order.
     *
     * Loads on demand if missing. Returns a reference owned by the cache.
     */
    CrossCorrelation &get(int order) override;

    /**
     * @brief Convenience accessor for the Left matrix of a given order.
     *
     * Triggers lazy load if needed, then returns a const reference.
     */
    const Eigen::MatrixXd &getLMatrix(int order);

    /**
     * @brief Convenience accessor for the Right matrix of a given order.
     *
     * Triggers lazy load if needed, then returns a const reference.
     */
    const Eigen::MatrixXd &getRMatrix(int order);

    /**
     * @brief Filter family/type code associated with this cache.
     *
     * Set in the private constructor based on the template parameter @p T.
     * (E.g., Interpol or Legendre.)
     */
    int getType() const { return this->type; }

protected:
    /**
     * @brief Filter family/type code (matches template parameter @p T).
     */
    int type;

    /**
     * @brief Base path to filter/correlation library on disk.
     *
     * Reserved for potential use by loaders. Actual path resolution is
     * currently handled inside CrossCorrelation (see details::find_filters()).
     */
    std::string libPath; ///< Base path to filter library

private:
    /**
     * @brief Private constructor enforces the singleton pattern.
     *
     * Initializes @ref type based on T; see .cpp for validation.
     */
    CrossCorrelationCache();

    // Non-copyable / non-assignable — keeps the singleton unique.
    CrossCorrelationCache(CrossCorrelationCache<T> const &ccc) = delete;
    CrossCorrelationCache<T> &operator=(CrossCorrelationCache<T> const &ccc) = delete;
};

} // namespace mrcpp