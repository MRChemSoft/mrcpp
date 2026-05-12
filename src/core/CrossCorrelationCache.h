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
 * @brief Bind a local reference @p X to the singleton CrossCorrelationCache<T>
 *
 * @details Expands to <tt>CrossCorrelationCache<T> &X = CrossCorrelationCache<T>::getInstance()</tt>
 * Example: <tt>getCrossCorrelationCache(Interpol, ccc); const auto& L = ccc.getLMatrix(order);</tt>
 */
#define getCrossCorrelationCache(T, X) CrossCorrelationCache<T> &X = CrossCorrelationCache<T>::getInstance()

/**
 * @class CrossCorrelationCache
 * @brief Thread-safe cache for @ref CrossCorrelation objects, keyed by order
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
 * @tparam T Filter family tag (int constant), e.g. Interpol or Legendre
 */
template <int T> class CrossCorrelationCache final : public ObjectCache<CrossCorrelation> {
public:
    /**
     * @brief Access the unique cache instance for the template family @p T
     *
     * Uses a function-local static (Meyers singleton). Thread-safe in C++11+.
     */
    static CrossCorrelationCache<T> &getInstance() {
        static CrossCorrelationCache<T> theCrossCorrelationCache;
        return theCrossCorrelationCache;
    }

    /**
     * @brief Ensure that the entry for @p order is present in the cache
     *
     * If absent, constructs a new @ref CrossCorrelation(order, type) and
     * inserts it. See .cpp for locking and memory accounting.
     */
    void load(int order) override;

    /**
     * @brief Retrieve the cached @ref CrossCorrelation for @p order
     *
     * Loads on demand if missing. Returns a reference owned by the cache.
     */
    CrossCorrelation &get(int order) override;

    /**
     * @brief Convenience accessor returning the Left cross-correlation matrix for a given order
     * @param order Polynomial order k
     * @return Const reference to the left matrix (lazy-loads the entry if absent)
     */
    const Eigen::MatrixXd &getLMatrix(int order);

    /**
     * @brief Convenience accessor returning the Right cross-correlation matrix for a given order
     * @param order Polynomial order \f$ k \f$
     * @return Const reference to the right cross-correlation coefficient matrix (lazy-loads if absent)
     *
     * @details The right matrix collects the cross-correlation coefficients
     * \f[
     *   C^{(+)}_{ijp} = \int_0^1 dz \int_0^1 dx\, \phi_i(x)\, \phi_j(x-z)\, \phi_p(z)
     * \f]
     * with \f$ i,j = 0,\ldots,k \f$ and \f$ p = 0,\ldots,2k+1 \f$, arranged row-wise with the
     * flattened \f$(i,j)\f$ index as the row and \f$ p \f$ as the column index
     */
    const Eigen::MatrixXd &getRMatrix(int order);

    /**
     * @brief Filter family/type code associated with this cache
     *
     * Set in the private constructor based on the template parameter @p T.
     * (E.g., Interpol or Legendre.)
     */
    int getType() const { return this->type; }

protected:
    int type;           ///< Filter family/type code matching template parameter @p T (e.g., Interpol or Legendre)
    std::string libPath; ///< Reserved base path for the filter library (path resolution is currently delegated to CrossCorrelation)

private:
    /**
     * @brief Private constructor enforces the singleton pattern
     *
     * Initializes @ref type based on T; see .cpp for validation.
     */
    CrossCorrelationCache();

    // Non-copyable / non-assignable — keeps the singleton unique.
    CrossCorrelationCache(CrossCorrelationCache<T> const &ccc) = delete;
    CrossCorrelationCache<T> &operator=(CrossCorrelationCache<T> const &ccc) = delete;
};

} // namespace mrcpp