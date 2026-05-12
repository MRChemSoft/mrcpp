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

#include <Eigen/Core>

#include "GaussQuadrature.h"
#include "ObjectCache.h"

namespace mrcpp {

/**
 * @def getQuadratureCache(X)
 * @brief Convenience macro to bind a local reference @p X to the global
 *        (singleton) QuadratureCache instance
 *
 * Expands to:
 *   QuadratureCache &X = QuadratureCache::getInstance()
 *
 * Example:
 *   getQuadratureCache(qc);
 *   const auto& w = qc.getWeights(order);
 */
#define getQuadratureCache(X) QuadratureCache &X = QuadratureCache::getInstance()

/**
 * @class QuadratureCache
 * @brief Process-wide singleton cache for Gauss–Legendre quadrature rules, keyed by order
 *
 * @details
 * Constructing GaussQuadrature objects (root finding via Newton iterations) can be costly if repeated.
 * This cache stores one instance per polynomial order and hands out references on demand. Domain
 * settings — integration bounds \f$ [A,B] \f$ and the number of sub-intervals — are kept in the cache
 * and applied when a new entry is loaded via load(). Inherits from ObjectCache<GaussQuadrature> for
 * integer-indexed storage and memory accounting. The singleton is accessed via getInstance(); first-time
 * loads are protected by OpenMP locks in the .cpp while subsequent reads are lock-free.
 *
 * @see GaussQuadrature for the cached object type
 */
class QuadratureCache final : public ObjectCache<GaussQuadrature> {
public:
    /**
     * @brief Access the singleton instance
     */
    static QuadratureCache &getInstance() {
        static QuadratureCache theQuadratureCache;
        return theQuadratureCache;
    }

    /**
     * @brief Ensure the quadrature of a given @p order is loaded
     *
     * Implemented in the .cpp: constructs/initializes a GaussQuadrature that
     * reflects the current @ref A, @ref B, and @ref intervals settings and
     * inserts it into the underlying ObjectCache if absent
     */
    void load(int order);

    /**
     * @brief Retrieve the cached quadrature for @p order (lazy-loads if needed)
     * @return Reference to the GaussQuadrature object owned by the cache
     */
    GaussQuadrature &get(int order);

    /**
     * @name Convenience accessors (fetch vectors directly)
     * @{
     * @brief Get the vector of abscissas (roots) for a given order
     */
    const Eigen::VectorXd &getRoots(int i) { return get(i).getRoots(); }

    /**
     * @brief Get the vector of weights for a given order
     */
    const Eigen::VectorXd &getWeights(int i) { return get(i).getWeights(); }
    /** @} */

    /**
     * @brief Set the number of sub-intervals for subsequently loaded quadrature rules
     * @param i Number of equal sub-intervals (\f$ \geq 1 \f$)
     *
     * @details Stores the value in #intervals so that new calls to load() will construct
     * GaussQuadrature objects with the updated partitioning. Existing cached entries are not affected.
     */
    void setIntervals(int i);

    /**
     * @brief Set the integration bounds used by subsequently loaded quadrature rules
     * @param a Lower bound \f$ A \f$
     * @param b Upper bound \f$ B \f$ (must satisfy \f$ a < b \f$)
     *
     * @details Stores the new bounds in #A and #B. Existing cached entries are unaffected until
     * explicitly unloaded and reloaded.
     */
    void setBounds(double a, double b);

    /** @return Current number of sub-intervals recorded in the cache */
    int getIntervals() const { return this->intervals; }

    /** @return Current upper integration bound B */
    double getUpperBound() const { return this->B; }

    /** @return Current lower integration bound A */
    double getLowerBound() const { return this->A; }

private:
    double A;      ///< Lower integration bound (default set in private constructor)
    double B;      ///< Upper integration bound (default set in private constructor)
    int intervals; ///< Number of equal sub-intervals used when constructing new GaussQuadrature entries

    /**
     * @brief Private constructor; initializes default bounds and interval count
     *
     * @details Enforces the singleton pattern — use getInstance() to obtain the cache
     */
    QuadratureCache();

    /// Private destructor; cache cleans up its owned objects via ObjectCache
    ~QuadratureCache();

    // Non-copyable / non-assignable to maintain singleton semantics.
    QuadratureCache(QuadratureCache const &qc) = delete;
    QuadratureCache &operator=(QuadratureCache const &qc) = delete;
};

} // namespace mrcpp