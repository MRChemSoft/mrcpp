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
 *        (singleton) QuadratureCache instance.
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
 * @brief Process-wide cache for Gaussian quadrature rules (roots & weights).
 *
 * High-level
 * ----------
 * Gaussian quadrature (Gauss-Legendre in MRCPP) is parameterized by:
 *   • order (number of nodes/weights),
 *   • integration domain [A, B],
 *   • optional replication over multiple equal sub-intervals ("intervals").
 *
 * Constructing GaussQuadrature objects repeatedly can be costly; this cache
 * stores one instance per order (and current domain/interval settings) and
 * hands out references on demand.
 *
 * Design
 * ------
 * - Singleton per process (Meyers' singleton via getInstance()).
 * - Inherits from ObjectCache<GaussQuadrature> which provides basic
 *   load/unload/get plumbing indexed by an integer id (here: order).
 * - Domain control:
 *     setBounds(a,b)   → set global integration bounds [A,B]
 *     setIntervals(i)  → split [A,B] into @p i equal sub-intervals (if used)
 *   These settings influence how GaussQuadrature is created in load(order).
 *
 * Thread-safety
 * -------------
 * The base ObjectCache does not synchronize by itself; specialized caches
 * typically guard first-time loads with OpenMP locks in the .cpp. Users should
 * assume the cache is safe to read concurrently after an entry is present.
 *
 * Typical usage
 * -------------
 *   auto& qc = QuadratureCache::getInstance();
 *   qc.setBounds(-1.0, 1.0);
 *   const Eigen::VectorXd& x = qc.getRoots(quad_order);
 *   const Eigen::VectorXd& w = qc.getWeights(quad_order);
 */
class QuadratureCache final : public ObjectCache<GaussQuadrature> {
public:
    /**
     * @brief Access the singleton instance.
     */
    static QuadratureCache &getInstance() {
        static QuadratureCache theQuadratureCache;
        return theQuadratureCache;
    }

    /**
     * @brief Ensure the quadrature of a given @p order is loaded.
     *
     * Implemented in the .cpp: constructs/initializes a GaussQuadrature that
     * reflects the current @ref A, @ref B, and @ref intervals settings and
     * inserts it into the underlying ObjectCache if absent.
     */
    void load(int order);

    /**
     * @brief Retrieve the cached quadrature for @p order (lazy-loads if needed).
     * @return Reference to the GaussQuadrature object owned by the cache.
     */
    GaussQuadrature &get(int order);

    /**
     * @name Convenience accessors (fetch vectors directly)
     * @{
     * @brief Get the vector of abscissas (roots) for a given order.
     */
    const Eigen::VectorXd &getRoots(int i) { return get(i).getRoots(); }

    /**
     * @brief Get the vector of weights for a given order.
     */
    const Eigen::VectorXd &getWeights(int i) { return get(i).getWeights(); }
    /** @} */

    /**
     * @brief Set the number of equal sub-intervals for composite quadrature.
     *
     * Interpretation:
     *  - If intervals > 1, the base interval [A,B] can be partitioned into
     *    `intervals` equal pieces and the quadrature replicated/shifted.
     *  - Exact semantics depend on GaussQuadrature; this cache records the
     *    value so that new loads honor it.
     */
    void setIntervals(int i);

    /**
     * @brief Set the integration bounds used by subsequently loaded rules.
     * @param a Lower bound A
     * @param b Upper bound B
     *
     * Newly created GaussQuadrature objects will target [A,B]. Existing
     * cached entries are unaffected until explicitly unloaded/reloaded.
     */
    void setBounds(double a, double b);

    /** @return Current number of sub-intervals recorded in the cache. */
    int getIntervals() const { return this->intervals; }

    /** @return Current upper integration bound B. */
    double getUpperBound() const { return this->B; }

    /** @return Current lower integration bound A. */
    double getLowerBound() const { return this->A; }

private:
    /**
     * @brief Lower and upper bounds of the integration domain.
     *
     * Defaults are set in the private constructor (see .cpp). Changing these
     * affects only future loads; existing cached rules remain as created.
     */
    double A;
    double B;

    /**
     * @brief Number of equal sub-intervals used to tile [A,B].
     *
     * When >1, the cache can generate composite quadrature by replicating the
     * base rule on each sub-interval (implementation in .cpp / GaussQuadrature).
     */
    int intervals;

    /**
     * @brief Private constructor initializes default bounds/intervals.
     *
     * Enforces the singleton pattern; use getInstance() to access the cache.
     */
    QuadratureCache();

    /// Private destructor; cache cleans up its owned objects via ObjectCache.
    ~QuadratureCache();

    // Non-copyable / non-assignable to maintain singleton semantics.
    QuadratureCache(QuadratureCache const &qc) = delete;
    QuadratureCache &operator=(QuadratureCache const &qc) = delete;
};

} // namespace mrcpp