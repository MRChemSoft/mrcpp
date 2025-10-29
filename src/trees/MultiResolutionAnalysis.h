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

#include "BoundingBox.h"
#include "core/MWFilter.h"
#include "core/ScalingBasis.h"

#include "MRCPP/mrcpp_declarations.h"

namespace mrcpp {

/**
 * @file MultiResolutionAnalysis.h
 * @brief Declaration of the MultiResolutionAnalysis class template.
 *
 * @details
 * A **MultiResolutionAnalysis (MRA)** bundles the information that must be
 * shared by compatible functions and operators:
 * - the computational domain (see @ref BoundingBox),
 * - the multiresolution scaling basis (see @ref ScalingBasis), and
 * - a maximum refinement depth.
 *
 * Two objects (functions/operators) can only be combined if their MRAs are
 * equal, i.e. identical domain, basis order/type, and depth.
 *
 * @par Example
 * @code{.cpp}
 * using MRA3 = mrcpp::MultiResolutionAnalysis<3>;
 *
 * // Domain: [-4, 4]^3 with automatically chosen root scale
 * mrcpp::BoundingBox<3> world({-4, 4});
 *
 * // Build a 3D MRA with Legendre order=7 and maxDepth=12
 * mrcpp::ScalingBasis basis(Legendre, /*order=*/7);
 * MRA3 mra(world, basis, /*depth=*/12);
 *
 * // Query information
 * int order     = mra.getOrder();
 * int maxScale  = mra.getMaxScale();
 * auto &box     = mra.getWorldBox();
 * auto &sbasis  = mra.getScalingBasis();
 * @endcode
 */

/**
 * @class MultiResolutionAnalysis
 * @tparam D Spatial dimension (1, 2, or 3).
 *
 * @brief Collects the computational domain and multiresolution basis.
 *
 * @details
 * The MRA fixes:
 * - the **world box** (domain tiling and scaling),
 * - the **scaling basis** (type and polynomial order), and
 * - the **maximum depth** of refinement relative to the world’s root scale.
 *
 * The combination of these parameters determines the finest admissible scale
 * via @ref getMaxScale.
 */
template <int D>
class MultiResolutionAnalysis final {
public:
    /**
     * @brief Construct from a symmetric domain and a basis order.
     *
     * @param[in] bb    Domain bounds as either [0,L] or [-L,L] (L>0).
     * @param[in] order Polynomial order of the scaling basis.
     * @param[in] depth Maximum refinement depth (relative to root scale).
     *
     * @details
     * The scaling basis type is chosen by MRCPP defaults for the given @p order.
     * The root scale is inferred from @p bb to keep the per-dimension scaling
     * factor in (1, 2).
     */
    MultiResolutionAnalysis(std::array<int, 2> bb, int order, int depth = MaxDepth);

    /**
     * @brief Construct from a preconfigured @ref BoundingBox and basis order.
     *
     * @param[in] bb    Computational domain (possibly periodic).
     * @param[in] order Polynomial order of the scaling basis.
     * @param[in] depth Maximum refinement depth.
     */
    MultiResolutionAnalysis(const BoundingBox<D> &bb, int order, int depth = MaxDepth);

    /**
     * @brief Construct from a @ref BoundingBox and a fully specified @ref ScalingBasis.
     *
     * @param[in] bb    Computational domain.
     * @param[in] sb    Scaling basis (type and order).
     * @param[in] depth Maximum refinement depth.
     */
    MultiResolutionAnalysis(const BoundingBox<D> &bb, const ScalingBasis &sb, int depth = MaxDepth);

    /** @brief Copy constructor. */
    MultiResolutionAnalysis(const MultiResolutionAnalysis<D> &mra);

    /** @brief Deleted assignment (MRAs are intended to be immutable after construction). */
    MultiResolutionAnalysis &operator=(const MultiResolutionAnalysis &mra) = delete;

    /** @brief Return polynomial order of the scaling basis. */
    int getOrder() const { return this->basis.getScalingOrder(); }

    /** @brief Maximum refinement depth relative to the world’s root scale. */
    int getMaxDepth() const { return this->maxDepth; }

    /**
     * @brief Absolute finest scale index.
     *
     * @details
     * This is the sum of the world root scale and @ref getMaxDepth, i.e.
     * the maximum scale the MRA allows trees to reach.
     */
    int getMaxScale() const { return this->world.getScale() + this->maxDepth; }

    /** @brief Low-level filter associated with the current basis. */
    const MWFilter &getFilter() const { return *this->filter; }

    /** @brief The scaling basis specification (type and order). */
    const ScalingBasis &getScalingBasis() const { return this->basis; }

    /** @brief The computational domain (world box). */
    const BoundingBox<D> &getWorldBox() const { return this->world; }

    /**
     * @brief Convenience: compute a minimal length scale from a tolerance.
     *
     * @param[in] epsilon Target tolerance.
     * @return A distance proportional to \f$\sqrt{\epsilon\,2^{-\mathrm{maxScale}}}\f$.
     */
    double calcMinDistance(double epsilon) const { return std::sqrt(epsilon * std::pow(2.0, -getMaxScale())); }

    /**
     * @brief Convenience: compute a maximal relevant distance.
     *
     * @details The exact definition is basis-dependent and implemented in
     * the corresponding source file.
     */
    double calcMaxDistance() const;

    /** @brief Root (coarsest) scale index of the world. */
    int getRootScale() const { return this->world.getScale(); }

    /**
     * @brief Equality: same world, same basis (type & order), same depth.
     *
     * @note Two MRAs must compare equal to allow mixing functions/operators.
     */
    bool operator==(const MultiResolutionAnalysis<D> &mra) const;

    /** @brief Inequality. */
    bool operator!=(const MultiResolutionAnalysis<D> &mra) const;

    /** @brief Human-readable diagnostics to stdout. */
    void print() const;

protected:
    const int maxDepth;           ///< Maximum refinement depth permitted by this MRA.
    const ScalingBasis basis;     ///< Scaling basis (type and polynomial order).
    const BoundingBox<D> world;   ///< Computational domain description.
    MWFilter *filter;             ///< Low-level filter derived from @ref basis.

    /** @brief Internal helper to instantiate @ref filter based on @ref basis. */
    void setupFilter();
};

} // namespace mrcpp