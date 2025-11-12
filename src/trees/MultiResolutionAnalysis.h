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
 * @class MultiResolutionAnalysis
 * @tparam D Spatial dimension (1, 2, or 3)
 *
 * @brief Class for MultiResolutionAnalysis templates
 *
 * @details
 * The MultiResolutionAnalysis (MRA) objects bundles information that must be shared for 
 * compatible functions and operators:
 * - Computational domain (see @ref BoundingBox)
 * - MultiResolution scaling basis, as a polynomial order (see @ref ScalingBasis)
 * - Maximum refinement depth, relative to the world’s root scale (= @ref maxDepth by default)
 * 
 * Class also contains useful functions to compare MRA objects, 
 * find max and min box sizes and print a human readable diagnostic for the MRA.
 */
template <int D>
class MultiResolutionAnalysis final {
public:
    /**
     * @brief Construct from a symmetric domain and a basis order
     *
     * @param[in] bb    2-element integer array defining domain bounds
     * @param[in] order Polynomial order of the multiwavelet basis
     * @param[in] depth Maximum refinement depth (relative to root scale). Default is \ref MaxDepth
     *
     * @details
     * Constructor of the MultiResolutionAnalysis class from scratch.
     * The scaling basis type is chosen by MRCPP defaults for the given @p order.
     * The root scale is inferred from @p bb to keep the per-dimension scaling factor in (1, 2).
     */
    MultiResolutionAnalysis(std::array<int, 2> bb, int order, int depth = MaxDepth);

    /**
     * @brief Constructs MultiResolutionAnalysis object from a pre-existing @ref BoundingBox object
     *
     * @param[in] bb    BoundingBox object representing the computational domain
     * @param[in] order Polynomial order of the multiwavelet basis
     * @param[in] depth Maximum refinement depth (relative to root scale). Default is \ref MaxDepth
     *
     * @details
     * Creates a MRA object from pre-existing BoundingBox, @p bb, object with a polynomial, @ref p, order to set the basis 
     * and the maximum amount of allowed refinement in a node, @p depth.
     */
    MultiResolutionAnalysis(const BoundingBox<D> &bb, int order, int depth = MaxDepth);

    /**
     * @brief Construct from a @ref BoundingBox and a fully specified @ref ScalingBasis.
     *
     * @param[in] bb    BoundingBox object representing the computational domain
     * @param[in] sb    Polynomial basis (MW) as a ScalingBasis object
     * @param[in] depth Maximum refinement depth (relative to root scale). Default is \ref MaxDepth
     
     * @details
     * Creates a MRA object from pre-existing BoundingBox, @p bb, and ScalingBasis, @p sb, objects 
     * and the maximum amount of allowed refinement in a node, @p depth.
     */
    MultiResolutionAnalysis(const BoundingBox<D> &bb, const ScalingBasis &sb, int depth = MaxDepth);

    /**
     * @brief Copy constructor for a MultiResolutionAnalysis object composed of computational domain (world) and a polynomial basis (Multiwavelets)
     * @param[in] mra Pre-existing MRA object
     * @details Copy a MultiResolutionAnalysis object without modifying the original
     */
    MultiResolutionAnalysis(const MultiResolutionAnalysis<D> &mra);

    /** @brief Deleted assignment (MRAs are intended to be immutable after construction). */
    MultiResolutionAnalysis &operator=(const MultiResolutionAnalysis &mra) = delete;

    /*
    * Getters
    */

    int getOrder() const { return this->basis.getScalingOrder(); }              ///< @return Polynomial order of the scaling basis
    int getMaxDepth() const { return this->maxDepth; }                          ///< @return Maximum refinement depth relative to the world’s root scale
    int getMaxScale() const { return this->world.getScale() + this->maxDepth; } ///< @return Sum of world root scale and maximum refinement depth, @ref getMaxDepth
    int getRootScale() const { return this->world.getScale(); }                 ///< @return World root scale

    const MWFilter &getFilter() const { return *this->filter; }                 ///< @return Low-level filter associated with the current basis
    const ScalingBasis &getScalingBasis() const { return this->basis; }         ///< @return Scaling basis type and order
    const BoundingBox<D> &getWorldBox() const { return this->world; }           ///< @return Computational domain (world box)

    /**
     * @brief Convenience: compute a minimal length scale from a tolerance
     * @param[in] epsilon Target tolerance
     * @return A distance proportional to \f$\sqrt{\epsilon\,2^{-\mathrm{maxScale}}}\f$.
     */
    double calcMinDistance(double epsilon) const { return std::sqrt(epsilon * std::pow(2.0, -getMaxScale())); }

    /**
     * @brief Convenience: compute a maximal relevant distance
     * @return Maximum distance of computational (world) domain
     * @note The exact definition is basis-dependent
     */
    double calcMaxDistance() const;

    /**
     * @brief Equality operator for the MultiResolutionAnalysis class (basis, domain, depth)
     *
     * @param[in] mra: MRA object, taken by constant reference
     * @returns True if both MRAs have the same polynomial basis, computational domain and maximum depth
     *
     * @note Two MRAs must be equal to allow mixing functions/operators
     */
    bool operator==(const MultiResolutionAnalysis<D> &mra) const;

    /**
     * @brief Inequality operator for the MultiResolutionAnalysis class (basis, domain, depth)
     * @param[in] mra: MRA object, taken by constant reference
     * @returns True if MRAs have different polynomial basis, computational domain or maximum depth
     */
    bool operator!=(const MultiResolutionAnalysis<D> &mra) const;

    void print() const;           ///< @brief Displays human-readable diagnostics of MRA to outputfile

protected:
    const int maxDepth;           ///< Maximum refinement depth permitted by this MRA
    const ScalingBasis basis;     ///< Scaling basis (type and polynomial order)
    const BoundingBox<D> world;   ///< Computational domain description
    MWFilter *filter;             ///< Low-level filter derived from @ref basis

    /**
     * @brief Internal helper to instantiate @ref filter based on @ref basis
     */
    void setupFilter();
};

} // namespace mrcpp