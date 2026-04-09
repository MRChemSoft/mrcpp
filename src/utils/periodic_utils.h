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
/**
 * @file
 * @brief Periodic boundary utilities for node indices and real-space coordinates.
 *
 * This header declares helpers for enforcing periodic boundary conditions (PBC)
 * on both discrete tree indices (@ref mrcpp::NodeIndex) and continuous
 * coordinates (@ref mrcpp::Coord). The utilities are templated on dimension
 * @p D and support selectively periodic directions via a boolean mask.
 *
 * Typical use cases:
 * - Normalizing a node index to the canonical unit cell before lookup.
 * - Wrapping real-space coordinates into the primary cell when sampling or
 *   exporting data.
 * - Applying periodicity per axis (e.g., 2D slab periodic in x/y but not z).
 */

#include "MRCPP/mrcpp_declarations.h"

namespace mrcpp {
/**
 * @namespace mrcpp::periodic
 * @brief Helpers for periodic index/coordinate manipulation.
 *
 * The functions here assume MRCPP’s convention where the canonical cell is the
 * unit hypercube, and indices/coordinates are normalized accordingly:
 * - Discrete indices: @ref NodeIndex<D> logically cover tiles of the unit cell
 *   at a given resolution (scale). These helpers re-map out-of-range indices
 *   back into the unit cell modulo the periodic axes.
 * - Continuous coordinates: @ref Coord<D> (double-valued) are wrapped by
 *   subtracting/adding integer lattice vectors along periodic axes so that
 *   the result lies in the half-open interval [0, 1) per periodic dimension.
 */
namespace periodic {

/**
 * @brief Check whether a node index lies inside the unit cell.
 *
 * @tparam D Spatial dimension.
 * @param idx Node index to test (scale and per-dimension integer indices).
 * @return @c true if @p idx is within the canonical unit cell bounds in all
 *         dimensions; @c false if any component is outside.
 *
 * @details “Inside” means the discrete index components fall in the valid range
 * for the node’s scale with no modular wrap required. This does not modify
 * @p idx and performs a pure check. Use @ref index_manipulation to fold an
 * index back into the unit cell when periodicity is intended.
 */
template <int D>
bool in_unit_cell(NodeIndex<D> idx);

/**
 * @brief Fold a node index into the unit cell under per-axis periodicity.
 *
 * @tparam D Spatial dimension.
 * @param[in,out] idx Node index to normalize; on return, the per-axis integer
 *                    index components are mapped into the unit-cell range for
 *                    the node’s scale when the corresponding axis is periodic.
 * @param periodic Boolean mask of length @p D; @c true marks an axis as
 *                 periodic, @c false leaves that axis unchanged (no wrapping).
 *
 * @details
 * For each periodic axis, the index component is reduced modulo the extent at
 * the node’s scale so that the resulting index is in-range. For non-periodic
 * axes, the index is left as-is (and may remain out-of-bounds if provided so).
 *
 * @note This function is idempotent for already in-range indices on periodic
 * axes and is a no-op for non-periodic axes.
 */
template <int D>
void index_manipulation(NodeIndex<D> &idx, const std::array<bool, D> &periodic);

/**
 * @brief Wrap a coordinate into the unit cell under per-axis periodicity.
 *
 * @tparam D Spatial dimension.
 * @param[in,out] r Coordinate to normalize; each periodic component is wrapped
 *                  into the half-open interval [0, 1).
 * @param periodic Boolean mask of length @p D; @c true marks an axis as
 *                 periodic, @c false leaves that axis unchanged (no wrapping).
 *
 * @details
 * For each periodic axis, the component @f$r_d@f$ is replaced by
 * @f$r_d - \lfloor r_d \rfloor@f$, producing a value in [0, 1). Non-periodic
 * axes are not modified. This is equivalent to applying @c std::floor-based
 * fractional reduction to each periodic component.
 *
 * @warning If your simulation cell is scaled or shifted relative to the unit
 * cube, convert to reduced coordinates before calling this function, or adjust
 * the values accordingly (e.g., divide by box length) and convert back.
 */
template <int D>
void coord_manipulation(Coord<D> &r, const std::array<bool, D> &periodic);

} // namespace periodic
} // namespace mrcpp