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

#include "CornerBandWidth.h"

namespace mrcpp {


/** @brief Checks if the distance to diagonal is lesser than the operator band width.
 *
 * @param[in] oTransl: distance to diagonal
 * @param[in] o_depth: scaling order
 * @param[in] idx: index corresponding to one of the matrices \f$ A, B, C \f$ or \f$ T \f$.
 * 
 * @returns True if \b oTransl is outside of the corner band (close to diagonal) and False otherwise. 
 * 
 */ 
bool CornerBandWidth::isOutsideBand(int oTransl, int o_depth, int idx) const
{
    return abs(oTransl) < getWidth(o_depth, idx);
}


} // namespace mrcpp
