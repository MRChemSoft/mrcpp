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

/**
 * @file GaussExp.cpp (lightweight connector)
 *
 * @brief Ties together Gaussian primitives and exponential utilities, and
 *        exposes (via forward declaration) the templated overlap routine
 *        without pulling in heavier headers that could cause cycles.
 *
 * What this TU does
 * -----------------
 * - Includes:
 *     - "GaussExp.h": utilities for Gaussian/exponential expressions used
 *       elsewhere in MRCPP (e.g., Boys integrals, screened interactions).
 *     - "Gaussian.h": the definition of `GaussFunc<D>`, i.e., a Cartesian
 *       Gaussian primitive storing powers, center, exponent(s), and a coefficient.
 * - Declares (but does not define) the templated function
 *   `function_utils::calc_overlap<D>(const GaussFunc<D>&, const GaussFunc<D>&)`.
 *   The definition lives in the function-utils implementation unit
 *   (see `function_utils.cpp`), which provides the Obara–Saika-based 1D core.
 *
 * Why only a forward declaration here?
 * ------------------------------------
 * - To avoid including a potentially heavy implementation header (and risking
 *   circular dependencies), we forward-declare the template in the *same*
 *   namespace `mrcpp::function_utils`. This enables use sites that only need
 *   the signature to compile quickly, while the actual template definition
 *   will be instantiated by the linker when the corresponding .cpp is linked.
 *
 * Notes on templates and linkage
 * ------------------------------
 * - Because this is only a declaration, any translation unit that actually
 *   *uses* `calc_overlap<D>` must see the template **definition** (e.g., by
 *   including the proper header or by relying on explicit instantiations
 *   provided in the implementation TU). MRCPP provides common explicit
 *   instantiations (e.g., D = 1, 2, 3) in `function_utils.cpp`.
 *
 * Example usage
 * -------------
 * @code
 * #include "Gaussian.h"
 * // (this file is included transitively somewhere)
 * using mrcpp::GaussFunc;
 * using mrcpp::function_utils::calc_overlap;
 *
 * GaussFunc<3> gA(...), gB(...);
 * double S = calc_overlap<3>(gA, gB); // calls Obara–Saika-backed routine
 * @endcode
 */

#include "GaussExp.h"
#include "Gaussian.h"

namespace mrcpp {

// Forward declaration only: definition is provided in function_utils.cpp.
// Keeping this here avoids heavy includes and potential include cycles.
namespace function_utils {
template <int D>
double calc_overlap(const GaussFunc<D> &a, const GaussFunc<D> &b);
} // namespace function_utils

} // namespace mrcpp