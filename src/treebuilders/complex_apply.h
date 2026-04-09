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

#include "trees/FunctionTreeVector.h"

namespace mrcpp {

/**
 * @file
 * @brief Complex wrapper utilities for multiwavelet trees and operators.
 *
 * @details
 * This header declares a lightweight wrapper, @ref ComplexObject, that groups
 * pointers to the real and imaginary parts of an object (e.g., a function
 * tree or a convolution operator). It also declares an `apply` routine that
 * applies a (possibly complex) convolution operator to a (possibly complex)
 * function, writing the result to a (possibly complex) output.
 *
 * The pattern keeps real and imaginary parts as separate objects for memory
 * locality and to reuse existing real-valued kernels, while allowing users to
 * orchestrate complex arithmetic at a higher level.
 */

/**
 * @brief Aggregates pointers to the real and imaginary parts of an object.
 *
 * @tparam MWClass Underlying class type of the wrapped objects
 *         (e.g., `FunctionTree<D>` or `ConvolutionOperator<D>`).
 *
 * @details
 * The struct is a non-owning pair of pointers. It does **not** manage
 * lifetime—callers must ensure both referenced objects outlive the wrapper.
 *
 * @note
 * The members are intentionally public for ergonomic access in kernels.
 */
template <typename MWClass>
struct ComplexObject {
    /** @brief Pointer to the real component (non-owning). */
    MWClass* real;
    /** @brief Pointer to the imaginary component (non-owning). */
    MWClass* imaginary;

    /**
     * @brief Construct from lvalue references to the real and imaginary parts.
     * @param realPart      Reference to the real component.
     * @param imaginaryPart Reference to the imaginary component.
     */
    ComplexObject(MWClass& realPart, MWClass& imaginaryPart)
        : real(&realPart)
        , imaginary(&imaginaryPart) {}
};

// clang-format off
//template <int D> class FunctionTree;
//template <int D> class ConvolutionOperator;

/**
 * @brief Apply a (complex) convolution operator to a (complex) function.
 *
 * @tparam D Spatial dimensionality of the multiwavelet representation.
 *
 * @param prec     Target accuracy. If `absPrec == false`, this is interpreted
 *                 as a **relative** tolerance; otherwise as an **absolute** tolerance.
 * @param out      Destination complex function trees (real/imag). On return,
 *                 contains \f$ \text{oper} \{\text{inp}\} \f$ within the requested
 *                 accuracy.
 * @param oper     Complex convolution operator (real/imag components).
 * @param inp      Input complex function trees to be transformed.
 * @param maxIter  Optional cap on internal refinement/iteration steps.
 *                 Use `-1` (default) for the implementation’s automatic choice.
 * @param absPrec  When `true`, treat `prec` as absolute; when `false`, as relative.
 *
 * @pre
 * - `out.real`, `out.imaginary`, `inp.real`, `inp.imaginary`,
 *   `oper.real`, and `oper.imaginary` are non-null and represent
 *   consistent discretizations (same MRA/order/domain).
 *
 * @post
 * - `out` holds the complex result. Implementations typically compute:
 *   \f[
 *     \Re(\text{out}) = \Re(\text{oper})\Re(\text{inp})
 *                      - \Im(\text{oper})\Im(\text{inp}),
 *     \qquad
 *     \Im(\text{out}) = \Re(\text{oper})\Im(\text{inp})
 *                      + \Im(\text{oper})\Re(\text{inp}),
 *   \f]
 *   with adaptive refinement to honor `prec`.
 *
 * @note
 * The exact refinement strategy and stopping criteria are backend-dependent.
 * For reproducibility across runs/nodes, set the relevant MPI/OpenMP controls
 * prior to calling.
 */
template <int D>
void apply
(
    double prec, ComplexObject< FunctionTree<D> >& out,
    ComplexObject< ConvolutionOperator<D> >& oper, ComplexObject< FunctionTree<D> >& inp,
    int maxIter = -1, bool absPrec = false
);
// clang-format on

} // namespace mrcpp