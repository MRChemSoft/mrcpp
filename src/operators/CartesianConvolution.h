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

#include "ConvolutionOperator.h"

namespace mrcpp {

// Forward declaration to avoid pulling the full header into users of this file.
template <int D> class GaussExp;

/**
 * @class CartesianConvolution
 * @brief 3D separable convolution operator assembled from a 1D Gaussian expansion.
 *
 * This operator represents a Cartesian, rank-R, separable convolution in 3D,
 * where the separation rank R equals the number of terms in a provided 1D
 * Gaussian expansion, `GaussExp<1>`.
 *
 * ### How it is constructed (see .cpp)
 * The implementation builds three *blocks* of 1D operator trees from the
 * same Gaussian expansion, corresponding to monomial prefactors of degree
 * 0, 1 and 2 (i.e. powers `{0}`, `{1}`, `{2}`), and stores them
 * contiguously. These blocks can then be assigned independently to the
 * x/y/z axes, enabling vector/tensor kernels that differ only by the
 * Cartesian polynomial factor.
 *
 * After construction, the total number of internally stored operator trees is
 * `3 * sep_rank` (three monomial blocks, each of size `sep_rank`).
 *
 * ### Choosing the Cartesian components
 * Use setCartesianComponents(x, y, z) to select which monomial block
 * (0 → degree 0, 1 → degree 1, 2 → degree 2) is used along each axis. This
 * *rewires* the already built 1D factors—no rebuilding occurs.
 *
 * ### Precision
 * The constructor accepts a single build precision `prec`. The .cpp implementation
 * employs a slightly stricter precision for fitting the 1D kernel terms so that
 * the final composed 3D operator meets the requested tolerance.
 *
 * ### Ownership / lifetime
 * The class does not take ownership of the input `GaussExp<1>`; it only reads it
 * during construction. Internally created operator trees are owned by this object.
 *
 * ### Copy semantics
 * Copying is disabled (non-copyable) because the underlying operator trees are
 * heavy and managed resources. Move is not provided.
 *
 * ### Example
 * @code
 * MultiResolutionAnalysis<3> mra(...);
 * GaussExp<1> kernel = ...;           // ∑_{r=1}^R α_r e^{-β_r (x - x_r)^2}
 * double prec = 1e-8;
 *
 * CartesianConvolution conv(mra, kernel, prec);
 * // Use degree-1 along x, degree-0 along y, degree-2 along z:
 * conv.setCartesianComponents(/* x = * / 1, /* y = * / 0, /* z = * / 2);
 *
 * // conv can now be applied as a separable 3D convolution operator.
 * @endcode
 */
class CartesianConvolution : public ConvolutionOperator<3> {
public:
    /**
     * @brief Construct a 3D separable convolution operator from a 1D Gaussian expansion.
     *
     * @param mra     Multiresolution analysis defining the 3D basis/domain.
     * @param kernel  1D Gaussian expansion; its length sets the separation rank R.
     *                The implementation temporarily adjusts the monomial power of each
     *                Gaussian term to build three internal blocks (degrees 0, 1, 2),
     *                but does not take ownership of @p kernel.
     * @param prec    Target build precision for the assembled operator.
     *
     * @details
     * Internally, three batches of operator trees are built (for polynomial degrees
     * 0/1/2), each of size R, and stored contiguously. The final separable operator
     * exposes rank R; per-axis assignment of the blocks is deferred to
     * setCartesianComponents().
     */
    CartesianConvolution(const MultiResolutionAnalysis<3> &mra, GaussExp<1> &kernel, double prec);

    CartesianConvolution(const CartesianConvolution &oper) = delete;
    CartesianConvolution &operator=(const CartesianConvolution &oper) = delete;
    virtual ~CartesianConvolution() = default;

    /**
     * @brief Select which monomial block is used on each Cartesian axis.
     *
     * @param x  Block index for x-axis (0 → degree 0, 1 → degree 1, 2 → degree 2).
     * @param y  Block index for y-axis (same convention).
     * @param z  Block index for z-axis (same convention).
     *
     * @details
     * - This operation is O(R) for each axis and **does not rebuild** the operator;
     *   it remaps the already constructed 1D operator trees into the separable slots.
     * - Valid indices are {0,1,2}. Using the same block on multiple axes is allowed.
     */
    void setCartesianComponents(int x, int y, int z);

protected:
    /**
     * @brief Separation rank R of the operator (number of terms in the input 1D kernel).
     *
     * @details
     * The internal storage contains 3·R operator trees (three monomial blocks),
     * but the exposed separable rank for downstream composition is R.
     */
    int sep_rank;
};

} // namespace mrcpp