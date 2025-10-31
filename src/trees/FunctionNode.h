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

#include "FunctionTree.h"
#include "MWNode.h"

namespace mrcpp {

/**
 * @file FunctionNode.h
 * @brief Leaf/branch node type that stores function coefficients on a
 *        multiresolution tree.
 *
 * @details
 * A FunctionNode is a concrete MWNode specialized for function representations.
 * It holds scaling and wavelet coefficients, provides allocation and refinement
 * helpers, and exposes utilities for evaluation, coefficient access and
 * basic per-node operations such as integration and local dot products.
 *
 * Template parameters:
 * - D: spatial dimension (1, 2 or 3)
 * - T: scalar type (double or ComplexDouble)
 */

/**
 * @class FunctionNode
 * @tparam D Spatial dimension.
 * @tparam T Scalar type.
 * @brief Node of a FunctionTree that stores coefficients and implements
 *        function-specific operations.
 *
 * @note Construction is managed by FunctionTree and NodeAllocator. Users do not
 *       construct FunctionNode directly.
 */
template <int D, typename T> class FunctionNode final : public MWNode<D, T> {
public:
    /** @name Typed accessors */
    ///@{

    /** @brief Return the owning FunctionTree (non-const). */
    FunctionTree<D, T> &getFuncTree() { return static_cast<FunctionTree<D, T> &>(*this->tree); }

    /** @brief Return the parent node cast to FunctionNode (non-const). */
    FunctionNode<D, T> &getFuncParent() { return static_cast<FunctionNode<D, T> &>(*this->parent); }

    /** @brief Return the i-th child cast to FunctionNode (non-const). */
    FunctionNode<D, T> &getFuncChild(int i) { return static_cast<FunctionNode<D, T> &>(*this->children[i]); }

    /** @brief Return the owning FunctionTree (const). */
    const FunctionTree<D, T> &getFuncTree() const { return static_cast<const FunctionTree<D, T> &>(*this->tree); }

    /** @brief Return the parent node cast to FunctionNode (const). */
    const FunctionNode<D, T> &getFuncParent() const { return static_cast<const FunctionNode<D, T> &>(*this->parent); }

    /** @brief Return the i-th child cast to FunctionNode (const). */
    const FunctionNode<D, T> &getFuncChild(int i) const { return static_cast<const FunctionNode<D, T> &>(*this->children[i]); }

    ///@}

    /** @name Tree-structure overrides */
    ///@{

    /**
     * @brief Create children of this node.
     * @param coefs If true, initialize children by transferring coefficients
     *              from this node as appropriate for the basis.
     *
     * @details Allocates child nodes through the node allocator and updates
     * the internal topology. When coefs is true, scaling/wavelet blocks are
     * propagated so that the represented function is unchanged by the split.
     */
    void createChildren(bool coefs) override;

    /**
     * @brief Generate (allocate) children if absent.
     * @details Convenience wrapper that creates children without coefficient
     * transfer. Intended for topology building when coefficients are filled
     * later by a calculator.
     */
    void genChildren() override;

    /**
     * @brief Ensure a parent exists and is allocated.
     * @details Creates the parent node if missing and links this node into the
     * parent children array.
     */
    void genParent() override;

    /**
     * @brief Delete children of this node.
     * @details Deallocates child nodes and updates internal state. Coefficients
     * in this node remain untouched.
     */
    void deleteChildren() override;

    ///@}

    /**
     * @brief Integrate the node contribution over its spatial support.
     * @return The integral of the locally represented function on this node.
     *
     * @details Uses the current scaling basis to compute the exact contribution
     * from scaling and wavelet parts confined to this node. For orthonormal
     * wavelets the integral often reduces to the scaling block.
     */
    T integrate() const;

    /** @name Value and coefficient access */
    ///@{

    /**
     * @brief Set nodal values from a vector.
     * @param vec Column vector of size getNCoefs().
     *
     * @details The vector is interpreted in the node's value layout used by
     * the scaling basis. Typical use is for interpolating bases, where values
     * correspond to quadrature or interpolation points. Internally, node
     * coefficients are updated accordingly.
     */
    void setValues(const Eigen::Matrix<T, Eigen::Dynamic, 1> &vec);

    /**
     * @brief Extract nodal values into a vector.
     * @param[out] vec Column vector resized to getNCoefs().
     *
     * @details The returned values correspond to the basis-specific value
     * layout for this node (e.g. interpolation/expanded points).
     */
    void getValues(Eigen::Matrix<T, Eigen::Dynamic, 1> &vec);

    /**
     * @brief Write absolute values of coefficients into a raw buffer.
     * @param[out] absCoefs Pointer to memory of length getNCoefs().
     *
     * @details Fills absCoefs[i] = abs(coef[i]). For complex T this is the
     * magnitude; for real T this is std::abs. The ordering matches the node's
     * internal coefficient layout.
     */
    void getAbsCoefs(T *absCoefs);

    ///@}

    friend class FunctionTree<D, T>;
    friend class NodeAllocator<D, T>;

protected:
    /** @name Constructors and assignment (managed by the tree) */
    ///@{

    FunctionNode()
            : MWNode<D, T>() {}

    explicit FunctionNode(MWTree<D, T> *tree, int rIdx)
            : MWNode<D, T>(tree, rIdx) {}

    FunctionNode(MWNode<D, T> *parent, int cIdx)
            : MWNode<D, T>(parent, cIdx) {}

    FunctionNode(MWTree<D, T> *tree, const NodeIndex<D> &idx)
            : MWNode<D, T>(tree, idx) {}

    FunctionNode(const FunctionNode<D, T> &node) = delete;
    FunctionNode<D, T> &operator=(const FunctionNode<D, T> &node) = delete;
    ~FunctionNode() = default;

    ///@}

    /** @brief Evaluate the reconstructed function at r (using this node only). */
    T evalf(Coord<D> r);

    /** @brief Evaluate the scaling part at r. */
    T evalScaling(const Coord<D> &r) const;

    /** @brief Deallocate node-owned memory and reset local state. */
    void dealloc() override;

    /** @brief Recompress local coefficients after updates. */
    void reCompress() override;

    /** @brief Integration helper for Legendre scaling basis. */
    T integrateLegendre() const;

    /** @brief Integration helper for interpolating scaling basis. */
    T integrateInterpolating() const;

    /** @brief Integration helper when values representation is active. */
    T integrateValues() const;
};

/** @name Per-node local dot-product helpers (double) */
///@{

/**
 * @brief Dot product of scaling parts on matching nodes (double).
 * @return Sum over matching scaling blocks within the two nodes.
 */
template <int D> double dot_scaling(const FunctionNode<D, double> &bra, const FunctionNode<D, double> &ket);

/**
 * @brief Dot product of wavelet parts on matching nodes (double).
 * @return Sum over matching wavelet blocks within the two nodes.
 */
template <int D> double dot_wavelet(const FunctionNode<D, double> &bra, const FunctionNode<D, double> &ket);

///@}

/** @name Per-node local dot-product helpers (complex-complex) */
///@{

template <int D> ComplexDouble dot_scaling(const FunctionNode<D, ComplexDouble> &bra, const FunctionNode<D, ComplexDouble> &ket);
template <int D> ComplexDouble dot_wavelet(const FunctionNode<D, ComplexDouble> &bra, const FunctionNode<D, ComplexDouble> &ket);

///@}

/** @name Per-node local dot-product helpers (complex-real and real-complex) */
///@{

template <int D> ComplexDouble dot_scaling(const FunctionNode<D, ComplexDouble> &bra, const FunctionNode<D, double> &ket);
template <int D> ComplexDouble dot_wavelet(const FunctionNode<D, ComplexDouble> &bra, const FunctionNode<D, double> &ket);

template <int D> ComplexDouble dot_scaling(const FunctionNode<D, double> &bra, const FunctionNode<D, ComplexDouble> &ket);
template <int D> ComplexDouble dot_wavelet(const FunctionNode<D, double> &bra, const FunctionNode<D, ComplexDouble> &ket);

///@}

} // namespace mrcpp