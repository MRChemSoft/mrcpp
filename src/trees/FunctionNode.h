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
 * @class FunctionNode
 * @tparam D Spatial dimension (1, 2, or 3)
 * @tparam T Coefficient type (e.g. double, ComplexDouble)
 *
 * @brief Node of a @ref FunctionTree that stores coefficients and implements
 * function-specific operations
 *
 * @details A FunctionNode is a concrete @ref MWNode specialized for function
 * representations.  It holds scaling and wavelet coefficients, provides
 * allocation and refinement helpers, and exposes utilities for evaluation,
 * coefficient access and basic per-node operations such as integration and
 * local dot products.
 *
 * @note FunctionNodes are managed by @ref FunctionTree and @ref NodeAllocator.
 * Most users should not construct FunctionNode directly.
 */
template <int D, typename T> class FunctionNode final : public MWNode<D, T> {
public:
    FunctionTree<D, T> &getFuncTree() { return static_cast<FunctionTree<D, T> &>(*this->tree); }        ///< @return A reference to the tree this node belongs to, cast to a non-const @ref FunctionTree
    FunctionNode<D, T> &getFuncParent() { return static_cast<FunctionNode<D, T> &>(*this->parent); }    ///< @return A reference to the parent of this node, cast to a non-const @ref FunctionNode

    /**
     * @param i The index of the child
     * @return A reference to the child at the given index, cast to a non-const @ref FunctionNode
     */
    FunctionNode<D, T> &getFuncChild(int i) { return static_cast<FunctionNode<D, T> &>(*this->children[i]); }

    const FunctionTree<D, T> &getFuncTree() const { return static_cast<const FunctionTree<D, T> &>(*this->tree); }      ///< @return A reference to the tree this node belongs to, cast to a const @ref FunctionTree
    const FunctionNode<D, T> &getFuncParent() const { return static_cast<const FunctionNode<D, T> &>(*this->parent); }  ///< @return A reference to the parent of this node, cast to a const @ref FunctionNode

    /**
     * @param i The index of the child
     * @return A reference to the child at the given index, cast to a const @ref FunctionNode
     */
    const FunctionNode<D, T> &getFuncChild(int i) const { return static_cast<const FunctionNode<D, T> &>(*this->children[i]); }

    /**
     * @brief Create child nodes for this node and abort if it is already a
     * branch node (see @ref FlagBranchNode)
     * @param coefs If `true`, allocate coefficient chunk for child nodes
     *
     * @details This routine allocates child nodes via the tree's @ref NodeAllocator.
     * The tree's node counter is incremented by @ref MWTree<D, T>::incrementNodeCount.
     * Finally, this node is marked as both a branch node (see @ref FlagBranchNode)
     * and a non-end node (see @ref FlagEndNode).
     */
    void createChildren(bool coefs) override;

    /**
     * @brief Generates child nodes with the @ref FlagGenNode bit flag set, and
     * abort if this node is already a branch node (see @ref FlagBranchNode)
     *
     * @details This routine creates general or redundant child nodes for
     * temporary use. As a result, the tree's node counter remains unchanged,
     * and this node is marked only as a branch node (see @ref FlagBranchNode).
     */
    void genChildren() override;

    /**
     * @brief Generate a parent for this node and abort if it already has one
     *
     * @details This routine allocates the parent node via the tree's @ref
     * NodeAllocator and links this node into the parent children array. The
     * tree's node counter is incremented by @ref MWTree<D, T>::incrementNodeCount.
     */
    void genParent() override;

    /**
     * @brief Recursive deallocation of children and all their descendants.
     *
     * @details This routine uses base class function @ref MWTree<D, T>::deleteChildren
     * for the deallocation. Finally, this node is marked as an end node (see
     * @ref FlagEndNode).
     */
    void deleteChildren() override;

    /**
     * @brief Function integration
     * @return The integral of type @p T
     *
     * @details Wrapper for function integration, that requires different
     * methods depending on scaling type @ref FuncType. Integrates the function
     * represented on the node on the full support of the node. This routine
     * will return zero if the node does not have coeffcients, and abort if the
     * node has invalid type of scaling basis (Legendre or Interpol; see
     * MRCPP/constants.h).
     */
    T integrate() const;

    /**
     * @brief Set values from a vector to the node's coefficients, and update
     * metadata of the node
     * @param vec Column vector
     *
     * @details This routine calls @ref MWTree<D, T>::setCoefBlock to set
     * values from the vector, and update metadata of the node by caling
     * @ref MWTree<D, T>::cvTransform, and @ref MWTree<D, T>::mwTransform. The
     * node is marked as having coefficients, and its square norm and component
     * norms are also computed by @ref MWTree<D, T>::calcNorms.
     */
    void setValues(const Eigen::Matrix<T, Eigen::Dynamic, 1> &vec);

    /**
     * @brief Extract the node's coefficients into a vector
     * @param[out] vec Column vector resized to the number of coefficients of
     * the node, see @ref MWTree<D, T>::getNCoefs
     */
    void getValues(Eigen::Matrix<T, Eigen::Dynamic, 1> &vec);

    /**
     * @brief Get coefficients corresponding to absolute value of function
     * @param[out] absCoefs Coefficients of type @p T
     *
     * @note The absolute value of function is computed using std::norm().
     */
    void getAbsCoefs(T *absCoefs);

    friend class FunctionTree<D, T>;
    friend class NodeAllocator<D, T>;

protected:
    /**
     * @brief FunctionNode constructor
     *
     * @note This routine uses @ref MWNode default constructor.
     */
    FunctionNode()
            : MWNode<D, T>() {}

    /**
     * @brief FunctionNode constructor
     * @param[in] tree The MWTree the root node belongs to
     * @param[in] rIdx The integer specifying the corresponding root node
     *
     * @details Constructor for root nodes. It actually calls @ref MWNode
     * contructor MWNode<D, T>(tree, rIdx).
     */
    explicit FunctionNode(MWTree<D, T> *tree, int rIdx)
            : MWNode<D, T>(tree, rIdx) {}

    /**
     * @brief FunctionNode constructor
     * @param[in] tree The MWTree the root node belongs to
     * @param[in] idx The NodeIndex defining scale and translation of the node
     *
     * @details Constructor for an empty node, which calls @ref MWNode
     * contructor MWNode<D, T>(tree, idx).
     */
    FunctionNode(MWTree<D, T> *tree, const NodeIndex<D> &idx)
            : MWNode<D, T>(tree, idx) {}

    /**
     * @brief FunctionNode constructor
     * @param[in] parent Parent node
     * @param[in] cIdx Child index of the current node
     *
     * @details Constructor for leaf nodes. It invokes @ref MWNode constructor
     * MWNode<D, T>(parent, cIdx).
     */
    FunctionNode(MWNode<D, T> *parent, int cIdx)
            : MWNode<D, T>(parent, cIdx) {}

    FunctionNode(const FunctionNode<D, T> &node) = delete;
    FunctionNode<D, T> &operator=(const FunctionNode<D, T> &node) = delete;

    /// @brief Default destructor of FunctionNode
    ~FunctionNode() = default;

    /**
     * @brief Function evaluation
     * @param[in,out] r The sought after node through the coordinates of a point in space
     * @return The evaluated result of type @p T
     *
     * @details Evaluate all polynomials defined on the child node found by
     * @ref MWTree<D, T>::getChildIndex. Trigger an error if the node does not
     * have coefficients, or failed to find the child node. For periodic
     * systems, the coordinate r will be mapped to the [-1, 1] periodic cell if
     * it is outside the unit cell, see @ref periodic::coord_manipulation.
     */
    //FIXME I guess the evaluation is performed for the child node, not this
    //      node, because this routine calls getFuncChild(cIdx).evalScaling(r),
    //      where cIdx = this->getChildIndex(r)
    T evalf(Coord<D> r);

    /**
     * @brief Function evaluation
     * @param[in] r Coordinate where the evaluation is performed at
     * @return The evaluated result of type @p T
     */
    T evalScaling(const Coord<D> &r) const;

    /// @brief Deallocate the node and detach it from the tree it belongs to
    //FIXME this routine calls dealloc of NodeAllocator, which frees the memory
    //      of the node by calling ~MWNode() and it seems that coefficients are freed
    //      only for LooseNode.
    void dealloc() override;

    /**
     * @brief Update the coefficients of the node by an MW transform of the
     * scaling coefficients of the children
     * @note There is a specialization for @p D = 3,
     * see @ref FunctionNode<3>::reCompress.
     */
    //FIXME It is written in FunctionNode.cpp that, "Option to overwrite or add
    //      up existing coefficients". Not sure what it means, in particular "add up
    //      existing coefficients".
    //FIXME not sure if @ref FunctionNode<3>::reCompress works
    void reCompress() override;

    /**
     * @brief Function integration, Legendre basis
     * @return The integral of type @p T
     *
     * @details Integrate the function represented on the node on the full
     * support of the node. The Legendre basis is particularly easy to
     * integrate, as the work is already done when calculating its
     * coefficients. The coefficients of the node is defined as the projection
     * integral \f$ s_i = \int f(x)\phi_i(x)\mathrm{d}x \f$ and since the first
     * Legendre function is the constant 1, the first coefficient is simply the
     * integral of \f$ f(x) \f$.
     */
    T integrateLegendre() const;

    /**
     * @brief Function integration, Interpolating basis
     * @return The integral of type @p T
     *
     * @details Integrate the function represented on the node on the full
     * support of the node. A bit more involved than in the Legendre basis, as
     * is requires some coupling of quadrature weights.
     */
    T integrateInterpolating() const;

    /**
     * @brief Function integration, Interpolating basis
     * @return The integral of type @p T
     *
     * @details Integrate the function represented on the node on the full
     * support of the node. A bit more involved than in the Legendre basis, as
     * is requires some coupling of quadrature weights.
     */
    //FIXME This routine has exactly the same documentation comment as
    //      integrateInterpolating() in FunctionNode.cpp.
    T integrateValues() const;
};

//FIXME All comments of dot_scaling() are exactly the same in FunctionNode.cpp,
//      a bit boilerplate. But another important thing is the conjugate is actually
//      taken for bra and ket in case of complex values, instead of only bra.
//FIXME Comments of dot_wavelet() have the same problem.

/**
 * @brief Inner product of the functions represented by the scaling basis of
 * the nodes
 * @param[in] bra FunctionNode on bra
 * @param[in] ket FunctionNode on ket
 * @return The computed inner product
 *
 * @details Integrates the product of the functions represented by the scaling
 * basis on the node on the full support of the nodes. The scaling basis is
 * fully orthonormal, and the inner product is simply the dot product of the
 * coefficient vectors. Assumes the nodes have identical support.
 */
template <int D> double dot_scaling(const FunctionNode<D, double> &bra, const FunctionNode<D, double> &ket);

/**
 * @brief Inner product of the functions represented by the scaling basis of
 * the nodes
 * @param[in] bra FunctionNode on bra
 * @param[in] ket FunctionNode on ket
 * @return The computed inner product
 *
 * @details Integrates the product of the functions represented by the scaling
 * basis on the node on the full support of the nodes. The scaling basis is
 * fully orthonormal, and the inner product is simply the dot product of the
 * coefficient vectors. Assumes the nodes have identical support.
 * @note Conjugates of bra and ket will be taken.
 */
template <int D> ComplexDouble dot_scaling(const FunctionNode<D, ComplexDouble> &bra, const FunctionNode<D, ComplexDouble> &ket);

/**
 * @brief Inner product of the functions represented by the scaling basis of
 * the nodes
 * @param[in] bra FunctionNode on bra
 * @param[in] ket FunctionNode on ket
 * @return The computed inner product
 *
 * @details Integrates the product of the functions represented by the scaling
 * basis on the node on the full support of the nodes. The scaling basis is
 * fully orthonormal, and the inner product is simply the dot product of the
 * coefficient vectors. Assumes the nodes have identical support.
 * @note Conjugate of bra will be taken.
 */
template <int D> ComplexDouble dot_scaling(const FunctionNode<D, ComplexDouble> &bra, const FunctionNode<D, double> &ket);

/**
 * @brief Inner product of the functions represented by the scaling basis of
 * the nodes
 * @param[in] bra FunctionNode on bra
 * @param[in] ket FunctionNode on ket
 * @return The computed inner product
 *
 * @details Integrates the product of the functions represented by the scaling
 * basis on the node on the full support of the nodes. The scaling basis is
 * fully orthonormal, and the inner product is simply the dot product of the
 * coefficient vectors. Assumes the nodes have identical support.
 * @note Conjugate of ket will be taken.
 */
template <int D> ComplexDouble dot_scaling(const FunctionNode<D, double> &bra, const FunctionNode<D, ComplexDouble> &ket);


/**
 * @brief Inner product of the functions represented by the wavelet basis of
 * the nodes
 * @param[in] bra FunctionNode on bra
 * @param[in] ket FunctionNode on ket
 * @return The computed inner product
 *
 * @details Integrates the product of the functions represented by the wavelet
 * basis on the node on the full support of the nodes. The wavelet basis is
 * fully orthonormal, and the inner product is simply the dot product of the
 * coefficient vectors. Assumes the nodes have identical support.
 */
template <int D> double dot_wavelet(const FunctionNode<D, double> &bra, const FunctionNode<D, double> &ket);

/**
 * @brief Inner product of the functions represented by the wavelet basis of
 * the nodes
 * @param[in] bra FunctionNode on bra
 * @param[in] ket FunctionNode on ket
 * @return The computed inner product
 *
 * @details Integrates the product of the functions represented by the wavelet
 * basis on the node on the full support of the nodes. The wavelet basis is
 * fully orthonormal, and the inner product is simply the dot product of the
 * coefficient vectors. Assumes the nodes have identical support.
 * @note Conjugates of bra and ket will be taken.
 */
template <int D> ComplexDouble dot_wavelet(const FunctionNode<D, ComplexDouble> &bra, const FunctionNode<D, ComplexDouble> &ket);

/**
 * @brief Inner product of the functions represented by the wavelet basis of
 * the nodes
 * @param[in] bra FunctionNode on bra
 * @param[in] ket FunctionNode on ket
 * @return The computed inner product
 *
 * @details Integrates the product of the functions represented by the wavelet
 * basis on the node on the full support of the nodes. The wavelet basis is
 * fully orthonormal, and the inner product is simply the dot product of the
 * coefficient vectors. Assumes the nodes have identical support.
 * @note Conjugate of bra will be taken.
 */
template <int D> ComplexDouble dot_wavelet(const FunctionNode<D, ComplexDouble> &bra, const FunctionNode<D, double> &ket);

/**
 * @brief Inner product of the functions represented by the wavelet basis of
 * the nodes
 * @param[in] bra FunctionNode on bra
 * @param[in] ket FunctionNode on ket
 * @return The computed inner product
 *
 * @details Integrates the product of the functions represented by the wavelet
 * basis on the node on the full support of the nodes. The wavelet basis is
 * fully orthonormal, and the inner product is simply the dot product of the
 * coefficient vectors. Assumes the nodes have identical support.
 * @note Conjugate of ket will be taken.
 */
template <int D> ComplexDouble dot_wavelet(const FunctionNode<D, double> &bra, const FunctionNode<D, ComplexDouble> &ket);

//FIXME There are template specializations in FunctionNode.cpp, do we
//      need to document them as well?

} // namespace mrcpp
