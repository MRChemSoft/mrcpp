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

#include "MRCPP/macros.h"
#include "utils/math_utils.h"
#include "utils/omp_utils.h"

#include "HilbertPath.h"
#include "MWTree.h"
#include "NodeIndex.h"

namespace mrcpp {

/**
 * @file MWNode.h
 * @brief Base node for multiresolution (multiwavelet) trees.
 *
 * @details
 * A node stores scaling/wavelet coefficients for one cell at scale `n` and
 * translation `l` in `D` spatial dimensions. It also keeps structural
 * information (parent/children, Hilbert path, status flags) and provides
 * utilities to:
 *
 * - allocate/attach coefficient buffers,
 * - compute and cache norms (total, per-component, maximum scaled norms),
 * - perform CV/MW transforms on the node,
 * - navigate and generate parts of the tree (parents/children),
 * - fetch geometry (bounds, center) and quadrature/child evaluation points.
 *
 * This class is templated on spatial dimension `D` (1, 2, or 3) and on the
 * scalar type `T` (e.g., `double` or `ComplexDouble`).
 */

/**
 * @class MWNode
 * @tparam D Spatial dimension (1, 2, or 3).
 * @tparam T Scalar type of coefficients (e.g., double, ComplexDouble).
 *
 * @brief Base class for multiwavelet tree nodes.
 *
 * @note
 * Nodes are created and managed by @ref MWTree and specialized trees
 * (e.g., @ref FunctionTree). Most users should not instantiate nodes
 * directly; instead, operate at the tree level.
 */
template <int D, typename T>
class MWNode {
public:
    /**
     * @brief Copy-construct a node.
     * @param node      Source node.
     * @param allocCoef If true, allocate a new coefficient buffer.
     * @param SetCoef   If true and @p allocCoef is true, copy coefficients.
     */
    MWNode(const MWNode<D, T> &node, bool allocCoef = true, bool SetCoef = true);

    MWNode<D, T> &operator=(const MWNode<D, T> &node) = delete;
    virtual ~MWNode();

    /// @name Basis/order and topology queries
    ///@{
    int getKp1() const { return getMWTree().getKp1(); }
    int getKp1_d() const { return getMWTree().getKp1_d(); }
    int getOrder() const { return getMWTree().getOrder(); }
    int getScalingType() const { return getMWTree().getMRA().getScalingBasis().getScalingType(); }
    int getTDim() const { return (1 << D); }
    int getDepth() const { return getNodeIndex().getScale() - getMWTree().getRootScale(); }
    int getScale() const { return getNodeIndex().getScale(); }
    int getNChildren() const { return (isBranchNode()) ? getTDim() : 0; }
    int getSerialIx() const { return this->serialIx; }
    void setSerialIx(int Ix) { this->serialIx = Ix; }

    const NodeIndex<D> &getNodeIndex() const { return this->nodeIndex; }
    const HilbertPath<D> &getHilbertPath() const { return this->hilbertPath; }
    ///@}

    /// @name Geometry
    ///@{
    Coord<D> getCenter() const;
    Coord<D> getUpperBounds() const;
    Coord<D> getLowerBounds() const;

    bool hasCoord(const Coord<D> &r) const;
    ///@}

    /// @name Structural relations
    ///@{
    bool isCompatible(const MWNode<D, T> &node);
    bool isAncestor(const NodeIndex<D> &idx) const;
    bool isDecendant(const NodeIndex<D> &idx) const;
    ///@}

    /// @name Norms
    ///@{
    double getSquareNorm() const { return this->squareNorm; }
    double getMaxSquareNorm() const { return (maxSquareNorm > 0.0) ? maxSquareNorm : calcScaledSquareNorm(); }
    double getMaxWSquareNorm() const { return (maxWSquareNorm > 0.0) ? maxWSquareNorm : calcScaledWSquareNorm(); }

    double getScalingNorm() const;
    virtual double getWaveletNorm() const;
    double getComponentNorm(int i) const { return this->componentNorms[i]; }
    ///@}

    /// @name Coefficients access
    ///@{
    int getNCoefs() const { return this->n_coefs; }
    void getCoefs(Eigen::Matrix<T, Eigen::Dynamic, 1> &c) const;
    void printCoefs() const;

    T *getCoefs() { return this->coefs; }
    const T *getCoefs() const { return this->coefs; }
    ///@}

    /// @name Evaluation points (quadrature / children)
    ///@{
    void getPrimitiveQuadPts(Eigen::MatrixXd &pts) const;
    void getPrimitiveChildPts(Eigen::MatrixXd &pts) const;
    void getExpandedQuadPts(Eigen::MatrixXd &pts) const;
    void getExpandedChildPts(Eigen::MatrixXd &pts) const;
    ///@}

    /// @name Tree navigation (typed accessors)
    ///@{
    MWTree<D, T> &getMWTree() { return static_cast<MWTree<D, T> &>(*this->tree); }
    MWNode<D, T> &getMWParent() { return static_cast<MWNode<D, T> &>(*this->parent); }
    MWNode<D, T> &getMWChild(int i) { return static_cast<MWNode<D, T> &>(*this->children[i]); }

    const MWTree<D, T> &getMWTree() const { return static_cast<const MWTree<D, T> &>(*this->tree); }
    const MWNode<D, T> &getMWParent() const { return static_cast<const MWNode<D, T> &>(*this->parent); }
    const MWNode<D, T> &getMWChild(int i) const { return static_cast<const MWNode<D, T> &>(*this->children[i]); }
    ///@}

    /// @name Coefficients editing (block-wise)
    ///@{
    void zeroCoefs();
    void setCoefBlock(int block, int block_size, const T *c);
    void addCoefBlock(int block, int block_size, const T *c);
    void zeroCoefBlock(int block, int block_size);
    void attachCoefs(T *coefs);
    ///@}

    /// @name Norm bookkeeping
    ///@{
    void calcNorms();
    void zeroNorms();
    void clearNorms();
    ///@}

    /// @name Topology modification
    ///@{
    virtual void createChildren(bool coefs);
    virtual void genChildren();
    virtual void genParent();
    virtual void deleteChildren();
    virtual void deleteParent();
    ///@}

    /// @name Local transforms
    ///@{
    virtual void cvTransform(int kind, bool firstchild = false);
    virtual void mwTransform(int kind);
    ///@}

    /**
     * @brief Node-norm at an arbitrary index.
     * @param idx Target index (may be at a finer scale).
     * @return A node-wise norm consistent with the basis and scale.
     */
    double getNodeNorm(const NodeIndex<D> &idx) const;

    /// @name Status flags
    ///@{
    bool hasParent() const { return (parent != nullptr) ? true : false; }
    bool hasCoefs() const { return (this->status & FlagHasCoefs); }
    bool isEndNode() const { return (this->status & FlagEndNode); }
    bool isGenNode() const { return (this->status & FlagGenNode); }
    bool isRootNode() const { return (this->status & FlagRootNode); }
    bool isLeafNode() const { return not(this->status & FlagBranchNode); }
    bool isAllocated() const { return (this->status & FlagAllocated); }
    bool isBranchNode() const { return (this->status & FlagBranchNode); }
    bool isLooseNode() const { return (this->status & FlagLooseNode); }
    bool checkStatus(unsigned char mask) const { return (mask == (this->status & mask)); }

    void setHasCoefs() { SET_BITS(status, FlagHasCoefs | FlagAllocated); }
    void setIsEndNode() { SET_BITS(status, FlagEndNode); }
    void setIsGenNode() { SET_BITS(status, FlagGenNode); }
    void setIsRootNode() { SET_BITS(status, FlagRootNode); }
    void setIsLeafNode() { CLEAR_BITS(status, FlagBranchNode); }
    void setIsAllocated() { SET_BITS(status, FlagAllocated); }
    void setIsBranchNode() { SET_BITS(status, FlagBranchNode); }
    void setIsLooseNode() { SET_BITS(status, FlagLooseNode); }
    void clearHasCoefs() { CLEAR_BITS(status, FlagHasCoefs); }
    void clearIsEndNode() { CLEAR_BITS(status, FlagEndNode); }
    void clearIsGenNode() { CLEAR_BITS(status, FlagGenNode); }
    void clearIsRootNode() { CLEAR_BITS(status, FlagRootNode); }
    void clearIsAllocated() { CLEAR_BITS(status, FlagAllocated); }
    ///@}

    friend std::ostream &operator<<(std::ostream &o, const MWNode<D, T> &nd) { return nd.print(o); }

    // Friend classes that are allowed to operate on internals.
    friend class TreeBuilder<D, T>;
    friend class MultiplicationCalculator<D, T>;
    friend class NodeAllocator<D, T>;
    friend class MWTree<D, T>;
    friend class FunctionTree<D, T>;
    friend class OperatorTree;
    friend class FunctionNode<D, T>;
    friend class OperatorNode;
    friend class DerivativeCalculator<D, T>;
    bool isComplex = false;               ///< Helper flag for mixed-real/complex workflows.
    friend class FunctionTree<D, double>; ///< Allows complex trees to access real nodes when needed.
    friend class FunctionTree<D, ComplexDouble>;
    int childSerialIx{-1};                ///< Index of first child in a serialized view, or -1 for leaves.

protected:
    // -------- Ownership and hierarchy --------
    MWTree<D, T> *tree{nullptr};          ///< Tree the node belongs to.
    MWNode<D, T> *parent{nullptr};        ///< Parent node (nullptr for roots).
    MWNode<D, T> *children[1 << D];       ///< Array of 2^D children (valid if branch node).

    // -------- Norms (cached) --------
    double squareNorm{-1.0};              ///< Squared norm of all 2^D (k+1)^D coefficients.
    double componentNorms[1 << D];        ///< Squared norms of the 2^D components.
    double maxSquareNorm{-1.0};           ///< Maximum scaled squared norm among node and descendants.
    double maxWSquareNorm{-1.0};          ///< Maximum scaled wavelet squared norm among node and descendants.

    // -------- Coefficients --------
    T *coefs{nullptr};                    ///< Buffer of size 2^D (k+1)^D with MW coefficients.
    int n_coefs{0};                       ///< Number of coefficients in @ref coefs.

    // -------- Serialization helpers --------
    int serialIx{-1};                     ///< Index in a serialized traversal.
    int parentSerialIx{-1};               ///< Serialized index of parent, or -1 for roots.

    // -------- Indexing and space-filling path --------
    NodeIndex<D> nodeIndex;               ///< Scale and translation of this node.
    HilbertPath<D> hilbertPath;           ///< Current Hilbert path state for child ordering.

    // -------- Construction helpers --------
    MWNode();
    MWNode(MWTree<D, T> *tree, int rIdx);
    MWNode(MWTree<D, T> *tree, const NodeIndex<D> &idx);
    MWNode(MWNode<D, T> *parent, int cIdx);

    /// Free coefficient buffer and reset counters.
    virtual void dealloc();

    /// Crop node based on precision; may trigger refinement.
    bool crop(double prec, double splitFac, bool absPrec);

    /// Initialize thread lock (when OpenMP is enabled).
    void initNodeLock() { MRCPP_INIT_OMP_LOCK(); }

    /// Allocate coefficient buffer as `n_blocks * block_size`.
    virtual void allocCoefs(int n_blocks, int block_size);

    /// Release coefficient buffer.
    virtual void freeCoefs();

    /// Update cached maxima from descendants.
    void setMaxSquareNorm();

    /// Invalidate cached maxima for this branch.
    void resetMaxSquareNorm();

    /// Scaled total norm \f$ 2^{D n}\|c\|^2 \f$ (lazy).
    double calcScaledSquareNorm() const { return std::pow(2.0, D * getScale()) * getSquareNorm(); }

    /// Scaled wavelet norm \f$ 2^{D n}\|d\|^2 \f$ (lazy).
    double calcScaledWSquareNorm() const { return std::pow(2.0, D * getScale()) * getWaveletNorm(); }

    /// Component-wise norm computation hook.
    virtual double calcComponentNorm(int i) const;

    /// Recompress local representation after edits.
    virtual void reCompress();

    /// Push coefficients from parent to all children.
    virtual void giveChildrenCoefs(bool overwrite = true);

    /// Push coefficients from parent to a specific child.
    virtual void giveChildCoefs(int cIdx, bool overwrite = true);

    /// Pull coefficients from children to parent.
    virtual void giveParentCoefs(bool overwrite = true);

    /// Rebuild local buffer from children (inverse of giveChildrenCoefs).
    virtual void copyCoefsFromChildren();

    /// Child index for a target node index (same scale or finer).
    int getChildIndex(const NodeIndex<D> &nIdx) const;

    /// Child index for a spatial coordinate.
    int getChildIndex(const Coord<D> &r) const;

    /// Whether two nodes lie in different branches (fast check).
    bool diffBranch(const MWNode<D, T> &rhs) const;

    /// Retrieve node owning coordinate @p r at given depth (may create).
    MWNode<D, T> *retrieveNode(const Coord<D> &r, int depth);

    /// Retrieve node at index @p idx (may create).
    MWNode<D, T> *retrieveNode(const NodeIndex<D> &idx, bool create = false);

    /// Retrieve parent node for index @p idx (may create ancestors).
    MWNode<D, T> *retrieveParent(const NodeIndex<D> &idx);

    /// Lookup without generation (const).
    const MWNode<D, T> *retrieveNodeNoGen(const NodeIndex<D> &idx) const;

    /// Lookup without generation (mutable).
    MWNode<D, T> *retrieveNodeNoGen(const NodeIndex<D> &idx);

    /// Find node or end node by coordinate (const).
    const MWNode<D, T> *retrieveNodeOrEndNode(const Coord<D> &r, int depth) const;

    /// Find node or end node by coordinate (mutable).
    MWNode<D, T> *retrieveNodeOrEndNode(const Coord<D> &r, int depth);

    /// Find node or end node by index (const).
    const MWNode<D, T> *retrieveNodeOrEndNode(const NodeIndex<D> &idx) const;

    /// Find node or end node by index (mutable).
    MWNode<D, T> *retrieveNodeOrEndNode(const NodeIndex<D> &idx);

    /// Thread-safe child creation.
    void threadSafeCreateChildren();

    /// Thread-safe generation of children.
    void threadSafeGenChildren();

    /// Remove nodes generated during adaptive build.
    void deleteGenerated();

    /// Printable diagnostics for a node.
    virtual std::ostream &print(std::ostream &o) const;

    // --- Bit flags describing node state (see status member) ---
    static const unsigned char FlagBranchNode = B8(00000001);
    static const unsigned char FlagGenNode    = B8(00000010);
    static const unsigned char FlagHasCoefs   = B8(00000100);
    static const unsigned char FlagAllocated  = B8(00001000);
    static const unsigned char FlagEndNode    = B8(00010000);
    static const unsigned char FlagRootNode   = B8(00100000);
    static const unsigned char FlagLooseNode  = B8(01000000);

private:
    unsigned char status{0};  ///< Bit mask of @ref FlagBranchNode, @ref FlagGenNode, etc.

#ifdef MRCPP_HAS_OMP
    omp_lock_t omp_lock;      ///< Per-node lock for thread-safe edits (OpenMP).
#endif
};

} // namespace mrcpp