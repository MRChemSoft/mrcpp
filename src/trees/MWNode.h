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
 * @class MWNode
 * @tparam D Spatial dimension (1, 2, or 3).
 * @tparam T Coefficient type (e.g. double, ComplexDouble).
 *
 * @brief Base class for Multiwavelet nodes
 *
 * @details A MWNode will contain the scaling and wavelet coefficients
 * to represent functions or operators within a Multiwavelet
 * framework. The nodes are multidimensional. The dimensionality is
 * set through the template parameter D=1,2,3. In addition to the
 * coefficients, the node contains metadata such as the scale, the
 * translation index, the norm, pointers to parent node and child
 * nodes, pointer to the corresponding MWTree etc... See member and
 * data descriptions for details.
 *
 * @note
 * Nodes are created and managed by MWTree and specialized trees
 * (e.g., FunctionTree). Most users should not instantiate nodes
 * directly; instead, operate at the tree level.
 */
template <int D, typename T>
class MWNode {
public:
    /**
     * @brief MWNode copy constructor
     * @param[in] node  The original node
     * @param allocCoef If true, allocate MW coefficients and copy from the original node
     * @param SetCoef   If true and @p allocCoef is true, copy coefficients
     *
     * @details Creates loose nodes and optionally copy coefs. The node
     * does not "belong" to the tree: It cannot be accessed by traversing
     * the tree.
     */
    MWNode(const MWNode<D, T> &node, bool allocCoef = true, bool SetCoef = true);

    MWNode<D, T> &operator=(const MWNode<D, T> &node) = delete;

    /// @brief Recursive deallocation of a node and all its decendants
    virtual ~MWNode();

    /*
     * Getters and setters
     */
    int getOrder() const { return getMWTree().getOrder(); }                                         ///< @return Polynomial order k
    int getKp1() const { return getMWTree().getKp1(); }                                             ///< @return k+1
    int getKp1_d() const { return getMWTree().getKp1_d(); }                                         ///< @return (k+1)^D
    int getScalingType() const { return getMWTree().getMRA().getScalingBasis().getScalingType(); }  ///< @return The type of scaling basis (Legendre or Interpol; see MRCPP/constants.h)
    int getTDim() const { return (1 << D); }                                                        ///< @return 2^D (number of children per internal node)
    int getDepth() const { return getNodeIndex().getScale() - getMWTree().getRootScale(); }         ///< @return The depth of this node
    int getScale() const { return getNodeIndex().getScale(); }                                      ///< @return The scale of this node
    int getNChildren() const { return (isBranchNode()) ? getTDim() : 0; }                           ///< @return The number of children of this node
    int getSerialIx() const { return this->serialIx; }                                              ///< @return The index in the serial tree
    void setSerialIx(int Ix) { this->serialIx = Ix; }                                               ///< @param Ix The index in the serial tree

    const NodeIndex<D> &getNodeIndex() const { return this->nodeIndex; }                            ///< @return The index (scale and translation) for this node
    const HilbertPath<D> &getHilbertPath() const { return this->hilbertPath; }                      // TODO document this

    Coord<D> getCenter() const;         ///< @return The coordinates of the centre of the node
    Coord<D> getUpperBounds() const;    ///< @return The upper bounds of the D-interval defining the node
    Coord<D> getLowerBounds() const;    ///< @return The lower bounds of the D-interval defining the node

    /**
     * @brief Test if a given coordinate is within the boundaries of the node
     * @param[in] r Point coordinates
     */
    bool hasCoord(const Coord<D> &r) const;

    /// @warning This method is currently not implemented.
    bool isCompatible(const MWNode<D, T> &node);

    /**
     * @brief Test if the node is decending from a given NodeIndex, that is, if they have
     * overlapping support.
     * @param[in] idx the NodeIndex of the requested node
     */
    bool isAncestor(const NodeIndex<D> &idx) const;

    /// @warning This method is currently not implemented.
    bool isDecendant(const NodeIndex<D> &idx) const;

    double getSquareNorm() const { return this->squareNorm; }                                                       ///< @return Squared norm of all 2^D (k+1)^D coefficients
    double getMaxSquareNorm() const { return (maxSquareNorm > 0.0) ? maxSquareNorm : calcScaledSquareNorm(); }      ///< @return Largest squared norm among itself and descendants.
    double getMaxWSquareNorm() const { return (maxWSquareNorm > 0.0) ? maxWSquareNorm : calcScaledWSquareNorm(); }  ///< @return Largest wavelet squared norm among itself and descendants.

    /**
     * @brief Calculate and return the squared scaling norm
     * @return The scaling norm
    */
    double getScalingNorm() const;
    /**
     * @brief Calculate and return the squared wavelet norm
     * @return The squared wavelet norm
     */
    virtual double getWaveletNorm() const;
    /**
     * @param i The component index
     * @return The squared norm of the component at the given index
     */
    double getComponentNorm(int i) const { return this->componentNorms[i]; }

    int getNCoefs() const { return this->n_coefs; }                 ///< @return The number of coefficients
    /**
     * @brief Wraps the MW coefficients into an Eigen vector object
     * @param[out] c The coefficient matrix
     */
    void getCoefs(Eigen::Matrix<T, Eigen::Dynamic, 1> &c) const;

    void printCoefs() const; ///< @brief Printout of node coefficients

    T *getCoefs() { return this->coefs; }               ///< @return The 2^D (k+1)^D MW coefficients
    const T *getCoefs() const { return this->coefs; }   ///< @return The 2^D (k+1)^D MW coefficients

    /**
     * @brief Returns the quadrature points of this node
     *
     * @param[out] pts Quadrature points in a \f$ d \times (k+1) \f$ matrix form.
     *
     * @details The original quadrature points are fetched and then
     * dilated and translated. For each cartesian direction \f$ \alpha =
     * x,y,z... \f$ the set of quadrature points becomes \f$ x^\alpha_i =
     * 2^{-n} (x_i + l^\alpha \f$. By taking all possible
     * \f$(k+1)^d\f$ combinations, they will then define a d-dimensional
     * grid of quadrature points.
     */
    void getPrimitiveQuadPts(Eigen::MatrixXd &pts) const;

    /**
     * @brief Returns the quadrature points of this node
     *
     * @param[out] pts Quadrature points in a \f$ d \times (k+1) \f$ matrix form.
     *
     * @details The original quadrature points are fetched and then
     * dilated and translated to match the quadrature points in the
     * children of this node. For each cartesian direction \f$ \alpha = x,y,z... \f$
     * the set of quadrature points becomes \f$ x^\alpha_i = 2^{-n-1} (x_i + 2 l^\alpha + t^\alpha) \f$, where \f$ t^\alpha =
     * 0,1 \f$. By taking all possible \f$(k+1)^d\f$ combinations, they will
     * then define a d-dimensional grid of quadrature points for the child
     * nodes.
     */
    void getPrimitiveChildPts(Eigen::MatrixXd &pts) const;

    /**
     * @brief Returns the quadrature points of this node
     *
     * @param[out] pts Expanded quadrature points in a \f$ d \times
     * (k+1)^d \f$ matrix form.
     *
     * @details The primitive quadrature points are used to obtain a
     * tensor-product representation collecting all \f$ (k+1)^d \f$
     * vectors of quadrature points.
     */
    void getExpandedQuadPts(Eigen::MatrixXd &pts) const;

    /**
     * @brief Returns the quadrature points of this node
     *
     * @param[out] pts Expanded quadrature points in a \f$ d \times
     * 2^d(k+1)^d \f$ matrix form.
     *
     * @details The primitive quadrature points of the children are used to obtain a
     * tensor-product representation collecting all \f$ 2^d (k+1)^d \f$
     * vectors of quadrature points.
     */
    void getExpandedChildPts(Eigen::MatrixXd &pts) const;

    MWTree<D, T> &getMWTree() { return static_cast<MWTree<D, T> &>(*this->tree); }              ///< @return The tree this node belongs to
    MWNode<D, T> &getMWParent() { return static_cast<MWNode<D, T> &>(*this->parent); }          ///< @return The parent of this node

    /**
     * @param i The index of the child
     * @return The child at the given index
     */
    MWNode<D, T> &getMWChild(int i) { return static_cast<MWNode<D, T> &>(*this->children[i]); }

    const MWTree<D, T> &getMWTree() const { return static_cast<const MWTree<D, T> &>(*this->tree); }                ///< @return The tree this node belongs to
    const MWNode<D, T> &getMWParent() const { return static_cast<const MWNode<D, T> &>(*this->parent); }            ///< @return The parent of this node

    /**
     * @param i The index of the child
     * @return The child at the given index
     */
    const MWNode<D, T> &getMWChild(int i) const { return static_cast<const MWNode<D, T> &>(*this->children[i]); }

    /// @brief Sets all MW coefficients and the norms to zero
    void zeroCoefs();

    /**
     * @brief Assigns values to a block of coefficients
     * @param block The block index
     * @param block_size Size of the block
     * @param[in] c The input coefficients
     *
     * @details A block is typically containing one kind of coefficients
     * (given scaling/wavelet in each direction). Its size is then \f$
     * (k+1)^D \f$ and the index is between 0 and \f$ 2^D-1 \f$.
     */
    void setCoefBlock(int block, int block_size, const T *c);

    /**
     * @brief Adds values to a block of coefficients
     * @param block The block index
     * @param block_size Size of the block
     * @param[in] c The input coefficients
     *
     * @details A block is typically containing one kind of coefficients
     * (given scaling/wavelet in each direction). Its size is then \f$
     * (k+1)^D \f$ and the index is between 0 and \f$ 2^D-1 \f$.
     */
    void addCoefBlock(int block, int block_size, const T *c);

    /**
     * @brief Sets values of a block of coefficients to zero
     * @param[in] block The block index
     * @param[in] block_size Size of the block
     *
     * @details A block is typically containing one kind of coefficients
     * (given scaling/wavelet in each direction). Its size is then \f$
     * (k+1)^D \f$ and the index is between 0 and \f$ 2^D-1 \f$.
     */
    void zeroCoefBlock(int block, int block_size);

    /**
     * @brief Attach a set of coefficients to this node. Only used locally (the tree is not aware of this).
     * @param[in] coefs The coefficients to attach
     *
     * @note The number of coefficients must remain the same.
     */
    void attachCoefs(T *coefs);

    void calcNorms();   ///< @brief Calculate and store square norm and component norms, if allocated.
    void zeroNorms();   ///< @brief Set all norms to zero.
    void clearNorms();  ///< @brief Set all norms to Undefined.

    /*
     * Implemented in child classes
     */
    virtual void createChildren(bool coefs);
    virtual void genChildren();
    virtual void genParent();

    /**
     * @brief Recursive deallocation of children and all their descendants.
     *
     * @details Leaves node as LeafNode and children[] as null pointer.
     */
    virtual void deleteChildren();

    /// @brief Recursive deallocation of parent and all their forefathers.
    virtual void deleteParent();

    /**
     * @brief Coefficient-Value transform
     * @param operation Forward (coef->value) or backward (value->coef).
     *
     * @details This routine transforms the scaling coefficients of the node to the
     * function values in the corresponding quadrature roots (of its children).
     *
     * @note This routine assumes a 0/1 (scaling on child 0 and 1)
     *       representation, instead of s/d (scaling and wavelet).
     */
    virtual void cvTransform(int operation, bool firstchild = false); // TODO document firstchild parameter

    /**
     * @brief Multiwavelet transform
     * @param operation compression (s0,s1->s,d) or reconstruction (s,d->s0,s1).
     *
     * @details Application of the filters on one node to pass from a 0/1 (scaling
     * on child 0 and 1) representation to an s/d (scaling and
     * wavelet) representation. Bit manipulation is used in order to
     * determine the correct filters and whether to apply them or just
     * pass to the next couple of indexes. The starting coefficients are
     * preserved until the application is terminated, then they are
     * overwritten. With minor modifications this code can also be used
     * for the inverse mw transform (just use the transpose filters) or
     * for the application of an operator (using A, B, C and T parts of an
     * operator instead of G1, G0, H1, H0). This is the version where the
     * three directions are operated one after the other. Although this
     * is formally faster than the other algorithm, the separation of the
     * three dimensions prevent the possibility to use the norm of the
     * operator in order to discard a priori negligible contributions.
     *
     */
    virtual void mwTransform(int operation);

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
    int serialIx{-1};                     ///< Index in the serial tree
    int parentSerialIx{-1};               ///< Index of parent in the serial tree, or -1 for roots

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