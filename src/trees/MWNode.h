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
 * @tparam D Spatial dimension (1, 2, or 3)
 * @tparam T Coefficient type (e.g. double, ComplexDouble)
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
 * @note Nodes are created and managed by MWTree and specialized trees
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
     * overlapping support
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
     * @param[out] pts Quadrature points in a \f$ d \times (k+1) \f$ matrix form
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
     * @param[out] pts Quadrature points in a \f$ d \times (k+1) \f$ matrix form
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
     * (k+1)^d \f$ matrix form
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
     * 2^d(k+1)^d \f$ matrix form
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
     * @brief Recursive deallocation of children and all their descendants
     *
     * @details Leaves node as LeafNode and children[] as null pointer
     */
    virtual void deleteChildren();

    /// @brief Recursive deallocation of parent and all their forefathers.
    virtual void deleteParent();

    /**
     * @brief Coefficient-Value transform
     * @param operation Forward (coef->value) or backward (value->coef)
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
     * @param operation Compression (s0,s1->s,d) or reconstruction (s,d->s0,s1)
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
     */
    virtual void mwTransform(int operation);

    /**
     * @brief Gives the norm (absolute value) of the node at the given NodeIndex
     * @param[in] idx the NodeIndex of the requested node
     *
     * @details
     * Recursive routine to find the node with a given NodeIndex. When an EndNode is
     * found, do not generate any new node, but rather give the value of the norm
     * assuming the function is uniformly distributed within the node.
     */
    double getNodeNorm(const NodeIndex<D> &idx) const;

    /*
     * Getters and setters
     */
    bool hasParent() const { return (parent != nullptr) ? true : false; }   ///< @return Whether the node hsa a parent
    bool hasCoefs() const { return (this->status & FlagHasCoefs); }         ///< @return Whether the node has coefficients
    bool isEndNode() const { return (this->status & FlagEndNode); }         ///< @return Whether the node is an end node
    bool isGenNode() const { return (this->status & FlagGenNode); }         ///< @return Whether the node is a generated node
    bool isRootNode() const { return (this->status & FlagRootNode); }       ///< @return Whether the node is a root node
    bool isLeafNode() const { return not(this->status & FlagBranchNode); }  ///< @return Whether the node is a leaf node
    bool isAllocated() const { return (this->status & FlagAllocated); }     ///< @return Whether the node is fully allocated
    bool isBranchNode() const { return (this->status & FlagBranchNode); }   ///< @return Whether the node is a leaf node
    bool isLooseNode() const { return (this->status & FlagLooseNode); }     ///< @return Whether the node is a loose node

    /**
     * @brief Allows checking the state of a node against a state mask
     * @param mask The status mask to compare against
     * @return Whether the state of the node matches the given mask
     */
    bool checkStatus(unsigned char mask) const { return (mask == (this->status & mask)); }

    void setHasCoefs() { SET_BITS(status, FlagHasCoefs | FlagAllocated); }  ///< @brief Marks the node as having coefficients
    void setIsEndNode() { SET_BITS(status, FlagEndNode); }                  ///< @brief Marks the node as an end node
    void setIsGenNode() { SET_BITS(status, FlagGenNode); }                  ///< @brief Marks the node as a generated node
    void setIsRootNode() { SET_BITS(status, FlagRootNode); }                ///< @brief Marks the node as a root node
    void setIsLeafNode() { CLEAR_BITS(status, FlagBranchNode); }            ///< @brief Marks the node as a leaf node
    void setIsAllocated() { SET_BITS(status, FlagAllocated); }              ///< @brief Marks the node as allocated
    void setIsBranchNode() { SET_BITS(status, FlagBranchNode); }            ///< @brief Marks the node as a leaf node
    void setIsLooseNode() { SET_BITS(status, FlagLooseNode); }              ///< @brief Marks the node as a loose node
    void clearHasCoefs() { CLEAR_BITS(status, FlagHasCoefs); }              ///< @brief Clears the mark for having coefficients
    void clearIsEndNode() { CLEAR_BITS(status, FlagEndNode); }              ///< @brief Clears the mark for being an end node
    void clearIsGenNode() { CLEAR_BITS(status, FlagGenNode); }              ///< @brief Clears the mark for being a generated node
    void clearIsRootNode() { CLEAR_BITS(status, FlagRootNode); }            ///< @brief Clears the mark for being a root node
    void clearIsAllocated() { CLEAR_BITS(status, FlagAllocated); }          ///< @brief Clears the mark for being allocated

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
    bool isComplex = false;                 // TODO put as one of the flags
    friend class FunctionTree<D, double>;   // required if a ComplexDouble tree access a double node from another tree!
    friend class FunctionTree<D, ComplexDouble>;
    int childSerialIx{-1};                  ///< index of first child in a serial tree, or -1 for leaf nodes/end nodes

protected:
    MWTree<D, T> *tree{nullptr};            ///< Tree the node belongs to
    MWNode<D, T> *parent{nullptr};          ///< Parent node (nullptr for root nodes)
    MWNode<D, T> *children[1 << D];         ///< Array of 2^D children (valid if branch node)

    double squareNorm{-1.0};                ///< Squared norm of all 2^D (k+1)^D coefficients
    double componentNorms[1 << D];          ///< Squared norms of the separated 2^D components
    double maxSquareNorm{-1.0};             ///< Maximum squared norm among the node and descendants
    double maxWSquareNorm{-1.0};            ///< Maximum wavelet squared norm among the node and descendants
                                            ///< NB: must be set before used.
    T *coefs{nullptr};                      ///< The 2^D (k+1)^D MW coefficients
                                            ///< For example, in case of a one dimensional function \f$ f \f$
                                            ///< this array equals \f$ s_0, \ldots, s_k, d_0, \ldots, d_k \f$,
                                            ///< where scaling coefficients \f$ s_j = s_{jl}^n(f) \f$
                                            ///< and wavelet coefficients \f$ d_j = d_{jl}^n(f) \f$.
                                            ///< Here \f$ n, l \f$ are unique for every node.
    int n_coefs{0};                         ///< Number of coefficients in @ref coefs.

    int serialIx{-1};                       ///< Index in the serial tree
    int parentSerialIx{-1};                 ///< Index of the parent in the serial tree, or -1 for root nodes

    NodeIndex<D> nodeIndex;                 ///< Scale and translation of this node.
    HilbertPath<D> hilbertPath;             ///< Current Hilbert path state for child ordering.

    /**
     * @brief MWNode default constructor
     *
     * @details Should be used only by NodeAllocator to obtain
     *  virtual table pointers for the derived classes
     */
    MWNode();

    /**
     * @brief MWNode constructor
     * @param[in] tree The MWTree the root node belongs to
     * @param[in] rIdx The integer specifying the corresponding root node
     *
     * @details Constructor for root nodes. It requires the corresponding
     * MWTree and an integer to fetch the right NodeIndex.
     */
    MWNode(MWTree<D, T> *tree, int rIdx);

    /**
     * @brief MWNode constructor
     * @param[in] tree The MWTree the root node belongs to
     * @param[in] idx The NodeIndex defining scale and translation of the node
     *
     * @details Constructor for an empty node, given the corresponding MWTree and NodeIndex
     */
    MWNode(MWTree<D, T> *tree, const NodeIndex<D> &idx);

    /**
     * @brief MWNode constructor
     * @param[in] parent Parent node
     * @param[in] cIdx Child index of the current node
     *
     * @details Constructor for leaf nodes. It requires the corresponding
     * parent and an integer to identify the correct child.
     */
    MWNode(MWNode<D, T> *parent, int cIdx);

    // Implemented in child classes
    virtual void dealloc();

    /**
     * @brief Recurse down until an EndNode is found, and then crop children below the given precision threshold
     * @param prec The required precision
     * @param splitFac Factor used in the split check (larger factor means tighter threshold for finer nodes)
     * @param absPrec Flag to switch from relative (false) to absolute (true) precision.
     * @return Whether the crop was successful
     */
    bool crop(double prec, double splitFac, bool absPrec);

    /// @brief Initialize thread lock (when OpenMP is enabled).
    void initNodeLock() { MRCPP_INIT_OMP_LOCK(); }

    /**
     * @brief Allocate the coefs vector
     * @param n_blocks The number of blocks
     * @param block_size The size of a block
     *
     * @details This is only used by loose nodes, because the loose nodes
     * are not treated by the NodeAllocator class.
     */
    virtual void allocCoefs(int n_blocks, int block_size);

    /**
     * @brief Deallocate the coefs vector
     *
     * @details This is only used by loose nodes, because the loose nodes
     * are not treated by the NodeAllocator class.
     */
    virtual void freeCoefs();

    /**
     * @brief recursively set maxSquaredNorm and maxWSquareNorm of parent and descendants
     *
     * @details
     * normalization is such that a constant function gives constant value,
     * i.e. *not* same normalization as a squareNorm
     */
    void setMaxSquareNorm();

    /// @brief Recursively reset maxSquaredNorm and maxWSquareNorm of parent and descendants to value -1
    void resetMaxSquareNorm();

    /// @return The scaled square norm.
    double calcScaledSquareNorm() const { return std::pow(2.0, D * getScale()) * getSquareNorm(); }

    /// @return The scaled wavelet square norm.
    double calcScaledWSquareNorm() const { return std::pow(2.0, D * getScale()) * getWaveletNorm(); }

    /**
     * @brief Calculate the norm of one component (NOT the squared norm!)
     * @param i The component index
     * @return The single component norm
     */
    virtual double calcComponentNorm(int i) const;

    /**
     * @brief Update the coefficients of the node by a MW transform of the scaling
     * coefficients of the children.
     */
    virtual void reCompress();

    /**
     * @brief Forward MW transform from this node to its children
     * @param overwrite If true, the coefficients of the children are
     * overwritten. If false, the values are summed to the already present
     * ones.
     *
     * @details It performs forward MW transform inserting the result
     * directly in the right place for each child node. The children must
     * already be present and its memory allocated for this to work
     * properly.
     */
    virtual void giveChildrenCoefs(bool overwrite = true);

    /**
     * @brief Forward MW transform to compute scaling coefficients of a single child
     * @param[in] cIdx The child index
     * @param[in] overwrite If true, the coefficients of the children are
     * overwritten. If false, the values are summed to the already present
     * ones.
     *
     * @details It performs forward MW transform in place on a loose
     * node. The scaling coefficients of the selected child are then
     * copied/summed in the correct child node.
     */
    virtual void giveChildCoefs(int cIdx, bool overwrite = true);

    /** @brief Backward MW transform to compute scaling/wavelet coefficients of a parent
     *
     * @details Takes a MWParent and generates coefficients, reverse operation from
     * giveChildrenCoefs.
     *
     * @note This routine is only used in connection with Periodic Boundary Conditions
     */
    virtual void giveParentCoefs(bool overwrite = true);

    /**
     * @brief Copy scaling coefficients from children to parent
     *
     * @details Takes the scaling coefficients of the children and stores
     * them consecutively in the corresponding block of the parent,
     * following the usual bitwise notation.
     */
    virtual void copyCoefsFromChildren();

    /**
     * @brief Routine to find the path along the tree
     * @param[in] nIdx The sought after node through its NodeIndex
     *
     * @details Given the translation indices at the final scale, computes the child m
     * to be followed at the current scale in oder to get to the requested
     * node at the final scale. The result is the index of the child needed.
     * The index is obtained by bit manipulation of of the translation indices.
     */
    int getChildIndex(const NodeIndex<D> &nIdx) const;

    /**
     * @brief Routine to find the path along the tree
     * @param[in] r The sought after node through the coordinates of a point in space
     *
     * @details Given a point in space, determines which child should be followed
     * to get to the corresponding terminal node.
     */
    int getChildIndex(const Coord<D> &r) const;

    /**
     * @brief Fast check whether two nodes lie in different branches
     * @param rhs The node to compare against
     * @return true if two nodes lie in different branches
     */
    bool diffBranch(const MWNode<D, T> &rhs) const;

    /**
     * @brief Node retriever that ALWAYS returns the requested node
     *
     * @param[in] r The coordinates of a point in the node
     * @param depth The depth to descend
     * @return The node at the given coordinates
     *
     * @details Recursive routine to find and return the node with a given NodeIndex.
     * This routine always returns the appropriate node, and will generate nodes
     * that does not exist. Recursion starts at this node and ASSUMES the
     * requested node is in fact decending from this node.
     */
    MWNode<D, T> *retrieveNode(const Coord<D> &r, int depth);

    /**
     * @brief Node retriever that ALWAYS returns the requested node, possibly without coefs
     * @param[in] idx The NodeIndex of the requested node
     * @return The node at the given node index
     *
     * @details Recursive routine to find and return the node with a given NodeIndex. This
     * routine always returns the appropriate node, and will generate nodes that
     * does not exist. Recursion starts at this node and ASSUMES the requested
     * node is in fact descending from this node.
     * If create = true, the nodes are permanently added to the tree.
     */
    MWNode<D, T> *retrieveNode(const NodeIndex<D> &idx, bool create = false);

    /**
     * @brief Node retriever that ALWAYS returns the requested node
     * @param[in] idx The NodeIndex of the requested node
     * @return The node at the given node index
     *
     * @details Recursive routine to find and return the node with a given NodeIndex. This
     * routine always returns the appropriate node, and will generate nodes that
     * does not exist. Recursion starts at this node and ASSUMES the requested
     * node is in fact related to this node.
     *
     * @warning This routine is NOT thread safe! Must be used within omp critical.
     */
    MWNode<D, T> *retrieveParent(const NodeIndex<D> &idx);

    /**
     * @brief Const version of node retriever that NEVER generates
     * @param[in] idx The requested NodeIndex
     * @returns The requested node
     *
     * @details Recursive routine to find and return the node with a given NodeIndex.
     * This routine returns the appropriate Node, or a NULL pointer if
     * the node does not exist, or if it is a GenNode. Recursion starts at at this
     * node and ASSUMES the requested node is in fact decending from this node.
     */
    const MWNode<D, T> *retrieveNodeNoGen(const NodeIndex<D> &idx) const;

    /**
     * @brief Node retriever that NEVER generates.
     * @param[in] idx The requested NodeIndex
     * @returns The requested node
     *
     * @details Recursive routine to find and return the node with a given NodeIndex.
     * This routine returns the appropriate Node, or a NULL pointer if
     * the node does not exist, or if it is a GenNode. Recursion starts at at this
     * node and ASSUMES the requested node is in fact decending from this node.
     */
    MWNode<D, T> *retrieveNodeNoGen(const NodeIndex<D> &idx);

    /**
     * @brief Node retriever that returns requested Node or EndNode (const version)
     * @param[in] r The coordinates of a point in the node
     * @param depth The depth to descend
     * @return The node at the given coordinates
     *
     * @details Recursive routine to find and return the node given the
     * coordinates of a point in space.  This routine returns the
     * appropriate Node, or the EndNode on the path to the requested node,
     * and will never create or return GenNodes.  Recursion starts at at
     * this node and ASSUMES the requested node is in fact decending from
     * this node.
     */
    const MWNode<D, T> *retrieveNodeOrEndNode(const Coord<D> &r, int depth) const;

    /**
     * @brief Node retriever that returns requested Node or EndNode
     * @param[in] r The coordinates of a point in the node
     * @param depth The depth to descend
     * @return The node at the given coordinates
     *
     * @details Recursive routine to find and return the node given the
     * coordinates of a point in space.  This routine returns the
     * appropriate Node, or the EndNode on the path to the requested node,
     * and will never create or return GenNodes.  Recursion starts at at
     * this node and ASSUMES the requested node is in fact decending from
     * this node.
     */
    MWNode<D, T> *retrieveNodeOrEndNode(const Coord<D> &r, int depth);

    /**
     * @brief Node retriever that returns requested Node or EndNode (const version)
     * @param[in] idx The NodeIndex of the requested node
     * @return The requested node
     *
     * @details Recursive routine to find and return the node given the
     * coordinates of a point in space.  This routine returns the
     * appropriate Node, or the EndNode on the path to the requested node,
     * and will never create or return GenNodes.  Recursion starts at at
     * this node and ASSUMES the requested node is in fact decending from
     * this node.
     */
    const MWNode<D, T> *retrieveNodeOrEndNode(const NodeIndex<D> &idx) const;

    /**
     * @brief Node retriever that returns requested Node or EndNode
     * @param[in] idx The NodeIndex of the requested node
     * @return The requested node
     *
     * @details Recursive routine to find and return the node given the
     * coordinates of a point in space.  This routine returns the
     * appropriate Node, or the EndNode on the path to the requested node,
     * and will never create or return GenNodes.  Recursion starts at at
     * this node and ASSUMES the requested node is in fact decending from
     * this node.
     */
    MWNode<D, T> *retrieveNodeOrEndNode(const NodeIndex<D> &idx);

    /**
     * @brief Creates scaling coefficients of children
     *
     * @details If the node is a leaf node, it takes the scaling&wavelet
     * coefficients of the parent and it generates the scaling
     * coefficients for the children and stores
     * them consecutively in the corresponding block of the parent,
     * following the usual bitwise notation. The new node is permanently added to the tree.
     */
    void threadSafeCreateChildren();

    /**
     * @brief Generates scaling coefficients of children
     *
     * @details If the node is a leaf node, it takes the scaling&wavelet
     * coefficients of the parent and it generates the scaling
     * coefficients for the children and stores
     * them consecutively in the corresponding block of the parent,
     * following the usual bitwise notation.
     */
    void threadSafeGenChildren();


    /// @brief Deallocation of all generated nodes
    void deleteGenerated();

    /**
     * @brief Prints of the node content
     * @param[in,out] o The output stream
     */
    virtual std::ostream &print(std::ostream &o) const;

    // Bit flags describing node state
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