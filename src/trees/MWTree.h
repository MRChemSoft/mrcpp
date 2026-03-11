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
#include <map>
#include <memory>

#include "MRCPP/mrcpp_declarations.h"
#include "utils/omp_utils.h"

#include "MultiResolutionAnalysis.h"
#include "NodeAllocator.h"
#include "NodeBox.h"

namespace mrcpp {

class BankAccount;

/**
 * @class MWTree
 * @tparam D Spatial dimension (1, 2, or 3)
 * @tparam T Coefficient type (e.g. double, ComplexDouble)
 *
 * @brief Base class for Multiwavelet tree structures, such as FunctionTree and OperatorTree
 *
 * @details The MWTree class is the base class for all tree structures
 * needed for Multiwavelet calculations. The MWTree is a D-dimensional
 * tree structure of MWNodes. The tree starts from a set of root nodes
 * at a common given scale, defining the world box. The most common
 * settings are either a single root node or \f$ 2^D \f$ root
 * nodes. Other configurations are however allowed. For example, in 3D
 * one could have a total of 12 root nodes (a 2x2x3 set of root
 * nodes). Once the tree structure is generated, each node will have a
 * parent node (except for the root nodes) and \f$ 2^D \f$ child nodes
 * (except for leaf nodes). Most of the methods deal with traversing
 * the tree structure in different ways to fetch specific nodes. Some
 * of them will return a node present in the tree; some other methods
 * will generate the required node on the fly using the MW transform;
 * some methods will return an empty pointer if the node is not
 * present. See specific methods for details.
 */
template <int D, typename T> class MWTree {
public:
    /**
      * @brief MWTree constructor.
      *
      * @param[in] mra The multiresolution analysis object
      * @param[in] n The name of the tree (only for printing purposes)
      *
      * @details Creates an empty tree object, containing only the set of
      * root nodes. The information for the root node configuration to use
      * is in the mra object which is passed to the constructor.
      */
    MWTree(const MultiResolutionAnalysis<D> &mra, const std::string &n);

    MWTree(const MWTree<D, T> &tree) = delete;
    MWTree<D, T> &operator=(const MWTree<D, T> &tree) = delete;

    /// @brief MWTree destructor.
    virtual ~MWTree();

    /**
     * @brief Set the MW coefficients to zero, keeping the same tree structure
     * 
     * @details Keeps the node structure of the tree, even though the zero
     * function is representable at depth zero. One should then use \ref cropTree to remove
     * unnecessary nodes.
     */
    void setZero();

    /** @brief Remove all nodes in the tree
     *
     * @details Leaves the tree in the same state as after construction,
     * i.e. undefined tree structure containing only root nodes without
     * coefficients. The assigned memory, including branch and leaf
     * nodes, (nodeChunks in NodeAllocator) is NOT released, but is
     * immediately available to the new function.
     */
    void clear();

    double getSquareNorm() const { return this->squareNorm; } ///< @return The squared L2 norm of the function

    /** @brief Calculate the squared norm \f$ ||f||^2_{\ldots} \f$ of a function represented as a tree.
     *
     * @details The norm is calculated using endNodes only. The specific
     * type of norm which is computed will depend on the derived class.
     */
    void calcSquareNorm(bool deep = false);

    void clearSquareNorm() { this->squareNorm = -1.0; } //< @brief Mark the norm as undefined (sets it to -1)

    int getOrder() const { return this->order; }                            ///< @return Polynomial order k                
    int getKp1() const { return this->order + 1; }                          ///< @return k+1     
    int getKp1_d() const { return this->kp1_d; }                            ///< @return (k+1)^D
    int getDim() const { return D; }                                        ///< @return The spatial dimension D
    int getTDim() const { return (1 << D); }                                ///< @return 2^D (number of children per internal node)
    int getNNodes() const { return getNodeAllocator().getNNodes(); }        ///< @return The total number of nodes in this tree
    int getNNegScales() const { return this->nodesAtNegativeDepth.size(); } ///< @return The number of negative scales in this tree
    int getRootScale() const { return this->rootBox.getScale(); }           ///< @return The root scale of this tree
    int getDepth() const { return this->nodesAtDepth.size(); }              ///< @return The maximum depth of this tree
    int getSizeNodes() const;                                               ///< @return The size of all MW coefficients in the tree (in kB)
    /**
     * @brief Returns the total number of nodes in the tree, at given depth (not in use)
     * @param i Tree depth (0 depth is the coarsest scale) to count
     * @return Number of nodes at depth i
     */
    int getNNodesAtDepth(int i) const;

    NodeBox<D, T> &getRootBox() { return this->rootBox; }                   ///< @return The container of nodes
    const NodeBox<D, T> &getRootBox() const { return this->rootBox; }       ///< @return The container of nodes
    const MultiResolutionAnalysis<D> &getMRA() const { return this->MRA; }  ///< @return The MRA object used by this tree

    /** 
     * @brief Full Multiwavelet transform of the tree in either directions
     *
     * @param type TopDown (from roots to leaves) or BottomUp (from
     * leaves to roots) which specifies the direction of the MW transform
     * @param overwrite If true, the result will overwrite preexisting coefficients.
     *
     * @details It performs a Multiwavlet transform of the whole tree. The
     * input parameters will specify the direction (upwards or downwards)
     * and whether the result is added to the coefficients or it
     * overwrites them. See the documentation for the #mwTransformUp
     * and #mwTransformDown for details.
     * \f[
     * \pmatrix{
     * s_{nl}\\
     * d_{nl}
     * }
     * \rightleftarrows \pmatrix{
     * s_{n+1,2l}\\
     * s_{n+1,2l+1}
     * }
     * \f]
     */
    void mwTransform(int type, bool overwrite = true);

    /**
     * @brief Set the name of the tree
     * @param n The new name
     */
    void setName(const std::string &n) { this->name = n; }
    const std::string &getName() const { return this->name; }   ///< @return The name of the tree

    /**
     * @param r Spatial coordinates
     * @return The index of the root box containng r
     */
    int getRootIndex(Coord<D> r) const { return this->rootBox.getBoxIndex(r); }
    /**
     * @param nIdx Index of a node
     * @return The index of the root box containng nIdx
     */
    int getRootIndex(NodeIndex<D> nIdx) const { return this->rootBox.getBoxIndex(nIdx); }

    /**
     * @brief Finds and returns the node pointer with the given NodeIndex
     * @param nIdx The NodeIndex to search for
     * 
     * @details Recursive routine to find and return the node with a given
     * NodeIndex. This routine returns the appropriate Node, or a NULL
     * pointer if the node does not exist, or if it is a
     * GenNode. Recursion starts at the appropriate rootNode.
     * 
     * @return Pointer to the required node.
     */
    MWNode<D, T> *findNode(NodeIndex<D> nIdx);
    /**
     * @brief Finds and returns the node pointer with the given NodeIndex
     * @param nIdx The NodeIndex to search for
     * 
     * @details Recursive routine to find and return the node with a given
     * NodeIndex. This routine returns the appropriate Node, or a NULL
     * pointer if the node does not exist, or if it is a
     * GenNode. Recursion starts at the appropriate rootNode.
     * 
     * @return Pointer to the required node.
     */
    const MWNode<D, T> *findNode(NodeIndex<D> nIdx) const;

    /**
     * @brief Finds and returns the node reference with the given NodeIndex.
     * @param nIdx The NodeIndex to search for
     * @param create If true, previously non-existing nodes will be stored permanently in the tree
     * 
     * @details This routine ALWAYS returns the node you ask for. If the
     * node does not exist, it will be generated by MW
     * transform. Recursion starts at the appropriate rootNode and descends
     * from this.
     * 
     * @return Reference to the required node.
     * @note The nodes are permanently added to the tree if create = true.
     */
    MWNode<D, T> &getNode(NodeIndex<D> nIdx, bool create = false);

    /**
     * @brief Finds and returns the node (or EndNode) for the given NodeIndex.
     * @param nIdx The NodeIndex to search for
     *
     * @details This routine returns the Node you ask for, or the EndNode
     * on the path to the requested node, if the requested one is deeper
     * than the leaf node ancestor. It will never create or return
     * GenNodes.  Recursion starts at the appropriate rootNode and decends
     * from this.
     * 
     * @return Reference to the required node or EndNode.
     */
    MWNode<D, T> &getNodeOrEndNode(NodeIndex<D> nIdx);
    /**
     * @brief Finds and returns the node (or EndNode) for the given NodeIndex.
     * @param nIdx The NodeIndex to search for
     *
     * @details This routine returns the Node you ask for, or the EndNode
     * on the path to the requested node, if the requested one is deeper
     * than the leaf node ancestor. It will never create or return
     * GenNodes.  Recursion starts at the appropriate rootNode and decends
     * from this.
     * 
     * @return Reference to the required node or EndNode.
     */
    const MWNode<D, T> &getNodeOrEndNode(NodeIndex<D> nIdx) const;

    /** 
     * @brief Finds and returns the node at a given depth that contains a given coordinate.
     *
     * @param r Coordinates of an arbitrary point in space
     * @param depth Requested node depth from root scale
     *
     * @details This routine ALWAYS returns the node you ask for, and will
     * generate nodes that do not exist. Recursion starts at the
     * appropriate rootNode and decends from this.
     * 
     * @return Reference to the required node.
     */
    MWNode<D, T> &getNode(Coord<D> r, int depth = -1);

    /** 
     * @brief Finds and returns the node at a given depth that contains a given coordinate.
     *
     * @param r Coordinates of an arbitrary point in space
     * @param depth Requested node depth from root scale.
     *
     * @details This routine returns the Node you ask for, or the EndNode on
     * the path to the requested node, and will never create or return GenNodes.
     * Recursion starts at the appropriate rootNode and decends from this.
     * 
     * @return Reference to the required node or EndNode.
     */
    MWNode<D, T> &getNodeOrEndNode(Coord<D> r, int depth = -1);
    /** 
     * @brief Finds and returns the node at a given depth that contains a given coordinate.
     *
     * @param r Coordinates of an arbitrary point in space
     * @param depth Requested node depth from root scale.
     *
     * @details This routine returns the Node you ask for, or the EndNode on
     * the path to the requested node, and will never create or return GenNodes.
     * Recursion starts at the appropriate rootNode and decends from this.
     * 
     * @return Reference to the required node or EndNode.
     */
    const MWNode<D, T> &getNodeOrEndNode(Coord<D> r, int depth = -1) const;

    int getNEndNodes() const { return this->endNodeTable.size(); }  ///< @return The number of end nodes
    int getNRootNodes() const { return this->rootBox.size(); }      ///< @return The number of root nodes

    /**
     * @param i Index of the end node
     * @return Reference to the i-th end node
     */
    MWNode<D, T> &getEndMWNode(int i) { return *this->endNodeTable[i]; }
    /**
     * @param i Index of the root node
     * @return Reference to the i-th root node
     */
    MWNode<D, T> &getRootMWNode(int i) { return this->rootBox.getNode(i); }
    /**
     * @param i Index of the end node
     * @return Reference to the i-th end node
     */
    const MWNode<D, T> &getEndMWNode(int i) const { return *this->endNodeTable[i]; }
    /**
     * @param i Index of the root node
     * @return Reference to the i-th root node
     */
    const MWNode<D, T> &getRootMWNode(int i) const { return this->rootBox.getNode(i); }

    bool isPeriodic() const { return this->MRA.getWorldBox().isPeriodic(); } ///< @return Whether the world is periodic

    /** 
     * @brief Returns the list of all EndNodes
     * @details Copies the list of all EndNode pointers into a new vector and returns it.
     * @return The copied end-node table.
     */
    MWNodeVector<D, T> *copyEndNodeTable();
    MWNodeVector<D, T> *getEndNodeTable() { return &this->endNodeTable; } ///< @return The end-node table

    /** 
     * @brief Deletes all the nodes in the tree
     * @details This method will recursively delete all the nodes,
     * including the root nodes. Derived classes will call this method
     * when the object is deleted.
     */
    void deleteRootNodes();
    /** 
     * @brief Recreate the endNodeTable
     * @details the endNodeTable is first deleted and then rebuilt from
     * scratch. It makes use of the TreeIterator to traverse the tree.
     */
    void resetEndNodeTable();
    /// @brief Clear the end-node table
    void clearEndNodeTable() { this->endNodeTable.clear(); }

    //// @warning This method is currently not implemented.
    int countBranchNodes(int depth = -1);
    //// @warning This method is currently not implemented.
    int countLeafNodes(int depth = -1);
    //// @warning This method is currently not implemented.
    int countAllocNodes(int depth = -1);
    //// @warning This method is currently not implemented.
    int countNodes(int depth = -1);
    /// @brief Whether the tree coefficients are stored in the Bank
    bool isLocal = false; 

    /** 
     * @brief Gives serialIx of a node from its NodeIndex
     * @param nIdx The NodeIndex of the node
     *
     * @details Gives a unique integer for each nodes corresponding to the position
     * of the node in the serialized representation. Only works if isLocal == true.
     * 
     * @return The serial index of the node.
     */
    int getIx(NodeIndex<D> nIdx);

    /** 
     * @brief Sets values for maxSquareNorm and maxWSquaredNorm in all nodes
     * @details It defines the upper bound of the squared norm \f$
     * ||f||^2_{\ldots} \f$ in this node or its descendents.
     */
    void makeMaxSquareNorms();

    NodeAllocator<D, T> &getNodeAllocator() { return *this->nodeAllocator_p; }              ///< @return Reference to the node allocator.
    const NodeAllocator<D, T> &getNodeAllocator() const { return *this->nodeAllocator_p; }  ///< @return Reference to the node allocator.

    MWNodeVector<D, T> endNodeTable; ///< @brief Final projected nodes

    /**
     * @brief Fetch coefficients of a specific node stored in Bank
     * @param nIdx Node index
     * @param[out] data The node coefficients are copied into this array
     */
    void getNodeCoeff(NodeIndex<D> nIdx, T *data);

    bool conjugate() const { return this->conj; }           ///< @return Whether the tree is conjugated
    void setConjugate(bool conjug) { this->conj = conjug; } ///< @param conjug Set whether the tree is conjugated
    
    /**
     * @brief Print a formatted summary of the tree
     * @param o The output stream
     * @param tree The tree to print
     * @return The output stream
     */
    friend std::ostream &operator<<(std::ostream &o, const MWTree<D, T> &tree) { return tree.print(o); }

    // Friend classes
    friend class MWNode<D, T>;
    friend class FunctionNode<D, T>;
    friend class OperatorNode;
    friend class TreeBuilder<D, T>;
    friend class NodeAllocator<D, T>;

protected:
    // Parameters that are set in construction and should never change
    const MultiResolutionAnalysis<D> MRA; ///< Domain and basis

    // Constant parameters that are derived internally
    const int order;   ///< Polynomial order k
    const int kp1_d;   ///< (k+1)^D

    std::map<NodeIndex<D>, int> NodeIndex2serialIx; ///< To store nodes serialIx 

    // Parameters that are dynamic and can be set by user
    std::string name; ///< Name of this tree
    std::unique_ptr<NodeAllocator<D, T>> nodeAllocator_p{nullptr}; ///< Node allocator

    // Tree data
    double squareNorm;                      ///< Global squared L2 norm (-1 if undefined).
    NodeBox<D, T> rootBox;                  ///< The actual container of nodes
    std::vector<int> nodesAtDepth;          ///< Node counter
    std::vector<int> nodesAtNegativeDepth;  ///< Node counter

    /** 
     * @brief Regenerates all scaling coeffs by MW transformation of existing s/w-coeffs
     * on coarser scales
     * @param overwrite If true, the preexisting coefficients are overwritten
     *
     * @details The transformation starts at the rootNodes and proceeds
     * recursively all the way to the leaf nodes. The existing scaling
     * coefficeints will either be overwritten or added to. The latter
     * operation is generally used after the operator application.
     */
    virtual void mwTransformDown(bool overwrite);

    /** 
     * @brief Regenerates all s/d-coeffs by backtransformation
     *
     * @details It starts at the bottom of the tree (scaling coefficients
     * of the leaf nodes) and it generates the scaling and wavelet
     * coefficients of the parent node. It then proceeds recursively all the
     * way up to the root nodes. This is generally used after a function
     * projection to purify the coefficients obtained by quadrature at
     * coarser scales which are therefore not precise enough.
     */
    virtual void mwTransformUp();

    /** 
     * @brief Increments node counter by one for non-GenNodes
     * @param scale Scale of the node
     * @warning: This routine is not thread safe, and must NEVER be called
     * outside a critical region in parallel. It's way, way too expensive to
     * lock the tree, so don't even think about it.
     */
    void incrementNodeCount(int scale);

    /** 
     * @brief Decrements node counter by one for non-GenNodes
     * @param scale Scale of the node
     * @warning: This routine is not thread safe, and must NEVER be called
     * outside a critical region in parallel. It's way, way too expensive to
     * lock the tree, so don't even think about it.
     */
    void decrementNodeCount(int scale);

    BankAccount *NodesCoeff = nullptr; ///< Bank account for node coefficients

    bool conj{false}; ///< Whether the tree is conjugated

    /** 
     * @brief Prints a summary of the tree structure on the output file
     * @param o The output stream
     * @return The formatted output stream
     */
    virtual std::ostream &print(std::ostream &o) const;
};

} // namespace mrcpp