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
 * @file MWTree.h
 * @brief Base template for multiwavelet (MW) tree data structures.
 *
 * @details
 * An MW tree stores a hierarchical collection of @ref MWNode "MWNode<D,T>"
 * objects arranged as a 2^D-ary tree over a @ref MultiResolutionAnalysis
 * (computational domain + basis). It provides:
 *   - ownership and construction of the root nodes (via @ref NodeBox),
 *   - navigation and on-demand generation of nodes,
 *   - bookkeeping of per-depth node counts,
 *   - utilities for MW transforms, norms, and end-node tables, and
 *   - access to a @ref NodeAllocator for memory management.
 *
 * This class is a base for both **function** and **operator** trees.
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
 * @tparam D Spatial dimension (1, 2, or 3).
 * @tparam T Coefficient scalar type (e.g. double, ComplexDouble).
 *
 * @brief Base class for MW tree structures (e.g., FunctionTree, OperatorTree).
 *
 * @details
 * A tree is defined over a @ref MultiResolutionAnalysis (MRA). The set of root
 * nodes is determined by the MRA world box. Each node has up to 2^D children.
 * Some accessors only *find* existing nodes, while others may *create* nodes
 * lazily (e.g., by splitting and transferring coefficients).
 *
 * ### Node retrieval semantics
 * - @ref findNode returns a pointer or `nullptr` if the node is missing.
 * - @ref getNode and @ref getNodeOrEndNode return a reference and can
 *   create intermediate nodes if requested (see parameters).
 *
 * ### Norms
 * @ref calcSquareNorm computes the global L2 norm (squared) either from
 * the existing coefficients only, or by visiting descendants when `deep=true`.
 */
template <int D, typename T> class MWTree {
public:
    /**
     * @brief Construct an empty tree bound to an MRA.
     * @param mra Multi-resolution analysis (domain + basis).
     * @param n   A short name for logging/printing.
     *
     * @post Root nodes are created according to the MRA world box.
     *       No coefficients are computed; the tree is “undefined”.
     */
    MWTree(const MultiResolutionAnalysis<D> &mra, const std::string &n);

    /// Non-copyable.
    MWTree(const MWTree<D, T> &tree) = delete;
    /// Non-assignable.
    MWTree<D, T> &operator=(const MWTree<D, T> &tree) = delete;

    /// Virtual destructor.
    virtual ~MWTree();

    /**
     * @brief Set all existing node coefficients to zero (structure unchanged).
     * @note Does not refine/coarsen the tree, only zeroes values.
     */
    void setZero();

    /**
     * @brief Remove all nodes and reset to a freshly constructed state.
     * @post Root nodes are recreated; end-node table and counters cleared.
     */
    void clear();

    /** @name Norms */
    ///@{
    /// @return Global squared L2 norm of the representation (negative if undefined).
    double getSquareNorm() const { return this->squareNorm; }

    /**
     * @brief Recompute the global squared L2 norm.
     * @param deep If `true`, may traverse deeper to ensure accuracy.
     */
    void calcSquareNorm(bool deep = false);

    /// @brief Mark the norm as undefined (sets it to -1).
    void clearSquareNorm() { this->squareNorm = -1.0; }
    ///@}

    /** @name Basis/structure parameters */
    ///@{
    /// @return Polynomial order k.
    int getOrder() const { return this->order; }
    /// @return k+1.
    int getKp1() const { return this->order + 1; }
    /// @return (k+1)^D.
    int getKp1_d() const { return this->kp1_d; }
    /// @return Spatial dimension D.
    int getDim() const { return D; }
    /// @return 2^D (number of children per internal node).
    int getTDim() const { return (1 << D); }
    /// @return Total number of nodes currently allocated in the tree.
    int getNNodes() const { return getNodeAllocator().getNNodes(); }
    /// @return Number of records kept for negative-depth counts.
    int getNNegScales() const { return this->nodesAtNegativeDepth.size(); }
    /// @return Root scale (MRA world scale).
    int getRootScale() const { return this->rootBox.getScale(); }
    /// @return Number of depth levels for which counters are stored.
    int getDepth() const { return this->nodesAtDepth.size(); }
    /// @return Number of nodes counted at a given depth.
    int getNNodesAtDepth(int i) const;
    /// @return Approximate memory footprint of nodes (kB).
    int getSizeNodes() const;
    ///@}

    /** @name MRA / root access */
    ///@{
    /// @return Mutable root-node container.
    NodeBox<D, T> &getRootBox() { return this->rootBox; }
    /// @return Const root-node container.
    const NodeBox<D, T> &getRootBox() const { return this->rootBox; }
    /// @return MRA bound to this tree.
    const MultiResolutionAnalysis<D> &getMRA() const { return this->MRA; }
    ///@}

    /**
     * @brief Perform a multiresolution transform.
     * @param type Transform kind (implementation-defined selector).
     * @param overwrite If `true`, may reuse buffers for speed.
     * @note Typical directions are “top-down” and “bottom-up”; see implementation.
     */
    void mwTransform(int type, bool overwrite = true);

    /** @name Naming */
    ///@{
    /// Set a short descriptive name (used in logs).
    void setName(const std::string &n) { this->name = n; }
    /// Get the current name.
    const std::string &getName() const { return this->name; }
    ///@}

    /** @name Root-index helpers */
    ///@{
    /// @return Root-box index containing coordinate @p r, or -1 if out-of-bounds (non-periodic).
    int getRootIndex(Coord<D> r) const { return this->rootBox.getBoxIndex(r); }
    /// @return Root-box index containing node @p nIdx, or -1 if out-of-bounds (non-periodic).
    int getRootIndex(NodeIndex<D> nIdx) const { return this->rootBox.getBoxIndex(nIdx); }
    ///@}

    /** @name Node lookup / retrieval */
    ///@{
    /**
     * @brief Find an existing node.
     * @param nIdx Target node index.
     * @return Pointer to the node if present, otherwise `nullptr`.
     * @warning Does not create nodes.
     */
    MWNode<D, T> *findNode(NodeIndex<D> nIdx);

    /// Const overload of @ref findNode.
    const MWNode<D, T> *findNode(NodeIndex<D> nIdx) const;

    /**
     * @brief Get a node by index, optionally creating it.
     * @param nIdx Target node index.
     * @param create If `true`, missing nodes may be generated on demand.
     * @return Reference to the node.
     */
    MWNode<D, T> &getNode(NodeIndex<D> nIdx, bool create = false);

    /**
     * @brief Get a node or the “closest” end node containing it.
     * @param nIdx Target node index.
     * @return Reference to an existing node; may be an end node if exact match is absent.
     * @note Never creates new nodes.
     */
    MWNode<D, T> &getNodeOrEndNode(NodeIndex<D> nIdx);

    /// Const overload of @ref getNodeOrEndNode(NodeIndex).
    const MWNode<D, T> &getNodeOrEndNode(NodeIndex<D> nIdx) const;

    /**
     * @brief Get a node by spatial coordinate.
     * @param r     Spatial coordinate.
     * @param depth Desired depth; if negative, use current deepest.
     * @return Reference to the node; may create if required by the implementation.
     */
    MWNode<D, T> &getNode(Coord<D> r, int depth = -1);

    /**
     * @brief Get a node or containing end node by coordinate.
     * @param r     Spatial coordinate.
     * @param depth Desired depth; if negative, use current deepest.
     * @return Reference to an existing node; may be an end node.
     */
    MWNode<D, T> &getNodeOrEndNode(Coord<D> r, int depth = -1);

    /// Const overload of @ref getNodeOrEndNode(Coord,int).
    const MWNode<D, T> &getNodeOrEndNode(Coord<D> r, int depth = -1) const;
    ///@}

    /** @name End-node table */
    ///@{
    /// @return Number of nodes currently listed as “end nodes”.
    int getNEndNodes() const { return this->endNodeTable.size(); }
    /// @return Number of root nodes.
    int getNRootNodes() const { return this->rootBox.size(); }

    /// @return Mutable reference to i-th end node.
    MWNode<D, T> &getEndMWNode(int i) { return *this->endNodeTable[i]; }
    /// @return Mutable reference to i-th root node.
    MWNode<D, T> &getRootMWNode(int i) { return this->rootBox.getNode(i); }

    /// @return Const reference to i-th end node.
    const MWNode<D, T> &getEndMWNode(int i) const { return *this->endNodeTable[i]; }
    /// @return Const reference to i-th root node.
    const MWNode<D, T> &getRootMWNode(int i) const { return this->rootBox.getNode(i); }
    ///@}

    /// @return `true` if the underlying world box has any periodic directions.
    bool isPeriodic() const { return this->MRA.getWorldBox().isPeriodic(); }

    /**
     * @brief Copy the current end-node table.
     * @return New heap-allocated vector; caller takes ownership.
     */
    MWNodeVector<D, T> *copyEndNodeTable();

    /// @return Direct pointer to the internal end-node table.
    MWNodeVector<D, T> *getEndNodeTable() { return &this->endNodeTable; }

    /** @name Tree maintenance */
    ///@{
    /// Delete all root nodes and reset root structures.
    void deleteRootNodes();
    /// Rebuild the end-node table by traversing the tree.
    void resetEndNodeTable();
    /// Clear the end-node table without traversing.
    void clearEndNodeTable() { this->endNodeTable.clear(); }
    ///@}

    /** @name Node statistics (current tree) */
    ///@{
    /// Count branch (non-leaf) nodes; if depth < 0, count all depths.
    int countBranchNodes(int depth = -1);
    /// Count leaf nodes; if depth < 0, count all depths.
    int countLeafNodes(int depth = -1);
    /// Count allocated nodes; if depth < 0, count all depths.
    int countAllocNodes(int depth = -1);
    /// Count nodes; if depth < 0, count all depths.
    int countNodes(int depth = -1);
    ///@}

    /// If `true`, coefficients are stored externally (Bank); used by serialization tools.
    bool isLocal = false;

    /**
     * @brief Map a node index to its serial index (when stored locally).
     * @param nIdx Node index.
     * @return Serial index, or a negative value if not present.
     */
    int getIx(NodeIndex<D> nIdx);

    /**
     * @brief Precompute per-node maxima used by some adaptive algorithms.
     * @details Fills `maxSquareNorm` and `maxWSquareNorm` for all nodes.
     */
    void makeMaxSquareNorms();

    /** @name Allocator access */
    ///@{
    /// @return Mutable reference to the node allocator.
    NodeAllocator<D, T> &getNodeAllocator() { return *this->nodeAllocator_p; }
    /// @return Const reference to the node allocator.
    const NodeAllocator<D, T> &getNodeAllocator() const { return *this->nodeAllocator_p; }
    ///@}

    /// Vector of final projected nodes (end nodes).
    MWNodeVector<D, T> endNodeTable;

    /**
     * @brief Fetch coefficients of a specific node (when using a Bank).
     * @param nIdx Node index.
     * @param data Destination buffer of size (k+1)^D * 2^D.
     */
    void getNodeCoeff(NodeIndex<D> nIdx, T *data);

    /// @return Whether the tree is marked as conjugated (used by some ops).
    bool conjugate() const { return this->conj; }
    /// Set or clear the conjugation flag.
    void setConjugate(bool conjug) { this->conj = conjug; }

    /// Print tree summary to a stream.
    friend std::ostream &operator<<(std::ostream &o, const MWTree<D, T> &tree) { return tree.print(o); }

    // Friends that require access to internals
    friend class MWNode<D, T>;
    friend class FunctionNode<D, T>;
    friend class OperatorNode;
    friend class TreeBuilder<D, T>;
    friend class NodeAllocator<D, T>;

protected:
    /** @name Immutable construction-time state */
    ///@{
    const MultiResolutionAnalysis<D> MRA; ///< Domain and basis.

    const int order;   ///< Polynomial order k.
    const int kp1_d;   ///< (k+1)^D.
    ///@}

    /// Map node index -> serial index (used by local/banked storage).
    std::map<NodeIndex<D>, int> NodeIndex2serialIx;

    /** @name User-settable metadata */
    ///@{
    std::string name; ///< Short name for diagnostics.
    ///@}

    /// Node memory allocator.
    std::unique_ptr<NodeAllocator<D, T>> nodeAllocator_p{nullptr};

    /** @name Tree data & counters */
    ///@{
    double squareNorm;                 ///< Global squared L2 norm (-1 if undefined).
    NodeBox<D, T> rootBox;             ///< Container of root nodes.
    std::vector<int> nodesAtDepth;     ///< Per-depth node counts (depth >= 0).
    std::vector<int> nodesAtNegativeDepth; ///< For negative-depth bookkeeping.
    ///@}

    /** @name MW transforms (internals) */
    ///@{
    virtual void mwTransformDown(bool overwrite);
    virtual void mwTransformUp();
    ///@}

    /// Increment per-depth counters for a node at the given scale.
    void incrementNodeCount(int scale);
    /// Decrement per-depth counters for a node at the given scale.
    void decrementNodeCount(int scale);

    /// Optional external storage of coefficients.
    BankAccount *NodesCoeff = nullptr;

    /// Conjugation flag for algorithms that need it.
    bool conj{false};

    /// Print a formatted summary (override in derived classes if needed).
    virtual std::ostream &print(std::ostream &o) const;
};

} // namespace mrcpp