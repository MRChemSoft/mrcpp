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
 * @file NodeAllocator.h
 * @brief Chunked allocator for MWNode objects and their coefficient storage.
 *
 * @details
 * The allocator handles:
 *  - contiguous chunk allocation for **nodes** and **coefficients**,
 *  - a simple stack-like free list for fast allocation/deallocation,
 *  - optional backing via a shared memory block (@ref SharedMemory),
 *  - utility routines for compaction (@ref compress) and reassembly
 *    after structural edits, and
 *  - query helpers for chunk sizes and usage.
 *
 * It is used by both @ref FunctionTree and @ref OperatorTree.
 */

#pragma once

#include <vector>

#include "MRCPP/mrcpp_declarations.h"
#include "utils/omp_utils.h"

namespace mrcpp {

/**
 * @class NodeAllocator
 * @tparam D Spatial dimension (1, 2, or 3).
 * @tparam T Scalar coefficient type (e.g., double, ComplexDouble).
 *
 * @brief Chunked memory manager for @ref MWNode objects and their coefficients.
 *
 * @details
 * Nodes and their coefficient arrays are organized in **chunks** to reduce
 * allocation overhead and improve spatial locality. Indices into this storage
 * are referred to as *serial indices* (`sIdx`), which are stable within a
 * given tree instance until compaction or reassembly occurs.
 *
 * ### Thread-safety
 * When MRCPP is built with OpenMP support (`MRCPP_HAS_OMP`), critical regions
 * in the allocator use locks to avoid races during allocation and pointer
 * retrieval. Callers are still responsible for higher-level synchronization
 * of tree edits.
 */
template <int D, typename T> class NodeAllocator final {
public:
    /**
     * @brief Construct an allocator bound to an operator tree.
     * @param tree           Owning @ref OperatorTree instance.
     * @param mem            Optional shared-memory provider for coefficients (may be `nullptr`).
     * @param coefsPerNode   Number of coefficients per node.
     * @param nodesPerChunk  Maximum number of nodes per chunk.
     */
    NodeAllocator(OperatorTree *tree, SharedMemory<T> *mem, int coefsPerNode, int nodesPerChunk);

    /**
     * @brief Construct an allocator bound to a function tree.
     * @param tree           Owning @ref FunctionTree instance.
     * @param mem            Optional shared-memory provider for coefficients (may be `nullptr`).
     * @param coefsPerNode   Number of coefficients per node.
     * @param nodesPerChunk  Maximum number of nodes per chunk.
     */
    NodeAllocator(FunctionTree<D, T> *tree, SharedMemory<T> *mem, int coefsPerNode, int nodesPerChunk);

    /// Non-copyable.
    NodeAllocator(const NodeAllocator<D, T> &tree) = delete;
    /// Non-assignable.
    NodeAllocator<D, T> &operator=(const NodeAllocator<D, T> &tree) = delete;

    /// Destructor; releases all owned chunks (nodes and coefficients).
    ~NodeAllocator();

    /**
     * @brief Allocate a consecutive block of nodes.
     * @param nNodes Number of nodes to allocate.
     * @param coefs  If `true`, also ensure coefficient storage is available.
     * @return Serial index (`sIdx`) of the first newly allocated node.
     *
     * @note May grow the underlying chunk arrays if space is exhausted.
     */
    int alloc(int nNodes, bool coefs = true);

    /**
     * @brief Deallocate a node at serial index.
     * @param sIdx Serial index of the node to free.
     *
     * @warning Does not shrink chunks; it only marks the slot as free.
     */
    void dealloc(int sIdx);

    /**
     * @brief Deallocate coefficient arrays for all nodes.
     * @note Node objects remain allocated; only their coefficient buffers are freed.
     */
    void deallocAllCoeff();

    /**
     * @brief Pre-allocate a number of chunks.
     * @param nChunks Number of chunks to append.
     * @param coefs   If `true`, allocate coefficient chunks as well.
     *
     * @details Useful to avoid repeated growth when the final size is known.
     */
    void init(int nChunks, bool coefs = true);

    /**
     * @brief Compact allocated nodes to reduce fragmentation.
     * @return Number of nodes moved during compaction.
     *
     * @details After compaction, serial indices may change internally; users
     * should refresh any external mappings that depend on `sIdx`.
     */
    int compress();

    /**
     * @brief Rebuild internal pointers after external moves/shuffling.
     * @details Typically invoked after operations that reorder nodes without
     * using @ref compress.
     */
    void reassemble();

    /**
     * @brief Drop trailing unused chunks to release memory.
     * @return Number of chunks deleted.
     */
    int deleteUnusedChunks();

    /** @name Introspection */
    ///@{
    /// @return Number of nodes currently in use (allocated and not freed).
    int getNNodes() const { return this->nNodes; }
    /// @return Number of coefficients per node.
    int getNCoefs() const { return this->coefsPerNode; }
    /// @return Total number of allocated node chunks.
    int getNChunks() const { return this->nodeChunks.size(); }
    /// @return Number of chunks currently used by active nodes.
    int getNChunksUsed() const { return (this->topStack + this->maxNodesPerChunk - 1) / this->maxNodesPerChunk; }
    /// @return Size in bytes of one node chunk (nodes only).
    int getNodeChunkSize() const { return this->maxNodesPerChunk * this->sizeOfNode; }
    /// @return Size in bytes of one coefficient chunk.
    int getCoefChunkSize() const { return this->maxNodesPerChunk * this->coefsPerNode * sizeof(T); }
    /// @return Maximum number of nodes that fit in a single chunk.
    int getMaxNodesPerChunk() const { return this->maxNodesPerChunk; }
    ///@}

    /**
     * @brief Get pointer to the coefficient array for a node.
     * @param sIdx Serial index of the node.
     * @return Pointer to `T[coefsPerNode]` or `nullptr` if unavailable.
     */
    T *getCoef_p(int sIdx);

    /**
     * @brief Get pointer to a node object by serial index.
     * @param sIdx Serial index of the node.
     * @return Pointer to the @ref MWNode instance.
     */
    MWNode<D, T> *getNode_p(int sIdx);

    /// @return Pointer to the i-th coefficient chunk (contiguous block).
    T *getCoefChunk(int i) { return this->coefChunks[i]; }
    /// @return Pointer to the i-th node chunk (contiguous block).
    MWNode<D, T> *getNodeChunk(int i) { return this->nodeChunks[i]; }

    /// Print allocator status (chunks, usage, sizes) to stdout.
    void print() const;

protected:
    int nNodes{0};           ///< Number of nodes actually in use.
    int topStack{0};         ///< Index of the next free slot (stack top).
    int sizeOfNode{0};       ///< `sizeof(NodeType)` used in chunks.
    int coefsPerNode{0};     ///< Number of coefficients per node.
    int maxNodesPerChunk{0}; ///< Capacity (in nodes) of each chunk.

    std::vector<int> stackStatus{};          ///< Slot state (occupied/free).
    std::vector<T *> coefChunks{};           ///< Coefficient chunk base pointers.
    std::vector<MWNode<D, T> *> nodeChunks{};///< Node chunk base pointers.

    char *cvptr{nullptr};              ///< Vtable cookie to initialize node objects.
    MWNode<D, T> *last_p{nullptr};     ///< Pointer just past the last active node.
    MWTree<D, T> *tree_p{nullptr};     ///< Back-pointer to owning tree.
    SharedMemory<T> *shmem_p{nullptr}; ///< Optional shared-memory backend.

    /// @return Whether coefficients are backed by @ref SharedMemory.
    bool isShared() const { return (this->shmem_p != nullptr); }
    /// @return Owning tree (non-const).
    MWTree<D, T> &getTree() { return *this->tree_p; }
    /// @return Shared-memory provider (non-const).
    SharedMemory<T> &getMemory() { return *this->shmem_p; }

    /// Internal: get coefficient pointer w/o locking (caller must synchronize).
    T *getCoefNoLock(int sIdx);
    /// Internal: get node pointer w/o locking (caller must synchronize).
    MWNode<D, T> *getNodeNoLock(int sIdx);

    /// Move a block of nodes within chunks (used by @ref compress).
    void moveNodes(int nNodes, int srcIdx, int dstIdx);
    /// Append a new chunk; if `coefs` is true, also append a coefficient chunk.
    void appendChunk(bool coefs);

    /// Find next contiguous range of free slots starting at or after `sIdx`.
    int findNextAvailable(int sIdx, int nNodes) const;
    /// Find next occupied slot at or after `sIdx`.
    int findNextOccupied(int sIdx) const;

#ifdef MRCPP_HAS_OMP
    omp_lock_t omp_lock; ///< OpenMP lock for critical sections.
#endif
};

} // namespace mrcpp
