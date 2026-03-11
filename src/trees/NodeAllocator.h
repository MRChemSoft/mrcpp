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

#include <vector>

#include "MRCPP/mrcpp_declarations.h"
#include "utils/omp_utils.h"

namespace mrcpp {

/**
 * @class NodeAllocator
 * @tparam D Spatial dimension (1, 2, or 3)
 * @tparam T Coefficient type (e.g. double, ComplexDouble)
 *
 * @brief Chunked memory manager for @ref MWNode objects and their coefficients
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
     * @brief Construct an allocator bound to a function tree
     * @param[in] tree           Owning @ref FunctionTree instance
     * @param[in] mem            Optional shared-memory provider for coefficients (may be `nullptr`)
     * @param[in] coefsPerNode   Number of coefficients per node
     * @param[in] nodesPerChunk  Maximum number of nodes per chunk
     * 
     * @details Reserves space for chunk pointers to avoid excessive reallocation,
     * but does not allocate any chunks until needed.
     */
    NodeAllocator(FunctionTree<D, T> *tree, SharedMemory<T> *mem, int coefsPerNode, int nodesPerChunk);


    /**
     * @brief Construct an allocator bound to an operator tree
     * @param[in] tree           Owning @ref OperatorTree instance
     * @param[in] mem            Optional shared-memory provider for coefficients (may be `nullptr`)
     * @param[in] coefsPerNode   Number of coefficients per node
     * @param[in] nodesPerChunk  Maximum number of nodes per chunk
     * 
     * @details Reserves space for chunk pointers to avoid excessive reallocation,
     * but does not allocate any chunks until needed.
     */
    NodeAllocator(OperatorTree *tree, SharedMemory<T> *mem, int coefsPerNode, int nodesPerChunk);


    /// Non-copyable.
    NodeAllocator(const NodeAllocator<D, T> &tree) = delete;
    /// Non-assignable.
    NodeAllocator<D, T> &operator=(const NodeAllocator<D, T> &tree) = delete;

    /// Destructor; releases all owned chunks (nodes and coefficients).
    ~NodeAllocator();


    /**
     * @brief Get pointer to a node object by serial index
     * @param[in] sIdx Serial index of the node
     * @return Pointer to the @ref MWNode instance.
     */
    MWNode<D, T> *getNode_p(int sIdx);


    /**
     * @brief Get pointer to the coefficient array for a node
     * @param[in] sIdx Serial index of the node
     * @return Pointer to 'T[coefsPerNode]' or 'nullptr' if unavailable.
     */
    T *getCoef_p(int sIdx);



    /**
     * @brief Allocate a consecutive block of nodes
     * @param[in] nNodes Number of nodes to allocate
     * @param[in] coefs  If 'true', also ensure coefficient storage is available
     * @return Serial index ('sIdx') of the first newly allocated node (the top of the stack)
     * 
     * @details Allocates a block of @p nNodes consecutive nodes, returning 
     * the serial index of the first node in the block. If `coefs` is true,
     * coefficient arrays are also allocated for each node. If there is not
     * enough space in existing chunks, new chunks are allocated to satisfy
     * the request
     * 
     * @warning Does not initialize the node objects; caller is responsible
     * for placement-new or similar
     * 
     * @warning If insufficient space is available, and allocation of new
     * chunks fails, an exception is thrown and no nodes are allocated.
     * 
     * @throw std::bad_alloc if memory allocation fails.
     * @throw std::runtime_error if insufficient space is available after
     *         attempting to allocate new chunks.
     * @note May grow the underlying chunk arrays if space is exhausted.
     */
    int alloc(int nNodes, bool coefs = true);

    /**
     * @brief Deallocate a node at serial index
     * @param[in] sIdx Serial index of the node to free
     * @details Marks the node at serial index @p sIdx as free for future
     * allocations. Does not destroy the node object or its coefficient array.
     * It also updates the number of allocated nodes.
     *
     * @warning Does not shrink chunks; it only marks the slot as free.
     *
     * @throw std::out_of_range if @p sIdx is invalid.
     */
    void dealloc(int sIdx);

    /**
     * @brief Deallocate coefficient arrays for all nodes
     * @note Node objects remain allocated; only their coefficient buffers are freed.
     */
    void deallocAllCoeff();

    /**
     * @brief Pre-allocate a number of chunks
     * @param[in] nChunks Number of chunks to append
     * @param[in] coefs   If 'true', allocate coefficient chunks as well
     *
     * @details It reinitializes the allocator, allocating @p nChunks chunks
     * (both nodes and coefficients, if @p coefs is true). It resized the
     * stackStatus vectors with the new total capacity, and resets the
     * allocation stack.
     * 
     * @note This method clears any previously allocated nodes and
     * their coefficient buffers.
     *
     * @throw If nChunks <= 0
     */
    void init(int nChunks, bool coefs = true);

    /**
     * @brief Fill all holes in the chunks with occupied nodes, then remove all empty chunks
     * @return Number of nodes deleted during compaction
     *
     * @details After compaction, serial indices may change internally; users
     * should refresh any external mappings that depend on 'sIdx'.
     */
    int compress();

     /**
     * @brief Drop trailing unused chunks to release memory.
     * @return Number of chunks deleted
     * 
     * @details Scans chunks from the end towards the beginning, deleting any
     * chunks that are completely unused. Stops when a chunk with at least
     * one occupied node is found.
     */
    int deleteUnusedChunks();

    /**
     * @brief Traverse tree and redefine pointer, counter and tables
     * @details Typically invoked after operations that reorder nodes without
     * using @ref compress 
     */
    void reassemble();



    
    int getNNodes() const { return this->nNodes; } ///< @return Number of nodes currently in use (allocated and not freed).
    int getNCoefs() const { return this->coefsPerNode; } ///< @return Number of coefficients per node.
    int getNChunks() const { return this->nodeChunks.size(); } ///< @return Total number of allocated node chunks.
    int getNChunksUsed() const { return (this->topStack + this->maxNodesPerChunk - 1) / this->maxNodesPerChunk; } ///< @return Number of chunks currently used by active nodes.
    int getNodeChunkSize() const { return this->maxNodesPerChunk * this->sizeOfNode; } ///< @return Size in bytes of one node chunk (nodes only).
    int getCoefChunkSize() const { return this->maxNodesPerChunk * this->coefsPerNode * sizeof(T); } ///< @return Size in bytes of one coefficient chunk.
    int getMaxNodesPerChunk() const { return this->maxNodesPerChunk; } ///< @return Maximum number of nodes that fit in a single chunk.



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
