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
/**
 * @file
 * @brief MPI-facing declarations and a lightweight shared-memory helper for MRCPP.
 *
 * This header provides:
 * - Portable aliases for MPI types (working in non-MPI builds as no-ops).
 * - Public globals describing the current MPI topology/roles used by MRCPP.
 * - A templated @ref mrcpp::SharedMemory class to allocate a shared-memory
 *   window among ranks that share a physical node (MPI-3 RMA).
 * - Prototypes for shipping trees between ranks: @ref send_tree, @ref recv_tree, @ref share_tree.
 *
 * @note All MPI symbols are guarded by `MRCPP_HAS_MPI`. In non-MPI builds, dummy
 *       typedefs are supplied so that client code can still compile.
 */

#ifdef MRCPP_HAS_MPI
  #include <mpi.h>
#else
  /// Fallback alias so non-MPI builds can compile client code.
  using MPI_Comm    = int;
  /// Fallback alias so non-MPI builds can compile client code.
  using MPI_Win     = int;
  /// Fallback alias so non-MPI builds can compile client code.
  using MPI_Request = int;
#endif

namespace mrcpp {

/// Alias for MPI communicator (portable across MPI/non-MPI builds).
using mpi_comm    = MPI_Comm;
/// Alias for MPI window used by RMA (portable across MPI/non-MPI builds).
using mpi_win     = MPI_Win;
/// Alias for MPI request handle (portable across MPI/non-MPI builds).
using mpi_request = MPI_Request;

/**
 * @namespace mrcpp::mpi
 * @brief Runtime MPI topology, role flags, and communicators used internally by MRCPP.
 *
 * These externs are set during MRCPP's MPI initialization (see implementation)
 * and describe how the current process participates in computation and data
 * distribution. They are intentionally kept as simple PODs for easy broadcasting
 * and logging.
 */
namespace mpi {

/// If true, the code may choose numerically exact variants of some algorithms.
extern bool numerically_exact;
/// Requested per-node shared-memory window size (in MB) for shared allocations.
extern int  shared_memory_size;

/// Rank of this process in `MPI_COMM_WORLD`.
extern int world_rank;
/// Size of `MPI_COMM_WORLD`.
extern int world_size;

/// Rank within the MRCPP "worker" communicator.
extern int wrk_rank;
/// Size of the MRCPP "worker" communicator.
extern int wrk_size;

/// Rank within the node-local shared-memory communicator.
extern int share_rank;
/// Size of the node-local shared-memory communicator.
extern int share_size;

/// Rank inside the group communicator that clusters ranks by shared-memory groups.
extern int sh_group_rank;

/// True iff this rank belongs to the bank (data service) group.
extern int is_bank;
/// True iff this rank is a worker (i.e., bank client).
extern int is_bankclient;

/// Number of ranks dedicated to the bank (data service).
extern int bank_size;
/// Desired number of bank ranks per node (if configured).
extern int bank_per_node;

/// User/auto-configured OpenMP thread count hint for workers.
extern int omp_threads;
/// If non-zero, honor the environment's OMP thread count for sizing decisions.
extern int use_omp_num_threads;

/// Total number of bank ranks (including any special managers).
extern int tot_bank_size;

/// Upper bound for usable MPI tags (implementation specific).
extern int max_tag;

/// World-rank of the special task-manager bank (if any).
extern int task_bank;

/// Communicator for workers (orbital/function computations).
extern MPI_Comm comm_wrk;
/// Communicator that groups ranks which share physical memory on the same node.
extern MPI_Comm comm_share;
/// Communicator that orders ranks within a shared-memory group.
extern MPI_Comm comm_sh_group;
/// Communicator that includes all bank ranks (and possibly clients for RPC).
extern MPI_Comm comm_bank;

} // namespace mpi
} // namespace mrcpp

namespace mrcpp {

/**
 * @class SharedMemory
 * @brief Thin RAII wrapper around an MPI-3 shared-memory window (per node).
 *
 * A `SharedMemory<T>` instance allocates a node-local window using
 * `MPI_Win_allocate_shared` (only when compiled with MPI). The window can be used
 * to place data structures (e.g., coefficient chunks of a @ref FunctionTree)
 * accessible by all ranks on the same physical node without explicit messaging.
 *
 * @tparam T Element type of the memory window.
 *
 * @par Usage
 * - Construct on one or more ranks of `mpi::comm_share` to allocate a window.
 * - Use `sh_start_ptr`/`sh_end_ptr`/`sh_max_ptr` to manage a simple bump allocator.
 * - Call @ref clear to reset the bump pointer without freeing the window.
 *
 * @note In non-MPI builds, this class becomes a trivial holder and does not
 *       allocate any real shared memory.
 */
template <typename T>
class SharedMemory {
public:
    /**
     * @brief Create (or attach to) a node-local shared-memory window.
     * @param comm   Node-local communicator (typically @ref mpi::comm_share).
     * @param sh_size Window size in megabytes (MB). Only rank 0 in @p comm
     *                dictates the size; other ranks attach to it.
     *
     * @details
     * When `MRCPP_HAS_MPI` is enabled, this constructor calls:
     * - `MPI_Win_allocate_shared` on @p comm with the requested size on rank 0
     *   (and size 0 on others),
     * - `MPI_Win_shared_query` so that every rank obtains a base pointer
     *   into the same window,
     * - Initializes `sh_start_ptr`, `sh_end_ptr`, and `sh_max_ptr`.
     */
    SharedMemory(mrcpp::mpi_comm comm, int sh_size);

    /// Deleted copy constructor to avoid double-free of the MPI window.
    SharedMemory(const SharedMemory &mem) = delete;
    /// Deleted copy assignment.
    SharedMemory<T> &operator=(const SharedMemory<T> &mem) = delete;

    /**
     * @brief Destroy the shared window and release resources.
     * Calls `MPI_Win_free` when built with MPI.
     */
    ~SharedMemory();

    /**
     * @brief Reset the bump pointer so the whole window appears free.
     * Does not deallocate or shrink the MPI window.
     */
    void clear();

    /// Pointer to the beginning of the shared window.
    T *sh_start_ptr{nullptr};
    /// Pointer to one past the last used element (bump pointer).
    T *sh_end_ptr{nullptr};
    /// Pointer to one past the last available element (capacity end).
    T *sh_max_ptr{nullptr};
    /// Underlying MPI window handle.
    mrcpp::mpi_win sh_win{};
    /// Rank of this process within the shared-memory communicator.
    int rank{0};
};

template <int D, typename T> class FunctionTree;

/**
 * @brief Send a @ref FunctionTree to another rank (blocking).
 *
 * Transfers node/structure (and optionally coefficient) chunks to @p dst
 * using point-to-point MPI. If @p nChunks is negative, a small header with the
 * number of chunks is sent first.
 *
 * @tparam D Dimensionality of the function.
 * @tparam T Scalar type (`double` or @ref ComplexDouble).
 * @param tree   FunctionTree to send.
 * @param dst    Destination rank (in @p comm).
 * @param tag    Base MPI tag (chunk indices are offset from this).
 * @param comm   Communicator over which to send.
 * @param nChunks Number of chunks to send; if `<0`, the count is sent first.
 * @param coeff  If true, also send coefficient chunks; otherwise only structure.
 */
template <int D, typename T>
void send_tree(FunctionTree<D, T> &tree, int dst, int tag, mrcpp::mpi_comm comm, int nChunks = -1, bool coeff = true);

/**
 * @brief Receive a @ref FunctionTree from another rank (blocking).
 *
 * Reconstructs tree structure (and optionally coefficients) by receiving the
 * same chunk layout produced by @ref send_tree.
 *
 * @tparam D Dimensionality of the function.
 * @tparam T Scalar type (`double` or @ref ComplexDouble).
 * @param tree   Destination FunctionTree (reinitialized internally).
 * @param src    Source rank (in @p comm).
 * @param tag    Base MPI tag (must match sender).
 * @param comm   Communicator over which to receive.
 * @param nChunks Number of chunks to receive; if `<0`, read the header first.
 * @param coeff  If true, receive coefficient chunks; otherwise only structure.
 */
template <int D, typename T>
void recv_tree(FunctionTree<D, T> &tree, int src, int tag, mrcpp::mpi_comm comm, int nChunks = -1, bool coeff = true);

/**
 * @brief Share a @ref FunctionTree to all ranks in a node-local communicator.
 *
 * Used to mirror the latest version of a shared function across ranks that
 * participate in a shared-memory group, without reconstructing the tree from scratch.
 *
 * @tparam D Dimensionality of the function.
 * @tparam T Scalar type (`double` or @ref ComplexDouble).
 * @param tree FunctionTree to disseminate.
 * @param src  Rank that owns the up-to-date copy (in @p comm).
 * @param tag  Base tag used to coordinate transfers.
 * @param comm Communicator comprising the sharing ranks (e.g., @ref mpi::comm_share).
 */
template <int D, typename T>
void share_tree(FunctionTree<D, T> &tree, int src, int tag, mrcpp::mpi_comm comm);

} // namespace mrcpp
