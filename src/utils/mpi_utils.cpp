/*
 * MRCPP, a numerical library based on multiresolution analysis and
 * the multiwavelet basis which provide low-scaling algorithms as well as
 * rigorous error control in numerical computations.
 * Copyright (C) 2020 Stig Rune Jensen, Jonas Juselius, Luca Frediani and contributors.
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

#include "mpi_utils.h"
#include "Printer.h"
#include "Timer.h"
#include "trees/FunctionTree.h"
#include "trees/ProjectedNode.h"
#include "trees/ProjectedNodeAllocator.h"

namespace mrcpp {

/** @brief SharedMemory constructor
 *
 *  @param[in] comm: Communicator sharing resources
 *  @param[in] sh_size: Memory size, in MB
 */
SharedMemory::SharedMemory(mrcpp::mpi_comm comm, int sh_size)
        : sh_start_ptr(nullptr)
        , sh_end_ptr(nullptr)
        , sh_max_ptr(nullptr)
        , sh_win(0)
        , rank(0) {

#ifdef MRCPP_HAS_MPI
    MPI_Comm_rank(comm, &this->rank);
    // MPI_Aint types are used for adresses (can be larger than int)
    MPI_Aint size = (rank == 0) ? sh_size : 0; // rank 0 defines length of segment
    size *= 1024 * 1024;                       // ->MB NB: Must be multiplied separately, for not overflow int type.

    int disp_unit = 16; // in order for the compiler to keep aligned
    // size is in bytes
    MPI_Win_allocate_shared(size, disp_unit, MPI_INFO_NULL, comm, &this->sh_start_ptr, &this->sh_win);
    MPI_Win_fence(0, this->sh_win); // wait until finished
    MPI_Aint qsize = 0;
    int qdisp = 0;
    MPI_Win_shared_query(this->sh_win, 0, &qsize, &qdisp, &this->sh_start_ptr);
    MPI_Win_fence(0, this->sh_win);
    this->sh_max_ptr = this->sh_start_ptr + qsize / sizeof(double);
    this->sh_end_ptr = this->sh_start_ptr;
#endif
}

void SharedMemory::clear() {
#ifdef MRCPP_HAS_MPI
    this->sh_end_ptr = this->sh_start_ptr;
#endif
}

SharedMemory::~SharedMemory() {
#ifdef MRCPP_HAS_MPI
    // deallocates the memory block
    MPI_Win_free(&this->sh_win);
#endif
}

/** @brief Send FunctionTree to a given MPI rank using blocking communication
 *
 *  @param[in] tree: FunctionTree to send
 *  @param[in] dst: MPI rank to send to
 *  @param[in] tag: unique identifier
 *  @param[in] comm: Communicator that defines ranks
 *  @param[in] nChunks: Number of memory chunks to send
 *
 *  @details The number of memory chunks must be known before we can send the
 *  tree. This can be specified in the last argument if known a priori, in order
 *  to speed up communication, otherwise it will be communicated in a separate
 *  step before the main communication.
 */
template <int D> void send_tree(FunctionTree<D> &tree, int dst, int tag, mrcpp::mpi_comm comm, int nChunks, bool coeff) {
#ifdef MRCPP_HAS_MPI
    auto &allocator = tree.getProjectedNodeAllocator();

    if (nChunks < 0) {
        nChunks = allocator.getNChunksUsed();
        MPI_Send(&nChunks, sizeof(int), MPI_BYTE, dst, tag, comm);
        println(10, " Sending " << nChunks << " chunks");
    }

    Timer t1;
    for (int iChunk = 0; iChunk < nChunks; iChunk++) {
        MPI_Send(allocator.getNodeChunk(iChunk), allocator.getNodeChunkSize(), MPI_BYTE, dst, tag + iChunk + 1, comm);
        if (coeff)
            MPI_Send(allocator.getCoeffChunk(iChunk), allocator.getCoeffChunkSize(), MPI_BYTE, dst, tag + iChunk + 1001, comm);
    }
    println(10, " Time send                   " << std::setw(30) << t1.elapsed());
#endif
}

/** @brief Receive FunctionTree from a given MPI rank using blocking communication
 *
 *  @param[in] tree: FunctionTree to write into
 *  @param[in] src: MPI rank to receive from
 *  @param[in] tag: unique identifier
 *  @param[in] comm: Communicator that defines ranks
 *  @param[in] nChunks: Number of memory chunks to receive
 *
 *  @details The number of memory chunks must be known before we can receive the
 *  tree. This can be specified in the last argument if known a priori, in order
 *  to speed up communication, otherwise it will be communicated in a separate
 *  step before the main communication.
 */
template <int D> void recv_tree(FunctionTree<D> &tree, int src, int tag, mrcpp::mpi_comm comm, int nChunks, bool coeff) {
#ifdef MRCPP_HAS_MPI
    MPI_Status status;
    auto &allocator = tree.getProjectedNodeAllocator();

    if (nChunks < 0) {
        MPI_Recv(&nChunks, sizeof(int), MPI_BYTE, src, tag, comm, &status);
        println(10, " Receiving " << nChunks << " chunks");
    }

    Timer t1;
    for (int iChunk = 0; iChunk < nChunks; iChunk++) {
        allocator.initChunk(iChunk, coeff);
        MPI_Recv(allocator.getNodeChunk(iChunk), allocator.getNodeChunkSize(), MPI_BYTE, src, tag + iChunk + 1, comm, &status);
        if (coeff)
            MPI_Recv(allocator.getCoeffChunk(iChunk), allocator.getCoeffChunkSize(), MPI_BYTE, src, tag + iChunk + 1001, comm, &status);
    }
    println(10, " Time receive                " << std::setw(30) << t1.elapsed());

    Timer t2;
    allocator.rewritePointers(coeff);
    println(10, " Time rewrite pointers       " << std::setw(30) << t2.elapsed());
#endif
}

/** @brief Share a FunctionTree among MPI processes that share the same physical memory
 *
 *  @param[in] tree: FunctionTree to write into
 *  @param[in] src: MPI rank that last updated the function
 *  @param[in] tag: unique identifier
 *  @param[in] comm: Communicator that defines ranks
 *
 *  @details This function should be called every time a shared function is
 *  updated, in order to update the local memory of each MPI process.
 */
template <int D> void share_tree(FunctionTree<D> &tree, int src, int tag, mrcpp::mpi_comm comm) {
#ifdef MRCPP_HAS_MPI
    Timer t1;
    auto &allocator = tree.getProjectedNodeAllocator();

    int size, rank;
    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);

    for (int dst = 0; dst < size; dst++) {
        if (dst == src) continue;
        int dst_tag = tag * (dst + 1);
        if (rank == src) {
            int nChunks = allocator.getNChunks();
            println(10, " Sending " << nChunks << " chunks");
            MPI_Send(&nChunks, sizeof(int), MPI_BYTE, dst, dst_tag, comm);
            for (int iChunk = 0; iChunk < nChunks; iChunk++) {
                println(10, " Sending chunk " << iChunk);
                MPI_Send(allocator.getNodeChunk(iChunk), allocator.getNodeChunkSize(), MPI_BYTE, dst, dst_tag + iChunk + 1, comm);
            }
        }
        if (rank == dst) {
            MPI_Status status;

            int nChunks;
            MPI_Recv(&nChunks, sizeof(int), MPI_BYTE, src, dst_tag, comm, &status);
            println(10, " Received " << nChunks << " chunks");

            for (int iChunk = 0; iChunk < nChunks; iChunk++) {
                allocator.initChunk(iChunk);
                println(10, " Receiving chunk " << iChunk);
                MPI_Recv(allocator.getNodeChunk(iChunk), allocator.getNodeChunkSize(), MPI_BYTE, src, dst_tag + iChunk + 1, comm, &status);
            }
            allocator.rewritePointers();
        }
    }
    println(10, " Time share                  " << std::setw(30) << t1.elapsed());
#endif
}

template void send_tree<1>(FunctionTree<1> &tree, int dst, int tag, mrcpp::mpi_comm comm, int nChunks, bool coeff);
template void send_tree<2>(FunctionTree<2> &tree, int dst, int tag, mrcpp::mpi_comm comm, int nChunks, bool coeff);
template void send_tree<3>(FunctionTree<3> &tree, int dst, int tag, mrcpp::mpi_comm comm, int nChunks, bool coeff);
template void recv_tree<1>(FunctionTree<1> &tree, int src, int tag, mrcpp::mpi_comm comm, int nChunks, bool coeff);
template void recv_tree<2>(FunctionTree<2> &tree, int src, int tag, mrcpp::mpi_comm comm, int nChunks, bool coeff);
template void recv_tree<3>(FunctionTree<3> &tree, int src, int tag, mrcpp::mpi_comm comm, int nChunks, bool coeff);
template void share_tree<1>(FunctionTree<1> &tree, int src, int tag, mrcpp::mpi_comm comm);
template void share_tree<2>(FunctionTree<2> &tree, int src, int tag, mrcpp::mpi_comm comm);
template void share_tree<3>(FunctionTree<3> &tree, int src, int tag, mrcpp::mpi_comm comm);

} // namespace mrcpp
