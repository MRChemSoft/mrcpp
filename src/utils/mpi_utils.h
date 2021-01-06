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

#pragma once

#ifdef MRCPP_HAS_MPI
#include <mpi.h>
namespace mrcpp {
using mpi_comm = MPI_Comm;
using mpi_win = MPI_Win;
using mpi_request = MPI_Request;
} // namespace mrcpp
#else
namespace mrcpp {
using mpi_comm = int;
using mpi_win = int;
using mpi_request = int;
} // namespace mrcpp
#endif

namespace mrcpp {

/** @class SharedMemory
 *
 *  @brief Shared memory block within a compute node
 *
 *  @details This class defines a shared memory window in a shared MPI
 *  communicator. In order to allocate a FunctionTree in shared memory,
 *  simply pass a SharedMemory object to the FunctionTree constructor.
 */
class SharedMemory {
public:
    SharedMemory(mrcpp::mpi_comm comm, int sh_size);
    SharedMemory(const SharedMemory &mem) = delete;
    SharedMemory &operator=(const SharedMemory &mem) = delete;
    ~SharedMemory();

    void clear(); // show shared memory as entirely available

    double *sh_start_ptr;  // start of shared block
    double *sh_end_ptr;    // end of used part
    double *sh_max_ptr;    // end of shared block
    mrcpp::mpi_win sh_win; // MPI window object
    int rank;              // rank among shared group
};

template <int D> class FunctionTree;

template <int D> void send_tree(FunctionTree<D> &tree, int dst, int tag, mrcpp::mpi_comm comm, int nChunks = -1, bool coeff = true);
template <int D> void recv_tree(FunctionTree<D> &tree, int src, int tag, mrcpp::mpi_comm comm, int nChunks = -1, bool coeff = true);
template <int D> void share_tree(FunctionTree<D> &tree, int src, int tag, mrcpp::mpi_comm comm);

} // namespace mrcpp
