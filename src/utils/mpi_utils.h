/*
 * MRCPP, a numerical library based on multiresolution analysis and
 * the multiwavelet basis which provide low-scaling algorithms as well as
 * rigorous error control in numerical computations.
 * Copyright (C) 2019 Stig Rune Jensen, Jonas Juselius, Luca Frediani and contributors.
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

#ifdef HAVE_MPI
#include <mpi.h>
#else
using MPI_Comm = int;
using MPI_Win = int;
using MPI_Request = int;
#endif

namespace mrcpp {

/** Share memory within a compute node
 */
class SharedMemory {
public:
    SharedMemory(MPI_Comm comm, int sh_size);
    SharedMemory(const SharedMemory &mem) = delete;
    SharedMemory &operator=(const SharedMemory &mem) = delete;
    ~SharedMemory();

    double *sh_start_ptr; // start of shared block
    double *sh_end_ptr;   // end of used part
    double *sh_max_ptr;   // end of shared block
    MPI_Win sh_win;       // MPI window object
    int rank;             // rank among shared group
};

template <int D> class FunctionTree;

template <int D>
void isend_tree(FunctionTree<D> &tree, int dst, int tag, MPI_Comm comm, MPI_Request *req, int nChunks = -1);
template <int D> void send_tree(FunctionTree<D> &tree, int dst, int tag, MPI_Comm comm, int nChunks = -1);
template <int D> void recv_tree(FunctionTree<D> &tree, int src, int tag, MPI_Comm comm, int nChunks = -1);
template <int D> void share_tree(FunctionTree<D> &tree, int src, int tag, MPI_Comm comm);

} // namespace mrcpp
