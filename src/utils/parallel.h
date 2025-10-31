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
 * @brief MPI/OpenMP orchestration and collectives for MRCPP.
 *
 * This header declares the process/thread orchestration utilities and the
 * common collective/point-to-point helpers used by MRCPP to distribute
 * multiresolution data structures across MPI ranks (and optionally coordinate
 * with OpenMP threads). It provides:
 *
 * - Initialization/finalization of the MRCPP MPI environment.
 * - Rank/topology helpers (e.g., “grand master”, ownership checks).
 * - Typed send/recv/broadcast for @ref mrcpp::CompFunction and trees.
 * - Element-wise allreduce helpers for Eigen vectors/matrices.
 *
 * All MPI symbols are no-ops in non-MPI builds (compiled without
 * `MRCPP_HAS_MPI`), allowing the same interface to work in serial.
 */

#include <Eigen/Core>

#include "CompFunction.h"
#include "mpi_utils.h"
#include "trees/MultiResolutionAnalysis.h"
#include <map>
#include <vector>

// define a class for things that can be sent with MPI

using namespace Eigen;

using IntVector    = Eigen::VectorXi;
using DoubleVector = Eigen::VectorXd;
using ComplexVector= Eigen::VectorXcd;

using IntMatrix    = Eigen::MatrixXi;
using DoubleMatrix = Eigen::MatrixXd;
using ComplexMatrix= Eigen::MatrixXcd;

namespace mrcpp {

/**
 * @namespace mrcpp::omp
 * @brief OpenMP runtime hints used by the parallel layer.
 */
namespace omp {
extern int n_threads; ///< Number of OpenMP threads MRCPP intends to use.
} // namespace omp

class Bank;              ///< Forward declaration of the in-memory data bank.
extern Bank dataBank;    ///< Global bank instance used by bank ranks.

/**
 * @namespace mrcpp::mpi
 * @brief MPI utilities, communicators, and collectives.
 *
 * Functions in this namespace act as thin wrappers around MPI and encode
 * MRCPP’s distribution policy for component functions and trees.
 */
namespace mpi {

/** @brief World ranks assigned to bank masters (control/data services). */
extern std::vector<int> bankmaster;

/**
 * @brief Initialize MRCPP’s MPI environment and process topology.
 *
 * Sets up communicators (workers, shared-memory groups, bank group),
 * partitions ranks into worker/bank roles, and configures OpenMP thread
 * counts per rank. Safe to call exactly once at program start.
 */
void initialize();

/**
 * @brief Finalize MRCPP’s MPI environment.
 *
 * Performs a global barrier, closes the global data bank (if present), and
 * calls `MPI_Finalize()` in MPI builds. Safe to call once at program exit.
 */
void finalize();

/**
 * @brief Rank barrier on a given communicator.
 * @param comm MPI communicator to synchronize.
 *
 * In non-MPI builds this is a no-op.
 */
void barrier(MPI_Comm comm);

/**
 * @brief Whether this rank is the global worker “grand master”.
 * @return @c true iff world rank is 0 and the rank is a worker (not a bank).
 */
bool grand_master();

/**
 * @brief Whether this rank is the master of its shared-memory group.
 * @return @c true iff rank is 0 within @ref mpi::comm_share.
 */
bool share_master();

/**
 * @name Ownership helpers
 * @brief Determine whether an object/function is owned by this rank.
 * @{
 */

/**
 * @brief Ownership test for an index.
 * @param j Global function index.
 * @return @c true if @c j maps to this rank under MRCPP’s block-cyclic policy.
 */
bool my_func(int j);

/**
 * @brief Ownership test for a component function (const ref).
 * @param func Component function to test.
 * @return @c true if @p func belongs to this rank (by @c func.rank()).
 */
bool my_func(const CompFunction<3> &func);

/**
 * @brief Ownership test for a component function (pointer).
 * @param func Pointer to component function.
 * @return @c true if @p func belongs to this rank (by @c func->rank()).
 */
bool my_func(CompFunction<3> *func);
/** @} */

/**
 * @brief Free memory held by functions not owned by this rank.
 * @param Phi Vector of component functions; foreign entries are freed in place.
 */
void free_foreign(CompFunctionVector &Phi);

/**
 * @name Point-to-point transfers for component functions
 * @brief Send/receive/share a @ref mrcpp::CompFunction across ranks.
 * @{
 */

/**
 * @brief Send a component function to a destination rank.
 * @param func Function to send.
 * @param dst Destination world rank.
 * @param tag Message tag base (submessages will offset from this).
 * @param comm Communicator (default: worker communicator).
 *
 * Sends the function header followed by its component trees. Assumes the
 * receiver uses @ref recv_function with the same @p tag and @p comm.
 */
void send_function(const CompFunction<3> &func, int dst, int tag, MPI_Comm comm = mpi::comm_wrk);

/**
 * @brief Receive a component function from a source rank.
 * @param func Function to receive into (resized as needed).
 * @param src Source world rank.
 * @param tag Message tag base (must match sender).
 * @param comm Communicator (default: worker communicator).
 */
void recv_function(CompFunction<3> &func, int src, int tag, MPI_Comm comm = mpi::comm_wrk);

/**
 * @brief Update shared-memory replicas of a function after modification.
 * @param func Function to share (must be marked shared).
 * @param src Rank that produced the update.
 * @param tag Base tag for the transfer.
 * @param comm Communicator that defines the sharing group.
 *
 * Only has effect if the function was allocated in shared memory.
 */
void share_function(CompFunction<3> &func, int src, int tag, MPI_Comm comm);
/** @} */

/**
 * @brief Reduce (sum/accumulate) a function onto rank 0 of @p comm.
 * @param prec Cropping precision applied after each accumulation.
 * @param func Function buffer holding the local contribution; on rank 0 it
 *             becomes the global sum; on other ranks it may be left unchanged.
 * @param comm Communicator over which to reduce.
 *
 * Uses a binary-tree pattern to send odd ranks to preceding even ranks; the
 * receiver adds and crops to control growth.
 */
void reduce_function(double prec, CompFunction<3> &func, MPI_Comm comm);

/**
 * @brief Broadcast a function from rank 0 to all ranks in @p comm.
 * @param func Buffer to receive (or hold, on root) the broadcasted function.
 * @param comm Communicator to broadcast over.
 *
 * Implements a reverse of the binary-tree pattern used by
 * @ref reduce_function.
 */
void broadcast_function(CompFunction<3> &func, MPI_Comm comm);

/**
 * @name Tree collectives (no coefficient payload)
 * @brief Perform collectives on @ref mrcpp::FunctionTree without coefficients.
 * @{
 */

/**
 * @brief Reduce (union) grids from all ranks to rank 0, excluding coeffs.
 * @tparam T Coefficient scalar type of the tree.
 * @param tree Output/input tree on each rank; on rank 0 it becomes the union.
 * @param comm Communicator.
 */
template <typename T>
void reduce_Tree_noCoeff(mrcpp::FunctionTree<3, T> &tree, MPI_Comm comm);

/**
 * @brief Build local union grid, reduce to rank 0, then broadcast to all.
 * @tparam T Coefficient scalar type.
 * @param tree Target tree to hold the global union grid (no coeffs).
 * @param Phi  Vector of trees whose grids contribute to the union.
 * @param comm Communicator.
 */
template <typename T>
void allreduce_Tree_noCoeff(mrcpp::FunctionTree<3, T> &tree,
                            std::vector<FunctionTree<3, T>> &Phi,
                            MPI_Comm comm);

/**
 * @brief Broadcast a no-coeff tree from rank 0 to all ranks.
 * @tparam T Coefficient scalar type.
 * @param tree Tree to broadcast/receive.
 * @param comm Communicator.
 */
template <typename T>
void broadcast_Tree_noCoeff(mrcpp::FunctionTree<3, T> &tree, MPI_Comm comm);

/**
 * @brief Build union grid from owned components in @p Phi, allreduce to all.
 * @tparam T Coefficient scalar type.
 * @param tree Output tree receiving the global union grid.
 * @param Phi  Vector of component functions contributing their grids.
 * @param comm Communicator.
 */
template <typename T>
void allreduce_Tree_noCoeff(mrcpp::FunctionTree<3, T> &tree,
                            std::vector<CompFunction<3>> &Phi,
                            MPI_Comm comm);
/** @} */

/**
 * @name Element-wise allreduce (sum) helpers
 * @brief Sum across ranks into every rank for Eigen containers.
 * @{
 */

/** @brief In-place element-wise sum allreduce for integer vectors. */
void allreduce_vector(IntVector &vec, MPI_Comm comm);
/** @brief In-place element-wise sum allreduce for double vectors. */
void allreduce_vector(DoubleVector &vec, MPI_Comm comm);
/** @brief In-place element-wise sum allreduce for complex vectors. */
void allreduce_vector(ComplexVector &vec, MPI_Comm comm);

/** @brief In-place element-wise sum allreduce for integer matrices. */
void allreduce_matrix(IntMatrix &vec, MPI_Comm comm);
/** @brief In-place element-wise sum allreduce for double matrices. */
void allreduce_matrix(DoubleMatrix &mat, MPI_Comm comm);
/** @brief In-place element-wise sum allreduce for complex matrices. */
void allreduce_matrix(ComplexMatrix &mat, MPI_Comm comm);
/** @} */

} // namespace mpi

} // namespace mrcpp