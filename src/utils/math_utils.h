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

/* \file math_utils.h
 *
 * \breif Collection of misc math funcs.
 */

#pragma once
/**
 * @file
 * @brief Linear algebra and small numerical helpers built on top of Eigen.
 *
 * This header exposes a compact set of utilities frequently used across MRCPP:
 * - Common Eigen-based type aliases (`Vector`, `Matrix`, complex scalars).
 * - Small combinatorial helpers (factorial, binomial, integer powers).
 * - Tensor products and self–outer-products.
 * - Matrix norms (1, 2, and ∞).
 * - Multi-index/tensor expansion helpers for separable bases.
 * - Hermitian eigendecompositions and block diagonalization.
 * - Cartesian products of small sets/vectors.
 * - Euclidean distance for MRCPP coordinates.
 */

#include <Eigen/Core>
#include <Eigen/Eigenvalues>

#include "MRCPP/mrcpp_declarations.h"

/** @name Eigen type aliases
 *  @brief Short, explicit aliases for Eigen vectors/matrices used in MRCPP.
 *  @{ */
using IntVector    = Eigen::VectorXi;   ///< Column vector of integers.
using DoubleVector = Eigen::VectorXd;   ///< Column vector of doubles.
using ComplexVector= Eigen::VectorXcd;  ///< Column vector of complex doubles.

using IntMatrix    = Eigen::MatrixXi;   ///< Integer matrix.
using DoubleMatrix = Eigen::MatrixXd;   ///< Double-precision matrix.
using ComplexMatrix= Eigen::MatrixXcd;  ///< Complex double-precision matrix.

using ComplexDouble= std::complex<double>; ///< Convenience alias for complex<double>.
/** @} */

namespace mrcpp {
/**
 * @namespace mrcpp::math_utils
 * @brief Numerical utilities layered on Eigen; header-only declarations.
 */
namespace math_utils {

/**
 * @brief Binomial coefficient \f$\binom{n}{j}\f$.
 * @param n Non-negative integer.
 * @param j Non-negative integer, \f$0 \le j \le n\f$.
 * @return \f$\frac{n!}{(n-j)!\,j!}\f$ as a double.
 * @note For out-of-domain inputs, the implementation may log an error.
 */
double binomial_coeff(int n, int j);

/**
 * @brief Pascal row of binomial coefficients.
 * @param order Row index \f$n\f$.
 * @return Vector \f$[\binom{n}{0}, \ldots, \binom{n}{n}]\f$.
 */
Eigen::VectorXd get_binomial_coefs(unsigned int order);

/**
 * @brief Factorial for non-negative integers.
 * @param n \f$n \ge 0\f$.
 * @return \f$n!\f$ as a double.
 */
double factorial(int n);

/**
 * @brief Integer power \f$m^e\f$ for \f$e \ge 0\f$ (loop-based; exact for small ranges).
 * @param m Base (integer).
 * @param e Exponent (integer, \f$e \ge 0\f$).
 * @return \f$m^e\f$ as an int.
 */
int ipow(int m, int e);

/** @name Tensor/Kronecker products
 *  @brief Kronecker products and outer products for vectors/matrices.
 *  @{
 */

/**
 * @brief Kronecker product \f$A \otimes B\f$.
 */
Eigen::MatrixXd tensor_product(const Eigen::MatrixXd &A, const Eigen::MatrixXd &B);

/**
 * @brief Kronecker product of a matrix and a column vector (treated as \f$B\f$).
 */
Eigen::MatrixXd tensor_product(const Eigen::MatrixXd &A, const Eigen::VectorXd &B);

/**
 * @brief Kronecker product of a column vector and a matrix (treated as \f$B\f$).
 */
Eigen::MatrixXd tensor_product(const Eigen::VectorXd &A, const Eigen::MatrixXd &B);

/**
 * @brief Outer product \f$A B^\top\f$ of two column vectors.
 */
Eigen::MatrixXd tensor_product(const Eigen::VectorXd &A, const Eigen::VectorXd &B);

/**
 * @brief Self outer-product \f$A \otimes A\f$ into a flat vector.
 * @param A Input column vector.
 * @param B Output vector (size must be \f$\mathrm{size}(A)^2\f$).
 */
void tensor_self_product(const Eigen::VectorXd &A, Eigen::VectorXd &B);

/**
 * @brief Self outer-product \f$A A^\top\f$ into a matrix.
 * @param A Input column vector.
 * @param B Output matrix (square, same dimension as \f$A\f$).
 */
void tensor_self_product(const Eigen::VectorXd &A, Eigen::MatrixXd &B);
/** @} */

/** @name Matrix norms
 *  @brief Induced matrix norms consistent with Eigen semantics.
 *  @{
 */
/**
 * @brief Infinity norm \f$\|M\|_\infty\f$ (max row 1-norm).
 */
double matrix_norm_inf(const Eigen::MatrixXd &M);

/**
 * @brief 1-norm \f$\|M\|_1\f$ (max column 1-norm).
 */
double matrix_norm_1(const Eigen::MatrixXd &M);

/**
 * @brief Spectral norm \f$\|M\|_2\f$ (largest singular value).
 */
double matrix_norm_2(const Eigen::MatrixXd &M);
/** @} */

/**
 * @brief Apply a linear filter to a coefficient block (templated on scalar type).
 *
 * Conceptually computes and accumulates a matrix product of the form
 * \f$G \leftarrow G + F^\top \cdot \mathrm{filter}\f$, where \f$F\f$ and \f$G\f$
 * are views over `in` and `out` with shapes derived from \p kp1 and \p kp1_dm1.
 *
 * @tparam T Scalar type (`double` or `ComplexDouble`).
 * @param[out] out Output buffer (accumulation destination).
 * @param[in]  in  Input buffer (interpreted as a matrix view).
 * @param[in]  filter Dense filter matrix to apply.
 * @param[in]  kp1 Leading polynomial order + 1 (per MR basis).
 * @param[in]  kp1_dm1 \f$\text{kp1}^{D-1}\f$ helper (stride in the mapped view).
 * @param[in]  fac If zero, overwrite the destination; otherwise accumulate.
 * @warning Buffers must contain at least the required number of elements for the
 *          implicit matrix views.
 */
template <typename T>
void apply_filter(T *out, T *in, const Eigen::MatrixXd &filter, int kp1, int kp1_dm1, double fac);

/**
 * @brief Expand separable 1D coefficient blocks into an \f$ \text{dim}\f$-D tensor layout.
 *
 * Recursively multiplies along dimensions using the columns of \p primitive.
 *
 * @param dim   Spatial dimensionality (1–3).
 * @param dir   Current recursion direction (starting at 0).
 * @param kp1   Polynomial order + 1.
 * @param kp1_d \f$\text{kp1}^d\f$ where \f$d = \text{dim}\f$ (total coefficients).
 * @param primitive Matrix with primitive 1D basis values per dimension.
 * @param[in,out] expanded Buffer holding intermediate input and final expanded output.
 */
void tensor_expand_coefs(int dim, int dir, int kp1, int kp1_d, const Eigen::MatrixXd &primitive, Eigen::VectorXd &expanded);

/**
 * @brief Generate 2D sampling coordinates on a tensor grid spanned by primitive 1D points.
 * @param kp1 Points per axis.
 * @param primitive Matrix whose rows are the per-axis primitive coordinates.
 * @param[out] expanded Output matrix of size \f$(\text{kp1}^2) \times 2\f$.
 */
void tensor_expand_coords_2D(int kp1, const Eigen::MatrixXd &primitive, Eigen::MatrixXd &expanded);

/**
 * @brief Generate 3D sampling coordinates on a tensor grid spanned by primitive 1D points.
 * @param kp1 Points per axis.
 * @param primitive Matrix whose rows are the per-axis primitive coordinates.
 * @param[out] expanded Output matrix of size \f$(\text{kp1}^3) \times 3\f$.
 */
void tensor_expand_coords_3D(int kp1, const Eigen::MatrixXd &primitive, Eigen::MatrixXd &expanded);

/**
 * @brief Hermitian matrix power \f$A^b\f$ via eigendecomposition.
 * @param A Hermitian (self-adjoint) complex matrix.
 * @param b Real exponent.
 * @return \f$U\,\mathrm{diag}(\lambda_i^b)\,U^\dagger\f$ where \f$A = U\,\mathrm{diag}(\lambda_i)\,U^\dagger\f$.
 * @note Eigenvalues with magnitude near zero are guarded to avoid blow-ups for negative \p b.
 */
ComplexMatrix hermitian_matrix_pow(const ComplexMatrix &A, double b);

/**
 * @brief Diagonalize a Hermitian matrix.
 * @param A Input Hermitian matrix (not modified).
 * @param[out] diag Real vector of eigenvalues (ascending).
 * @return Matrix of eigenvectors as columns (unitary).
 */
ComplexMatrix diagonalize_hermitian_matrix(const ComplexMatrix &A, DoubleVector &diag);

/**
 * @brief In-place diagonalization of a Hermitian sub-block.
 *
 * Replaces the \f$n_\text{size}\times n_\text{size}\f$ block of @p M at
 * \f$(n_\text{start}, n_\text{start})\f$ with its eigenvalues on the diagonal
 * and writes the corresponding eigenvectors into the same block of @p U.
 *
 * @param[in,out] M Matrix containing the Hermitian sub-block to diagonalize.
 * @param[out] U Matrix receiving the block eigenvectors.
 * @param nstart Upper-left index of the block.
 * @param nsize  Size of the (square) block.
 */
void diagonalize_block(ComplexMatrix &M, ComplexMatrix &U, int nstart, int nsize);

/** @name Cartesian products
 *  @brief Simple, small-container cartesian products (for enumeration tasks).
 *  @{
 */
/**
 * @brief Cartesian product \f$A \times B\f$.
 * @tparam T Element type.
 * @param A First list.
 * @param B Second list.
 * @return Vector of pairs `[a, b]`.
 */
template <class T>
std::vector<std::vector<T>> cartesian_product(std::vector<T> A, std::vector<T> B);

/**
 * @brief Cartesian product \f{(l\_A) \times B\f}, where each element of @p l_A is itself a tuple.
 * @tparam T Element type.
 * @param l_A List of partial tuples.
 * @param B   Second list.
 * @return Concatenated tuples.
 */
template <class T>
std::vector<std::vector<T>> cartesian_product(std::vector<std::vector<T>> l_A, std::vector<T> B);

/**
 * @brief Repeated cartesian power \f$A^{\times \text{dim}}\f$.
 * @tparam T Element type.
 * @param a   Base list.
 * @param dim Number of repeats (\f$\ge 1\f$).
 * @return All length-\p dim tuples with elements from \p a.
 */
template <class T>
std::vector<std::vector<T>> cartesian_product(std::vector<T> a, int dim);
/** @} */

/**
 * @brief Euclidean distance between two D-dimensional coordinates.
 * @tparam D Dimension (compile-time).
 * @param a First point.
 * @param b Second point.
 * @return \f$\sqrt{\sum_{i=1}^D (a_i-b_i)^2}\f$.
 */
template <int D>
double calc_distance(const Coord<D> &a, const Coord<D> &b);

} // namespace math_utils
} // namespace mrcpp