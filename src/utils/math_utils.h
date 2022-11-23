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

#include <Eigen/Core>
#include <Eigen/Eigenvalues>

#include "MRCPP/mrcpp_declarations.h"

using IntVector = Eigen::VectorXi;
using DoubleVector = Eigen::VectorXd;
using ComplexVector = Eigen::VectorXcd;

using IntMatrix = Eigen::MatrixXi;
using DoubleMatrix = Eigen::MatrixXd;
using ComplexMatrix = Eigen::MatrixXcd;

using ComplexDouble = std::complex<double>;

namespace mrcpp {
namespace math_utils {

double binomial_coeff(int n, int j);
Eigen::VectorXd get_binomial_coefs(unsigned int order);

double factorial(int n);
int ipow(int m, int e);

Eigen::MatrixXd tensor_product(const Eigen::MatrixXd &A, const Eigen::MatrixXd &B);
Eigen::MatrixXd tensor_product(const Eigen::MatrixXd &A, const Eigen::VectorXd &B);
Eigen::MatrixXd tensor_product(const Eigen::VectorXd &A, const Eigen::MatrixXd &B);
Eigen::MatrixXd tensor_product(const Eigen::VectorXd &A, const Eigen::VectorXd &B);

void tensor_self_product(const Eigen::VectorXd &A, Eigen::VectorXd &B);
void tensor_self_product(const Eigen::VectorXd &A, Eigen::MatrixXd &B);

double matrix_norm_inf(const Eigen::MatrixXd &M);
double matrix_norm_1(const Eigen::MatrixXd &M);
double matrix_norm_2(const Eigen::MatrixXd &M);

void apply_filter(double *out, double *in, const Eigen::MatrixXd &filter, int kp1, int kp1_dm1, double fac);

void tensor_expand_coefs(int dim, int dir, int kp1, int kp1_d, const Eigen::MatrixXd &primitive, Eigen::VectorXd &expanded);

void tensor_expand_coords_2D(int kp1, const Eigen::MatrixXd &primitive, Eigen::MatrixXd &expanded);
void tensor_expand_coords_3D(int kp1, const Eigen::MatrixXd &primitive, Eigen::MatrixXd &expanded);

ComplexMatrix hermitian_matrix_pow(const ComplexMatrix &A, double b);
ComplexMatrix diagonalize_hermitian_matrix(const ComplexMatrix &A, DoubleVector &diag);
void diagonalize_block(ComplexMatrix &M, ComplexMatrix &U, int nstart, int nsize);

template <class T> std::vector<std::vector<T>> cartesian_product(std::vector<T> A, std::vector<T> B);
template <class T> std::vector<std::vector<T>> cartesian_product(std::vector<std::vector<T>> l_A, std::vector<T> B);
template <class T> std::vector<std::vector<T>> cartesian_product(std::vector<T> a, int dim);

template <int D> double calc_distance(const Coord<D> &a, const Coord<D> &b);

} // namespace math_utils
} // namespace mrcpp
