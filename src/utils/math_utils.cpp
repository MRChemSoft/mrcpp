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

#include <cmath>
#include <cstdio>
#include <fstream>

#include "MRCPP/constants.h"
#include "Printer.h"
#include "math_utils.h"

#ifdef HAVE_BLAS
extern "C" {
#include BLAS_H
}
#endif

using namespace Eigen;

namespace mrcpp {

int math_utils::ipow(int m, int e) {
    if (e < 0) MSG_ABORT("Exponent cannot be negative: " << e)
    int result = 1;
    for (int i = 0; i < e; i++) { result *= m; }
    return result;
}

double math_utils::matrix_norm_2(const MatrixXd &M) {
    return M.lpNorm<2>();
}

double math_utils::matrix_norm_1(const MatrixXd &M) {
    return M.colwise().lpNorm<1>().maxCoeff();
}

double math_utils::matrix_norm_inf(const MatrixXd &M) {
    return M.rowwise().lpNorm<1>().maxCoeff();
}

double math_utils::binomial_coeff(int n, int j) {
    double binomial_n_j = 1.0;
    if (n < 0 || j < 0 || j > n) {
        MSG_ERROR("Negative argument or j > n is not defined.");
    } else {
        int k = 0;
        while (k < j) {
            binomial_n_j *= (double)(n - k);
            k += 1;
        }
        binomial_n_j /= factorial(j);
    }
    return binomial_n_j;
}

VectorXd math_utils::get_binomial_coefs(unsigned int order) {
    VectorXd coefs = VectorXd::Ones(order + 1);
    for (int k = 0; k <= (int)order; k++) { coefs[k] = math_utils::binomial_coeff(order, k); }
    return coefs;
}

double math_utils::factorial(int n) {
    int k = 1;
    double fac_n = 1.0;

    if (n < 0) {
        MSG_ABORT("Negative argument is not defined.");
    } else if (n == 0 || n == 1)
        return 1.0;
    else if (n > 1) {
        while (k <= n) {
            fac_n *= (double)k;
            k += 1;
        }
    }
    return fac_n;
}

MatrixXd math_utils::tensor_product(const MatrixXd &A, const MatrixXd &B) {
    int Ar = A.rows();
    int Ac = A.cols();
    int Br = B.rows();
    int Bc = B.cols();
    MatrixXd tprod(Ar * Br, Ac * Bc);
    for (int i = 0; i < Ar; i++) {
        for (int j = 0; j < Ac; j++) { tprod.block(i * Br, j * Bc, Br, Bc) = A(i, j) * B; }
    }
    return tprod;
}

MatrixXd math_utils::tensor_product(const MatrixXd &A, const VectorXd &B) {
    int Ar = A.rows();
    int Ac = A.cols();
    int Br = B.rows();
    MatrixXd tprod(Ar * Br, Ac);
    for (int i = 0; i < Br; i++) { tprod.block(i * Br, 0, Ar, Ac) = A * B(i); }
    return tprod;
}

MatrixXd math_utils::tensor_product(const VectorXd &A, const MatrixXd &B) {
    int Ar = A.rows();
    int Br = B.rows();
    int Bc = B.cols();
    MatrixXd tprod(Ar * Br, Ar);
    for (int i = 0; i < Ar; i++) { tprod.block(i * Br, 0, Br, Bc) = A(i) * B; }
    return tprod;
}

MatrixXd math_utils::tensor_product(const VectorXd &A, const VectorXd &B) {
    int Ar = A.rows();
    int Br = B.rows();
    MatrixXd tprod(Ar, Br);
    for (int i = 0; i < Ar; i++) { tprod.block(i, 0, 1, Br) = A(i) * B.transpose(); }
    return tprod;
}

void math_utils::tensor_self_product(const VectorXd &A, VectorXd &tprod) {
    int Ar = A.rows();
    for (int i = 0; i < Ar; i++) { tprod.segment(i * Ar, Ar) = A(i) * A; }
}

void math_utils::tensor_self_product(const VectorXd &A, MatrixXd &tprod) {
    int Ar = A.rows();
    for (int i = 0; i < Ar; i++) { tprod.block(i, 0, 1, Ar) = A(i) * A; }
}

template <typename T>
void math_utils::apply_filter(T *out, T *in, const MatrixXd &filter, int kp1, int kp1_dm1, double fac) {
    if constexpr (std::is_same<T, double>::value) {
        Map<MatrixXd> f(in, kp1, kp1_dm1);
        Map<MatrixXd> g(out, kp1_dm1, kp1);
        if (fac < MachineZero) {
            g.noalias() = f.transpose() * filter;
        } else {
            g.noalias() += f.transpose() * filter;
        }
    } else if constexpr (std::is_same<T, ComplexDouble>::value) {
        Map<MatrixXcd> f(in, kp1, kp1_dm1);
        Map<MatrixXcd> g(out, kp1_dm1, kp1);
        if (fac < MachineZero) {
            g.noalias() = f.transpose() * filter;
        } else {
            g.noalias() += f.transpose() * filter;
        }
    } else
        NOT_IMPLEMENTED_ABORT;
}

void math_utils::tensor_expand_coefs(int dim, int dir, int kp1, int kp1_d, const MatrixXd &primitive, VectorXd &expanded) {
    if (dir < dim - 1) {
        int idx = math_utils::ipow(kp1, dir + 1);
        int nelem = idx * kp1;
        int pos = kp1_d - nelem;
        int inpos = kp1_d - idx;
        for (int i = 0; i < kp1; i++) expanded.segment(pos + i * idx, idx) = expanded.segment(inpos, idx) * primitive.col(dir + 1)(i);
        tensor_expand_coefs(dim, dir + 1, kp1, kp1_d, primitive, expanded);
    }
}

void math_utils::tensor_expand_coords_2D(int kp1, const MatrixXd &primitive, MatrixXd &expanded) {
    int n = 0;
    for (int i = 0; i < kp1; i++) {
        for (int j = 0; j < kp1; j++) {
            expanded(0, n) = primitive(0, j);
            expanded(1, n) = primitive(1, i);
            n++;
        }
    }
}

void math_utils::tensor_expand_coords_3D(int kp1, const MatrixXd &primitive, MatrixXd &expanded) {
    int n = 0;
    for (int i = 0; i < kp1; i++) {
        for (int j = 0; j < kp1; j++) {
            for (int k = 0; k < kp1; k++) {
                expanded(0, n) = primitive(0, k);
                expanded(1, n) = primitive(1, j);
                expanded(2, n) = primitive(2, i);
                n++;
            }
        }
    }
}

ComplexMatrix math_utils::diagonalize_hermitian_matrix(const ComplexMatrix &A, DoubleVector &diag) {
    Eigen::SelfAdjointEigenSolver<ComplexMatrix> es(A.cols());
    es.compute(A);
    diag = es.eigenvalues();
    return es.eigenvectors();
}

ComplexMatrix math_utils::hermitian_matrix_pow(const ComplexMatrix &A, double b) {
    DoubleVector diag;
    ComplexMatrix U = diagonalize_hermitian_matrix(A, diag);

    DoubleMatrix B = DoubleMatrix::Zero(A.rows(), A.cols());
    for (int i = 0; i < diag.size(); i++) {
        if (std::abs(diag(i)) < mrcpp::MachineZero) {
            B(i, i) = 0.0;
        } else {
            B(i, i) = std::pow(diag(i), b);
        }
    }
    return U * B * U.adjoint();
}

void math_utils::diagonalize_block(ComplexMatrix &A, ComplexMatrix &U, int nstart, int nsize) {
    Eigen::SelfAdjointEigenSolver<ComplexMatrix> es(nsize);
    es.compute(A.block(nstart, nstart, nsize, nsize));
    ComplexMatrix ei_vec = es.eigenvectors();
    ComplexVector ei_val = es.eigenvalues().cast<ComplexDouble>();
    U.block(nstart, nstart, nsize, nsize) = ei_vec;
    A.block(nstart, nstart, nsize, nsize) = ei_val.asDiagonal();
}

template <int D>
double math_utils::calc_distance(const Coord<D> &a, const Coord<D> &b) {
    double r = 0.0;
    for (int i = 0; i < D; i++) { r += std::pow(a[i] - b[i], 2.0); }
    return std::sqrt(r);
}

template <class T>
std::vector<std::vector<T>> math_utils::cartesian_product(std::vector<T> A, std::vector<T> B) {
    std::vector<std::vector<T>> output;
    for (auto &a : A) {
        for (auto &b : B) output.push_back(std::vector<T>{a, b});
    }
    return output;
}

template <class T>
std::vector<std::vector<T>> math_utils::cartesian_product(std::vector<std::vector<T>> l_A, std::vector<T> B) {
    std::vector<std::vector<T>> output;
    for (auto A : l_A) {
        for (auto &b : B) {
            A.push_back(b);
            output.push_back(A);
            A.pop_back();
        }
    }
    return output;
}

template <class T>
std::vector<std::vector<T>> math_utils::cartesian_product(std::vector<T> A, int dim) {
    std::vector<std::vector<T>> output;
    if (dim < 0) MSG_ABORT("Dimension has to be 1 or greater")
    if (dim == 1) {
        for (auto v : A) output.push_back(std::vector<T>{v});
    } else {
        output = cartesian_product(A, A);
        for (auto i = 0; i < dim - 2; i++) output = cartesian_product(output, A);
    }
    return output;
}

template void math_utils::apply_filter<double>(double *out, double *in, const Eigen::MatrixXd &filter, int kp1, int kp1_dm1, double fac);
template void math_utils::apply_filter<ComplexDouble>(ComplexDouble *out, ComplexDouble *in, const Eigen::MatrixXd &filter, int kp1, int kp1_dm1, double fac);

template double math_utils::calc_distance<1>(const Coord<1> &a, const Coord<1> &b);
template double math_utils::calc_distance<2>(const Coord<2> &a, const Coord<2> &b);
template double math_utils::calc_distance<3>(const Coord<3> &a, const Coord<3> &b);

template std::vector<std::vector<int>> math_utils::cartesian_product(std::vector<int> A, std::vector<int> B);
template std::vector<std::vector<int>> math_utils::cartesian_product(std::vector<std::vector<int>> l_A, std::vector<int> B);
template std::vector<std::vector<int>> math_utils::cartesian_product(std::vector<int> A, int dim);

template std::vector<std::vector<double>> math_utils::cartesian_product(std::vector<double> A, std::vector<double> B);
template std::vector<std::vector<double>> math_utils::cartesian_product(std::vector<std::vector<double>> l_A, std::vector<double> B);
template std::vector<std::vector<double>> math_utils::cartesian_product(std::vector<double> A, int dim);

} // namespace mrcpp