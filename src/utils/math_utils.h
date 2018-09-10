/* \file math_utils.h
 *
 * \breif Collection of misc math funcs.
 */

#pragma once

#include <Eigen/Core>

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

void apply_filter(double *out, double *in, const Eigen::MatrixXd &filter,
                  int kp1, int kp1_dm1, double fac);

void tensor_expand_coefs(int dim, int dir, int kp1, int kp1_d,
                         const Eigen::MatrixXd &primitive,
                         Eigen::VectorXd &expanded);

void tensor_expand_coords_2D(int kp1, const Eigen::MatrixXd &primitive, Eigen::MatrixXd &expanded);
void tensor_expand_coords_3D(int kp1, const Eigen::MatrixXd &primitive, Eigen::MatrixXd &expanded);

double calc_distance(int D, const double *a, const double *b);

} // namespace math_utils
} // namespace mrcpp
