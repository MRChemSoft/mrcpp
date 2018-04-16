/* \file mwmath.h
 *
 * \breif Collection of misc math funcs.
 */

#pragma once

#include <Eigen/Core>

namespace mrcpp {
namespace mwmath {

double binomialCoeff(int n, int j);
Eigen::VectorXd getBinomialCoefs(unsigned int order);

double factorial(int n);
int ipow(int m, int e);

Eigen::MatrixXd tensorproduct(const Eigen::MatrixXd &A, const Eigen::MatrixXd &B);
Eigen::MatrixXd tensorproduct(const Eigen::MatrixXd &A, const Eigen::VectorXd &B);
Eigen::MatrixXd tensorproduct(const Eigen::VectorXd &A, const Eigen::MatrixXd &B);
Eigen::MatrixXd tensorproduct(const Eigen::VectorXd &A, const Eigen::VectorXd &B);

void tensorSelfProduct(const Eigen::VectorXd &A, Eigen::VectorXd &B);
void tensorSelfProduct(const Eigen::VectorXd &A, Eigen::MatrixXd &B);

double matrixNormInfinity(const Eigen::MatrixXd &M);
double matrixNorm2(const Eigen::MatrixXd &M);
double matrixNorm1(const Eigen::MatrixXd &M);

double matrixNormInfinity(const Eigen::VectorXd &M);
double matrixNorm2(const Eigen::VectorXd &M);
double matrixNorm1(const Eigen::VectorXd &M);

void applyFilter(double *out, double *in, const Eigen::MatrixXd &filter,
                 int kp1, int kp1_dm1, double fac);

void tensorExpandCoefs(int dim, int dir, int kp1, int kp1_d,
                       const Eigen::MatrixXd &primitive,
                       Eigen::VectorXd &expanded);

void tensorExpandCoords_2D(int kp1, const Eigen::MatrixXd &primitive, Eigen::MatrixXd &expanded);
void tensorExpandCoords_3D(int kp1, const Eigen::MatrixXd &primitive, Eigen::MatrixXd &expanded);

double calcDistance(int D, const double *a, const double *b);

} //namespace mwmath
} //namespace mrcpp
