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

/**
 * @file ABGVCalculator.cpp
 * @brief Local block assembly for the Alpert–Beylkin–Gines–Vozovoi derivative operator.
 *
 * @details
 * This module implements the calculator that fills the per-node matrix blocks for the
 * ABGV derivative operator used in multiresolution form. It is consumed by a
 * TreeBuilder to populate an OperatorTree with local stencil entries expressed in
 * the chosen scaling basis (interpolating or Legendre).
 *
 * The assembly depends on:
 * - the basis type and quadrature order,
 * - precomputed endpoint values of basis functions on the reference interval [0,1],
 * - a basis-dependent local derivative matrix K,
 * - two boundary weights A and B that select central, forward, backward, or
 *   semi-local differences.
 *
 * For each operator node the calculator determines the relative logical offset
 * between interacting cells. Only three cases produce non-zero local couplings:
 * left neighbor, same cell, and right neighbor. The four component blocks of the
 * 2-by-2 cell coupling are then filled accordingly, rescaled to the current level,
 * compressed to multiwavelet form, and cached with per-node norms.
 */

#include "ABGVCalculator.h"
#include "core/InterpolatingBasis.h"
#include "core/LegendreBasis.h"
#include "core/QuadratureCache.h"
#include "trees/MWNode.h"
#include "utils/Printer.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

namespace mrcpp {

/**
 * @brief Construct an ABGVCalculator and precompute basis-dependent tables.
 *
 * @param basis Scaling basis that defines quadrature order and function family.
 * @param a     Left boundary weight that controls semi-local coupling.
 * @param b     Right boundary weight that controls semi-local coupling.
 *
 * @details
 * The constructor allocates and fills:
 * - K: a kp1-by-kp1 local derivative matrix assembled on the reference cell,
 * - valueZero: endpoint values phi_i(0) for all basis indices,
 * - valueOne:  endpoint values phi_i(1) for all basis indices.
 *
 * The exact formulas are basis dependent and computed in calcKMatrix and
 * calcValueVectors respectively.
 */
ABGVCalculator::ABGVCalculator(const ScalingBasis &basis, double a, double b)
        : A(a)
        , B(b) {
    int kp1 = basis.getQuadratureOrder();
    this->K = MatrixXd::Zero(kp1, kp1);
    this->valueZero = VectorXd::Zero(kp1);
    this->valueOne = VectorXd::Zero(kp1);
    calcKMatrix(basis);
    calcValueVectors(basis);
}

/**
 * @brief Precompute endpoint values of scaling functions on [0, 1].
 *
 * @param basis Scaling basis.
 *
 * @details
 * - Interpolating basis: values are obtained by direct evaluation at 0 and 1.
 * - Legendre basis on [0, 1]: closed-form values are used.
 *   For index i we set
 *     valueOne(i)  = sqrt(2*i + 1),
 *     valueZero(i) = (-1)^i * sqrt(2*i + 1).
 */
void ABGVCalculator::calcValueVectors(const ScalingBasis &basis) {
    int kp1 = basis.getQuadratureOrder();
    double sqrtCoef[kp1];
    for (int i = 0; i < kp1; i++) { sqrtCoef[i] = std::sqrt(2.0 * i + 1.0); }

    switch (basis.getScalingType()) {
        case Interpol:
            for (int i = 0; i < kp1; i++) {
                this->valueZero(i) = basis.getFunc(i).evalf(0.0);
                this->valueOne(i) = basis.getFunc(i).evalf(1.0);
            }
            break;
        case Legendre:
            for (int i = 0; i < kp1; i++) {
                double val = sqrtCoef[i];
                this->valueOne(i) = val;
                if (IS_ODD(i)) val *= -1.0;
                this->valueZero(i) = val;
            }
            break;
        default:
            MSG_ERROR("Invalid scaling type");
            break;
    }
}

/**
 * @brief Assemble the local derivative matrix K on the reference cell.
 *
 * @param basis Scaling basis.
 *
 * @details
 * The construction of K depends on the basis family:
 * - Interpolating basis: K(i,j) = 2 * sqrt(w_j) * d(phi_i)/dx evaluated at x_j,
 *   where (x_j, w_j) are Gauss–Legendre quadrature nodes and weights provided
 *   by QuadratureCache. The factor 2 accounts for mapping from [-1,1] to [0,1].
 * - Legendre basis: a closed-form sparse pattern is used where K(j,i) is non-zero
 *   only if (i - j) is odd, in which case K(j,i) = 2 * sqrt(2i+1) * sqrt(2j+1).
 */
void ABGVCalculator::calcKMatrix(const ScalingBasis &basis) {
    int kp1 = basis.getQuadratureOrder();
    double sqrtCoef[kp1];
    for (int i = 0; i < kp1; i++) { sqrtCoef[i] = std::sqrt(2.0 * i + 1.0); }
    getQuadratureCache(qCache);
    const VectorXd &roots = qCache.getRoots(kp1);
    const VectorXd &weights = qCache.getWeights(kp1);
    VectorXd sqrtWeights = weights.array().sqrt();

    switch (basis.getScalingType()) {
        case Interpol:
            for (int i = 0; i < kp1; i++) {
                Polynomial poly = basis.getFunc(i).calcDerivative();
                for (int j = 0; j < kp1; j++) { this->K(i, j) = 2.0 * sqrtWeights(j) * poly.evalf(roots(j)); }
            }
            break;
        case Legendre:
            for (int i = 0; i < kp1; i++) {
                for (int j = i; j < kp1; j++) {
                    if (IS_ODD((i - j))) { this->K(j, i) = 2.0 * sqrtCoef[i] * sqrtCoef[j]; }
                }
            }
            break;
        default:
            MSG_ERROR("Invalid scaling type");
            break;
    }
}

/**
 * @brief Fill the local operator block for a given operator node and finalize it.
 *
 * @param node Operator node to be populated.
 *
 * @details
 * The node couples two 1D intervals at the same scale; its logical index encodes
 * which pair is assembled. Let l = idx[1] - idx[0]. Three cases are handled:
 *
 * - l =  0: intra-cell coupling. All four sub-blocks are filled using endpoint
 *   values and K, with boundary weights A and B selecting central or semi-local
 *   behavior.
 * - l = +1: right neighbor coupling. Only the block that mixes left and right
 *   components is filled, proportional to B.
 * - l = -1: left neighbor coupling. Only the symmetric block is filled,
 *   proportional to A.
 *
 * After filling, all entries are scaled by 2^(n+1) where n = idx.getScale() to
 * account for the derivative scaling at that level, then the node is transformed
 * with compression, marked as having coefficients, and its norms are computed.
 */
void ABGVCalculator::calcNode(MWNode<2> &node) {
    node.zeroCoefs();

    const auto &idx = node.getNodeIndex();
    int l = idx[1] - idx[0];
    int np1 = idx.getScale() + 1;
    int kp1 = node.getKp1();
    int kp1_d = node.getKp1_d();
    double two_np1 = std::pow(2.0, np1);
    double *coefs = node.getCoefs();

    double a = this->A;
    double b = this->B;
    switch (idx[1] - idx[0]) {
        case 1:
            if (b > MachineZero) {
                for (int i = 0; i < kp1; i++) {
                    double zero_i = this->valueZero(i);
                    for (int j = 0; j < kp1; j++) {
                        int idx = i * kp1 + j;
                        double one_j = this->valueOne(j);
                        coefs[1 * kp1_d + idx] = -b * zero_i * one_j;
                    }
                }
            }
            break;
        case 0:
            for (int i = 0; i < kp1; i++) {
                double zero_i = this->valueZero(i);
                double one_i = this->valueOne(i);
                for (int j = 0; j < kp1; j++) {
                    double K_ij = this->K(i, j);
                    double one_j = this->valueOne(j);
                    double one_ij = one_i * one_j;
                    double zero_j = this->valueZero(j);
                    double zero_ij = zero_i * zero_j;

                    int idx = i * kp1 + j;
                    coefs[0 * kp1_d + idx] = (1.0 - a) * one_ij - (1.0 - b) * zero_ij - K_ij;
                    coefs[1 * kp1_d + idx] = a * one_i * zero_j;
                    coefs[2 * kp1_d + idx] = -b * zero_i * one_j;
                    coefs[3 * kp1_d + idx] = (1.0 - a) * one_ij - (1.0 - b) * zero_ij - K_ij;
                }
            }
            break;
        case -1:
            if (a > MachineZero) {
                for (int i = 0; i < kp1; i++) {
                    double one_i = this->valueOne(i);
                    for (int j = 0; j < kp1; j++) {
                        int idx = i * kp1 + j;
                        double zero_j = this->valueZero(j);
                        coefs[2 * kp1_d + idx] = a * one_i * zero_j;
                    }
                }
            }
            break;
        default:
            // When periodic do nothing, else it should never end up here.
            break;
    }
    for (int i = 0; i < node.getNCoefs(); i++) { coefs[i] *= two_np1; }
    node.mwTransform(Compression);
    node.setHasCoefs();
    node.calcNorms();
}

} // namespace mrcpp
