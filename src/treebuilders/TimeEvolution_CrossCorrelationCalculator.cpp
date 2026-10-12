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

#include "TimeEvolution_CrossCorrelationCalculator.h"
#include "trees/FunctionTree.h"
#include "trees/MWNode.h"
#include "utils/Printer.h"

#include <array>
#include <cmath>
#include <vector>

using Eigen::MatrixXd;
using Eigen::VectorXcd;

namespace mrcpp {

namespace {

// Largest phase rate (|l| + 1) / 2a of the kernel, in radians per box, for which a
// block is computed from the power-integral series
constexpr double MaxSeriesFrequency = 1.0;

// n-point Gauss-Legendre rule on [0, 1]. GaussQuadrature stops at MaxGaussOrder = 42
// points, fewer than calcKernelCoefs uses for orders above 21.
void gauss_legendre_01(int n, std::vector<double> &x, std::vector<double> &w) {
    auto legendre = [n](double z, double &p_n, double &dp_n) {
        double p0 = 1.0, p1 = z;
        for (int j = 2; j <= n; j++) {
            double p2 = ((2.0 * j - 1.0) * z * p1 - (j - 1.0) * p0) / j;
            p0 = p1;
            p1 = p2;
        }
        p_n = p1;
        dp_n = n * (z * p1 - p0) / (z * z - 1.0);
    };
    x.resize(n);
    w.resize(n);
    for (int i = 0; i < (n + 1) / 2; i++) {
        double z = std::cos(pi * (i + 0.75) / (n + 0.5));
        double p_n, dp_n;
        for (int iter = 0; iter < 100; iter++) {
            legendre(z, p_n, dp_n);
            double dz = p_n / dp_n;
            z -= dz;
            if (std::abs(dz) < 1.0e-16) break;
        }
        legendre(z, p_n, dp_n);
        double weight = 1.0 / ((1.0 - z * z) * dp_n * dp_n); // 2 / ((1 - z^2) P_n'^2) on [-1, 1], halved
        x[i] = 0.5 * (1.0 - z);
        x[n - 1 - i] = 0.5 * (1.0 + z);
        w[i] = weight;
        w[n - 1 - i] = weight;
    }
}

} // namespace

TimeEvolution_CrossCorrelationCalculator::TimeEvolution_CrossCorrelationCalculator(double time, int max_Jpower, SchrodingerEvolution_CrossCorrelation *cross_correlation)
        : time(time)
        , max_Jpower(max_Jpower)
        , cross_correlation(cross_correlation) {
    // Loaded here, outside the parallel tree build
    const int order = cross_correlation->getOrder();
    getCrossCorrelationCache(Legendre, ccc);
    this->lMat = ccc.getLMatrix(order).cast<ComplexDouble>();
    this->rMat = ccc.getRMatrix(order).cast<ComplexDouble>();

    // order + 21 points per interval: calcKernelCoefs integrates Legendre polynomials up
    // to order 2k + 1 times at most 6 radians of the kernel's phase on each
    gauss_legendre_01(20 + order + 1, this->quadRoots, this->quadWeights);
}

/** @brief Legendre coefficients of the free-particle kernel on one box
 *
 * @param[in] m: translation of the box \f$ [m, m + 1] \f$, in units of the box length.
 * @param[in] a: dimensionless time \f$ t / h^2 \f$.
 *
 * @returns The coefficients of order \f$ 2k + 1 \f$ of
 * \f$ K_a(u) = (4 \pi i a)^{-1/2} \exp(i u^2 / 4a) \f$, normalized like the scaling
 * functions of a LegendreBasis, which is what the cross-correlation matrices expect.
 *
 * @details The phase of the kernel turns at most \f$ (|m| + 1) / 2a \f$ radians per unit
 * length on the box, so the box is split into intervals of at most 6 radians, each
 * integrated with the Gauss-Legendre rule set up in the constructor.
 */
VectorXcd TimeEvolution_CrossCorrelationCalculator::calcKernelCoefs(int m, double a) const {
    const int K = 2 * this->cross_correlation->getOrder() + 1;
    const ComplexDouble prefactor = 1.0 / std::sqrt(ComplexDouble(0.0, 4.0 * pi * a));
    const double max_frequency = (std::abs(m) + 1.0) / (2.0 * a);
    const int n_int = static_cast<int>(std::ceil(max_frequency / 6.0)) + 1;

    VectorXcd coefs = VectorXcd::Zero(K + 1);
    for (int s = 0; s < n_int; s++) {
        for (size_t q = 0; q < this->quadRoots.size(); q++) {
            const double x = (s + this->quadRoots[q]) / n_int;
            const double u = m + x;
            const ComplexDouble kernel = prefactor * std::exp(ComplexDouble(0.0, u * u / (4.0 * a))) * (this->quadWeights[q] / n_int);
            // Normalized Legendre polynomials sqrt(2i + 1) P_i(2x - 1) by recurrence
            const double y = 2.0 * x - 1.0;
            double p0 = 1.0, p1 = y;
            coefs(0) += kernel;
            coefs(1) += kernel * std::sqrt(3.0) * y;
            for (int i = 2; i <= K; i++) {
                double p2 = ((2.0 * i - 1.0) * y * p1 - (i - 1.0) * p0) / i;
                p0 = p1;
                p1 = p2;
                coefs(i) += kernel * std::sqrt(2.0 * i + 1.0) * p2;
            }
        }
    }
    return coefs;
}

/** @param[in] node: ...
 *  @details This will ... (work in progress)
 *
 *
 *
 *
 */
void TimeEvolution_CrossCorrelationCalculator::calcNode(MWNode<2, ComplexDouble> &node) {
    node.zeroCoefs();
    int type = node.getMWTree().getMRA().getScalingBasis().getScalingType();
    switch (type) {
        case Interpol: {
            MSG_ERROR("Invalid scaling type");
            break;
        }
        case Legendre: {
            applyCcc(node);
            break;
        }
        default:
            MSG_ERROR("Invalid scaling type");
            break;
    }
    node.mwTransform(Compression);
    node.setHasCoefs();
    node.calcNorms();
}

/** @param[in] node: ...
 *  @details This will ... (work in progress)
 *
 *
 *
 *
 */
// template <int T>
void TimeEvolution_CrossCorrelationCalculator::applyCcc(MWNode<2, ComplexDouble> &node) {
    // std::cout << node;
    //  The scale of J power integrals:
    // int scale = node.getScale() + 1;  //scale = n = (n - 1) + 1

    int t_dim = node.getTDim();  // t_dim = 4
    int kp1_d = node.getKp1_d(); // kp1_d = (k + 1)^2

    VectorXcd vec_o = VectorXcd::Zero(t_dim * kp1_d);
    const NodeIndex<2> &idx = node.getNodeIndex();

    // The power integrals live on the scale of the children, n + 1, and depend on the
    // dimensionless time a = t / h^2, where h = s 2^{-(n + 1)} is the length of a box
    // at that scale and s the scaling factor of the world box (a = t 4^{n + 1} for s = 1).
    const double sf = node.getMWTree().getMRA().getWorldBox().getScalingFactor(0);
    const double a = this->time * std::pow(4, node.getScale() + 1) / (sf * sf);

    // The four children only see the translations 2L - 1, 2L and 2L + 1, where L is
    // the translation of this node, so three matrix blocks cover them all.
    //
    // The power-integral series loses accuracy when the kernel turns fast across a box
    // and diverges in double precision once (|l| + 1) / 2a is above about 20. Above
    // MaxSeriesFrequency the block is instead built from the kernel's Legendre
    // coefficients on the two boxes it overlaps, as CrossCorrelationCalculator does
    // for Gaussian kernels.
    const int l_min = 2 * (idx[1] - idx[0]) - 1;
    const size_t n_matrices = this->cross_correlation->Matrix.size();
    std::array<VectorXcd, 3> blocks;
    std::array<VectorXcd, 4> kernel_coefs; // boxes l_min - 1, ..., l_min + 2, computed on demand
    for (int s = 0; s < 3; s++) {
        const int l_b = l_min + s;
        VectorXcd &block = blocks[s];
        block = VectorXcd::Zero(kp1_d);
        if ((std::abs(l_b) + 1.0) / (2.0 * a) <= MaxSeriesFrequency) {
            const auto J_l = JpowerIntegrals::calculate_J_power_integrals(l_b, a, this->max_Jpower);
            int block_index = 0;
            for (int p = 0; p <= node.getOrder(); p++)
                for (int j = 0; j <= node.getOrder(); j++) {
                    for (size_t k = 0; 2 * k + p + j < J_l.size() and k < n_matrices; k++) {
                        // The kernel is complex; the whole power integral is kept.
                        const ComplexDouble J = J_l[2 * k + p + j];
                        block(block_index) += J * cross_correlation->Matrix[k](p, j); // by default eigen library reads a transpose matrix from a file
                    }
                    block_index++;
                }
        } else {
            for (int m : {s, s + 1}) // boxes l_b - 1 and l_b
                if (kernel_coefs[m].size() == 0) kernel_coefs[m] = calcKernelCoefs(l_min - 1 + m, a);
            block = this->lMat * kernel_coefs[s] + this->rMat * kernel_coefs[s + 1];
        }
    }

    for (int i = 0; i < t_dim; i++) {
        NodeIndex<2> l = idx.child(i);
        int l_b = l[1] - l[0];
        vec_o.segment(i * kp1_d, kp1_d) = blocks[l_b - l_min];
    }

    ComplexDouble *coefs = node.getCoefs();
    for (int i = 0; i < t_dim * kp1_d; i++) {
        // auto scaling_factor = node.getMWTree().getMRA().getWorldBox().getScalingFactor(0);
        coefs[i] = vec_o(i);
        // std::cout<< "coefs[i] = " << coefs[i] << std::endl;
    }
}

} // namespace mrcpp
