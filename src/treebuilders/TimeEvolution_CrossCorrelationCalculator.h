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

#include "TreeCalculator.h"
#include "core/CrossCorrelationCache.h"
#include "core/SchrodingerEvolution_CrossCorrelation.h"
#include "functions/JpowerIntegrals.h"

#include <vector>

namespace mrcpp {

/** @class TimeEvolution_CrossCorrelationCalculator
 *
 * @brief An efficient way to calculate ... (work in progress)
 *
 * @details An efficient way to calculate ... having the form
 * \f$ \ldots = \ldots \f$
 *
 *
 *
 */
class TimeEvolution_CrossCorrelationCalculator final : public TreeCalculator<2, ComplexDouble> {
public:
    /// @param[in] time: the time step \f$ t \f$ of the semigroup \f$ \exp(i t \partial_x^2) \f$.
    /// @param[in] max_Jpower: highest power integral \f$ \widetilde J_m \f$ kept in the expansion.
    /// @details The power integrals are computed node by node, only for the translations
    /// each node needs.
    TimeEvolution_CrossCorrelationCalculator(double time, int max_Jpower, SchrodingerEvolution_CrossCorrelation *cross_correlation);
    // private:
    double time;
    int max_Jpower;
    SchrodingerEvolution_CrossCorrelation *cross_correlation;

    void calcNode(MWNode<2, ComplexDouble> &node) override;

    // template <int T>
    void applyCcc(MWNode<2, ComplexDouble> &node);

private:
    Eigen::MatrixXcd lMat;           ///< Cross-correlation matrices of the scaling basis,
    Eigen::MatrixXcd rMat;           ///< as used for Gaussian kernels in CrossCorrelationCalculator
    std::vector<double> quadRoots;   ///< Gauss-Legendre rule on [0, 1] used for
    std::vector<double> quadWeights; ///< the kernel's Legendre coefficients

    Eigen::VectorXcd calcKernelCoefs(int m, double a) const;
    // template <int T> void applyCcc(MWNode<2, ComplexDouble> &node, CrossCorrelationCache<T> &ccc);
};

} // namespace mrcpp
