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

#include <map>

#include "TreeCalculator.h"
#include "core/CrossCorrelationCache.h"
#include "core/SchrodingerEvolution_CrossCorrelation.h"
#include "functions/JpowerIntegrals.h"

namespace mrcpp {

/** @class TimeEvolution_CrossCorrelationCalculator
 *
 * @brief Calculator for time-evolution cross-correlation nodes (work in progress).
 */
class TimeEvolution_CrossCorrelationCalculator final : public TreeCalculator<2> {
public:
    TimeEvolution_CrossCorrelationCalculator(std::map<int, JpowerIntegrals *> &J,
                                             SchrodingerEvolution_CrossCorrelation *cross_correlation,
                                             bool imaginary)
            : J_power_inetgarls(J)
            , cross_correlation(cross_correlation)
            , imaginary(imaginary) {}

    void calcNode(MWNode<2> &node) override;
    void applyCcc(MWNode<2> &node);

    // Inputs
    std::map<int, JpowerIntegrals *> J_power_inetgarls;
    SchrodingerEvolution_CrossCorrelation *cross_correlation;

    /// If false → use Re(J); if true → use Im(J).
    bool imaginary;
};

/** @class DerivativeCrossCorrelationCalculator
 *
 * @brief Calculator for derivative-operator cross-correlation nodes (work in progress).
 *
 * Note: Uses DerivativePowerIntegrals (real-valued table).
 */
class DerivativeCrossCorrelationCalculator final : public TreeCalculator<2> {
public:
    DerivativeCrossCorrelationCalculator(std::map<int, DerivativePowerIntegrals *> &J,
                                         SchrodingerEvolution_CrossCorrelation *cross_correlation)
            : J_power_inetgarls(J)
            , cross_correlation(cross_correlation) {}

    void calcNode(MWNode<2> &node) override;
    void applyCcc(MWNode<2> &node);

    // Inputs
    std::map<int, DerivativePowerIntegrals *> J_power_inetgarls;
    SchrodingerEvolution_CrossCorrelation *cross_correlation;
};

} // namespace mrcpp