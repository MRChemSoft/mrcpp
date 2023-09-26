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
class TimeEvolution_CrossCorrelationCalculator final : public TreeCalculator<2>
{
public:
    TimeEvolution_CrossCorrelationCalculator
    (JpowerIntegrals &J, SchrodingerEvolution_CrossCorrelation *cross_correlation, bool imaginary)
        : J_power_inetgarls(&J), cross_correlation(cross_correlation), imaginary(imaginary){}
//private:
    JpowerIntegrals *J_power_inetgarls;
    SchrodingerEvolution_CrossCorrelation *cross_correlation;
    
    /// @brief If False then the calculator is using th real part of integrals, otherwise - the imaginary part.
    bool imaginary;

    void calcNode(MWNode<2> &node) override;

    //template <int T>
    void applyCcc(MWNode<2> &node);
    //template <int T> void applyCcc(MWNode<2> &node, CrossCorrelationCache<T> &ccc);
};

} // namespace mrcpp
