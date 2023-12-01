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

#include "MWOperator.h"
#include "ConvolutionOperator.h"
#include "IdentityConvolution.h"
#include "core/SchrodingerEvolution_CrossCorrelation.h"

namespace mrcpp {


/**
 * @brief time evolition operator
 * 
 * 
 * 
 * 
 */
template <int D> class TimeEvolutionOperator : public ConvolutionOperator<D>   //One can use ConvolutionOperator instead as well
{
public:
    TimeEvolutionOperator(const MultiResolutionAnalysis<D> &mra, double prec, double time, int finest_scale, bool imaginary, int max_Jpower = 20);
    TimeEvolutionOperator(const TimeEvolutionOperator &oper) = delete;
    TimeEvolutionOperator &operator=(const TimeEvolutionOperator &oper) = delete;
    virtual ~TimeEvolutionOperator() = default;

    double getBuildPrec() const { return this->build_prec; }

//protected:
/*
    TimeEvolutionOperator(const MultiResolutionAnalysis<D> &mra)
        : MWOperator<D>(mra, mra.getRootScale(), -10) {}
    TimeEvolutionOperator(const MultiResolutionAnalysis<D> &mra, int root, int reach)
        : MWOperator<D>(mra, root, reach) {}
*/

    void initialize(double time, int finest_scale, bool imaginary, int max_Jpower);
    void setBuildPrec(double prec) { this->build_prec = prec; }

    double build_prec{-1.0};
    SchrodingerEvolution_CrossCorrelation *cross_correlation{nullptr};
    
    //Difficult to impliment:
    //int crop(double prec, double splitFac = 1.0, bool absPrec = true);
};

} // namespace mrcpp
