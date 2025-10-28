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

#include "ConvolutionOperator.h"
#include "MWOperator.h"
#include "core/SchrodingerEvolution_CrossCorrelation.h"

namespace mrcpp {

/** @class TimeEvolutionOperator
 *
 * @brief Semigroup of the free-particle Schrodinger equation
 *
 * @details Represents the semigroup
 * \f$
 *      \exp \left( i t \partial_x^2 \right)
 *      .
 * \f$
 * Matrix elements (actual operator tree) of the operator can be obtained by calling getComponent(0, 0).
 *
 * @note So far implementation is done for Legendre scaling functions in 1d.
 *
 * \todo: Extend to D dimensinal on a general interval [a, b] in the future.
 *
 */
template <int D>
class TimeEvolutionOperator : public ConvolutionOperator<D> // One can use ConvolutionOperator instead as well
{
public:
    TimeEvolutionOperator(const MultiResolutionAnalysis<D> &mra, double prec, double time, int finest_scale, bool imaginary, int max_Jpower = 30);
    TimeEvolutionOperator(const MultiResolutionAnalysis<D> &mra, double prec, double time, bool imaginary, int max_Jpower = 30);
    TimeEvolutionOperator(const TimeEvolutionOperator &oper) = delete;
    TimeEvolutionOperator &operator=(const TimeEvolutionOperator &oper) = delete;
    virtual ~TimeEvolutionOperator() = default;

    double getBuildPrec() const { return this->build_prec; }

protected:
    void initialize(double time, int finest_scale, bool imaginary, int max_Jpower);
    void initialize(double time, bool imaginary, int max_Jpower);
    void initializeSemiUniformly(double time, bool imaginary, int max_Jpower);

    void setBuildPrec(double prec) { this->build_prec = prec; }

    double build_prec{-1.0};
    SchrodingerEvolution_CrossCorrelation *cross_correlation{nullptr};
};


/** @class SmoothDerivative
 *
 * @brief Derivative with cutt off high frequencies
 *
 * @details Represents the convolution with the kernel
 * \f$
 *      \mathcal F^{-1} \left( i \xi \chi_R(\xi) \right)
 *      .
 * \f$
 * Matrix elements (actual operator tree) of the operator can be obtained by calling getComponent(0, 0).
 *
 * @note So far implementation is done for Legendre scaling functions in 1d.
 *
 * \todo: Extend to D dimensinal on a general interval [a, b] in the future.
 *
 */
template <int D>
class SmoothDerivative : public ConvolutionOperator<D> // One can use ConvolutionOperator instead as well
{
public:
//    SmoothDerivative(const MultiResolutionAnalysis<D> &mra, double prec, double time, int finest_scale, bool imaginary, int max_Jpower = 30);
    SmoothDerivative(const MultiResolutionAnalysis<D> &mra, double prec, double cut_off, int max_Jpower = 30);
    SmoothDerivative(const SmoothDerivative &oper) = delete;
    SmoothDerivative &operator=(const SmoothDerivative &oper) = delete;
    virtual ~SmoothDerivative() = default;

    double getBuildPrec() const { return this->build_prec; }

protected:
    void initialize(double cut_off, int max_Jpower);
    
    void setBuildPrec(double prec) { this->build_prec = prec; }

    double build_prec{-1.0};
    SchrodingerEvolution_CrossCorrelation *cross_correlation{nullptr};
};

} // namespace mrcpp
