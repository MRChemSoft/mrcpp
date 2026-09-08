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
    /// @brief Pass as `finest_scale` to refine adaptively instead of uniformly.
    static constexpr int Adaptive = -1;

    /** @brief Semigroup \f$ \exp(i t \partial_x^2) \f$.
     *
     * @param[in] finest_scale: uniform refinement down to this scale;
     * `Adaptive` (the default) refines adaptively.
     * @param[in] max_Jpower: number of power integrals used in the expansion.
     * Construction cost is linear in this and there is no negligibility
     * truncation, so the series is always computed in full. Convergence is
     * slower for *small* time steps, since the recurrence coefficients scale
     * as 1/(time*4^n): measured at order 4 and prec 1e-7, six powers suffice
     * for time steps down to 1e-3, and the default leaves ample margin.
     *
     * @note Applying this to a real `FunctionTree` aborts; project the input as
     * `ComplexDouble` first, or let the `CompFunction` overload promote it.
     */
    TimeEvolutionOperator(const MultiResolutionAnalysis<D> &mra, double prec, double time, int finest_scale = Adaptive, int max_Jpower = 30);

    /// @brief Rejects the legacy signatures carrying the removed `imaginary` argument.
    TimeEvolutionOperator(const MultiResolutionAnalysis<D> &mra, double prec, double time, bool imaginary, int max_Jpower = 30) = delete;
    TimeEvolutionOperator(const MultiResolutionAnalysis<D> &mra, double prec, double time, int finest_scale, bool imaginary, int max_Jpower = 30) = delete;
    TimeEvolutionOperator(const TimeEvolutionOperator &oper) = delete;
    TimeEvolutionOperator &operator=(const TimeEvolutionOperator &oper) = delete;
    virtual ~TimeEvolutionOperator() = default;

    double getBuildPrec() const { return this->build_prec; }

protected:
    void initialize(double time, int finest_scale, int max_Jpower);
    void initialize(double time, int max_Jpower);
    void initializeSemiUniformly(double time, int max_Jpower);

    void setBuildPrec(double prec) { this->build_prec = prec; }

    double build_prec{-1.0};
    SchrodingerEvolution_CrossCorrelation *cross_correlation{nullptr};
};

} // namespace mrcpp
