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
#include <memory>

namespace mrcpp {

/**
 * @class TimeEvolutionOperator
 * @brief Semigroup of the free-particle Schrodinger equation
 */
template <int D>
class TimeEvolutionOperator : public ConvolutionOperator<D> {
public:
    TimeEvolutionOperator(const MultiResolutionAnalysis<D> &mra,
                          double prec, double time, int finest_scale,
                          bool imaginary, int max_Jpower = 30);

    TimeEvolutionOperator(const MultiResolutionAnalysis<D> &mra,
                          double prec, double time,
                          bool imaginary, int max_Jpower = 30);

    TimeEvolutionOperator(const TimeEvolutionOperator &) = delete;
    TimeEvolutionOperator &operator=(const TimeEvolutionOperator &) = delete;
    virtual ~TimeEvolutionOperator() = default;

    double getBuildPrec() const { return this->build_prec; }

protected:
    void initialize(double time, int finest_scale, bool imaginary, int max_Jpower);
    void initialize(double time, bool imaginary, int max_Jpower);
    void initializeSemiUniformly(double time, bool imaginary, int max_Jpower);

    void setBuildPrec(double prec) { this->build_prec = prec; }

    double build_prec{-1.0};

    // Own the correlation tables to avoid dangling pointer bugs
    std::unique_ptr<SchrodingerEvolution_CrossCorrelation> cross_correlation_;
};

} // namespace mrcpp