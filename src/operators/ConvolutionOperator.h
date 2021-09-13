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
#include "trees/FunctionTreeVector.h"

namespace mrcpp {

template <int D> class ConvolutionOperator : public MWOperator {
public:
    ConvolutionOperator(const MultiResolutionAnalysis<D> &mra, double bprec, double kprec);
    ConvolutionOperator(const MultiResolutionAnalysis<D> &mra, double bprec, double kprec, int root, int reach);
    ConvolutionOperator(const ConvolutionOperator &oper) = delete;
    ConvolutionOperator &operator=(const ConvolutionOperator &oper) = delete;
    ~ConvolutionOperator() override;

protected:
    double kern_prec;
    double build_prec;
    MultiResolutionAnalysis<1> kern_mra;
    FunctionTreeVector<1> kern_exp;

    void initializeOperator(GreensKernel &greens_kernel);
    void clearKernel();

    double calcMinDistance(const MultiResolutionAnalysis<D> &MRA, double epsilon) const;
    double calcMaxDistance(const MultiResolutionAnalysis<D> &MRA) const;
};

} // namespace mrcpp
