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

namespace mrcpp {

/** @class ConvolutionOperator
 *
 * @brief Convolution defined by a Gaussian expansion
 * 
 * @details Represents the operator
 * \f[
 * 	    T = \sum_{m=1}^M
 *  	\text{sign} (\alpha_m) \bigotimes \limits_{d = 1}^D T_d
 *      \left( \beta_m, \sqrt[D]{| \alpha_m |} \right)
 *      ,
 * \f]
 * where each
 * \f$ T_d \left( \beta, \alpha \right) \f$
 * is the convolution operator with one-dimensional Gaussian kernel
 * \f$ k(x_d) = \alpha \exp \left( - \beta x_d^2 \right) \f$.
 * Operator
 * \f$ T \f$
 * is obtained from the Gaussian expansion
 * \f[
 *      \sum_{m=1}^M \alpha_m \exp \left( - \beta_m |x|^2 \right)
 * \f]
 * which is passed as a parameter to the first two constructors.
 *  
 * @note Every \f$ T_d \left( \beta_m, \sqrt[D]{| \alpha_m |} \right) \f$ is the same
 * operator associated with the one-dimensional variable \f$ x_d \f$ for \f$ d = 1, \ldots, D \f$.
 * 
 * \todo: One may want to change the logic so that \f$ D \f$-root is evaluated on the previous step,
 * namely, when \f$ \alpha_m, \beta_m \f$ are calculated.
 * 
 */
template <int D> class ConvolutionOperator : public MWOperator<D> {
public:
    ConvolutionOperator(const MultiResolutionAnalysis<D> &mra, GaussExp<1> &kernel, double prec);
    ConvolutionOperator(const MultiResolutionAnalysis<D> &mra, GaussExp<1> &kernel, double prec, int root, int reach);
    ConvolutionOperator(const ConvolutionOperator &oper) = delete;
    ConvolutionOperator &operator=(const ConvolutionOperator &oper) = delete;
    virtual ~ConvolutionOperator() = default;

    double getBuildPrec() const { return this->build_prec; }

protected:
    ConvolutionOperator(const MultiResolutionAnalysis<D> &mra)
        : MWOperator<D>(mra, mra.getRootScale(), -10) {}
    ConvolutionOperator(const MultiResolutionAnalysis<D> &mra, int root, int reach)
        : MWOperator<D>(mra, root, reach) {}

    void initialize(GaussExp<1> &kernel, double k_prec, double o_prec);
    void setBuildPrec(double prec) { this->build_prec = prec; }

    MultiResolutionAnalysis<1> getKernelMRA() const;

    double build_prec{-1.0};
};

} // namespace mrcpp
