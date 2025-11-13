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

#include "TimeEvolutionOperator.h"

#include "core/InterpolatingBasis.h"
#include "core/LegendreBasis.h"

#include "functions/GaussExp.h"
#include "functions/Gaussian.h"

#include "treebuilders/CrossCorrelationCalculator.h"
#include "treebuilders/DefaultCalculator.h"
#include "treebuilders/OperatorAdaptor.h"
#include "treebuilders/SplitAdaptor.h"
#include "treebuilders/TreeBuilder.h"
#include "treebuilders/grid.h"
#include "treebuilders/project.h"
#include "treebuilders/TimeEvolution_CrossCorrelationCalculator.h"

#include "trees/BandWidth.h"
#include "trees/CornerOperatorTree.h"
#include "trees/FunctionTreeVector.h"
#include "trees/OperatorNode.h"

#include "utils/Printer.h"
#include "utils/Timer.h"
#include "utils/math_utils.h"

#include <cmath>
#include <map>
#include <memory>

namespace mrcpp {

/* ========================= TimeEvolutionOperator ========================= */

template <int D>
TimeEvolutionOperator<D>::TimeEvolutionOperator(const MultiResolutionAnalysis<D> &mra,
                                                double prec,
                                                double time,
                                                int finest_scale,
                                                bool imaginary,
                                                int max_Jpower)
    : ConvolutionOperator<D>(mra, mra.getRootScale(), -10) {
    int oldlevel = Printer::setPrintLevel(0);
    this->setBuildPrec(prec);

    // Allocate on heap and keep the pointer (header currently uses raw pointer)
    // NOTE: this will be freed when/if you change the header to use unique_ptr
    this->cross_correlation = new SchrodingerEvolution_CrossCorrelation(
        30, mra.getOrder(), mra.getScalingBasis().getScalingType()
    );

    initialize(time, finest_scale, imaginary, max_Jpower);

    this->initOperExp(1);
    Printer::setPrintLevel(oldlevel);
}

template <int D>
TimeEvolutionOperator<D>::TimeEvolutionOperator(const MultiResolutionAnalysis<D> &mra,
                                                double prec,
                                                double time,
                                                bool imaginary,
                                                int max_Jpower)
    : ConvolutionOperator<D>(mra, mra.getRootScale(), -10) {
    int oldlevel = Printer::setPrintLevel(0);
    this->setBuildPrec(prec);

    this->cross_correlation = new SchrodingerEvolution_CrossCorrelation(
        30, mra.getOrder(), mra.getScalingBasis().getScalingType()
    );

    initialize(time, imaginary, max_Jpower);

    this->initOperExp(1);
    Printer::setPrintLevel(oldlevel);
}

template <int D>
void TimeEvolutionOperator<D>::initialize(double time, bool imaginary, int max_Jpower) {
    const int N = 18; // adaptive depth cap

    double o_prec = this->build_prec;
    auto   o_mra  = this->getOperatorMRA();

    auto o_tree = std::make_unique<CornerOperatorTree>(o_mra, o_prec);

    std::map<int, JpowerIntegrals*> J;
    for (int n = 0; n <= N + 1; ++n)
        J[n] = new JpowerIntegrals(time * std::pow(4.0, n), n, max_Jpower);

    TimeEvolution_CrossCorrelationCalculator calculator(J, this->cross_correlation, imaginary);

    OperatorAdaptor adaptor(o_prec, o_mra.getMaxScale(), true);
    TreeBuilder<2> builder;
    builder.build(*o_tree, calculator, adaptor, N);

    // Postprocess
    Timer trans_t;
    o_tree->mwTransform(BottomUp);
    o_tree->removeRoughScaleNoise();
    o_tree->calcSquareNorm();
    o_tree->setupOperNodeCache();
    print::time(10, "Time transform", trans_t);
    print::separator(10, ' ');

    this->raw_exp.push_back(std::move(o_tree));

    for (int n = 0; n <= N + 1; ++n) delete J[n];
}

template <int D>
void TimeEvolutionOperator<D>::initialize(double time, int finest_scale, bool imaginary, int max_Jpower) {
    double o_prec = this->build_prec;
    auto   o_mra  = this->getOperatorMRA();

    TreeBuilder<2> builder;
    SplitAdaptor<2> uniform(o_mra.getMaxScale(), true);

    const int N = finest_scale;
    const double threshold = o_prec / 1000.0;

    std::map<int, JpowerIntegrals*> J;
    for (int n = 0; n <= N + 1; ++n)
        J[n] = new JpowerIntegrals(time * std::pow(4.0, n), n, max_Jpower, threshold);

    TimeEvolution_CrossCorrelationCalculator calculator(J, this->cross_correlation, imaginary);

    auto o_tree = std::make_unique<CornerOperatorTree>(o_mra, o_prec);
    builder.build(*o_tree, calculator, uniform, N);

    // Postprocess
    Timer trans_t;
    o_tree->mwTransform(BottomUp);
    o_tree->calcSquareNorm();
    o_tree->setupOperNodeCache();
    print::time(10, "Time transform", trans_t);
    print::separator(10, ' ');

    this->raw_exp.push_back(std::move(o_tree));

    for (int n = 0; n <= N + 1; ++n) delete J[n];
}

/* ========================= Optional: SmoothDerivative =========================
   NOTE:
   The derivative calculator (and its MWNode/indices) is written for 2D.
   If you use it, keep it 2D-only, or generalize the calculator to ND.
   The implementation below is kept as in your version but with the safe
   TreeBuilder usage and without forcing invalid template instantiations.
*/

template <int D>
SmoothDerivative<D>::SmoothDerivative(const MultiResolutionAnalysis<D> &mra,
                                      double prec,
                                      double cut_off,
                                      int max_Jpower)
    : ConvolutionOperator<D>(mra, mra.getRootScale(), -10) {
    int oldlevel = Printer::setPrintLevel(0);
    this->setBuildPrec(prec);

    this->cross_correlation = new SchrodingerEvolution_CrossCorrelation(
        30, mra.getOrder(), mra.getScalingBasis().getScalingType()
    );

    initialize(cut_off, max_Jpower);

    this->initOperExp(1);
    Printer::setPrintLevel(oldlevel);
}

template <int D>
void SmoothDerivative<D>::initialize(double cut_off, int max_Jpower) {
    const int N = 18;

    double o_prec = this->build_prec;
    auto   o_mra  = this->getOperatorMRA();

    auto o_tree = std::make_unique<OperatorTree>(o_mra, o_prec);

    std::map<int, DerivativePowerIntegrals*> J;
    for (int n = 0; n <= N + 1; ++n)
        J[n] = new DerivativePowerIntegrals(cut_off, n, max_Jpower);

    DerivativeCrossCorrelationCalculator calculator(J, this->cross_correlation);
    OperatorAdaptor adaptor(o_prec, o_mra.getMaxScale(), true);

    TreeBuilder<2> builder;                   // calculators are 2D
    builder.build(*o_tree, calculator, adaptor, N);

    // Postprocess
    Timer trans_t;
    o_tree->mwTransform(BottomUp);
    o_tree->removeRoughScaleNoise();
    o_tree->calcSquareNorm();
    o_tree->setupOperNodeCache();
    print::time(10, "Time transform", trans_t);
    print::separator(10, ' ');

    this->raw_exp.push_back(std::move(o_tree));

    for (int n = 0; n <= N + 1; ++n) delete J[n];
}

/* Explicit instantiations actually used */
template class TimeEvolutionOperator<1>;
template class TimeEvolutionOperator<2>;
template class TimeEvolutionOperator<3>;

/* If you only use the derivative in 2D, instantiate only <2>. */
template class SmoothDerivative<2>;

} // namespace mrcpp