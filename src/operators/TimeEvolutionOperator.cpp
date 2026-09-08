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
//#include "MRCPP/MWOperators"

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

#include "trees/BandWidth.h"
#include "trees/CornerOperatorTree.h"
#include "trees/FunctionTreeVector.h"

#include "utils/Printer.h"
#include "utils/Timer.h"
#include "utils/math_utils.h"

#include "treebuilders/TimeEvolution_CrossCorrelationCalculator.h"

#include <vector>

#include "trees/OperatorNode.h"

namespace mrcpp {

/** @brief An adaptive constructor for TimeEvolutionOperator class.
 *
 * @param[in] mra: MRA.
 * @param[in] prec: precision.
 * @param[in] time: the time moment (step).
 *
 * @details Constructs the complete complex Schrodinger semigroup at a given time moment.
 *
 * @note For technical reasons the operator tree is constructed no deeper than to scale \f$ n = 18 \f$.
 */
template <int D>
TimeEvolutionOperator<D>::TimeEvolutionOperator(const MultiResolutionAnalysis<D> &mra, double prec, double time, int finest_scale, int max_Jpower)
        : ConvolutionOperator<D>(mra, mra.getRootScale(), -10) {
    if (max_Jpower <= 0) MSG_ABORT("max_Jpower must be positive");
    if (finest_scale != Adaptive and finest_scale < mra.getRootScale()) MSG_ABORT("finest_scale is above the root scale");

    int oldlevel = Printer::setPrintLevel(0);
    this->setBuildPrec(prec);

    // The first argument counts cross-correlation matrices, not power integrals:
    // applyCcc consumes one matrix per two power-integral orders, so 30 covers
    // max_Jpower up to 60. It is deliberately not tied to max_Jpower.
    SchrodingerEvolution_CrossCorrelation cross_correlation(30, mra.getOrder(), mra.getScalingBasis().getScalingType());
    this->cross_correlation = &cross_correlation;

    // Adaptive is the sentinel for "no uniform scale given".
    if (finest_scale == Adaptive) {
        initialize(time, max_Jpower);
    } else {
        initialize(time, finest_scale, max_Jpower);
    }
    this->cross_correlation = nullptr; // the object above dies with this scope

    this->initOperExp(1); // one separable term
    Printer::setPrintLevel(oldlevel);
}
/** @brief Creates the complex operator
 *
 * @details Adaptive down to scale \f$ N = 18 \f$.
 * This scale limit bounds the amount of JpowerIntegrals
 * to be calculated.
 * @note In future work we plan to optimize calculation of JpowerIntegrals so that we calculate
 * only needed ones, while building the tree (in progress).
 *
 */
template <int D> void TimeEvolutionOperator<D>::initialize(double time, int max_Jpower) {
    int N = 18;

    double o_prec = this->build_prec;
    auto o_mra = this->getOperatorMRA();
    auto o_tree = std::make_unique<CornerOperatorTree<ComplexDouble>>(o_mra, o_prec);

    std::map<int, JpowerIntegrals *> J;
    for (int n = 0; n <= N + 1; n++) J[n] = new JpowerIntegrals(time * std::pow(4, n), n, max_Jpower);
    TimeEvolution_CrossCorrelationCalculator calculator(J, this->cross_correlation);

    OperatorAdaptor<ComplexDouble> adaptor(o_prec, o_mra.getMaxScale(), true);

    mrcpp::TreeBuilder<2, ComplexDouble> builder;
    builder.build(*o_tree, calculator, adaptor, N);

    // Postprocess to make the operator functional
    Timer trans_t;
    o_tree->mwTransform(BottomUp);
    o_tree->removeRoughScaleNoise();
    // o_tree->clearSquareNorm(); //does not affect printing
    o_tree->calcSquareNorm();
    o_tree->setupOperNodeCache();

    print::time(10, "Time transform", trans_t);
    print::separator(10, ' ');

    this->raw_exp_cplx.push_back(std::move(o_tree));

    for (int n = 0; n <= N + 1; n++) delete J[n];
}

/** @brief Creates the complex operator
 *
 * @details Uniform down to finest scale.
 *
 */
template <int D> void TimeEvolutionOperator<D>::initialize(double time, int finest_scale, int max_Jpower) {
    double o_prec = this->build_prec;
    auto o_mra = this->getOperatorMRA();

    // Setup uniform tree builder
    TreeBuilder<2, ComplexDouble> builder;
    SplitAdaptor<2, ComplexDouble> uniform(o_mra.getMaxScale(), true);

    int N = finest_scale;
    double threshold = o_prec / 1000.0;
    std::map<int, JpowerIntegrals *> J;
    for (int n = 0; n <= N + 1; n++) J[n] = new JpowerIntegrals(time * std::pow(4, n), n, max_Jpower, threshold);
    TimeEvolution_CrossCorrelationCalculator calculator(J, this->cross_correlation);

    auto o_tree = std::make_unique<CornerOperatorTree<ComplexDouble>>(o_mra, o_prec);
    builder.build(*o_tree, calculator, uniform, N); // Expand 1D kernel into 2D operator

    // Postprocess to make the operator functional
    Timer trans_t;
    o_tree->mwTransform(BottomUp);
    o_tree->calcSquareNorm();
    o_tree->setupOperNodeCache();
    print::time(10, "Time transform", trans_t);
    print::separator(10, ' ');

    this->raw_exp_cplx.push_back(std::move(o_tree));

    for (int n = 0; n <= N + 1; n++) delete J[n];
}

/** @brief Creates the operator (in progress)
 *
 * @details Tree construction starts uniformly and then continues adaptively down to scale \f$ N = 18 \f$.
 * This scale limit bounds the amount of JpowerIntegrals
 * to be calculated.
 * @note This method is not ready for use and should not be used (in progress).
 *
 */
template <int D> void TimeEvolutionOperator<D>::initializeSemiUniformly(double time, int max_Jpower) {
    MSG_ABORT("Not implemented");

    double o_prec = this->build_prec;
    auto o_mra = this->getOperatorMRA();

    mrcpp::TreeBuilder<2, ComplexDouble> builder;
    mrcpp::SplitAdaptor<2, ComplexDouble> uniform(o_mra.getMaxScale(), true);

    int N = 18;

    auto o_tree = std::make_unique<CornerOperatorTree<ComplexDouble>>(o_mra, o_prec);
    DefaultCalculator<2, ComplexDouble> intitial_calculator;
    for (auto n = 0; n < 4; n++) builder.build(*o_tree, intitial_calculator, uniform, 1);

    double threshold = o_prec / 1000.0;
    std::map<int, mrcpp::JpowerIntegrals *> J;
    for (int n = 0; n <= N + 1; n++) J[n] = new mrcpp::JpowerIntegrals(time * std::pow(4, n), n, max_Jpower, threshold);
    mrcpp::TimeEvolution_CrossCorrelationCalculator calculator(J, this->cross_correlation);

    OperatorAdaptor<ComplexDouble> adaptor(o_prec, o_mra.getMaxScale());
    builder.build(*o_tree, calculator, adaptor, 13);

    // Postprocess to make the operator functional
    Timer trans_t;
    o_tree->mwTransform(mrcpp::BottomUp);
    o_tree->removeRoughScaleNoise();
    o_tree->calcSquareNorm();
    o_tree->setupOperNodeCache();
    print::time(10, "Time transform", trans_t);
    print::separator(10, ' ');

    this->raw_exp_cplx.push_back(std::move(o_tree));

    for (int n = 0; n <= N + 1; n++) delete J[n];
}

template class TimeEvolutionOperator<1>;
template class TimeEvolutionOperator<2>;
template class TimeEvolutionOperator<3>;

} // namespace mrcpp
