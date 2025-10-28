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

/**
 * @file TimeEvolutionOperator.cpp
 * @brief Construction of (real/imaginary) parts of the Schrödinger time-evolution
 *        operator in the multiwavelet framework.
 *
 * The implementation builds a separable, multi-resolution representation of the
 * free-particle time-evolution semigroup
 * \f[
 *   U(t) = e^{\, i t \Delta}
 * \f]
 * (or its real/imaginary part), using cross-correlation calculators and
 * precomputed power integrals \f$ \widetilde J_m \f$ (see @ref JpowerIntegrals).
 * Two build modes are provided:
 *  - **Adaptive** down to a fixed scale \f$N=18\f$, bounding the number of
 *    power integrals.
 *  - **Uniform** down to a user-specified finest scale.
 *
 * Assembly follows the standard operator pipeline:
 * projection/lifting → multiwavelet transform → cache/init of operator blocks.
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

/**
 * @brief Uniform constructor.
 *
 * @tparam D Spatial dimension (1, 2, or 3).
 * @param mra            Target @ref MultiResolutionAnalysis defining domain and basis.
 * @param prec           Build precision for assembly and pruning.
 * @param time           Time parameter \f$ t \f$ of the semigroup.
 * @param finest_scale   Uniform build depth (finest level) of the operator tree.
 * @param imaginary      If `true`, build the imaginary part; otherwise, the real part.
 * @param max_Jpower     Maximum number of power-integral terms \f$ \widetilde J_m \f$ to retain.
 *
 * @details
 * Builds a **uniform** operator down to @p finest_scale. Internally sets up a
 * @ref SchrodingerEvolution_CrossCorrelation calculator and calls
 * the uniform @ref initialize(double,int,bool,int) overload. The operator
 * expansion is finalized via @ref initOperExp(1).
 */
template <int D>
TimeEvolutionOperator<D>::TimeEvolutionOperator(const MultiResolutionAnalysis<D> &mra,
                                                double prec,
                                                double time,
                                                int finest_scale,
                                                bool imaginary,
                                                int max_Jpower)
        : ConvolutionOperator<D>(mra, mra.getRootScale(), -10) // One can use ConvolutionOperator instead as well
{
    int oldlevel = Printer::setPrintLevel(0);
    this->setBuildPrec(prec);

    SchrodingerEvolution_CrossCorrelation cross_correlation(30, mra.getOrder(), mra.getScalingBasis().getScalingType());
    this->cross_correlation = &cross_correlation;

    initialize(time, finest_scale, imaginary, max_Jpower); // will go outside of the constructor in future

    this->initOperExp(1); // Important to finalize component mapping
    Printer::setPrintLevel(oldlevel);
}

/**
 * @brief Adaptive constructor.
 *
 * @tparam D Spatial dimension (1, 2, or 3).
 * @param mra         Target @ref MultiResolutionAnalysis.
 * @param prec        Build precision.
 * @param time        Time parameter \f$ t \f$ of the semigroup.
 * @param imaginary   If `true`, build the imaginary part; otherwise, the real part.
 * @param max_Jpower  Maximum number of power-integral terms \f$ \widetilde J_m \f$ to retain.
 *
 * @details
 * Builds an **adaptive** operator down to a fixed scale \f$N=18\f$, which keeps the number
 * of necessary power integrals bounded. The assembly uses a
 * @ref TimeEvolution_CrossCorrelationCalculator fed by a per-scale
 * map of @ref JpowerIntegrals.
 */
template <int D>
TimeEvolutionOperator<D>::TimeEvolutionOperator(const MultiResolutionAnalysis<D> &mra,
                                                double prec,
                                                double time,
                                                bool imaginary,
                                                int max_Jpower)
        : ConvolutionOperator<D>(mra, mra.getRootScale(), -10) // One can use ConvolutionOperator instead as well
{
    int oldlevel = Printer::setPrintLevel(0);
    this->setBuildPrec(prec);

    SchrodingerEvolution_CrossCorrelation cross_correlation(30, mra.getOrder(), mra.getScalingBasis().getScalingType());
    this->cross_correlation = &cross_correlation;

    initialize(time, imaginary, max_Jpower); // will go outside of the constructor in future

    this->initOperExp(1); // Important to finalize component mapping
    Printer::setPrintLevel(oldlevel);
}

/**
 * @brief Adaptive build: create real or imaginary part of the operator.
 *
 * @param time        Time parameter \f$ t \f$.
 * @param imaginary   If `true`, build the imaginary part; otherwise, the real part.
 * @param max_Jpower  Maximum number of power-integral terms \f$ \widetilde J_m \f$ per scale.
 *
 * @details
 * Builds **adaptively** down to scale \f$ N = 18 \f$. For each scale
 * \f$ n=0,\dots,N+1 \f$ a corresponding @ref JpowerIntegrals object is created
 * with parameter \f$ a = t\,4^n \f$. The operator is assembled using a
 * @ref TimeEvolution_CrossCorrelationCalculator and finalized by multiwavelet
 * transform, rough-scale noise removal, square-norm evaluation, and cache setup.
 *
 * @note The fixed depth ensures a bounded number of power integrals while building.
 * Future work aims to compute only the power integrals actually needed during build.
 */
template <int D>
void TimeEvolutionOperator<D>::initialize(double time, bool imaginary, int max_Jpower) {
    int N = 18;

    double o_prec = this->build_prec;
    auto o_mra = this->getOperatorMRA();
    auto o_tree = std::make_unique<CornerOperatorTree>(o_mra, o_prec);

    std::map<int, JpowerIntegrals *> J;
    for (int n = 0; n <= N + 1; n++) J[n] = new JpowerIntegrals(time * std::pow(4, n), n, max_Jpower);
    TimeEvolution_CrossCorrelationCalculator calculator(J, this->cross_correlation, imaginary);

    OperatorAdaptor adaptor(o_prec, o_mra.getMaxScale(), true);

    mrcpp::TreeBuilder<2> builder;
    builder.build(*o_tree, calculator, adaptor, N);

    // Postprocess to make the operator functional
    Timer trans_t;
    o_tree->mwTransform(BottomUp);
    o_tree->removeRoughScaleNoise();
    // o_tree->clearSquareNorm(); // does not affect printing
    o_tree->calcSquareNorm();
    o_tree->setupOperNodeCache();

    print::time(10, "Time transform", trans_t);
    print::separator(10, ' ');

    this->raw_exp.push_back(std::move(o_tree));

    for (int n = 0; n <= N + 1; n++) delete J[n];
}

/**
 * @brief Uniform build: create real or imaginary part of the operator.
 *
 * @param time          Time parameter \f$ t \f$.
 * @param finest_scale  Finest (uniform) scale to which the operator tree is constructed.
 * @param imaginary     If `true`, build the imaginary part; otherwise, the real part.
 * @param max_Jpower    Maximum number of power-integral terms \f$ \widetilde J_m \f$ per scale.
 *
 * @details
 * Builds **uniformly** down to @p finest_scale using a @ref SplitAdaptor.
 * A threshold of \f$ \text{prec}/1000 \f$ is used while creating
 * @ref JpowerIntegrals for scales \f$ n=0,\dots,N+1 \f$ with
 * \f$ a = t\,4^n \f$. The resulting @ref CornerOperatorTree is then transformed,
 * squared-normed, and cached for later application.
 */
template <int D>
void TimeEvolutionOperator<D>::initialize(double time, int finest_scale, bool imaginary, int max_Jpower) {
    double o_prec = this->build_prec;
    auto o_mra = this->getOperatorMRA();

    // Setup uniform tree builder
    TreeBuilder<2> builder;
    SplitAdaptor<2> uniform(o_mra.getMaxScale(), true);

    int N = finest_scale;
    double threshold = o_prec / 1000.0;
    std::map<int, JpowerIntegrals *> J;
    for (int n = 0; n <= N + 1; n++) J[n] = new JpowerIntegrals(time * std::pow(4, n), n, max_Jpower, threshold);
    TimeEvolution_CrossCorrelationCalculator calculator(J, this->cross_correlation, imaginary);

    auto o_tree = std::make_unique<CornerOperatorTree>(o_mra, o_prec);
    builder.build(*o_tree, calculator, uniform, N); // Expand 1D kernel into 2D operator

    // Postprocess to make the operator functional
    Timer trans_t;
    o_tree->mwTransform(BottomUp);
    o_tree->calcSquareNorm();
    o_tree->setupOperNodeCache();
    print::time(10, "Time transform", trans_t);
    print::separator(10, ' ');

    this->raw_exp.push_back(std::move(o_tree));

    for (int n = 0; n <= N + 1; n++) delete J[n];
}

/**
 * @brief Semi-uniform build (prototype; not ready for production).
 *
 * @param time        Time parameter \f$ t \f$.
 * @param imaginary   If `true`, build the imaginary part; otherwise, the real part.
 * @param max_Jpower  Maximum number of power-integral terms \f$ \widetilde J_m \f$ per scale.
 *
 * @details
 * Starts with a small uniform prefix of the operator tree and continues adaptively
 * down to \f$ N = 18 \f$. **Not implemented**—kept as a placeholder for future work.
 *
 * @warning This method deliberately aborts at runtime.
 */
template <int D>
void TimeEvolutionOperator<D>::initializeSemiUniformly(double time, bool imaginary, int max_Jpower) {
    MSG_ERROR("Not implemented yet method.");

    double o_prec = this->build_prec;
    auto o_mra = this->getOperatorMRA();

    mrcpp::TreeBuilder<2> builder;
    mrcpp::SplitAdaptor<2> uniform(o_mra.getMaxScale(), true);

    int N = 18;

    auto o_tree = std::make_unique<CornerOperatorTree>(o_mra, o_prec);
    DefaultCalculator<2> intitial_calculator;
    for (auto n = 0; n < 4; n++) builder.build(*o_tree, intitial_calculator, uniform, 1);

    double threshold = o_prec / 1000.0;
    std::map<int, mrcpp::JpowerIntegrals *> J;
    for (int n = 0; n <= N + 1; n++) J[n] = new mrcpp::JpowerIntegrals(time * std::pow(4, n), n, max_Jpower, threshold);
    mrcpp::TimeEvolution_CrossCorrelationCalculator calculator(J, this->cross_correlation, imaginary);

    OperatorAdaptor adaptor(o_prec, o_mra.getMaxScale());
    builder.build(*o_tree, calculator, adaptor, 13);

    // Postprocess to make the operator functional
    Timer trans_t;
    o_tree->mwTransform(mrcpp::BottomUp);
    o_tree->removeRoughScaleNoise();
    o_tree->calcSquareNorm();
    o_tree->setupOperNodeCache();
    print::time(10, "Time transform", trans_t);
    print::separator(10, ' ');

    this->raw_exp.push_back(std::move(o_tree));

    for (int n = 0; n <= N + 1; n++) delete J[n];
}

/* Explicit template instantiations */
template class TimeEvolutionOperator<1>;
template class TimeEvolutionOperator<2>;
template class TimeEvolutionOperator<3>;

} // namespace mrcpp
