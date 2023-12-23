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

#include "functions/Gaussian.h"
#include "functions/GaussExp.h"

#include "treebuilders/CrossCorrelationCalculator.h"
#include "treebuilders/DefaultCalculator.h"
#include "treebuilders/OperatorAdaptor.h"
#include "treebuilders/SplitAdaptor.h"
#include "treebuilders/TreeBuilder.h"
#include "treebuilders/grid.h"
#include "treebuilders/project.h"

#include "trees/BandWidth.h"
#include "trees/FunctionTreeVector.h"
#include "trees/CornerOperatorTree.h"

#include "utils/Printer.h"
#include "utils/Timer.h"
#include "utils/math_utils.h"

#include "treebuilders/TimeEvolution_CrossCorrelationCalculator.h"

#include <vector>

#include "trees/OperatorNode.h"


namespace mrcpp {


/** @brief A constructor for TimeEvolutionOperator class.
 *
 * @param[in] mra: MRA.
 * @param[in] prec: precision.
 * @param[in] time: the time moment (step).
 * @param[in] finest_scale: the operator constructed down to this scale.
 * @param[in] imaginary: defines the real (faulse) or imaginary (true) part of the semigroup.
 * @param[in] max_Jpower: maximum amount of power integrals used.
 *
 * @details Constructs either real or imaginary part of the Schrodinger semigroup at a given time moment.
 *
 *
 *
 */
template <int D>
TimeEvolutionOperator<D>::TimeEvolutionOperator
(const MultiResolutionAnalysis<D> &mra, double prec, double time, int finest_scale, bool imaginary, int max_Jpower)
    : ConvolutionOperator<D>(mra, mra.getRootScale(), -10)   //One can use ConvolutionOperator instead as well
{
    int oldlevel = Printer::setPrintLevel(0);
    this->setBuildPrec(prec);

    SchrodingerEvolution_CrossCorrelation cross_correlation(30, mra.getOrder(), mra.getScalingBasis().getScalingType() );
    this->cross_correlation = &cross_correlation;

    initialize(time, finest_scale, imaginary, max_Jpower);     //will go outside of the constructor

    this->initOperExp(1);   //this turns out to be important
    Printer::setPrintLevel(oldlevel);
}


/** @brief A constructor for TimeEvolutionOperator class.
 *
 * @param[in] mra: MRA.
 * @param[in] prec: precision.
 * @param[in] time: the time moment (step).
 * @param[in] imaginary: defines the real (faulse) or imaginary (true) part of the semigroup.
 * @param[in] max_Jpower: maximum amount of power integrals used.
 *
 * @details Constructs either real or imaginary part of the Schrodinger semigroup at a given time moment.
 * 
 * 
 * 
 */
template <int D>
TimeEvolutionOperator<D>::TimeEvolutionOperator
(const MultiResolutionAnalysis<D> &mra, double prec, double time, bool imaginary, int max_Jpower)
    : ConvolutionOperator<D>(mra, mra.getRootScale(), -10)   //One can use ConvolutionOperator instead as well
{
    int oldlevel = Printer::setPrintLevel(0);
    this->setBuildPrec(prec);
    
    SchrodingerEvolution_CrossCorrelation cross_correlation(30, mra.getOrder(), mra.getScalingBasis().getScalingType() );
    this->cross_correlation = &cross_correlation;

    initialize(time, imaginary, max_Jpower);     //will go outside of the constructor

    this->initOperExp(1);   //this turns out to be important 
    Printer::setPrintLevel(oldlevel);
}



/** @brief Creates Re or Im of operator
 *
 * @details Uniform down to finest scale so far... (in progress)
 *
 *
 *
 *
 *
 *
 */
template <int D>
void TimeEvolutionOperator<D>::initialize(double time, bool imaginary, int max_Jpower)
{
    double o_prec = this->build_prec;
    auto o_mra = this->getOperatorMRA();

    // Setup uniform tree builder
    mrcpp::TreeBuilder<2> builder;
    mrcpp::SplitAdaptor<2> uniform(o_mra.getMaxScale(), true);

    int N = 19;

    auto o_tree = std::make_unique<CornerOperatorTree>(o_mra, o_prec);
    DefaultCalculator<2> intitial_calculator;
    //for (auto n = 0; n < 6; n++) builder.build(*o_tree, intitial_calculator, uniform, 1);
    //o_tree->setZero();

    double threshold = o_prec / 1000.0;
    std::map<int, mrcpp::JpowerIntegrals *> J;
    for( int n = 0; n <= N+1; n ++ )
        J[n] = new mrcpp::JpowerIntegrals(time * std::pow(4, n), n, max_Jpower, threshold);
    mrcpp::TimeEvolution_CrossCorrelationCalculator calculator(J, this->cross_correlation, imaginary);

    OperatorAdaptor adaptor(o_prec, o_mra.getMaxScale(), true); // Splits all nodes
    builder.build(*o_tree, calculator, adaptor, 12); // Expand 1D kernel into 2D operator
//    builder.build(*o_tree, calculator, adaptor, 1); // Expand 1D kernel into 2D operator
    // Postprocess to make the operator functional
    Timer trans_t;
    o_tree->mwTransform(BottomUp);

    std::cout << "Removing rubbish" << std::endl;
        
    o_tree->removeRubbish();

    o_tree->mwTransform(BottomUp);
    o_tree->clearSquareNorm();
    o_tree->calcSquareNorm();
    o_tree->setupOperNodeCache();

    print::time(10, "Time transform", trans_t);
    print::separator(10, ' ');

    this->raw_exp.push_back(std::move(o_tree));

    for( int n = 0; n <= N+1; n ++ )
        delete J[n];
}

/** @brief Creates Re or Im of operator
 *
 * @details Uniform down to finest scale so far... (in progress)
 *
 * 
 * 
 * 
 * 
 *  
 */
template <int D>
void TimeEvolutionOperator<D>::initialize(double time, int finest_scale, bool imaginary, int max_Jpower)
{
    double o_prec = this->build_prec;
    auto o_mra = this->getOperatorMRA();

    // Setup uniform tree builder
    mrcpp::TreeBuilder<2> builder;
    mrcpp::SplitAdaptor<2> uniform(o_mra.getMaxScale(), true);

    int N = finest_scale;
    double a = time * std::pow(4, N + 1);
    double threshold = o_prec / 1000.0;
    std::map<int, mrcpp::JpowerIntegrals *> J;
    for( int n = 0; n <= N+1; n ++ )
        J[n] = new mrcpp::JpowerIntegrals(time * std::pow(4, n), n, max_Jpower, threshold);
    mrcpp::TimeEvolution_CrossCorrelationCalculator calculator(J, this->cross_correlation, imaginary);

    auto o_tree = std::make_unique<OperatorTree>(o_mra, o_prec);
    DefaultCalculator<2> intitial_calculator;
    for (auto n = 0; n < 6; n++) builder.build(*o_tree, intitial_calculator, uniform, 1);

    double threshold = o_prec / 100.0;
    std::map<int, mrcpp::JpowerIntegrals *> J;
    for( int n = 0; n <= N+1; n ++ )
        J[n] = new mrcpp::JpowerIntegrals(time * std::pow(4, n), n, max_Jpower, threshold);
    mrcpp::TimeEvolution_CrossCorrelationCalculator calculator(J, this->cross_correlation, imaginary);
//    mrcpp::TimeEvolution_CrossCorrelationCalculator Im_calculator(J, this->cross_correlation, true);

    OperatorAdaptor adaptor(o_prec, o_mra.getMaxScale()); // Splits all nodes
    builder.build(*o_tree, calculator, adaptor, 12); // Expand 1D kernel into 2D operator

    // Postprocess to make the operator functional
    Timer trans_t;
    o_tree->mwTransform(mrcpp::BottomUp);
    o_tree->calcSquareNorm();
    o_tree->setupOperNodeCache();
    print::time(10, "Time transform", trans_t);
    print::separator(10, ' ');

    this->raw_exp.push_back(std::move(o_tree));

    for( int n = 0; n <= N+1; n ++ )
        delete J[n];
}

/** @brief SHOULD Reduce the precision of the tree by deleting nodes
 *
 * @param prec: New precision criterion
 * @param splitFac: Splitting factor: 1, 2 or 3
 * @param absPrec: Use absolute precision
 *
 * @details NOT implemented
 * This will run the tree building algorithm in "reverse", starting
 * from the leaf nodes, and perform split checks on each node based on the given
 * precision and the local wavelet norm.
 *
 * @note The splitting factor appears in the threshold for the wavelet norm as
 * \f$ ||w|| < 2^{-sn/2} ||f|| \epsilon \f$. In principal, `s` should be equal
 * to the dimension; in practice, it is set to `s=1`.
 *
 *
 */
/*
template <int D> int TimeEvolutionOperator<D>::crop(double prec, double splitFac, bool absPrec) {
    for (int i = 0; i < this->rootBox.size(); i++) {
        MWNode<D> &root = this->getRootMWNode(i);
        root.crop(prec, splitFac, absPrec);
    }
    int nChunks = this->getNodeAllocator().compress();
    this->resetEndNodeTable();
    this->calcSquareNorm();
    return nChunks;
}
*/




template class TimeEvolutionOperator<1>;
template class TimeEvolutionOperator<2>;
template class TimeEvolutionOperator<3>;

} // namespace mrcpp
