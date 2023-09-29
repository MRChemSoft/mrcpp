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
#include "treebuilders/OperatorAdaptor.h"
#include "treebuilders/SplitAdaptor.h"
#include "treebuilders/TreeBuilder.h"
#include "treebuilders/grid.h"
#include "treebuilders/project.h"

#include "trees/BandWidth.h"
#include "trees/FunctionTreeVector.h"
#include "trees/OperatorTree.h"

#include "utils/Printer.h"
#include "utils/Timer.h"
#include "utils/math_utils.h"

#include "treebuilders/TimeEvolution_CrossCorrelationCalculator.h"

#include <vector>


namespace mrcpp {

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





/** @brief Creates Re and Im of operator
 *
 * @details Uniform down to N = 3 so far... (in progress)
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
    double treshold = o_prec / 100.0;
    mrcpp::JpowerIntegrals J(a, N + 1, max_Jpower, treshold);
    mrcpp::TimeEvolution_CrossCorrelationCalculator calculator(J, this->cross_correlation, imaginary);
//    mrcpp::TimeEvolution_CrossCorrelationCalculator Im_calculator(J, this->cross_correlation, true);

    auto o_tree = std::make_unique<OperatorTree>(o_mra, o_prec);
    builder.build(*o_tree, calculator, uniform, N ); // Expand 1D kernel into 2D operator

    // Postprocess to make the operator functional
    Timer trans_t;
    o_tree->mwTransform(mrcpp::BottomUp);
    o_tree->calcSquareNorm();
    o_tree->setupOperNodeCache();
    //o_tree->cropTree(o_prec); //there is no method 'crop' in 'mrcpp::OperatorTree'
    print::time(10, "Time transform", trans_t);
    print::separator(10, ' ');

    this->raw_exp.push_back(std::move(o_tree));
}




template class TimeEvolutionOperator<1>;
template class TimeEvolutionOperator<2>;
template class TimeEvolutionOperator<3>;

} // namespace mrcpp
