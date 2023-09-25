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
TimeEvolutionOperator<D>::TimeEvolutionOperator(const MultiResolutionAnalysis<D> &mra, double prec)
    : MWOperator<D>(mra, mra.getRootScale(), -10)   //One can use ConvolutionOperator instead as well
{
    int oldlevel = Printer::setPrintLevel(0);

//    std::cout << (this->cross_correlation == nullptr) << std::endl;

    this->setBuildPrec(prec);
    auto o_prec = prec;
    //auto k_prec = prec / 10.0;
    initialize(o_prec);
    this->initOperExp(1);   //this turns out to be important 

    Printer::setPrintLevel(oldlevel);
}

/*
template <int D>
TimeEvolutionOperator<D>::TimeEvolutionOperator(const MultiResolutionAnalysis<D> &mra, GaussExp<1> &kernel, double prec, int root, int reach)
        : MWOperator<D>(mra, root, reach) {
    int oldlevel = Printer::setPrintLevel(0);

    this->setBuildPrec(prec);
    auto o_prec = prec;
    auto k_prec = prec / 100.0;
    //    initialize(kernel, k_prec, o_prec);
    //    this->initOperExp(kernel.size());

    Printer::setPrintLevel(oldlevel);
}
*/



/** @brief Creates Re and Im of operator
 *
 * @details Uniform down to N = 3 so far... (in progress)
 * 
 */
template <int D>
void TimeEvolutionOperator<D>::initialize(double o_prec)
{


    //auto k_mra = this->getKernelMRA();
    auto o_mra = this->getOperatorMRA();

    //TreeBuilder<2> builder;
    //OperatorAdaptor adaptor(o_prec, o_mra.getMaxScale());

    // Initialize world in the unit cube [0,1]
    auto order = o_mra.getOrder();
    auto basis = mrcpp::LegendreBasis(order);
    auto world = mrcpp::BoundingBox<1>(0);
    auto MRA = mrcpp::MultiResolutionAnalysis<1>(world, basis, o_mra.getMaxScale());

    // Create dummy operator with a single OperatorTree and clear it
    //mrcpp::IdentityConvolution<1> O(MRA, o_prec);
    //mrcpp::OperatorTree &oTree = O.getComponent(0,0);
    //oTree.clear();

    // Setup uniform tree builder
    mrcpp::TreeBuilder<2> builder;
    mrcpp::SplitAdaptor<2> uniform(o_mra.getMaxScale(), true);

    // Create dummy kernel K(x,y) = 0 and use the ProjectionCalculator to project it into
    // the nodes of the OperatorTree. You need to replace this with your new calculator.
    // Create 1D kernel world with order 2k + 1 and domain [-1,1]
    auto kBasis = mrcpp::LegendreBasis(2*order + 1);
    auto kWorld = mrcpp::BoundingBox<1>(0, {-1}, {2});
    auto kMRA = mrcpp::MultiResolutionAnalysis<1>(kWorld, kBasis, o_mra.getMaxScale());
    mrcpp::FunctionTree<1> kTree(kMRA);

    // Project dummy 1D kernel K(x-y) = 0 and pass it to the CC calculator
    auto kFunc = [] (const mrcpp::Coord<1> &r) { return 0.0; };
    mrcpp::project<1>(o_prec/10, kTree, kFunc);
    mrcpp::TimeEvolution_CrossCorrelationCalculator calculator(kTree);

    int N = 3, M = 10;
    double a = 0.5;
    double treshold = 1.0e-15;
    mrcpp::JpowerIntegrals J(a, N, M, treshold);
    mrcpp::TimeEvolution_CrossCorrelationCalculator Re_calculator(J, this->cross_correlation, false);
    mrcpp::TimeEvolution_CrossCorrelationCalculator Im_calculator(J, this->cross_correlation, true);

    
    // Builds a uniform output tree up to scale N using the given calculator
    //builder.build(oTree, calculator, uniform, N);
    //or like this:
    auto o_tree = std::make_unique<OperatorTree>(o_mra, o_prec);
    builder.build(*o_tree, calculator, uniform, N); // Expand 1D kernel into 2D operator


    // Postprocess to make the operator functional
    Timer trans_t;
    o_tree->mwTransform(mrcpp::BottomUp);
    o_tree->calcSquareNorm();
    o_tree->setupOperNodeCache();
    //o_tree->cropTree(o_prec); //there is no method 'crop' in 'mrcpp::OperatorTree'
    print::time(10, "Time transform", trans_t);
    print::separator(10, ' ');



    //This part needs to be replaced with the code written by Evgueni
    /*    
    for (int i = 0; i < kernel.size(); i++) {
        // Rescale Gaussian for D-dim application
        auto *k_func = kernel.getFunc(i).copy();
        k_func->setCoef(std::pow(k_func->getCoef(), 1.0/D));

        FunctionTree<1> k_tree(k_mra);
        mrcpp::build_grid(k_tree, *k_func);    // Generate empty grid to hold narrow Gaussian
        mrcpp::project(k_prec, k_tree, *k_func); // Project Gaussian starting from the empty grid
        delete k_func;

        CrossCorrelationCalculator calculator(k_tree);
        auto o_tree = std::make_unique<OperatorTree>(o_mra, o_prec);
        builder.build(*o_tree, calculator, adaptor, -1); // Expand 1D kernel into 2D operator

        Timer trans_t;
        o_tree->mwTransform(BottomUp);
        o_tree->calcSquareNorm();
        o_tree->setupOperNodeCache();
        print::time(10, "Time transform", trans_t);
        print::separator(10, ' ');

        
    }
    */
        //this->raw_exp.push_back(std::move(&oTree));
    this->raw_exp.push_back(std::move(o_tree));



}



/*
// This one is most likely not necessary for the TimeEvolutionOperator //
template <int D>
MultiResolutionAnalysis<1> TimeEvolutionOperator<D>::getKernelMRA() const {
    const BoundingBox<D> &box = this->MRA.getWorldBox();
    const ScalingBasis &basis = this->MRA.getScalingBasis();

    int type = basis.getScalingType();
    int kern_order = 2 * basis.getScalingOrder() + 1;

    ScalingBasis *kern_basis = nullptr;
    if (type == Interpol) {
        kern_basis = new InterpolatingBasis(kern_order);
    } else if (type == Legendre) {
        kern_basis = new LegendreBasis(kern_order);
    } else {
        MSG_ABORT("Invalid scaling type");
    }

    int root = this->oper_root;
    int reach = this->oper_reach + 1;
    if (reach < 0) {
        for (int i = 0; i < D; i++) {
            if (box.size(i) > reach) reach = box.size(i);
        }
    }
    auto start_l = std::array<int, 1>{-reach};
    auto tot_l = std::array<int, 1>{2 * reach};
    // Zero in argument since operators are only implemented
    // for uniform scaling factor
    auto sf = std::array<double, 1>{box.getScalingFactor(0)};
    BoundingBox<1> kern_box(root, start_l, tot_l, sf);
    MultiResolutionAnalysis<1> kern_mra(kern_box, *kern_basis);
    delete kern_basis;
    return kern_mra;
}

template <int D>
std::vector<std::complex<double>> TimeEvolutionOperator<D>::calculate_J_power_inetgarls(int l, double a, int M)
{}

template <int D>
void TimeEvolutionOperator<D>::set_J_power_inetgarls()
{}
*/

template class TimeEvolutionOperator<1>;
template class TimeEvolutionOperator<2>;
template class TimeEvolutionOperator<3>;

} // namespace mrcpp
