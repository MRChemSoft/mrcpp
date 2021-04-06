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

#include "ConvolutionOperator.h"
#include "GreensKernel.h"
#include "functions/Gaussian.h"
#include "treebuilders/CrossCorrelationCalculator.h"
#include "treebuilders/OperatorAdaptor.h"
#include "treebuilders/TreeBuilder.h"
#include "treebuilders/grid.h"
#include "treebuilders/project.h"
#include "trees/BandWidth.h"
#include "trees/FunctionTreeVector.h"
#include "trees/OperatorTree.h"
#include "utils/Printer.h"
#include "utils/Timer.h"
#include "utils/math_utils.h"

using namespace Eigen;

namespace mrcpp {

template <int D>
ConvolutionOperator<D>::ConvolutionOperator(const MultiResolutionAnalysis<D> &mra, double pr)
        : MWOperator(mra.getOperatorMRA())
        , kern_mra(mra.getKernelMRA())
        , prec(pr) {
    if (this->prec < 0.0) MSG_ERROR("Negative build prec");
}

template <int D> ConvolutionOperator<D>::~ConvolutionOperator() {
    this->clearKernel();
}

template <int D> void ConvolutionOperator<D>::initializeOperator(GreensKernel &greens_kernel) {
    int max_scale = this->oper_mra.getMaxScale();

    TreeBuilder<2> builder;
    OperatorAdaptor adaptor(this->prec, max_scale);

    for (int i = 0; i < greens_kernel.size(); i++) {
        Gaussian<1> &k_func = *greens_kernel[i];
        auto *k_tree = new FunctionTree<1>(this->kern_mra);
        mrcpp::build_grid(*k_tree, k_func);               // Generate empty grid to hold narrow Gaussian
        mrcpp::project(this->prec / 10, *k_tree, k_func); // Project Gaussian starting from the empty grid
        CrossCorrelationCalculator calculator(*k_tree);

        auto *o_tree = new OperatorTree(this->oper_mra, this->prec);
        builder.build(*o_tree, calculator, adaptor, -1); // Expand 1D kernel into 2D operator

        Timer trans_t;
        o_tree->mwTransform(BottomUp);
        o_tree->calcSquareNorm();
        o_tree->setupOperNodeCache();
        print::time(10, "Time transform", trans_t);
        print::separator(10, ' ');

        this->kern_exp.push_back(std::make_tuple(1.0, k_tree));
        this->oper_exp.push_back(o_tree);
    }
}

template <int D> void ConvolutionOperator<D>::clearKernel() {
    // namespace explicitly needed to disambiguate...
    mrcpp::clear(this->kern_exp, true);
}

template <int D>
double ConvolutionOperator<D>::calcMinDistance(const MultiResolutionAnalysis<D> &MRA, double epsilon) const {
    int maxScale = MRA.getMaxScale();
    return std::sqrt(epsilon * std::pow(2.0, -maxScale));
}

template <int D> double ConvolutionOperator<D>::calcMaxDistance(const MultiResolutionAnalysis<D> &MRA) const {
    const Coord<D> &lb = MRA.getWorldBox().getLowerBounds();
    const Coord<D> &ub = MRA.getWorldBox().getUpperBounds();
    return math_utils::calc_distance<D>(lb, ub);
}

template class ConvolutionOperator<1>;
template class ConvolutionOperator<2>;
template class ConvolutionOperator<3>;

} // namespace mrcpp
