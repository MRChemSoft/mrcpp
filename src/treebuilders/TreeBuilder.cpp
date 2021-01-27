/*
 * MRCPP, a numerical library based on multiresolution analysis and
 * the multiwavelet basis which provide low-scaling algorithms as well as
 * rigorous error control in numerical computations.
 * Copyright (C) 2020 Stig Rune Jensen, Jonas Juselius, Luca Frediani and contributors.
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

#include "TreeBuilder.h"
#include "TreeAdaptor.h"
#include "TreeCalculator.h"
#include "trees/MWNode.h"
#include "trees/MWTree.h"
#include "utils/Printer.h"
#include "utils/Timer.h"
#include "utils/tree_utils.h"

namespace mrcpp {

template <int D>
void TreeBuilder<D>::build(MWTree<D> &tree, TreeCalculator<D> &calculator, TreeAdaptor<D> &adaptor, int maxIter) const {
    Timer calc_t(false), split_t(false), norm_t(false);
    println(10, " == Building tree");

    MWNodeVector<D> *newVec = nullptr;
    MWNodeVector<D> *workVec = calculator.getInitialWorkVector(tree);

    double sNorm = 0.0;
    double wNorm = 0.0;

    int iter = 0;
    while (workVec->size() > 0) {
        printout(10, "  -- #" << std::setw(3) << iter << ": Calculated ");
        printout(10, std::setw(6) << workVec->size() << " nodes ");
        calc_t.resume();
        calculator.calcNodeVector(*workVec);
        calc_t.stop();

        norm_t.resume();
        if (iter == 0) { sNorm = calcScalingNorm(*workVec); }
        wNorm += calcWaveletNorm(*workVec);

        if (sNorm < 0.0 or wNorm < 0.0) {
            tree.squareNorm = -1.0;
        } else {
            // approximate norm for thresholding only
            // exact norm is recomputed after mwTransform
            tree.squareNorm = sNorm + wNorm;
        }
        println(10, std::setw(24) << tree.squareNorm);
        norm_t.stop();

        split_t.resume();
        newVec = new MWNodeVector<D>;
        if (iter >= maxIter and maxIter >= 0) workVec->clear();
        adaptor.splitNodeVector(*newVec, *workVec);
        split_t.stop();

        delete workVec;
        workVec = newVec;
        iter++;
    }
    tree.resetEndNodeTable();
    delete workVec;

    print::separator(10, ' ');
    print::time(10, "Time calc", calc_t);
    print::time(10, "Time norm", norm_t);
    print::time(10, "Time split", split_t);
}

template <int D> void TreeBuilder<D>::clear(MWTree<D> &tree, TreeCalculator<D> &calculator) const {
    println(10, " == Clearing tree");

    Timer clean_t;
    MWNodeVector<D> nodeVec;
    tree_utils::make_node_table(tree, nodeVec);
    calculator.calcNodeVector(nodeVec); // clear all coefficients
    clean_t.stop();

    tree.clearSquareNorm();

    println(10, "  -- #  1: Cleared      " << std::setw(6) << nodeVec.size() << " nodes");
    print::separator(10, ' ');
    print::time(10, "Time clean", clean_t);
    print::separator(10, ' ');
}

template <int D> int TreeBuilder<D>::split(MWTree<D> &tree, TreeAdaptor<D> &adaptor, bool passCoefs) const {
    println(10, " == Refining tree");

    Timer split_t;
    MWNodeVector<D> newVec;
    MWNodeVector<D> *workVec = tree.copyEndNodeTable();
    adaptor.splitNodeVector(newVec, *workVec);
    if (passCoefs) {
        for (int i = 0; i < workVec->size(); i++) {
            MWNode<D> &node = *(*workVec)[i];
            if (node.isBranchNode()) { node.giveChildrenCoefs(true); }
        }
    }
    delete workVec;
    tree.resetEndNodeTable();
    split_t.stop();

    printout(10, "  -- #  0: Split        ");
    printout(10, std::setw(6) << newVec.size() << " nodes\n");

    print::separator(10, ' ');
    print::time(10, "Time split", split_t);
    print::separator(10, ' ');

    return newVec.size();
}

template <int D> void TreeBuilder<D>::calc(MWTree<D> &tree, TreeCalculator<D> &calculator) const {
    println(10, " == Calculating tree");

    Timer calc_t;
    MWNodeVector<D> *workVec = calculator.getInitialWorkVector(tree);
    calculator.calcNodeVector(*workVec);
    printout(10, "  -- #" << std::setw(3) << 0 << ": Calculated ");
    printout(10, std::setw(6) << workVec->size() << " nodes ");
    delete workVec;
    calc_t.stop();

    tree.calcSquareNorm();

    print::separator(10, ' ');
    print::time(10, "Time calc", calc_t);
}

template <int D> double TreeBuilder<D>::calcScalingNorm(const MWNodeVector<D> &vec) const {
    double sNorm = 0.0;
    for (int i = 0; i < vec.size(); i++) {
        const MWNode<D> &node = *vec[i];
        sNorm += node.getScalingNorm();
    }
    return sNorm;
}

template <int D> double TreeBuilder<D>::calcWaveletNorm(const MWNodeVector<D> &vec) const {
    double wNorm = 0.0;
    for (int i = 0; i < vec.size(); i++) {
        const MWNode<D> &node = *vec[i];
        wNorm += node.getWaveletNorm();
    }
    return wNorm;
}

template class TreeBuilder<1>;
template class TreeBuilder<2>;
template class TreeBuilder<3>;

} // namespace mrcpp
