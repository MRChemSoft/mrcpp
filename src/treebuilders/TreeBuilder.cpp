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

#include "TreeBuilder.h"
#include "TreeAdaptor.h"
#include "TreeCalculator.h"
#include "trees/MWNode.h"
#include "trees/MWTree.h"
#include "trees/NodeAllocator.h"
#include "utils/Printer.h"
#include "utils/Timer.h"
#include "utils/tree_utils.h"

namespace mrcpp {

template <int D, typename T>
void TreeBuilder<D, T>::build(MWTree<D, T> &tree,
                              TreeCalculator<D, T> &calculator,
                              TreeAdaptor<D, T> &adaptor,
                              int maxIter) const {
    Timer calc_t(false), split_t(false), norm_t(false);
    println(10, " == Building tree");

    MWNodeVector<D, T> *newVec = nullptr;
    MWNodeVector<D, T> *workVec = calculator.getInitialWorkVector(tree);

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
        if (iter == 0) sNorm = calcScalingNorm(*workVec);
        wNorm += calcWaveletNorm(*workVec);

        if (sNorm < 0.0 || wNorm < 0.0) {
            tree.squareNorm = -1.0;
        } else {
            tree.squareNorm = sNorm + wNorm;
        }
        println(10, std::setw(24) << tree.squareNorm);
        norm_t.stop();

        split_t.resume();
        newVec = new MWNodeVector<D, T>;
        if (iter >= maxIter && maxIter >= 0) {
            workVec->clear();
        }
        adaptor.splitNodeVector(*newVec, *workVec);
        split_t.stop();

        delete workVec;
        workVec = newVec;
        iter++;
    }

    tree.resetEndNodeTable();
    delete workVec;

    print::separator(10, ' ');
    print::time(10, "Time calc",  calc_t);
    print::time(10, "Time norm",  norm_t);
    print::time(10, "Time split", split_t);
}

template <int D, typename T>
void TreeBuilder<D, T>::clear(MWTree<D, T> &tree, TreeCalculator<D, T> &calculator) const {
    println(10, " == Clearing tree");

    Timer clean_t;
    MWNodeVector<D, T> nodeVec;
    tree_utils::make_node_table(tree, nodeVec);
    calculator.calcNodeVector(nodeVec);
    clean_t.stop();

    tree.clearSquareNorm();

    println(10, "  -- #  1: Cleared      " << std::setw(6) << nodeVec.size() << " nodes");
    print::separator(10, ' ');
    print::time(10, "Time clean", clean_t);
    print::separator(10, ' ');
}

template <int D, typename T>
int TreeBuilder<D, T>::split(MWTree<D, T> &tree, TreeAdaptor<D, T> &adaptor, bool passCoefs) const {
    println(10, " == Refining tree");

    Timer split_t;
    MWNodeVector<D, T> newVec;
    MWNodeVector<D, T> *workVec = tree.copyEndNodeTable();

    adaptor.splitNodeVector(newVec, *workVec);

    if (passCoefs) {
        for (int i = 0; i < workVec->size(); i++) {
            MWNode<D, T> &node = *(*workVec)[i];
            if (node.isBranchNode()) {
                node.giveChildrenCoefs(true);
            }
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

template <int D, typename T>
void TreeBuilder<D, T>::calc(MWTree<D, T> &tree, TreeCalculator<D, T> &calculator) const {
    println(10, " == Calculating tree");

    Timer calc_t;
    MWNodeVector<D, T> *workVec = calculator.getInitialWorkVector(tree);
    calculator.calcNodeVector(*workVec);
    printout(10, "  -- #" << std::setw(3) << 0 << ": Calculated ");
    printout(10, std::setw(6) << workVec->size() << " nodes ");
    delete workVec;
    calc_t.stop();

    tree.calcSquareNorm();

    print::separator(10, ' ');
    print::time(10, "Time calc", calc_t);
}

template <int D, typename T>
double TreeBuilder<D, T>::calcScalingNorm(const MWNodeVector<D, T> &vec) const {
    double sNorm = 0.0;
    for (int i = 0; i < vec.size(); i++) {
        const MWNode<D, T> &node = *vec[i];
        if (node.getDepth() >= 0) sNorm += node.getScalingNorm();
    }
    return sNorm;
}

template <int D, typename T>
double TreeBuilder<D, T>::calcWaveletNorm(const MWNodeVector<D, T> &vec) const {
    double wNorm = 0.0;
    for (int i = 0; i < vec.size(); i++) {
        const MWNode<D, T> &node = *vec[i];
        if (node.getDepth() >= 0) wNorm += node.getWaveletNorm();
    }
    return wNorm;
}

template class TreeBuilder<1, double>;
template class TreeBuilder<2, double>;
template class TreeBuilder<3, double>;

template class TreeBuilder<1, ComplexDouble>;
template class TreeBuilder<2, ComplexDouble>;
template class TreeBuilder<3, ComplexDouble>;

} // namespace mrcpp
