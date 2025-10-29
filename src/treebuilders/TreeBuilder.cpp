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
 * @file TreeBuilder.cpp
 * @brief Generic driver that orchestrates adaptive construction, refinement,
 *        and coefficient (re)calculation of multiwavelet trees.
 *
 * @details
 * A TreeBuilder manages the high-level loop:
 *  1) pick a work set of nodes (via a TreeCalculator-provided policy),
 *  2) compute coefficients on those nodes (calculator),
 *  3) estimate norms to drive thresholding,
 *  4) ask the TreeAdaptor where to split next,
 *  5) iterate until the work set is empty or a maximum iteration is reached.
 *
 * The builder never changes numerical kernels; it delegates all math to
 * a TreeCalculator and all grid-refinement policy to a TreeAdaptor.
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

/**
 * @brief Adaptive build of a tree using a calculator/adaptor pair.
 *
 * @param[in,out] tree       Target tree to be populated/refined.
 * @param[in,out] calculator Computes node coefficients & provides initial work set.
 * @param[in,out] adaptor    Decides which nodes to split next (refinement policy).
 * @param[in]     maxIter    Maximum refinement iterations; negative => unbounded.
 *
 * @details
 * Loop invariant:
 *  - `workVec` holds the nodes to be (re)computed at the current iteration.
 *  - After computing, the builder updates an approximate squared norm
 *    (scaling + wavelet) to drive relative thresholding elsewhere.
 *  - The adaptor produces the next `workVec` by splitting according to
 *    its policy. If `maxIter >= 0` and `iter >= maxIter`, splitting is
 *    disabled and the loop terminates after coefficients are computed.
 *
 * @note
 * The approximate norm written into `tree.squareNorm` is for thresholding and
 * progress reporting only. A precise norm is expected to be recomputed later
 * (e.g., after a bottom-up transform).
 */
template <int D, typename T>
void TreeBuilder<D, T>::build(MWTree<D, T> &tree,
                              TreeCalculator<D, T> &calculator,
                              TreeAdaptor<D, T> &adaptor,
                              int maxIter) const {
    Timer calc_t(false), split_t(false), norm_t(false);
    println(10, " == Building tree");

    MWNodeVector<D, T> *newVec = nullptr;
    MWNodeVector<D, T> *workVec = calculator.getInitialWorkVector(tree);

    double sNorm = 0.0;  // accumulated scaling contribution (approx.)
    double wNorm = 0.0;  // accumulated wavelet contribution (approx.)

    int iter = 0;
    while (workVec->size() > 0) {
        printout(10, "  -- #" << std::setw(3) << iter << ": Calculated ");
        printout(10, std::setw(6) << workVec->size() << " nodes ");

        // 1) Compute coefficients on current work set
        calc_t.resume();
        calculator.calcNodeVector(*workVec);
        calc_t.stop();

        // 2) Update approximate norms used for thresholding/progress only
        norm_t.resume();
        if (iter == 0) sNorm = calcScalingNorm(*workVec);
        wNorm += calcWaveletNorm(*workVec);

        if (sNorm < 0.0 || wNorm < 0.0) {
            // Propagate "unknown" / invalid norm
            tree.squareNorm = -1.0;
        } else {
            // Approximate norm (exact one will be recomputed later)
            tree.squareNorm = sNorm + wNorm;
        }
        println(10, std::setw(24) << tree.squareNorm);
        norm_t.stop();

        // 3) Decide and perform refinement for the next iteration
        split_t.resume();
        newVec = new MWNodeVector<D, T>;
        if (iter >= maxIter && maxIter >= 0) {
            // Respect iteration cap: stop splitting
            workVec->clear();
        }
        adaptor.splitNodeVector(*newVec, *workVec);
        split_t.stop();

        delete workVec;
        workVec = newVec;
        iter++;
    }

    // Invalidate cached end-node table because the grid changed
    tree.resetEndNodeTable();
    delete workVec;

    print::separator(10, ' ');
    print::time(10, "Time calc",  calc_t);
    print::time(10, "Time norm",  norm_t);
    print::time(10, "Time split", split_t);
}

/**
 * @brief Remove all coefficients from the tree (fixed grid), using the calculator
 *        to "clear" node data.
 *
 * @param[in,out] tree       Target MW tree.
 * @param[in,out] calculator Calculator invoked to clear coefficients for nodes.
 *
 * @details
 * - The grid topology is preserved.
 * - `tree.squareNorm` is reset.
 */
template <int D, typename T>
void TreeBuilder<D, T>::clear(MWTree<D, T> &tree, TreeCalculator<D, T> &calculator) const {
    println(10, " == Clearing tree");

    Timer clean_t;
    MWNodeVector<D, T> nodeVec;
    tree_utils::make_node_table(tree, nodeVec);
    calculator.calcNodeVector(nodeVec); // calculator is responsible for zeroing/clearing
    clean_t.stop();

    tree.clearSquareNorm();

    println(10, "  -- #  1: Cleared      " << std::setw(6) << nodeVec.size() << " nodes");
    print::separator(10, ' ');
    print::time(10, "Time clean", clean_t);
    print::separator(10, ' ');
}

/**
 * @brief Split (refine) the current leaf nodes according to an adaptor policy.
 *
 * @param[in,out] tree     Target tree to refine.
 * @param[in,out] adaptor  Adaptor that decides which nodes to split.
 * @param[in]     passCoefs If true, transfer parent coefficients to children
 *                          (preserving function representation).
 *
 * @return Number of newly created child nodes (i.e., number of splits * children).
 *
 * @details
 * - The end-node table is reset after refinement.
 * - If `passCoefs == true` and a refined node remains a branch node, the parent
 *   distributes its coefficients to the children (e.g., via projection / exact transfer).
 */
template <int D, typename T>
int TreeBuilder<D, T>::split(MWTree<D, T> &tree, TreeAdaptor<D, T> &adaptor, bool passCoefs) const {
    println(10, " == Refining tree");

    Timer split_t;
    MWNodeVector<D, T> newVec;             // newly created nodes (unused beyond counting)
    MWNodeVector<D, T> *workVec = tree.copyEndNodeTable();  // current leaves

    adaptor.splitNodeVector(newVec, *workVec);

    if (passCoefs) {
        for (int i = 0; i < workVec->size(); i++) {
            MWNode<D, T> &node = *(*workVec)[i];
            if (node.isBranchNode()) {
                // Transfer coefficients from parent to children
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

/**
 * @brief Recalculate coefficients on the calculator-provided work set
 *        without refinement.
 *
 * @param[in,out] tree       Target tree.
 * @param[in,out] calculator Calculator used to compute node coefficients.
 *
 * @details
 * Computes on the initial work vector (as defined by the calculator) and then
 * recomputes the exact squared norm of the tree.
 */
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

/**
 * @brief Sum of scaling contributions (approximate) across a vector of nodes.
 *
 * @param[in] vec Node vector from the current iteration.
 * @return Approximate sum of scaling norms for nodes with depth >= 0.
 */
template <int D, typename T>
double TreeBuilder<D, T>::calcScalingNorm(const MWNodeVector<D, T> &vec) const {
    double sNorm = 0.0;
    for (int i = 0; i < vec.size(); i++) {
        const MWNode<D, T> &node = *vec[i];
        if (node.getDepth() >= 0) sNorm += node.getScalingNorm();
    }
    return sNorm;
}

/**
 * @brief Sum of wavelet contributions (approximate) across a vector of nodes.
 *
 * @param[in] vec Node vector from the current iteration.
 * @return Approximate sum of wavelet norms for nodes with depth >= 0.
 */
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
