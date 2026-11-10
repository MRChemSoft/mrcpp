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

#include "CornerOperatorTree.h"
#include "BandWidth.h"
#include "OperatorNode.h"
#include "utils/Printer.h"

using namespace Eigen;

namespace mrcpp {

/** @brief Calculates band widths of the non-standard form matrices.
 *
 * @param[in] prec: Precision used for thresholding
 *
 * @details For each depth and each of \f$ T, A, B, C \f$ two distances to the diagonal are
 * recorded. The largest one with a considerable value is stored as the band width, which
 * ConvolutionCalculator takes as an upper bound on \f$ |l| \f$ as for any OperatorTree.
 * The smallest one is the edge below which the corner matrices vanish.
 *
 */
template <typename T> void CornerOperatorTree<T>::calcBandWidth(double prec) {
    this->clearBandWidth();
    this->bandWidth = new BandWidth(this->getDepth());
    this->innerWidth = MatrixXi::Constant(this->getDepth(), 4, -1);

    VectorXi max_transl;
    this->getMaxTranslations(max_transl);

    if (prec < 0.0) prec = this->normPrec;
    double thrs = std::max(MachinePrec, prec / 10.0); // should be enough due to oscillating behaviour of corner matrix elements (it's affected by polynomial order)

    for (int depth = 0; depth < this->getDepth(); depth++) {
        for (int l = 0; l <= max_transl[depth]; l++) {
            const MWNode<2, T> &node = this->getNode(depth, l);
            for (int k = 0; k < 4; k++) {
                if (node.getComponentNorm(k) > thrs) {
                    if (this->innerWidth(depth, k) < 0) this->innerWidth(depth, k) = l;
                    this->bandWidth->setWidth(depth, k, l);
                }
            }
        }
    }
    println(100, "\nOperator BandWidth" << *this->bandWidth);
}

/** @brief Checks if the distance to diagonal lies outside the band of considerable values.
 *
 * @param[in] oTransl: distance to diagonal
 * @param[in] o_depth: scaling order
 * @param[in] idx: index corresponding to one of the matrices \f$ A, B, C \f$ or \f$ T \f$.
 *
 * @returns True if \b oTransl is beyond the band width, or, for \f$ A, B, C \f$, closer to the
 * diagonal than their first considerable value; False otherwise.
 *
 */
template <typename T> bool CornerOperatorTree<T>::isOutsideBand(int oTransl, int o_depth, int idx) {
    int l = abs(oTransl);
    if (l > this->bandWidth->getWidth(o_depth, idx)) return true;
    // T is not screened near the diagonal: for D > 1 it is part of the mixed terms,
    // and the pure T term is already restricted to the root in ConvolutionCalculator::calcNode.
    if (idx == 0) return false;
    return l < this->innerWidth(o_depth, idx);
}

template class CornerOperatorTree<double>;
template class CornerOperatorTree<ComplexDouble>;

} // namespace mrcpp
