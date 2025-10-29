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

void CornerOperatorTree::calcBandWidth(double prec) {
    if (this->bandWidth == nullptr) clearBandWidth();
    this->bandWidth = new BandWidth(getDepth());

    VectorXi max_transl;
    getMaxTranslations(max_transl);

    if (prec < 0.0) prec = this->normPrec;
    double thrs = std::max(MachinePrec, prec / 10.0);

    for (int depth = 0; depth < this->getDepth(); depth++) {
        int l = (1 << depth) - 1;
        this->bandWidth->setWidth(depth, 0, l);
        bool done = false;

        while (not done) {
            done = true;
            MWNode<2> *node = findNode(NodeIndex<2>(depth, {l, 0}));
            for (int k = 1; k < 4; k++) {
                if ((node != nullptr) && (node->getComponentNorm(k) > thrs)) {
                    this->bandWidth->setWidth(depth, k, l);
                    done = false;
                }
            }
            if (--l < 0) break;
        }
    }
    println(100, "\nOperator BandWidth" << *this->bandWidth);
}

bool CornerOperatorTree::isOutsideBand(int oTransl, int o_depth, int idx) {
    return std::abs(oTransl) < this->bandWidth->getWidth(o_depth, idx);
}

} // namespace mrcpp