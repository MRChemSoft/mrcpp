/*
 * MRCPP, a numerical library based on multiresolution analysis and
 * the multiwavelet basis which provide low-scaling algorithms as well as
 * rigorous error control in numerical computations.
 * Copyright (C) 2019 Stig Rune Jensen, Jonas Juselius, Luca Frediani and contributors.
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

#include "MWOperator.h"
#include "trees/BandWidth.h"
#include "utils/Timer.h"
#include "utils/Printer.h"

using namespace Eigen;

namespace mrcpp {

void MWOperator::clear(bool dealloc) {
    if (dealloc) {
        for (int i = 0; i < this->oper_exp.size(); i++) {
            if (this->oper_exp[i] != 0) delete this->oper_exp[i];
        }
    }
    this->oper_exp.clear();
}

OperatorTree& MWOperator::getComponent(int i) {
    if (this->oper_exp[i] == 0) MSG_ERROR("Invalid component");
    if (i < 0 or i >= this->oper_exp.size()) MSG_ERROR("Out of bounds");
    return *this->oper_exp[i];
}

const OperatorTree& MWOperator::getComponent(int i) const {
    if (this->oper_exp[i] == 0) MSG_ERROR("Invalid component");
    if (i < 0 or i >= this->oper_exp.size()) MSG_ERROR("Out of bounds");
    return *this->oper_exp[i];
}

int MWOperator::getMaxBandWidth(int depth) const {
    int maxWidth = -1;
    if (depth < 0 ) {
        maxWidth = this->band_max.maxCoeff();
    } else if (depth < this->band_max.size() ) {
        maxWidth = this->band_max(depth);
    }
    return maxWidth;
}

void MWOperator::clearBandWidths() {
    for (unsigned int i = 0; i < this->oper_exp.size(); i++) {
        this->oper_exp[i]->clearBandWidth();
    }
}

void MWOperator::calcBandWidths(double prec) {
    int maxDepth = 0;
    // First compute BandWidths and find depth of the deepest component
    for (unsigned int i = 0; i < this->oper_exp.size(); i++) {
        OperatorTree &oTree = *this->oper_exp[i];
        oTree.calcBandWidth(prec);
        const BandWidth &bw = oTree.getBandWidth();
        int depth = bw.getDepth();
        if (depth > maxDepth) {
            maxDepth = depth;
        }
    }
    this->band_max = VectorXi(maxDepth + 1);
    this->band_max.setConstant(-1);
    // Find the largest effective bandwidth at each scale
    for (unsigned int i = 0; i < this->oper_exp.size(); i++) {
        const OperatorTree &oTree = *this->oper_exp[i];
        const BandWidth &bw = oTree.getBandWidth();
        for (int n = 0; n <= bw.getDepth(); n++) { // scale loop
            for (int j = 0; j < 4; j++) { //component loop
                int w = bw.getWidth(n, j);
                if (w > this->band_max(n)) {
                    this->band_max(n) = w;
                }
            }
        }
    }
    println(20, "  Maximum bandwidths:\n" << this->band_max << std::endl);
}

} // namespace mrcpp
