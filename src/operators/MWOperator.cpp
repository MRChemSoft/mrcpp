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

#include "MWOperator.h"
#include "trees/BandWidth.h"
#include "utils/Printer.h"
#include "utils/Timer.h"

using namespace Eigen;

namespace mrcpp {

template <int D>
void MWOperator<D>::initOperExp(int M) {
    if (this->raw_exp.size() < M) MSG_ABORT("Incompatible raw expansion");
    this->oper_exp.clear();
    for (int m = 0; m < M; m++) {
        std::array<OperatorTree *, D> otrees;
        otrees.fill(nullptr);
        this->oper_exp.push_back(otrees);
    }
    for (int i = 0; i < M; i++)
        for (int d = 0; d < D; d++) assign(i, d, this->raw_exp[i].get());
}

template <int D>
OperatorTree &MWOperator<D>::getComponent(int i, int d) {
    if (i < 0 or i >= this->oper_exp.size()) MSG_ERROR("Index out of bounds");
    if (d < 0 or d >= D) MSG_ERROR("Dimension out of bounds");
    if (this->oper_exp[i][d] == nullptr) MSG_ERROR("Invalid component");
    return *this->oper_exp[i][d];
}

template <int D>
const OperatorTree &MWOperator<D>::getComponent(int i, int d) const {
    if (i < 0 or i >= this->oper_exp.size()) MSG_ERROR("Index out of bounds");
    if (d < 0 or d >= D) MSG_ERROR("Dimension out of bounds");
    if (this->oper_exp[i][d] == nullptr) MSG_ERROR("Invalid component");
    return *this->oper_exp[i][d];
}

template <int D>
int MWOperator<D>::getMaxBandWidth(int depth) const {
    int maxWidth = -1;
    if (depth < 0) {
        maxWidth = *std::max_element(this->band_max.begin(), this->band_max.end());
    } else if (depth < this->band_max.size()) {
        maxWidth = this->band_max[depth];
    }
    return maxWidth;
}

template <int D>
void MWOperator<D>::clearBandWidths() {
    for (auto &i : this->oper_exp)
        for (int d = 0; d < D; d++) i[d]->clearBandWidth();
}

template <int D>
void MWOperator<D>::calcBandWidths(double prec) {
    int maxDepth = 0;
    for (auto &i : this->oper_exp) {
        for (int d = 0; d < D; d++) {
            OperatorTree &oTree = *i[d];
            oTree.calcBandWidth(prec);
            const BandWidth &bw = oTree.getBandWidth();
            int depth = bw.getDepth();
            if (depth > maxDepth) maxDepth = depth;
        }
    }
    this->band_max = std::vector<int>(maxDepth + 1, -1);

    for (auto &i : this->oper_exp) {
        for (int d = 0; d < D; d++) {
            const OperatorTree &oTree = *i[d];
            const BandWidth &bw = oTree.getBandWidth();
            for (int n = 0; n <= bw.getDepth(); n++) {
                for (int j = 0; j < 4; j++) {
                    int w = bw.getWidth(n, j);
                    if (w > this->band_max[n]) this->band_max[n] = w;
                }
            }
        }
    }
    println(20, "  Maximum bandwidths:");
    for (auto bw : this->band_max) println(20, bw);
    println(20, std::endl);
}

template <int D>
MultiResolutionAnalysis<2> MWOperator<D>::getOperatorMRA() const {
    const BoundingBox<D> &box = this->MRA.getWorldBox();
    const ScalingBasis &basis = this->MRA.getScalingBasis();

    int reach = this->oper_reach + 1;
    if (reach < 0) {
        for (int i = 0; i < D; i++) {
            if (box.size(i) > reach) reach = box.size(i);
        }
    }
    auto l = std::array<int, 2>{};
    auto nbox = std::array<int, 2>{reach, reach};
    auto sf = std::array<double, 2>{box.getScalingFactor(0), box.getScalingFactor(0)};

    BoundingBox<2> oper_box(this->oper_root, l, nbox, sf);
    auto oper_mra = MultiResolutionAnalysis<2>(oper_box, basis);
    return oper_mra;
}

template class MWOperator<1>;
template class MWOperator<2>;
template class MWOperator<3>;

} // namespace mrcpp