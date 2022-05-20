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

template<int D>
MWOperator<D>::MWOperator(const MultiResolutionAnalysis<D> &mra, int root, int reach)
        : oper_root(root)
        , oper_reach(reach)
        , MRA(mra) {
    const BoundingBox<D> &box = this->MRA.getWorldBox();
    if (this->oper_root > box.getScale()) MSG_ABORT("Operator root cannot be higher than world root");
    if (!box.isPeriodic() && this->oper_root != box.getScale()) MSG_ABORT("Operator root must equal world root for non-periodic");
    if (this->oper_reach <= 0) {
        // Set default reach based on world size
        for (int d = 0; d < D; d++) this->oper_reach = std::max(this->oper_reach, box.size(d));
    }
}

template <int D>
OperatorTree &MWOperator<D>::getComponent(int i) {
    if (this->oper_exp[i] == nullptr) MSG_ERROR("Invalid component");
    if (i < 0 or i >= this->oper_exp.size()) MSG_ERROR("Out of bounds");
    return *this->oper_exp[i];
}

template <int D>
const OperatorTree &MWOperator<D>::getComponent(int i) const {
    if (this->oper_exp[i] == nullptr) MSG_ERROR("Invalid component");
    if (i < 0 or i >= this->oper_exp.size()) MSG_ERROR("Out of bounds");
    return *this->oper_exp[i];
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
    for (auto &i : this->oper_exp) i->clearBandWidth();
}

template <int D>
void MWOperator<D>::calcBandWidths(double prec) {
    int maxDepth = 0;
    // First compute BandWidths and find depth of the deepest component
    for (auto &i : this->oper_exp) {
        OperatorTree &oTree = *i;
        oTree.calcBandWidth(prec);
        const BandWidth &bw = oTree.getBandWidth();
        int depth = bw.getDepth();
        if (depth > maxDepth) maxDepth = depth;
    }
    this->band_max = std::vector<int>(maxDepth + 1, -1);

    // Find the largest effective bandwidth at each scale
    for (auto &i : this->oper_exp) {
        const OperatorTree &oTree = *i;
        const BandWidth &bw = oTree.getBandWidth();
        for (int n = 0; n <= bw.getDepth(); n++) { // scale loop
            for (int j = 0; j < 4; j++) {          // component loop
                int w = bw.getWidth(n, j);
                if (w > this->band_max[n]) this->band_max[n] = w;
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

    int o_size = this->oper_reach + 1;

    auto l = std::array<int, 2>{};
    auto nbox = std::array<int, 2>{o_size, o_size};
    // Zero in argument since operators are only implemented
    // for uniform scaling factor
    auto sf = std::array<double, 2>{box.getScalingFactor(0), box.getScalingFactor(0)};

    BoundingBox<2> oper_box(this->oper_root, l, nbox, sf);
    auto oper_mra = MultiResolutionAnalysis<2>(oper_box, basis);
    return oper_mra;
}

template class MWOperator<1>;
template class MWOperator<2>;
template class MWOperator<3>;

} // namespace mrcpp
