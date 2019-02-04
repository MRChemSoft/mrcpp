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

#include "MultiResolutionAnalysis.h"
#include "MRCPP/core/FilterCache.h"
#include "MRCPP/core/InterpolatingBasis.h"
#include "MRCPP/core/LegendreBasis.h"
#include "MRCPP/utils/Printer.h"

namespace mrcpp {

template <int D>
MultiResolutionAnalysis<D>::MultiResolutionAnalysis(const MultiResolutionAnalysis<D> &mra)
        : maxDepth(mra.maxDepth)
        , basis(mra.basis)
        , world(mra.world) {
    if (getMaxDepth() > MaxDepth) MSG_FATAL("Beyond MaxDepth");
    if (getMaxScale() > MaxScale) MSG_FATAL("Beyond MaxScale");
    setupFilter();
}

template <int D>
MultiResolutionAnalysis<D>::MultiResolutionAnalysis(const BoundingBox<D> &bb, const ScalingBasis &sb, int depth)
        : maxDepth(depth)
        , basis(sb)
        , world(bb) {
    if (getMaxDepth() > MaxDepth) MSG_FATAL("Beyond MaxDepth");
    if (getMaxScale() > MaxScale) MSG_FATAL("Beyond MaxScale");
    setupFilter();
}

template <int D> MultiResolutionAnalysis<1> MultiResolutionAnalysis<D>::getKernelMRA() const {
    const BoundingBox<D> &box = getWorldBox();
    const ScalingBasis &basis = getScalingBasis();

    int type = basis.getScalingType();
    int kern_order = 2 * basis.getScalingOrder() + 1;

    ScalingBasis *kern_basis = nullptr;
    if (type == Interpol) {
        kern_basis = new InterpolatingBasis(kern_order);
    } else if (type == Legendre) {
        kern_basis = new LegendreBasis(kern_order);
    } else {
        MSG_FATAL("Invalid scaling type");
    }

    int max_l = (box.isPeriodic()) ? 10 : 0;
    for (int i = 0; i < D; i++) {
        if (box.size(i) > max_l) { max_l = box.size(i); }
    }
    auto start_l = std::array<int, 1>{-max_l};
    auto tot_l = std::array<int, 1>{2 * max_l};
    // Zero in argument since operators are only implemented
    // for uniform scaling factor
    auto sf = std::array<double, 1>{box.getScalingFactor(0)};
    BoundingBox<1> kern_box(box.getScale(), start_l, tot_l, sf);
    MultiResolutionAnalysis<1> mra(kern_box, *kern_basis);
    delete kern_basis;
    return mra;
}

template <int D> MultiResolutionAnalysis<2> MultiResolutionAnalysis<D>::getOperatorMRA() const {
    const BoundingBox<D> &box = getWorldBox();
    const ScalingBasis &basis = getScalingBasis();

    int maxn = (box.isPeriodic()) ? 10 : 0;
    for (int i = 0; i < D; i++) {
        if (box.size(i) > maxn) { maxn = box.size(i); }
    }
    auto l = std::array<int, 2>{};
    auto nbox = std::array<int, 2>{maxn, maxn};
    // Zero in argument since operators are only implemented
    // for uniform scaling factor
    auto sf = std::array<double, 2>{box.getScalingFactor(0), box.getScalingFactor(0)};
    NodeIndex<2> idx(box.getScale());
    BoundingBox<2> oper_box(box.getScale(), l, nbox, sf);
    return MultiResolutionAnalysis<2>(oper_box, basis);
}

template <int D> bool MultiResolutionAnalysis<D>::operator==(const MultiResolutionAnalysis<D> &mra) const {
    if (this->basis != mra.basis) return false;
    if (this->world != mra.world) return false;
    if (this->maxDepth != mra.maxDepth) return false;
    return true;
}

template <int D> bool MultiResolutionAnalysis<D>::operator!=(const MultiResolutionAnalysis<D> &mra) const {
    if (this->basis != mra.basis) return true;
    if (this->world != mra.world) return true;
    if (this->maxDepth != mra.maxDepth) return true;
    return false;
}

template <int D> void MultiResolutionAnalysis<D>::print() const {
    println(0, std::endl);
    println(0, "============================================================");
    println(0, "                  MultiResolution Analysis                  ");
    println(0, "------------------------------------------------------------");
    println(0, this->basis);
    println(0, "------------------------------------------------------------");
    println(0, this->world);
    println(0, "============================================================");
    println(0, std::endl);
}

template <int D> void MultiResolutionAnalysis<D>::setupFilter() {
    getLegendreFilterCache(lfilters);
    getInterpolatingFilterCache(ifilters);
    int k = this->basis.getScalingOrder();
    int type = this->basis.getScalingType();
    switch (type) {
        case Legendre:
            this->filter = &lfilters.get(k);
            break;
        case Interpol:
            this->filter = &ifilters.get(k);
            break;
        default:
            MSG_ERROR("Invalid scaling basis selected.")
    }
}

template class MultiResolutionAnalysis<1>;
template class MultiResolutionAnalysis<2>;
template class MultiResolutionAnalysis<3>;

} // namespace mrcpp
