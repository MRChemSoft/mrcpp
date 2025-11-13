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

#include "MultiResolutionAnalysis.h"
#include "core/FilterCache.h"
#include "core/InterpolatingBasis.h"
#include "core/LegendreBasis.h"
#include "utils/Printer.h"
#include "utils/math_utils.h"

namespace mrcpp {

template <int D>
MultiResolutionAnalysis<D>::MultiResolutionAnalysis(std::array<int, 2> bb, int order, int depth)
        : maxDepth(depth)
        , basis(InterpolatingBasis(order))
        , world(bb) {
    if (getMaxDepth() > MaxDepth) MSG_ABORT("Beyond MaxDepth");
    if (getMaxScale() > MaxScale) MSG_ABORT("Beyond MaxScale");
    setupFilter();
}

template <int D>
MultiResolutionAnalysis<D>::MultiResolutionAnalysis(const BoundingBox<D> &bb, int order, int depth)
        : maxDepth(depth)
        , basis(InterpolatingBasis(order))
        , world(bb) {
    if (getMaxDepth() > MaxDepth) MSG_ABORT("Beyond MaxDepth");
    if (getMaxScale() > MaxScale) MSG_ABORT("Beyond MaxScale");
    setupFilter();
}

template <int D>
MultiResolutionAnalysis<D>::MultiResolutionAnalysis(const MultiResolutionAnalysis<D> &mra)
        : maxDepth(mra.maxDepth)
        , basis(mra.basis)
        , world(mra.world) {
    if (getMaxDepth() > MaxDepth) MSG_ABORT("Beyond MaxDepth");
    if (getMaxScale() > MaxScale) MSG_ABORT("Beyond MaxScale");
    setupFilter();
}

template <int D>
MultiResolutionAnalysis<D>::MultiResolutionAnalysis(const BoundingBox<D> &bb, const ScalingBasis &sb, int depth)
        : maxDepth(depth)
        , basis(sb)
        , world(bb) {
    if (getMaxDepth() > MaxDepth) MSG_ABORT("Beyond MaxDepth");
    if (getMaxScale() > MaxScale) MSG_ABORT("Beyond MaxScale");
    setupFilter();
}

template <int D> bool MultiResolutionAnalysis<D>::operator==(const MultiResolutionAnalysis<D> &mra) const {
    if (this->basis != mra.basis) return false;
    if (this->world != mra.world) return false;
    if (this->maxDepth != mra.maxDepth) return false;
    return true;
}

template <int D> bool MultiResolutionAnalysis<D>::operator!=(const MultiResolutionAnalysis<D> &mra) const {
    if (this->basis != mra.basis) return true;
    if (this->world != mra.world)
        std::cout << "diff world " << this->world << std::endl
                  << "and  "
                  << " " << mra.world << std::endl;
    if (this->world != mra.world) return true;
    if (this->maxDepth != mra.maxDepth) return true;
    return false;
}

template <int D> void MultiResolutionAnalysis<D>::print() const {
    print::separator(0, ' ');
    print::header(0, "MultiResolution Analysis");
    println(0, this->basis);
    print::separator(0, '-');
    println(0, this->world);
    print::separator(0, '=', 2);
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

template <int D> double MultiResolutionAnalysis<D>::calcMaxDistance() const {
    const Coord<D> &lb = getWorldBox().getLowerBounds();
    const Coord<D> &ub = getWorldBox().getUpperBounds();
    return math_utils::calc_distance<D>(lb, ub);
}

template class MultiResolutionAnalysis<1>;
template class MultiResolutionAnalysis<2>;
template class MultiResolutionAnalysis<3>;

} // namespace mrcpp
