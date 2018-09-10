#include "MultiResolutionAnalysis.h"
#include "core/InterpolatingBasis.h"
#include "core/LegendreBasis.h"
#include "core/FilterCache.h"
#include "utils/Printer.h"

namespace mrcpp {

template<int D>
MultiResolutionAnalysis<D>::MultiResolutionAnalysis(const MultiResolutionAnalysis<D> &mra)
        : maxDepth(mra.maxDepth),
          basis(mra.basis),
          world(mra.world) {
    if (getMaxDepth() > MaxDepth) MSG_FATAL("Beyond MaxDepth");
    if (getMaxScale() > MaxScale) MSG_FATAL("Beyond MaxScale");
    setupFilter();
}

template<int D>
MultiResolutionAnalysis<D>::MultiResolutionAnalysis(const BoundingBox<D> &bb,
                        const ScalingBasis &sb, int depth)
        : maxDepth(depth),
          basis(sb),
          world(bb) {
    if (getMaxDepth() > MaxDepth) MSG_FATAL("Beyond MaxDepth");
    if (getMaxScale() > MaxScale) MSG_FATAL("Beyond MaxScale");
    setupFilter();
}

template<int D>
MultiResolutionAnalysis<1> MultiResolutionAnalysis<D>::getKernelMRA() const {
    const BoundingBox<D> &box = getWorldBox();
    const ScalingBasis &basis = getScalingBasis();

    int type = basis.getScalingType();
    int kern_order = 2*basis.getScalingOrder() + 1;

    ScalingBasis *kern_basis = 0;
    if (type == Interpol) {
        kern_basis = new InterpolatingBasis(kern_order);
    } else if (type == Legendre) {
        kern_basis = new LegendreBasis(kern_order);
    } else {
        MSG_FATAL("Invalid scaling type");
    }

    int max_l = 0;
    for (int i = 0; i < D; i++) {
        if (box.size(i) > max_l) {
            max_l = box.size(i);
        }
    }
    int start_l = -max_l;
    int tot_l = 2*max_l;
    NodeIndex<1> idx(box.getScale(), &start_l);
    BoundingBox<1> kern_box(idx, &tot_l);
    MultiResolutionAnalysis<1> mra(kern_box, *kern_basis);
    delete kern_basis;
    return mra;
}

template<int D>
MultiResolutionAnalysis<2> MultiResolutionAnalysis<D>::getOperatorMRA() const {
    const BoundingBox<D> &box = getWorldBox();
    const ScalingBasis &basis = getScalingBasis();

    int maxn = 0;
    for (int i = 0; i < D; i++) {
        if (box.size(i) > maxn) {
            maxn = box.size(i);
        }
    }
    int nbox[2] = { maxn, maxn};
    NodeIndex<2> idx(box.getScale());
    BoundingBox<2> oper_box(idx, nbox);
    return MultiResolutionAnalysis<2>(oper_box, basis);
}

template<int D>
bool MultiResolutionAnalysis<D>::operator==(const MultiResolutionAnalysis<D> &mra) const {
    if (this->basis != mra.basis) return false;
    if (this->world != mra.world) return false;
    if (this->maxDepth != mra.maxDepth) return false;
    return true;
}

template<int D>
bool MultiResolutionAnalysis<D>::operator!=(const MultiResolutionAnalysis<D> &mra) const {
    if (this->basis != mra.basis) return true;
    if (this->world != mra.world) return true;
    if (this->maxDepth != mra.maxDepth) return true;
    return false;
}

template<int D>
void MultiResolutionAnalysis<D>::print() const {
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

template<int D>
void MultiResolutionAnalysis<D>::setupFilter() {
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
