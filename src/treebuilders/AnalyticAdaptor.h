#pragma once

#include "TreeAdaptor.h"
#include "constants.h"

namespace mrcpp {

template<int D>
class AnalyticAdaptor final : public TreeAdaptor<D> {
public:
    AnalyticAdaptor(const RepresentableFunction<D> &f, int ms) : TreeAdaptor<D>(ms), func(&f) { }

private:
    const RepresentableFunction<D> *func;

    bool splitNode(const MWNode<D> &node) const {
        int scale = node.getScale();
        int nQuadPts = node.getKp1();
        if (this->func->isVisibleAtScale(scale, nQuadPts)) {
            return false;
        }
        double lb[3], ub[3];
        node.getBounds(lb, ub);
        if (this->func->isZeroOnInterval(lb, ub)) {
            return false;
        }
        return true;
    }
};

}
