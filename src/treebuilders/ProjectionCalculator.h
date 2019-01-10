#pragma once

#include "TreeCalculator.h"

namespace mrcpp {

template<int D>
class ProjectionCalculator final : public TreeCalculator<D> {
public:
    ProjectionCalculator(const RepresentableFunction<D> &inp_func,
                         const std::array<double, D> &sf)
                         : func(&inp_func),
                           scaling_factor(sf) { }

private:
    const RepresentableFunction<D> *func;
    const std::array<double, D> scaling_factor;
    void calcNode(MWNode<D> &node);
};

}
