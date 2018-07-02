#pragma once

#include "TreeCalculator.h"

namespace mrcpp {

template<int D>
class ProjectionCalculator final : public TreeCalculator<D> {
public:
    ProjectionCalculator(const RepresentableFunction<D> &inp_func) : func(&inp_func) { }

protected:
    const RepresentableFunction<D> *func;

    virtual void calcNode(MWNode<D> &node);
};

}
