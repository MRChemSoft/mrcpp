#ifndef ANALYTICPROJECTOR_H
#define ANALYTICPROJECTOR_H

#include "TreeCalculator.h"

template<int D>
class ProjectionCalculator : public TreeCalculator<D> {
public:
    ProjectionCalculator(const RepresentableFunction<D> &inp_func) : func(&inp_func) { }
    virtual ~ProjectionCalculator() { }

protected:
    const RepresentableFunction<D> *func;

    virtual void calcNode(MWNode<D> &node);
};

#endif // ANALYTICPROJECTOR_H
