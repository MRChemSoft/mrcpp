#ifndef ANALYTICPROJECTOR_H
#define ANALYTICPROJECTOR_H

#include "TreeCalculator.h"

template<int D>
class ProjectionCalculator : public TreeCalculator<D> {
public:
    ProjectionCalculator(RepresentableFunction<D> &inp_func)
            : func(&inp_func) { }
    virtual ~ProjectionCalculator() { }
    virtual bool computesCoefs() const { return true; }

protected:
    RepresentableFunction<D> *func;

    virtual void calcNode(MWNode<D> &node) const;
};

#endif // ANALYTICPROJECTOR_H
