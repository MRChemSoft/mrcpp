#pragma once

#include "constants.h"
#include "mrcpp_declarations.h"

template<int D>
class GridGenerator {
public:
    GridGenerator(int ms = MaxScale) : maxScale(ms) { }
    virtual ~GridGenerator() { }

    int getMaxScale() const { return this->maxScale; }
    void setMaxScale(int ms) { this->maxScale = ms; }

    void operator()(FunctionTree<D> &out,
                    const RepresentableFunction<D> &inp,
                    int maxIter = -1) const;
    void operator()(FunctionTree<D> &out,
                    FunctionTree<D> &inp,
                    int max_iter = -1) const;
    void operator()(FunctionTree<D> &out,
                    FunctionTreeVector<D> &inp,
                    int max_iter = -1) const;
protected:
    int maxScale;
};

