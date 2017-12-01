#pragma once

#include "constants.h"
#include "mrcpp_declarations.h"

template<int D>
class MWDerivative {
public:
    MWDerivative(int ms = MaxScale) : maxScale(ms) { }
    virtual ~MWDerivative() { }

    int getMaxScale() const { return this->maxScale; }
    void setMaxScale(int ms) { this->maxScale = ms; }

    void operator()(FunctionTree<D> &out,
                    DerivativeOperator<D> &oper,
                    FunctionTree<D> &inp,
                    int dir = -1) const;
protected:
    int maxScale;
};

