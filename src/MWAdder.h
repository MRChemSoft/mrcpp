#pragma once

#include "constants.h"
#include "mrcpp_declarations.h"

template<int D>
class MWAdder {
public:
    MWAdder(double pr = -1.0, int ms = MaxScale) : prec(pr), maxScale(ms) { }
    virtual ~MWAdder() { }

    double getPrecision() const { return this->prec; }
    int getMaxScale() const { return this->maxScale; }

    void setPrecision(double pr) { this->prec = pr; }
    void setMaxScale(int ms) { this->maxScale = ms; }

    void operator()(FunctionTree<D> &out,
                    double a, FunctionTree<D> &tree_a,
                    double b, FunctionTree<D> &tree_b,
                    int maxIter = -1) const;
    void operator()(FunctionTree<D> &out,
                    FunctionTreeVector<D> &inp,
                    int maxIter = -1) const;
protected:
    double prec;
    int maxScale;
};

