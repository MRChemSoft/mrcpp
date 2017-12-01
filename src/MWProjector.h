#pragma once

#include <functional>

#include "constants.h"
#include "mrcpp_declarations.h"

template<int D>
class MWProjector {
public:
    MWProjector(double pr = -1.0, int ms = MaxScale) : prec(pr), maxScale(ms) { }
    virtual ~MWProjector() { }

    double getPrecision() const { return this->prec; }
    int getMaxScale() const { return this->maxScale; }

    void setPrecision(double pr) { this->prec = pr; }
    void setMaxScale(int ms) { this->maxScale = ms; }

    void operator()(FunctionTree<D> &out,
                    std::function<double (const double *r)> func,
                    int maxIter = -1) const;
    void operator()(FunctionTree<D> &out,
                    RepresentableFunction<D> &inp,
                    int maxIter = -1) const;
protected:
    double prec;
    int maxScale;
};

