#pragma once

#include "constants.h"
#include "mrcpp_declarations.h"

template<int D>
class MWConvolution {
public:
    MWConvolution(double pr = -1.0, int ms = MaxScale, bool ap = false)
        : prec(pr), absPrec(ap), maxScale(ms) { }
    virtual ~MWConvolution() { }

    double getPrecision() const { return this->prec; }
    bool getAbsPrec() const { return this->absPrec; }
    int getMaxScale() const { return this->maxScale; }

    void setPrecision(double pr) { this->prec = pr; }
    void setAbsPrec(bool ap) { this->absPrec = ap; }
    void setMaxScale(int ms) { this->maxScale = ms; }

    void operator()(FunctionTree<D> &out,
                    ConvolutionOperator<D> &oper,
                    FunctionTree<D> &inp,
                    int max_iter = -1) const;
protected:
    double prec;
    bool absPrec;
    int maxScale;
};

