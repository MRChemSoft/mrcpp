#ifndef GRIDCLEANER_H
#define GRIDCLEANER_H

#include "WaveletAdaptor.h"
#include "DefaultCalculator.h"
#include "mrcpp_declarations.h"

template<int D>
class GridCleaner {
public:
    GridCleaner(double pr = -1.0, int max_scale = MaxScale)
        : prec(pr), maxScale(max_scale) { }
    virtual ~GridCleaner() { }

    double getPrecision() const { return this->prec; }
    void setPrecision(double pr) { this->prec = pr; }

    int operator()(MWTree<D> &out) {
        DefaultCalculator<D> calculator;
        WaveletAdaptor<D> adaptor(this->prec, this->maxScale);
        int nSplit = clean(out, calculator, adaptor);
        return nSplit;
    }

protected:
    double prec;
    int maxScale;

    int clean(MWTree<D> &tree,
              TreeCalculator<D> &calculator,
              TreeAdaptor<D> &adaptor) const;
};

#endif // GRIDCLEANER_H
