#ifndef GRIDCLEANER_H
#define GRIDCLEANER_H

#include "TreeBuilder.h"
#include "WaveletAdaptor.h"
#include "DefaultCalculator.h"
#include "mrcpp_declarations.h"

template<int D>
class GridCleaner : public TreeBuilder<D> {
public:
    GridCleaner(double pr = -1.0) : prec(pr) { }
    virtual ~GridCleaner() { }

    void setPrecision(double pr) { this->prec = pr; }
    void multPrecision(double fac) { this->prec *= fac; }

    int operator()(MWTree<D> &out) {
        this->adaptor = new WaveletAdaptor<D>(this->prec, MaxScale);
        this->calculator = new DefaultCalculator<D>();
        int nSplit = clean(out);
        this->clearCalculator();
        this->clearAdaptor();
        return nSplit;
    }

protected:
    double prec;

    int clean(MWTree<D> &tree) const;
};

#endif // GRIDCLEANER_H
