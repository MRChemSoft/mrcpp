#ifndef TREEBUILDER_H
#define TREEBUILDER_H

#include "mrcpp_declarations.h"

template<int D>
class TreeBuilder {
public:
    TreeBuilder(double pr, int max_scale)
        : prec(pr), maxScale(max_scale) { }
    virtual ~TreeBuilder() { }

    double getPrecision() const { return this->prec; }
    void setPrecision(double pr) { this->prec = pr; }

protected:
    double prec;
    int maxScale;

    double calcScalingNorm(const MWNodeVector &vec) const;
    double calcWaveletNorm(const MWNodeVector &vec) const;

    void build(MWTree<D> &tree,
               TreeCalculator<D> &calculator,
               TreeAdaptor<D> &adaptor,
               int maxIter) const;
};

#endif // TREEBUILDER_H
