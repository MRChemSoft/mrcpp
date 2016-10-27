#ifndef TREEBUILDER_H
#define TREEBUILDER_H

#include "mrcpp_declarations.h"

template<int D>
class TreeBuilder {
public:
    TreeBuilder();
    virtual ~TreeBuilder();

protected:
    TreeAdaptor<D> *adaptor;
    TreeCalculator<D> *calculator;

    void clearCalculator();
    void clearAdaptor();

    double calcScalingNorm(const MWNodeVector &vec) const;
    double calcWaveletNorm(const MWNodeVector &vec) const;

    void build(MWTree<D> &tree, int maxIter) const;
};

#endif // TREEBUILDER_H
