#pragma once

#include "mrcpp_declarations.h"

namespace mrcpp {

template<int D>
class TreeBuilder {
public:
    TreeBuilder() { }
    virtual ~TreeBuilder() { }

    void build(MWTree<D> &tree,
               TreeCalculator<D> &calculator,
               TreeAdaptor<D> &adaptor,
               int maxIter) const;

    int split(MWTree<D> &tree, TreeAdaptor<D> &adaptor, bool passCoefs) const;
    void calc(MWTree<D> &tree, TreeCalculator<D> &calculator) const;
    void clear(MWTree<D> &tree, TreeCalculator<D> &calculator) const;

protected:
    double calcScalingNorm(const MWNodeVector &vec) const;
    double calcWaveletNorm(const MWNodeVector &vec) const;
};

}
