#pragma once

#include "mrcpp_declarations.h"

namespace mrcpp {

template<int D>
class TreeBuilder final {
public:
    void build(MWTree<D> &tree, TreeCalculator<D> &calculator, TreeAdaptor<D> &adaptor, int maxIter) const;
    void clear(MWTree<D> &tree, TreeCalculator<D> &calculator) const;
    void calc(MWTree<D> &tree, TreeCalculator<D> &calculator) const;
    int split(MWTree<D> &tree, TreeAdaptor<D> &adaptor, bool passCoefs) const;

private:
    double calcScalingNorm(const MWNodeVector<D> &vec) const;
    double calcWaveletNorm(const MWNodeVector<D> &vec) const;
};

}
