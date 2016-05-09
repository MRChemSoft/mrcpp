#ifndef TREEBUILDER_H
#define TREEBUILDER_H

#include "mwrepr_declarations.h"

template<int D>
class TreeBuilder {
public:
    TreeBuilder(const MultiResolutionAnalysis<D> &mra);
    virtual ~TreeBuilder();

protected:
    TreeAdaptor<D> *adaptor;
    TreeCalculator<D> *calculator;
    const MultiResolutionAnalysis<D> MRA;

    void clearCalculator();
    void clearAdaptor();

    double calcScalingNorm(const MWNodeVector &vec) const;
    double calcWaveletNorm(const MWNodeVector &vec) const;

    void build(MWTree<D> &tree, int maxIter) const;
};

#endif // TREEBUILDER_H
