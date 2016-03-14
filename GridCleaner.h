#ifndef GRIDCLEANER_H
#define GRIDCLEANER_H

#include "TreeBuilder.h"

template<int D>
class GridCleaner : public TreeBuilder<D> {
public:
    GridCleaner(const MultiResolutionAnalysis<D> &mra);
    GridCleaner(const MultiResolutionAnalysis<D> &mra, const TreeAdaptor<D> &a);
    virtual ~GridCleaner();

    int operator()(MWTree<D> &out) const {
        return clean(out);
    }

protected:
    int clean(MWTree<D> &tree) const;
};

#endif // GRIDCLEANER_H
