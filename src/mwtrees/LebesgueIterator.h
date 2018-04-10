#pragma once

#include "mwtrees/MWTree.h"
#include "mwtrees/MWNode.h"
#include "mwtrees/TreeIterator.h"

namespace mrcpp {

template<int D>
class LebesgueIterator: public TreeIterator<D> {
public:
    LebesgueIterator(MWTree<D> *tree, int dir = TopDown):
        TreeIterator<D>(dir) {
        this->init(tree);
    }
    virtual ~LebesgueIterator() {}
protected:
    int getChildIndex(int i) const { return i; }
};

}
