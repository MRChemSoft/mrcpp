#pragma once

#include "MWTree.h"
#include "MWNode.h"
#include "TreeIterator.h"

namespace mrcpp {

template<int D>
class LebesgueIterator final : public TreeIterator<D> {
public:
    LebesgueIterator(MWTree<D> *tree, int dir = TopDown):
        TreeIterator<D>(dir) {
        this->init(tree);
    }

protected:
    int getChildIndex(int i) const { return i; }
};

}
