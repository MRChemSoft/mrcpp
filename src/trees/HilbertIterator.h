#pragma once

#include "MWTree.h"
#include "MWNode.h"
#include "TreeIterator.h"
#include "HilbertPath.h"

namespace mrcpp {

template<int D>
class HilbertIterator final : public TreeIterator<D> {
public:
    HilbertIterator(MWTree<D> *tree, int dir = TopDown)
            : TreeIterator<D>(dir) {
        this->init(tree);
    }

protected:
    int getChildIndex(int i) const {
        const MWNode<D> &node = *this->state->node;
        const HilbertPath<D> &h = node.getHilbertPath();
        return h.getZIndex(i);
    }
};

}
