#pragma once

#include "TreeAdaptor.h"

namespace mrcpp {

template<int D>
class SplitAdaptor final : public TreeAdaptor<D> {
public:
    SplitAdaptor(int ms, bool sp) : TreeAdaptor<D>(ms), split(sp) { }

protected:
    bool split;

    bool splitNode(const MWNode<D> &node) const { return this->split; }
};

}
