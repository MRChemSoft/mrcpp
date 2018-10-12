#pragma once

#include "trees/MWNode.h"
#include "mrcpp_declarations.h"

namespace mrcpp {

template<int D>
class TreeAdaptor {
public:
    TreeAdaptor(int ms) : maxScale(ms) { }
    virtual ~TreeAdaptor() = default;

    void setMaxScale(int ms) { this->maxScale = ms; }

    void splitNodeVector(MWNodeVector<D> &out, MWNodeVector<D> &inp) const {
        for (int n = 0; n < inp.size(); n++) {
            MWNode<D> &node = *inp[n];
            // Can be BranchNode in operator application
            if (node.isBranchNode()) continue;
            if (node.getScale() + 2 > this->maxScale) continue;
            if (splitNode(node)) {
                node.createChildren();
                for (int i = 0; i < node.getNChildren(); i++) {
                    out.push_back(&node.getMWChild(i));
                }
            }
        }
    }

protected:
    int maxScale;

    virtual bool splitNode(const MWNode<D> &node) const = 0;
};

}
