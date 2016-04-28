#ifndef TREEADAPTOR_H
#define TREEADAPTOR_H

#include "mwrepr_declarations.h"
#include "TelePrompter.h"

template<int D>
class TreeAdaptor {
public:
    TreeAdaptor() { }
    TreeAdaptor(const TreeAdaptor<D> &adap) { }
    virtual ~TreeAdaptor() { }
    virtual TreeAdaptor<D> *copy() const { return new TreeAdaptor<D>(*this); }

    void splitNodeVector(MWNodeVector &out, MWNodeVector &inp) const {
        for (int n = 0; n < inp.size(); n++) {
            MWNode<D> &node = *inp[n];
            // Can be BranchNode in operator application
            if (node.isBranchNode()) continue;
            if (splitNode(node)) {
                node.createChildren();
                for (int i = 0; i < node.getNChildren(); i++) {
                    out.push_back(&node.getMWChild(i));
                }
            }
        }
    }

protected:
    virtual bool splitNode(const MWNode<D> &node) const { return false; }
};

#endif // TREEADAPTOR_H
