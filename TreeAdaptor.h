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

    MWNodeVector* splitNodeVector(MWNodeVector &nodeVec,
                                  MWNodeVector *no_split = 0) const {
        MWNodeVector *split = new MWNodeVector;
        int nNodes = nodeVec.size();
        for (int n = 0; n < nNodes; n++) {
            MWNode<D> &node = static_cast<MWNode<D> &>(*nodeVec[n]);
            if (splitNode(node)) {
                split->push_back(&node);
            } else if (no_split != 0) {
                no_split->push_back(&node);
            }
        }
        return split;
    }

protected:
    virtual bool splitNode(MWNode<D> &node) const { return false; }
};

#endif // TREEADAPTOR_H
