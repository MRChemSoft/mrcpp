#ifndef COPYADAPTOR_H
#define COPYADAPTOR_H

#include "mwrepr_declarations.h"
#include "TreeAdaptor.h"

template<int D>
class CopyAdaptor : public TreeAdaptor<D> {
public:
    CopyAdaptor(const FunctionTree<D> &t)
        : tree(&t) {
    }

    CopyAdaptor(const CopyAdaptor<D> &a)
        : tree(a.tree) {
    }

    virtual ~CopyAdaptor() { }

    virtual TreeAdaptor<D> *copy() const { return new CopyAdaptor<D>(*this); }

protected:
    const FunctionTree<D> *tree;

    virtual bool splitNode(MWNode<D> &node) const {
        const NodeIndex<D> &idx = node.getNodeIndex();
        const MWNode<D> *testNode = tree->findNode(idx);
        if (testNode == 0) {
            return false;
        } else if (testNode->isEndNode()) {
            return false;
        }
        return true;
    }
};

#endif // COPYADAPTOR_H
