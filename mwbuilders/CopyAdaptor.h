#ifndef COPYADAPTOR_H
#define COPYADAPTOR_H

#include "TreeAdaptor.h"
#include "FunctionTreeVector.h"
#include "constants.h"

template<int D>
class CopyAdaptor : public TreeAdaptor<D> {
public:
    CopyAdaptor(FunctionTree<D> &t, int ms = MaxScale)
            : TreeAdaptor<D>(ms) {
        tree_vec.push_back(&t);
    }
    CopyAdaptor(FunctionTreeVector<D> &t, int ms = MaxScale)
            : TreeAdaptor<D>(ms),
              tree_vec(t) {
    }
    virtual ~CopyAdaptor() {
    }

protected:
    FunctionTreeVector<D> tree_vec;

    virtual bool splitNode(const MWNode<D> &node) const {
        const NodeIndex<D> &idx = node.getNodeIndex();
        for (int i = 0; i < this->tree_vec.size(); i++) {
            const FunctionTree<D> &func_i = tree_vec.getFunc(i);
            const MWNode<D> *testNode = func_i.findNode(idx);
            if (testNode == 0) continue;
            if (testNode->isBranchNode()) return true;
        }
        return false;
    }
};

#endif // COPYADAPTOR_H
