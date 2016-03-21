#ifndef COPYADAPTOR_H
#define COPYADAPTOR_H

#include "TreeAdaptor.h"
#include "FunctionTreeVector.h"

template<int D>
class CopyAdaptor : public TreeAdaptor<D> {
public:
    CopyAdaptor(FunctionTree<D> &t) { tree_vec.push_back(t); }
    CopyAdaptor(FunctionTreeVector<D> &t) : tree_vec(t) { }
    CopyAdaptor(const CopyAdaptor<D> &a) : tree_vec(a.tree_vec) { }
    virtual ~CopyAdaptor() { }
    virtual TreeAdaptor<D> *copy() const { return new CopyAdaptor<D>(*this); }

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
