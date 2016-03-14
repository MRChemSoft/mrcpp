#ifndef COPYADAPTOR_H
#define COPYADAPTOR_H

#include "TreeAdaptor.h"

template<int D>
class CopyAdaptor : public TreeAdaptor<D> {
public:
    CopyAdaptor(FunctionTree<D> &t) { trees.push_back(&t); }
    CopyAdaptor(std::vector<FunctionTree<D> *> &t) : trees(t) { }
    CopyAdaptor(const CopyAdaptor<D> &a) : trees(a.trees) { }
    virtual ~CopyAdaptor() { }
    virtual TreeAdaptor<D> *copy() const { return new CopyAdaptor<D>(*this); }

protected:
    std::vector<FunctionTree<D> *> trees;

    virtual bool splitNode(const MWNode<D> &node) const {
        const NodeIndex<D> &idx = node.getNodeIndex();
        for (int n = 0; n < this->trees.size(); n++) {
            const MWNode<D> *testNode = trees[n]->findNode(idx);
            if (testNode == 0) continue;
            if (testNode->isBranchNode()) return true;
        }
        return false;
    }
};

#endif // COPYADAPTOR_H
