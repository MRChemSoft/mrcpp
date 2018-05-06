#pragma once

#include "TreeAdaptor.h"
#include "trees/FunctionTreeVector.h"

namespace mrcpp {

template<int D>
class CopyAdaptor : public TreeAdaptor<D> {
public:
    CopyAdaptor(FunctionTree<D> &t, int ms, int *bw);
    CopyAdaptor(FunctionTreeVector<D> &t, int ms, int *bw);
    virtual ~CopyAdaptor() { }

protected:
    int bandWidth[D];
    FunctionTreeVector<D> tree_vec;

    void setBandWidth(int *bw);
    bool splitNode(const MWNode<D> &node) const;
};

}
