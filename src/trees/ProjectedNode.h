#pragma once

#include "FunctionNode.h"

namespace mrcpp {

template<int D>
class ProjectedNode final : public FunctionNode<D> {
public:
    void createChildren();
    void genChildren();
    void deleteChildren();

    friend class SerialFunctionTree<D>;

protected:
    ProjectedNode() : FunctionNode<D>() { }
    ~ProjectedNode() { assert(this->tree == 0); }

    void dealloc();
    void reCompress();
};

}
