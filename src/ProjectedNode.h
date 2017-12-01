#pragma once

#include "FunctionNode.h"

template<int D>
class ProjectedNode: public FunctionNode<D> {
public:
    void createChildren();
    void genChildren();
    void deleteChildren();

    friend class SerialFunctionTree<D>;

protected:
    ProjectedNode() : FunctionNode<D>() { }
    virtual ~ProjectedNode() { assert(this->tree == 0); }

    void dealloc();
    void reCompress();
};

