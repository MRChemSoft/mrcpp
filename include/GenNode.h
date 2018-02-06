
/*
 *
 */

#pragma once

#include "FunctionNode.h"

namespace mrcpp {

template<int D>
class GenNode: public FunctionNode<D> {
public:
    double getWaveletNorm() const { return 0.0; }

    virtual void createChildren();
    virtual void genChildren();
    virtual void cvTransform(int kind);
    virtual void mwTransform(int kind);

    virtual void setValues(const Eigen::VectorXd &vec);
    virtual void getValues(Eigen::VectorXd &vec);

    friend class SerialFunctionTree<D>;

protected:
    GenNode() : FunctionNode<D>() { }
    virtual ~GenNode() { assert(this->tree == 0); }

    double calcComponentNorm(int i) const;
    virtual void dealloc();
    void reCompress();
};

}
